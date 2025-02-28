import sys
import gzip
import argparse
import yaml
from concurrent.futures import ProcessPoolExecutor, wait, ALL_COMPLETED
from multiprocessing import Manager
from queue import Empty
from tqdm import tqdm
from pathlib import Path
import os

from loguru import logger

import metrics

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=str, required=True)
    parser.add_argument("--out_pattern", type=str, default="out")
    parser.add_argument("--config", type=str, default=None, help='Path to config file')
    parser.add_argument("--num_processes", type=int, default=4)
    args = parser.parse_args()
    
    # Determine config path
    if args.config is not None:
        # Use user-provided config path
        config_path = args.config
    else:
        # Get the directory where the script is located
        if getattr(sys, 'frozen', False):
            # If running as a bundled executable
            script_dir = Path(sys.executable).parent
        else:
            # If running as a script
            script_dir = Path(os.path.dirname(os.path.abspath(__file__)))
        
        # Use default config from script directory
        config_path = os.path.join(script_dir, 'config.yaml')
    
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Config file not found at {config_path}")
    
    # Read the config
    with open(config_path, "r") as f:
        args.config = yaml.safe_load(f)
    
    return args

def process_sequence(record, queue_out, **kwargs):
    """Process a single sequence and put result in queue"""
    try:
        result = metrics.calculate_mertics(record, **kwargs)
        try:
            queue_out.put(result, timeout=5)  # Add timeout to prevent blocking forever
            return True
        except Exception as e:
            logger.error(f"Queue put error: {e}")
            return False
    except Exception as e:
        logger.error(f"Process error: {e}")
        return False

def process_in_batches(input_file, batch_size=1000, **kwargs):
    records = []
    count = 0
    opener = gzip.open if str(input_file).endswith('.gz') else open
    with opener(input_file, "rt") as f:
        while True:
            try:
                header = f.readline().strip()
                if not header:
                    break
                    
                sequence = f.readline().strip()
                plus_line = f.readline().strip()
                quality = f.readline().strip()
                
                if not all([sequence, plus_line, quality]):
                    logger.error("Incomplete FASTQ record")
                    break
                    
                record = {
                    'id': header.split()[0].lstrip('@'),
                    'sequence': sequence,
                    'quality': quality
                }
                
                if len(sequence) > kwargs["window_size"] + 10:
                    records.append(record)
                    count += 1
                    
                    if len(records) >= batch_size:
                        yield records
                        records = []
                        
            except Exception as e:
                logger.error(f"Error reading FASTQ: {e}")
                break
                
    if records:  # yield remaining records
        yield records

def main():
    args = parse_args()

    # configure logger
    logger.remove()
    logger.add(sys.stdout, level="INFO")
    logger.add(args.out_pattern + ".log", level="TRACE")

    # info about the parameters
    logger.info(f"Input: {args.input}")
    logger.info(f"Out pattern: {args.out_pattern}")
    logger.info(f"Config: {args.config}")

    manager = Manager()
    try:
        queue_output = manager.Queue()
        total_records = 0
        
        with ProcessPoolExecutor(max_workers=args.num_processes) as executor:
            pbar_process = tqdm(desc="Processing sequences", unit=" reads")
            total_processed = 0
            active_futures = set()
            
            try:
                # Process and write results simultaneously
                logger.info(f"Processing and writing results")
                completed = 0
                output_bed, output_table = args.out_pattern + ".bed", args.out_pattern + ".tsv"
                
                with open(output_bed, 'w') as f_bed, open(output_table, 'w') as f_table:
                    pbar_write = tqdm(desc="Processing and writing", unit=" records", total=total_records)
                    f_table.write("\t".join(["id", "length", "failed_regions", "failed_regions_length", "wanted_regions", "wanted_regions_length"]))
                    f_table.write("\n")
                    
                    # Process batches and write results immediately
                    try:
                        for batch in process_in_batches(args.input, batch_size=100, **args.config["metrics"]):
                            total_records += len(batch)
                            
                            # Clean up completed futures
                            active_futures = {f for f in active_futures if not f.done()}
                            
                            # Process available results
                            while True:
                                try:
                                    result = queue_output.get_nowait()
                                    # Write result immediately
                                    seq_id = result['id']
                                    wanted_regions = metrics.get_wanted_regions(
                                        result["metrics"]["vector_kmers_avg_qual_change"].size,
                                        result["metrics"]["failed_regions"]
                                    )
                                    
                                    total_retained = sum([region[1] + args.config["metrics"]["kmer_length"] - region[0] for region in wanted_regions])
                                    total_failed = result['length'] - total_retained

                                    record_stats = [
                                        seq_id, 
                                        result['length'], 
                                        len(result["metrics"]["failed_regions"]), 
                                        total_failed, 
                                        len(wanted_regions), 
                                        total_retained
                                    ]

                                    bed_entries = []
                                    for kstart, kend in wanted_regions:
                                        bed_entries.append(f"{seq_id}\t{kstart}\t{kend + args.config['metrics']['kmer_length']}")
                                    
                                    if bed_entries:
                                        f_bed.write("\n".join(bed_entries) + "\n")
                                        f_bed.flush()
                                    
                                    f_table.write("\t".join(map(str, record_stats)) + "\n")
                                    f_table.flush()
                                    
                                    completed += 1
                                    pbar_write.update(1)
                                except Empty:
                                    break
                                
                            # Submit new tasks
                            for record in batch:
                                future = executor.submit(process_sequence, record, queue_output, **args.config["metrics"])
                                active_futures.add(future)
                                
                    except Exception as e:
                        logger.error(f"Error during processing and writing: {e}")
                        raise
                    
                    pbar_write.close()

            except BrokenPipeError as e:
                logger.error(f"Broken pipe during processing: {e}")
                for future in active_futures:
                    future.cancel()
            except Exception as e:
                logger.error(f"Unexpected error during processing: {e}")
                raise

            # Wait for remaining tasks
            if active_futures:
                done, not_done = wait(active_futures, timeout=60)
                if not_done:
                    logger.warning(f"Cancelling {len(not_done)} incomplete tasks")
                    for future in not_done:
                        future.cancel()

            # Process remaining results
            logger.info(f"Processing remaining results ({total_records - total_processed} remaining)")
            while total_processed < total_records:
                try:
                    result = queue_output.get(timeout=1)
                    total_processed += 1
                    pbar_process.update(1)
                except Empty:
                    break
                except Exception as e:
                    logger.error(f"Error processing remaining results: {e}")
                    break

            pbar_process.close()

    except Exception as e:
        logger.error(f"Critical error: {e}")
        sys.exit(1)
    finally:
        try:
            manager.shutdown()
        except Exception as e:
            logger.error(f"Error shutting down manager: {e}")

if __name__ == "__main__":
    main()
