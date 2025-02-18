import sys
import gzip
import argparse
import yaml
from concurrent.futures import ProcessPoolExecutor, wait, ALL_COMPLETED
from multiprocessing import Manager
from queue import Empty
from tqdm import tqdm

from loguru import logger

import metrics

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=str, required=True)
    parser.add_argument("--out_pattern", type=str, default="out")
    parser.add_argument("--config", type=str, default="config.yaml")
    parser.add_argument("--num_processes", type=int, default=4)
    return parser.parse_args()

def process_sequence(record, queue_out, **kwargs):
    """Process a single sequence and put result in queue"""
    try:
        result = metrics.calculate_mertics(record, **kwargs)
        queue_out.put(result)
        return True
    except Exception as e:
        logger.error(f"Process error: {e}")
        return False

def main():
    # parser args
    args = parse_args()

    # configure logger
    logger.remove()
    logger.add(sys.stdout, level="INFO")
    logger.add(args.out_pattern + ".log", level="TRACE")

    # info about the parameters
    logger.info(f"Input: {args.input}")
    logger.info(f"Out pattern: {args.out_pattern}")
    logger.info(f"Config: {args.config}")

    # load config
    with open(args.config, "r") as f:
        config = yaml.safe_load(f)

    # info about the config
    logger.info(f"Config: {config}")

    # Create progress bars
    pbar_read = tqdm(desc="Reading sequences", unit=" reads")

    records = []
    
    try:
        # Read input file
        opener = gzip.open if str(args.input).endswith('.gz') else open
        with opener(args.input, "rt") as f:
            while True:
                header = f.readline().strip()
                if not header:  # EOF
                    break
                sequence = f.readline().strip()
                plus = f.readline().strip()
                quality = f.readline().strip()
                
                record = {
                    'id': header.split()[0].lstrip('@'),
                    'sequence': sequence,
                    'quality': quality
                }
                records.append(record)
                pbar_read.update(1)
        pbar_read.close()
        
        logger.info(f"Read {len(records)} records")

        with Manager() as manager:
            queue_output = manager.Queue()
        
            with ProcessPoolExecutor(max_workers=args.num_processes) as executor:
                pbar_process = tqdm(desc="Processing sequences", unit=" reads", total=len(records))
                    
                futures = [executor.submit(process_sequence, record, queue_output, **config["metrics"]) for record in records]
                for future in futures:
                    future.result()
                    pbar_process.update(1)
                pbar_process.close()
                
                logger.info(f"Processed {len(futures)} records, writing results to {args.out_pattern}.txt")
        
                # Process results and write to output file
                completed = 0
                total_records = len(futures)
                output_bed, output_table = args.out_pattern + ".bed", args.out_pattern + ".tsv"
                
                with open(output_bed, 'w') as f_bed, open(output_table, 'w') as f_table:
                    pbar_write = tqdm(desc="Writing results", unit=" records", total=total_records)
                    f_table.write("\t".join(["id", "length", "failed_regions", "failed_regions_length", "wanted_regions", "wanted_regions_length"]))
                    f_table.write("\n")
                    while completed < total_records:
                        try:
                            result = queue_output.get(timeout=0.1)
                        except Empty:
                            # Check if all futures are done
                            done, _ = wait(futures, timeout=0, return_when=ALL_COMPLETED)
                            if len(done) == len(futures):
                                if queue_output.empty():
                                    break
                            continue
                        
                        try:
                            # Write result to file
                            seq_id = result['id']
                            wanted_regions = metrics.get_wanted_regions(result["metrics"]["vector_kmers_avg_qual_change"].size, result["metrics"]["failed_regions"])
                            logger.trace(f"{seq_id} wanted regions: {wanted_regions}")

                            total_retained = sum([region[1] + config["metrics"]["kmer_length"] - region[0] for region in wanted_regions])
                            total_failed = result['length'] - total_retained

                            record_stats = [
                                seq_id, 
                                result['length'], 
                                len(result["metrics"]["failed_regions"]), total_failed, 
                                len(wanted_regions), total_retained
                            ]

                            bed_entries = []
                            for kstart, kend in wanted_regions:
                                bed_entries.append(f"{seq_id}\t{kstart}\t{kend + config['metrics']['kmer_length']}")
                            f_bed.write("\n".join(bed_entries))
                            f_bed.write("\n")
                            f_bed.flush()
                            
                            f_table.write("\t".join(list(map(str, record_stats))))
                            f_table.write("\n")
                            f_table.flush()
                            
                            completed += 1
                            pbar_write.update(1)
                        except Exception as e:
                            logger.error(f"Error: {e}")

                # Check for any exceptions in futures
                for future in futures:
                    future.result()
            
    except Exception as e:
        logger.error(f"Error: {e}")


if __name__ == "__main__":
    main()
