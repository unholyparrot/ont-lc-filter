import sys
import gzip
import argparse
import yaml

from tqdm import tqdm
from plotly.subplots import make_subplots
import plotly.graph_objects as go

from loguru import logger

import metrics

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=str, required=True)
    parser.add_argument("--interest_id", type=str, required=False)
    parser.add_argument("--interest_file", type=str, required=False)
    parser.add_argument("--out_pattern", type=str, default="out")
    parser.add_argument("--config", type=str, default="config.yaml")
    return parser.parse_args()

def make_plot(record, out_name, desired_q_value, min_q_change, kmer_length, window_size):

    record_id = record['id']
    norm_kmers_difference = record["metrics"]['norm_kmers_difference']
    norm_kmer_complexity = record["metrics"]['norm_kmer_complexity']
    vector_kmers_avg_qual = record["metrics"]['vector_kmers_avg_qual']
    vector_kmers_avg_qual_change = record["metrics"]['vector_kmers_avg_qual_change']
    failed_regions = record["metrics"]['failed_regions']

    fig = make_subplots(rows=4, cols=1, shared_xaxes=True)
    
    fig.add_trace(
        go.Scatter(y=norm_kmers_difference),
        row=1, col=1
    )
    fig.add_trace(
        go.Scatter(y=norm_kmer_complexity, ),
        row=2, col=1, 
    )
    fig.add_trace(
        go.Scatter(y=vector_kmers_avg_qual),
        row=3, col=1, 
    )
    fig.add_hrect(y0=0, y1=desired_q_value, row=3, col=1,
                 line_width=0, fillcolor="orange", opacity=0.2)
    fig.add_trace(
        go.Scatter(y=vector_kmers_avg_qual_change),
        row=4, col=1
    )
    fig.add_hrect(y0=-min_q_change, y1=min_q_change,
                 row=4, col=1, line_width=0, fillcolor="orange", opacity=0.2)
    
    for region in failed_regions:
        for i in [1, 2, 3, 4]:
            fig.add_vrect(x0=region[0], x1=region[1],
                          row=i, col=1, line_width=0, fillcolor="black", opacity=0.2)
    
    fig.update_layout(height=600, width=1800, title_text=f"READ ID {record_id}, k-{kmer_length}, W={window_size}")
    fig.update_yaxes(title_text=r"Outer", range=[-0.1, 2.1], row=1, col=1)
    fig.update_yaxes(title_text=r"Inner", range=[-0.1, 1.1], row=2, col=1)
    fig.update_yaxes(title_text=r"Q", range=[-0.1, 52], row=3, col=1)
    fig.update_yaxes(title_text=r"Q ch", row=4, col=1)
    
    fig.write_html(out_name)

def main():
    # parser args
    args = parse_args()

    # configure logger
    logger.remove()
    logger.add(sys.stdout, level="INFO")

    logger.info(f"Input: {args.input}")
    logger.info(f"Out pattern: {args.out_pattern}")
    logger.info(f"Config: {args.config}")

    # load config
    with open(args.config, "r") as f:
        config = yaml.safe_load(f)

    if args.interest_id:
        interest = {args.interest_id}
    elif args.interest_file:
        with open(args.interest_file, "r") as f:
            interest = {*f.read().splitlines()}
    else:
        raise ValueError("Either --interest_id or --interest_file must be provided")
    
    logger.info(f"Interest size: {len(interest)}")

    records = []

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
            if record['id'] in interest:
                records.append(record)

    logger.info(f"Found {len(records)} records")

    for record in tqdm(records):
        read_metrics = metrics.calculate_mertics(record, **config["metrics"])
        
        out_name = f"{args.out_pattern}_{record['id']}.html"
        make_plot(read_metrics, out_name, config["metrics"]["desired_q_value"], config["metrics"]["min_q_change"], config["metrics"]["kmer_length"], config["metrics"]["window_size"])

if __name__ == "__main__":
    main()
