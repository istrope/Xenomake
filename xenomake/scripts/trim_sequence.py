import subprocess
import os
import argparse
import logging
import time

def run_cutadapt(input_ubam, output_ubam, cores):
    """Run cutadapt directly on a UBAM file."""
    cutadapt_cmd = [
        "cutadapt",
        "-q", "20",
        "--trim-n",
        "--poly-a",
        f"-j {cores}",
        "-o", output_ubam,
        input_ubam
    ]
    try:
        start_time = time.time()
        subprocess.run(cutadapt_cmd, check=True)
        elapsed_time = time.time() - start_time
        logging.info(f"Cutadapt processing completed in {elapsed_time:.2f} seconds")
    except subprocess.CalledProcessError as e:
        logging.error(f"Cutadapt failed: {e}")
        raise

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)  # Set up logging to output to console
    parser = argparse.ArgumentParser(description="Trim adapters and poly-A/C/G/T tails from a UBAM file.")
    parser.add_argument("input_ubam", type=str, help="Path to the input UBAM file with tagged reads.")
    parser.add_argument("output_ubam", type=str, help="Path to the output UBAM file with trimmed reads.")
    parser.add_argument("--threads", type=int, default=1, help="Number of CPU cores to use.")

    args = parser.parse_args()
    
    run_cutadapt(args.input_ubam, args.output_ubam, args.threads)
