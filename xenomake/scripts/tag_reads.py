import pysam
import argparse
import time
import logging

def parse_range(range_str):
    start, end = map(int, range_str.split('-'))
    return start, end

def tag_pair(read1, read2, barcode_start, barcode_end, umi_start, umi_end):
    """
    Tags both reads in a pair with the same cell barcodes and UMIs from read1.
    """
    # Extract barcode and UMI from read1 sequence
    barcode = read1.query_sequence[barcode_start - 1:barcode_end]
    umi = read1.query_sequence[umi_start - 1:umi_end]
    
    # Tag both reads with the extracted barcode and UMI
    for read in [read1, read2]:
        read.set_tag("CB", barcode, value_type="Z")  # Cell Barcode
        read.set_tag("MI", umi, value_type="Z")      # UMI
    
    return read1, read2

def process_and_tag_reads(input_bam, output_bam, barcode_range, umi_range, chunk_size=10000, log_interval=100000):
    """
    Processes and tags paired-end reads with cell barcodes and UMIs from read1.
    """
    barcode_start, barcode_end = barcode_range
    umi_start, umi_end = umi_range
    total_reads = 0
    last_logged_reads = -log_interval  # Ensures the first log occurs at 0
    start_time = time.time()
    
    with pysam.AlignmentFile(input_bam, "rb", check_sq=False) as input_file:
        with pysam.AlignmentFile(output_bam, "wb", header=input_file.header) as output_file:
            read_buffer = {}
            chunk = []

            for read in input_file:
                qname = read.query_name


                if qname in read_buffer:
                    # If a paired read is already in the buffer, this is read2
                    read1 = read_buffer.pop(qname)
                    read1, read = tag_pair(read1, read, barcode_start, barcode_end, umi_start, umi_end)
                    # Add both reads to the chunk
                    chunk.extend([read1, read])
                    total_reads += 2

                else:
                    # Store read1 in the buffer if not already paired
                    read_buffer[qname] = read

                # Process chunk if it reaches the specified size
                if len(chunk) >= chunk_size:
                    for tagged_read in chunk:
                        output_file.write(tagged_read)
                    chunk.clear()  # Clear chunk after writing

                if total_reads - last_logged_reads >= log_interval:
                    last_logged_reads = total_reads
                    elapsed_time = time.time() - start_time
                    logging.info(f"Processed {total_reads} reads in {elapsed_time:.2f} seconds")

            # Write any remaining reads in the last chunk
            for tagged_read in chunk:
                output_file.write(tagged_read)

            # Final log for total reads processed
            elapsed_time = time.time() - start_time
            logging.info(f"Finished processing {total_reads} reads in {elapsed_time:.2f} seconds")

# Command-line interface
if __name__ == "__main__":
    print("Script started")
    parser = argparse.ArgumentParser(description="Tag barcodes and UMIs in a BAM file for paired-end reads")
    parser.add_argument("--input_bam", required=True, help="Path to input BAM/UBAM file")
    parser.add_argument("--output_bam", required=True, help="Path to output BAM file")
    parser.add_argument("--barcode_range", required=True, help="Range for cell barcode (e.g., '1-16')")
    parser.add_argument("--umi_range", required=True, help="Range for UMI (e.g., '17-28')")
    parser.add_argument("--chunk_size", type=int, default=10000, help="Number of reads per chunk")
    parser.add_argument("--log_interval", type=int, default=1000000, help="Number of reads between log messages")

    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

    barcode_range = parse_range(args.barcode_range)
    umi_range = parse_range(args.umi_range)
    
    process_and_tag_reads(args.input_bam, args.output_bam, barcode_range, umi_range, args.chunk_size, args.log_interval)
