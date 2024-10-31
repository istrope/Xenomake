import pysam
import argparse
from multiprocessing import Pool

def parse_range(range_str):
    """
    Parse a range string (e.g., "1-16") into start and end positions.
    
    Parameters:
    - range_str: String representing the range (e.g., "1-16").
    
    Returns:
    - Tuple with start and end positions as integers.
    """
    start, end = map(int, range_str.split('-'))
    return start, end

def tag_single_read(read, barcode_start, barcode_end, umi_start, umi_end):
    """
    Tags a single read with cell barcodes and UMIs.
    
    Parameters:
    - read: The read object.
    - barcode_start: Start position of the cell barcode in the read sequence.
    - barcode_end: End position of the cell barcode in the read sequence.
    - umi_start: Start position of the UMI in the read sequence.
    - umi_end: End position of the UMI in the read sequence.
    
    Returns:
    - The modified read with added tags.
    """
    if read.is_read1:  # Assuming barcodes and UMIs are in R1
        # Extract barcode and UMI from read sequence
        barcode = read.query_sequence[barcode_start-1:barcode_end]
        umi = read.query_sequence[umi_start-1:umi_end]
        
        # Add tags to the read
        read.set_tag("CB", barcode, value_type="Z")  # Cell Barcode
        read.set_tag("UB", umi, value_type="Z")  # UMI
    
    return read

def process_reads_chunk(chunk, barcode_start, barcode_end, umi_start, umi_end):
    """
    Processes a chunk of reads, tagging them with barcodes and UMIs.
    """
    tagged_reads = []
    for read in chunk:
        tagged_read = tag_single_read(read, barcode_start, barcode_end, umi_start, umi_end)
        tagged_reads.append(tagged_read)
    return tagged_reads

def tag_barcodes_and_umis(input_bam, output_bam, barcode_range, umi_range, threads):
    """
    Tags reads with cell barcodes and UMIs based on input ranges.
    
    Parameters:
    - input_bam: Path to input BAM/UBAM file containing paired reads.
    - output_bam: Path to output BAM file where tagged reads will be saved.
    - barcode_range: Tuple representing the start and end positions for the cell barcode.
    - umi_range: Tuple representing the start and end positions for the UMI.
    - threads: Number of threads to use for multithreading.
    """
    barcode_start, barcode_end = barcode_range
    umi_start, umi_end = umi_range
    
    input_file = pysam.AlignmentFile(input_bam, "rb", threads=threads)
    output_file = pysam.AlignmentFile(output_bam, "wb", header=input_file.header)

    chunk_size = 10000  # Number of reads per chunk
    chunk = []
    
    with Pool(processes=threads) as pool:
        for read in input_file:
            chunk.append(read)
            if len(chunk) >= chunk_size:
                # Process the chunk of reads in parallel
                results = pool.apply_async(process_reads_chunk, (chunk, barcode_start, barcode_end, umi_start, umi_end))
                for tagged_read in results.get():
                    output_file.write(tagged_read)
                chunk = []  # Reset chunk

        # Process the remaining reads in the final chunk
        if chunk:
            results = pool.apply_async(process_reads_chunk, (chunk, barcode_start, barcode_end, umi_start, umi_end))
            for tagged_read in results.get():
                output_file.write(tagged_read)

    input_file.close()
    output_file.close()

# Command-line interface
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Tag barcodes and UMIs in a BAM file")
    parser.add_argument("--input_bam", required=True, help="Path to input BAM/UBAM file")
    parser.add_argument("--output_bam", required=True, help="Path to output BAM file")
    parser.add_argument("--barcode_range", required=True, help="Range for cell barcode (e.g., '1-16')")
    parser.add_argument("--umi_range", required=True, help="Range for UMI (e.g., '17-28')")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads for parallel processing")

    args = parser.parse_args()

    # Parse input ranges for barcode and UMI
    barcode_range = parse_range(args.barcode_range)
    umi_range = parse_range(args.umi_range)
    
    # Run the tagging function
    tag_barcodes_and_umis(args.input_bam, args.output_bam, barcode_range, umi_range, args.threads)
