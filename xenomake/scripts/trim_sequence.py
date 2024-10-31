import pysam
import subprocess
import os
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed

def process_chunk(chunk_fastq, trimmed_fastq):
    """Run cutadapt on a single chunk of reads."""
    cutadapt_cmd = [
        "cutadapt",
        "--detect-adapter-forwards",
        "-q", "20",
        "--trim-n",
        "--poly-a",
        "--poly-g",
        "--poly-t",
        "--poly-c",
        "-o", trimmed_fastq,
        chunk_fastq
    ]
    subprocess.run(cutadapt_cmd, check=True)

def trim_ubam_multithreaded(input_ubam, output_ubam, num_threads=4, tmp_dir="temp_chunks"):
    """
    Multithreaded trimming of adapters and poly-A/C/G/T tails from a UBAM file.
    """
    os.makedirs(tmp_dir, exist_ok=True)
    chunk_files = []
    trimmed_files = []

    # Step 1: Split the UBAM file into smaller FASTQ files (chunks)
    with pysam.AlignmentFile(input_ubam, "rb") as ubam:
        chunk_size = 10000  # Adjust chunk size as needed
        read_chunk = []
        chunk_count = 0

        for read in ubam:
            read_chunk.append(read)
            if len(read_chunk) >= chunk_size:
                chunk_fastq = os.path.join(tmp_dir, f"chunk_{chunk_count}.fastq")
                with open(chunk_fastq, "w") as fq_out:
                    for r in read_chunk:
                        fq_out.write(f"@{r.query_name}\n{r.query_sequence}\n+\n{''.join(chr(q + 33) for q in r.query_qualities)}\n")
                chunk_files.append(chunk_fastq)
                trimmed_files.append(os.path.join(tmp_dir, f"trimmed_{chunk_count}.fastq"))
                chunk_count += 1
                read_chunk = []

        # Last chunk
        if read_chunk:
            chunk_fastq = os.path.join(tmp_dir, f"chunk_{chunk_count}.fastq")
            with open(chunk_fastq, "w") as fq_out:
                for r in read_chunk:
                    fq_out.write(f"@{r.query_name}\n{r.query_sequence}\n+\n{''.join(chr(q + 33) for q in r.query_qualities)}\n")
            chunk_files.append(chunk_fastq)
            trimmed_files.append(os.path.join(tmp_dir, f"trimmed_{chunk_count}.fastq"))

    # Step 2: Run cutadapt in parallel on each chunk
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        futures = [executor.submit(process_chunk, chunk, trimmed) for chunk, trimmed in zip(chunk_files, trimmed_files)]
        for future in as_completed(futures):
            future.result()  # Wait for all trimming to complete

    # Step 3: Merge trimmed chunks back into the output UBAM file
    with pysam.AlignmentFile(output_ubam, "wb", header=pysam.AlignmentFile(input_ubam, "rb").header) as out_ubam:
        for trimmed_fastq, original_chunk in zip(trimmed_files, chunk_files):
            with pysam.FastxFile(trimmed_fastq) as trimmed_fq, pysam.AlignmentFile(input_ubam, "rb") as ubam:
                for original_read, trimmed_read in zip(ubam, trimmed_fq):
                    # Update the sequence and quality with trimmed values
                    original_read.query_sequence = trimmed_read.sequence
                    original_read.query_qualities = pysam.qualitystring_to_array(trimmed_read.quality)
                    
                    # Set the flag to unaligned
                    original_read.is_unmapped = True
                    
                    # Write the read to output UBAM with preserved tags
                    out_ubam.write(original_read)

    # Clean up temporary files
    for file in chunk_files + trimmed_files:
        os.remove(file)
    os.rmdir(tmp_dir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Multithreaded trimming of adapters and poly-A/C/G/T tails from a UBAM file.")
    parser.add_argument("input_ubam", type=str, help="Path to the input UBAM file with tagged reads.")
    parser.add_argument("output_ubam", type=str, help="Path to the output UBAM file with trimmed reads.")
    parser.add_argument("--num_threads", type=int, default=4, help="Number of threads to use for parallel processing (default: 4).")
    parser.add_argument("--tmp_dir", type=str, default="temp_chunks", help="Directory for temporary chunk files.")
    
    args = parser.parse_args()
    
    trim_ubam_multithreaded(args.input_ubam, args.output_ubam, args.num_threads, args.tmp_dir)
