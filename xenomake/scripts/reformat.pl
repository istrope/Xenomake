#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use PerlIO::gzip;
"""
Taken from: https://github.com/MingyuYang-Yale/DBiT-seq/blob/master/Pre-processing/reformat.pl
Author: MingyuYang-Yale
This is a script to reformat DBIT-seq reads to match specifications for downstream processing with this pipeline
"""

my $usage=<<USAGE;
Usage:
    perl reformat.pl --read1 <path_to_read1> --read2 <path_to_read2> --outdir <output_directory>
Options:
    --read1     Path to the read1 input FASTQ file (gzip compressed)
    --read2     Path to the read2 input FASTQ file (gzip compressed)
    --outdir    Output directory for reformatted data
USAGE

my ($read1, $read2, $outdir, $help);
GetOptions(
    "read1=s"  => \$read1,
    "read2=s"  => \$read2,
    "outdir=s" => \$outdir,
    "help"     => \$help,
);

die $usage if $help;
die $usage unless $read1 && $read2 && $outdir;

# Create output directory if it doesn't exist
mkdir $outdir unless (-d $outdir);

# Open log file
open my $log, ">", "$outdir/log" or die "Cannot open log file: $!";
print $log "Filter primer begin at: " . `date`;

# Open input FASTQ files and main output files for filtered reads
open my $in1, "<:gzip", $read1 or die "Cannot open $read1: $!";
open my $in2, "<:gzip", $read2 or die "Cannot open $read2: $!";

open my $out1, ">:gzip", "$outdir/reformatted.R1.fastq.gz" or die "Cannot open output file: $!";
open my $out2, ">:gzip", "$outdir/reformatted.R2.fastq.gz" or die "Cannot open output file: $!";

# Initialize counters
my $num = 0;

# Process each pair of reads in the input files
while (1) {
    # Read 4 lines for each read in the pair
    my $line1_1 = <$in1>;
    my $line1_2 = <$in1>;
    my $line1_3 = <$in1>;
    my $line1_4 = <$in1>;

    my $line2_1 = <$in2>;
    my $line2_2 = <$in2>;
    my $line2_3 = <$in2>;
    my $line2_4 = <$in2>;

    # Exit the loop if end of file is reached
    last unless (defined($line1_1) and defined($line2_1));
    chomp($line1_1, $line1_2, $line1_3, $line1_4, $line2_1, $line2_2, $line2_3, $line2_4);

    # Define the UMI and Barcode sequences (example positions)
    my $UMI = substr($line2_2, 22, 10);
    my $valueA = substr($line2_2, 70, 8);
    my $valueB = substr($line2_2, 32, 8);

    # Write filtered reads with concatenated barcode and UMI into the main output files
    my @header_split_1 = split / /, $line1_1;
    my @header_split_2 = split / /, $line2_1;
    print $out2 "$header_split_1[0]\n$line1_2\n$line1_3\n$line1_4\n";
    print $out1 "$header_split_2[0]\n$valueB$valueA$UMI\n$line2_3\n$line2_4\n";

    $num++;
}

# Log results
print $log "Total processed reads: $num\n";
print $log "Filter primer end at: " . `date`;

# Close all filehandles
close $in1;
close $in2;
close $out1;
close $out2;
close $log;
