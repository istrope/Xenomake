module load java/jdk-18.0.1.1
#java -jar ../human/picard.jar FilterSamReads I=../mouse/final.mapped.bam O=../mouse/remaining.bam READ_LIST_FILE=../human/human.mouse.common.readnames FILTER=excludeReadList
#java -jar ../human/picard.jar FilterSamReads I=../human/final.mapped.bam O=../human/remaining.bam READ_LIST_FILE=../human/human.mouse.common.readnames FILTER=excludeReadList
#java -jar picard.jar FilterSamReads I=../mouse/final.mapped.bam O=../mouse/overlapped.bam READ_LIST_FILE=human.mouse.common.readnames FILTER=includeReadList

cat overlap-ambiguous.fq |grep "@"|sed "s/@//g"|sort -u > overlap-ambiguous.readname
java -jar picard.jar FilterSamReads I=overlapped.bam O=overlap-ambiguous.bam READ_LIST_FILE=overlap-ambiguous.readname FILTER=includeReadList SORT_ORDER=queryname
java -jar picard.jar FilterSamReads I=../mouse/overlapped.bam O=../mouse/overlap-ambiguous.bam READ_LIST_FILE=overlap-ambiguous.readname FILTER=includeReadList SORT_ORDER=queryname

cat overlap-host.fq |grep "@"|sed "s/@//g"|sort -u > overlap-host.readname
java -jar picard.jar FilterSamReads I=../mouse/overlapped.bam O=../mouse/overlap-host.bam READ_LIST_FILE=overlap-host.readname FILTER=includeReadList SORT_ORDER=queryname

cat overlap-graft.fq |grep "@"|sed "s/@//g"|sort -u > overlap-graft.readname
java -jar picard.jar FilterSamReads I=overlapped.bam O=overlap-graft.bam READ_LIST_FILE=overlap-graft.readname FILTER=includeReadList SORT_ORDER=queryname

cat overlap-both.fq |grep "@"|sed "s/@//g"|sort -u > overlap-both.readname
java -jar picard.jar FilterSamReads I=overlapped.bam O=overlap-both.bam READ_LIST_FILE=overlap-both.readname FILTER=includeReadList SORT_ORDER=queryname
java -jar picard.jar FilterSamReads I=../mouse/overlapped.bam O=../mouse/overlap-both.bam READ_LIST_FILE=overlap-both.readname FILTER=includeReadList SORT_ORDER=queryname

eval "$(conda shell.bash hook)"
conda activate /project/qzhu/lab_conda/spacemake

python3 filter_mm_reads_both.py --in1 overlap-both.bam --in2 ../mouse/overlap-both.bam --out1 overlap-both.filtered.bam --out2 ../mouse/overlap-both.filtered.bam
python3 filter_mm_reads_both.py --in1 overlap-ambiguous.bam --in2 ../mouse/overlap-ambiguous.bam --out1 overlap-ambiguous.filtered.bam --out2 ../mouse/overlap-ambiguous.filtered.bam

python3 filter_mm_reads.py --in-bam overlap-graft.bam --out-bam overlap-graft.filtered.bam
python3 filter_mm_reads.py --in-bam ../mouse/overlap-host.bam --out-bam ../mouse/overlap-host.filtered.bam
