import pysam
import argparse

parser = argparse.ArgumentParser(description='find overlapping mapped reads between organisms')
parser.add_argument('--human',help='human mapped bam')
parser.add_argument('--mosue',help='mouse mapped bam')
parser.add_argument('--out',help='text file of overlapped reads')
args = parser.parse_args()

human = pysam.AlignmentFile(args.human,'rb')
mouse = pysam.AlignmentFile(args.mouse,'rb')


#extract readnames
human_readnames = human.query_name
mouse_readnames = mouse.query_name

#Find overlapping reads
overlapped = list(set(human_readnames) & set(mouse_readnames))

#write reads to file 
with open(args.out,'w') as fp:
    for item in overlapped:
        fp.write("%s\n" % item)
    print('Finished Writing Overlapped Reads')
