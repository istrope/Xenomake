import pysam
import argparse
import datetime

def find_overlaps(x,y):
    start_time = datetime.datetime.now()
    human = pysam.AlignmentFile(x,'rb')
    mouse = pysam.AlignmentFile(y,'rb')

    #extract readnames
    human_readnames = [x.qname for x in human]
    mouse_readnames = [x.qname for x in mouse]

    #Find overlapping reads
    overlapped = list(set(human_readnames) & set(mouse_readnames))
    end_time = datetime.datetime.now()
    print(end_time)
    return(overlapped)
    #write reads to file 



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='find overlapping mapped reads between organisms')
    
    parser.add_argument('--human',help='human mapped bam')
    parser.add_argument('--mosue',help='mouse mapped bam')
    parser.add_argument('--out',help='text file of overlapped reads')
    
    args = parser.parse_args()

    overlapped = find_overlaps(args.human,args.mouse)

    with open(args.out,'w') as fp:
        for item in overlapped:
            fp.write("%s\n" % item)









