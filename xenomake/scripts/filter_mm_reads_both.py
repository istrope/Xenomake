"""
Code taken from spacemake code under the same name filter_mm_reads.py
Parses through the tagged bam file and filters multi mapped reads and returns in bam format
see https://github.com/rajewsky-lab/spacemake/blob/master/spacemake/snakemake/scripts/filter_mm_reads.py for source code 
"""
import sys
import pysam
import datetime
import argparse
import numpy as np

counted_regions = ['UTR', 'CODING']

def select_alignment(alignments):
	read_names = [aln.query_name for aln in alignments]
	if read_names.count(read_names[0]) != len(read_names):
		print(read_names)
		raise Exception(f'input alignments do not come from the same read')

	def is_exonic(aln):
		if not aln.has_tag('XF'):
			return False
		#print(aln.get_tag("XF"))
		return aln.get_tag('XF') in counted_regions

	alignments_are_exonic = np.array([is_exonic(aln) for aln in alignments])
	exonic_ix = np.where(alignments_are_exonic  == True)[0]
	num_exonic = exonic_ix.shape[0]

	if num_exonic == 1:
		# if only one exonic reads from the group
		# return the exonic indices
		#print(alignments[exonic_ix[0]])
		return alignments[exonic_ix[0]]
	elif num_exonic > 1:
		#pass
		return None
		'''
		for aln in alignments:
			if aln.has_tag("gf") and aln.has_tag("gn"):
				tt = aln.get_tag("gf")
				tx = aln.get_tag("gn")
				print(aln.query_name, tt, tx, aln.reference_name, aln.reference_start, aln.get_tag("AS"))
		'''
	else:
		return None

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Filter out ambiguous multi-mapper reads')

	parser.add_argument('--in1', help='input bam 1')
	parser.add_argument('--out1', help='output bam 1')
	parser.add_argument('--in2', help='input bam 2')
	parser.add_argument('--out2', help='output bam 2')

	args = parser.parse_args()
	print(args)

	bam_in_1 = pysam.AlignmentFile(args.in1, "rb")
	bam_in_2 = pysam.AlignmentFile(args.in2, "rb")
	bam_out_1 = pysam.AlignmentFile(args.out1, 'wb', header= bam_in_1.header)
	bam_out_2 = pysam.AlignmentFile(args.out2, 'wb', header= bam_in_2.header)

	counter_1 = 0
	counter_2 = 0
	start_time = datetime.datetime.now()


	alns_1 = bam_in_1.fetch(until_eof=True)
	alns_2 = bam_in_2.fetch(until_eof=True)

	#print(next(alns_1))
	#sys.exit(0)

	ind_1 = 0
	ind_2 = 0
	multi_mappers_1 = []
	multi_mappers_2 = []

	to_write_1 = None
	to_write_2 = None

	counts = 0
	while True:
		ax_1 = next(alns_1, None)
		counts+=1
		if counts%1000000==0:
			print(counts)
		
		if ax_1==None:
			print("Stopped, because alns_1 reached end")
			ax_2 = next(alns_2, None)
			if ax_2==None:
				print("Stopped, because alns_2 reached end")
			break
		q_1 = ax_1.query_name
		mapped_number_1 = ax_1.get_tag("NH")
		if mapped_number_1 == 1:
			#bam_out_1.write(ax_1)
			to_write_1 = ax_1
		else:
			while True:
				if len(multi_mappers_1) < (mapped_number_1 - 1):
					multi_mappers_1.append(ax_1)
					ax_1 = next(alns_1)
				else:
					multi_mappers_1.append(ax_1)
					aln_to_keep_1 = select_alignment(multi_mappers_1)
					if aln_to_keep_1 is not None:
						aln_to_keep_1.flag = aln_to_keep_1.flag & ~(1<<8)
						#bam_out_1.write(aln_to_keep_1)
						to_write_1 = aln_to_keep_1
					else: #no alignment to keep
						to_write_1 = None
					multi_mappers_1 = []
					break

		ax_2 = next(alns_2, None)
		q_2 = ax_2.query_name

		#print(q_1, q_2)
		if q_1!=q_2:
			print("Error, query name differs!")
			sys.exit(1)

		mapped_number_2 = ax_2.get_tag("NH")
		if mapped_number_2 == 1:
			#bam_out_2.write(ax_2)
			to_write_2 = ax_2
		else:
			while True:
				if len(multi_mappers_2) < (mapped_number_2 - 1):
					multi_mappers_2.append(ax_2)
					ax_2 = next(alns_2)
				else:
					multi_mappers_2.append(ax_2)
					aln_to_keep_2 = select_alignment(multi_mappers_2)
					if aln_to_keep_2 is not None:
						aln_to_keep_2.flag = aln_to_keep_2.flag & ~(1<<8)
						#bam_out_2.write(aln_to_keep_2)
						to_write_2 = aln_to_keep_2
					else: # no alignment to keep
						to_write_2 = None
					multi_mappers_2 = []
					break

		if to_write_1==None and to_write_2!=None:
			bam_out_2.write(to_write_2)
		elif to_write_1!=None and to_write_2==None:
			bam_out_1.write(to_write_1)
		elif to_write_1!=None and to_write_2!=None:
			as1 = to_write_1.get_tag("AS")
			as2 = to_write_2.get_tag("AS")
			#print(as1, as2, as1-as2)
			if as1>as2:
				bam_out_1.write(to_write_1)
			else:
				bam_out_2.write(to_write_2)
				
		
	
	sys.exit(0)

	multi_mappers = []



	for aln in bam_in.fetch(until_eof=True):
		counter = counter + 1

		if counter % 1000000 == 0:
			finish_time = datetime.datetime.now()
			delta_seconds = (finish_time - start_time).seconds
			# restart time
			start_time = finish_time
			print(f'Processed 1 millon records in {delta_seconds} seconds, total records processed {counter}. current time: {finish_time}')

		mapped_number = aln.get_tag('NH')

		if mapped_number == 1:
			bam_out.write(aln)
		else:
			if len(multi_mappers) < (mapped_number - 1):
				# still some multimappers missing. we need to add the alignments
				# until the last one to the list
				multi_mappers.append(aln)
			else:
				# add the last alignment
				multi_mappers.append(aln)
				# decide which, if any, to keep
				aln_to_keep = select_alignment(multi_mappers)

				if aln_to_keep is not None:
					# set aln secondary flag to 0, so that it is flagged as primary
					# secondary flag is at 0x100, so 8th bit (starting from 0)
					aln_to_keep.flag = aln_to_keep.flag & ~(1<<8)
					bam_out.write(aln_to_keep)

				# reset multimapper list
				multi_mappers = []
