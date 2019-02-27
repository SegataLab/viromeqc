#!/usr/bin/env python3

import os
from Bio import SeqIO
import argparse as ap
import sys
import os
import numpy as np 

def read_params(args):
	p = ap.ArgumentParser(description = 'fastax_len_filter.py Parameters\n')
	p.add_argument('--min_len', required = True, default = None, type = int)
	p.add_argument('--min_qual', default = 0, type = int)
	p.add_argument('--no_anonim', action="store_true")
	p.add_argument('--count')
	p.add_argument('-i','--input', required = True, help="Input File (FASTQ)")
	p.add_argument('-o','--output', required = True, help="Output File (FASTQ)")


	return vars(p.parse_args())

screenQual=range(20,31)
qualCounter = dict((k,0) for k in screenQual)

counter = 0
allCounter = 0
rpl=[]



if __name__ == '__main__':


	args = read_params(sys.argv)

	if not os.path.isfile(args['input']):
		print("Error: file "+args['input']+' is not accessible!')
		sys.exit(1)

	if args['input'].endswith('.gz'):
		import gzip
		from functools import partial
		_open = partial(gzip.open, mode='rt')
	elif args['input'].endswith('.bz2'):
		import bz2
		from functools import partial

		_open = partial(bz2.open, mode='rt')
	else:
		_open = open



	min_len = args['min_len']
	with open(args['output'],'w') as outf:

		with _open(args['input']) as f:
			for r in SeqIO.parse(f, "fastq"):

				avQual= np.mean(r.letter_annotations['phred_quality'])
				
				if avQual >= args['min_qual'] and len(r) >= min_len:	 
					if not args['no_anonim']:
						r.id = r.id+'_'+str(allCounter)

					counter+=1
					rpl.append(r)

					for qu in screenQual:
						if avQual >= qu:
							qualCounter[qu]+=1

					
					if len(rpl) % 30000 == 0:
				
						SeqIO.write(rpl, outf, "fastq")
						rpl=[]

				allCounter+=1


			if len(rpl) > 0:
				SeqIO.write(rpl, outf, "fastq")
				rpl=[]
				 
			if args['count']:
				outCount = open(args['count'],'w')
				outCount.write(str(counter)+'\t'+str(allCounter)+'\t'+'\t'.join([str(ke)+':'+str(val) for ke,val in qualCounter.items()]))
				outCount.close()
