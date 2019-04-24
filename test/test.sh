if [ -f demo_SRR829034.fastq.bz2 ]; then
	bunzip2 demo_SRR829034.fastq.bz2
fi;
../viromeQC.py -i demo_SRR829034.fastq -o out.txt
