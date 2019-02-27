#!/usr/bin/env python3
import os
import sys
import argparse
import zipfile
import time
import tempfile
import subprocess
import pandas as pd


try:
	from Bio import SeqIO
	from Bio.Seq import Seq
	from Bio.SeqRecord import SeqRecord
	from Bio.Alphabet import IUPAC
except ImportError as e:
	metamlst_print("Failed in importing Biopython. Please check Biopython is installed properly on your system!",'FAIL',bcolors.FAIL)
	sys.exit(1)

__author__ = 'Moreno Zolfo (moreno.zolfo@unitn.it)'
__version__ = '1.0'
__date__ = '1 Feb 2019'



def byte_to_megabyte(byte):
    """
    Convert byte value to megabyte
    """

    return byte / (1024.0**2)


class ReportHook():
    def __init__(self):
        self.start_time = time.time()

    def report(self, blocknum, block_size, total_size):
        """
        Print download progress message
        """

        if blocknum == 0:
            self.start_time = time.time()
            if total_size > 0:
                sys.stderr.write("Downloading file of size: {:.2f} MB\n"
                                 .format(byte_to_megabyte(total_size)))
        else:
            total_downloaded = blocknum * block_size
            status = "{:3.2f} MB ".format(byte_to_megabyte(total_downloaded))

            if total_size > 0:
                percent_downloaded = total_downloaded * 100.0 / total_size
                # use carriage return plus sys.stderr to overwrite stderr
                download_rate = total_downloaded / (time.time() - self.start_time)
                estimated_time = (total_size - total_downloaded) / download_rate
                estimated_minutes = int(estimated_time / 60.0)
                estimated_seconds = estimated_time - estimated_minutes * 60.0
                status += ("{:3.2f} %  {:5.2f} MB/sec {:2.0f} min {:2.0f} sec "
                           .format(percent_downloaded,
                                   byte_to_megabyte(download_rate),
                                   estimated_minutes, estimated_seconds))

            status += "        \r"
            sys.stderr.write(status)


def download(url, download_file):
    """
    Download a file from a url
    """
    # try to import urllib.request.urlretrieve for python3
    try:
        from urllib.request import urlretrieve
    except ImportError:
        from urllib import urlretrieve

    if not os.path.isfile(download_file):
        try:
            sys.stderr.write("\nDownloading " + url + "\n")
            file, headers = urlretrieve(url, download_file,
                                        reporthook=ReportHook().report)
        except EnvironmentError:
            sys.stderr.write("\nWarning: Unable to download " + url + "\n")
    else:
        sys.stderr.write("\nFile {} already present!\n".format(download_file))



def print_version():
	print ("Version:\t"+__version__)
	print ("Author:\t\t"+__author__)
	print ("Software:\t"+'Virome QC')
	sys.exit(0)

class bcolors:
	HEADER = '\033[95m'
	OKBLUE = '\033[94m'
	OKGREEN = '\033[92m'
	WARNING = '\033[93m'
	FAIL = '\033[91m'
	ENDC = '\033[0m'
	OKGREEN2 = '\033[42m\033[30m'
	RED = '\033[1;91m'
	CYAN = '\033[0;37m'


def fancy_print(mesg,label,type,reline=False,newLine=False):
	opening = "\r" if reline else ''
	ending = "\r\n" if not reline or newLine else ''

	if len(mesg) < 65:
	
		sys.stdout.write(opening+mesg.ljust(66)+(type+'[ - '+label.center(5)+' - ]'+bcolors.ENDC).ljust(14)+ending)
	else: 
		c=0
		wds = []
		lines=[]
		for word in mesg.split(' '):

				if c + len(word)+2 > 65:
					print (' '.join(wds))
					c=0
					wds=[word]
					continue
				c = c+len(word)+2
				wds.append(word)
		sys.stdout.write(opening+(' '.join(wds)).ljust(66)+(type+'[ - '+label.center(5)+' - ]'+bcolors.ENDC).ljust(14)+ending)

	sys.stdout.flush()



CHECKER_PATH=os.path.abspath(os.path.dirname(__file__))
INDEX_PATH=CHECKER_PATH+"/index/"

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		description='Checks a virome FASTQ file for enrichment efficiency')

parser.add_argument("-i","--input", required=True, nargs="*", help="Raw Reads in FASTQ format. Supports multiple inputs (plain, gz o bz2)")
parser.add_argument("-o",'--output',required=True, help="output file")

parser.add_argument("--minlen", help="Minimum Read Length allowed",default='75')
parser.add_argument("--minqual", help="Minimum Read Average Phred quality",default='20')

parser.add_argument("--version", help="Prints version informations", action='store_true')
parser.add_argument("--bowtie2_threads", help="Number of Threads to use with Bowtie2",default='4') 
parser.add_argument("--diamond_threads", help="Number of Threads to use with Diamond",default='4') 

parser.add_argument("-w","--enrichment_preset", help="Calculate the enrichment basing on human or environmental metagenomes. Defualt: human-microbiome",default='human') 

parser.add_argument('--bowtie2_path', type=str, default='bowtie2',
        help="Full path to the bowtie2 command to use, deafult assumes "
             "that 'bowtie2 is present in the system path") 
parser.add_argument('--diamond_path', type=str, default='diamond',
        help="Full path to the diamond command to use, deafult assumes "
             "that 'diamond is present in the system path") 


args=parser.parse_args()

medians = pd.read_csv(CHECKER_PATH+'/medians.csv',sep='\t')

if args.version: print_version()
if args.enrichment_preset not in ['human','environmental']:
	fancy_print('-w should be human or environmental, no other preset are supported','FAIL',bcolors.FAIL)
	sys.exit(1)
	

#pre-flight check
for inputFile in args.input:
	if not os.path.isfile(inputFile):
		fancy_print("Error: file ",inputFile,'does not exist','ERROR',bcolors.FAIL)

commands = [['zcat', '-h'],[args.bowtie2_path, '-h'],[args.diamond_path, 'help']]


for sw in commands:
	try: 
		with open(os.devnull, 'w') as devnull:
			subprocess.check_call(sw, stdout=devnull, stderr=devnull)

	except Exception as e:
		fancy_print("Error, command not found: "+sw[0],'ERROR',bcolors.FAIL)

try:
	#download indexes if you don't have it
	to_download=[]
	fancy_print("Checking Database Files",'...',bcolors.OKBLUE,reline=True)

	if not os.path.isdir(INDEX_PATH):
		os.mkdir(INDEX_PATH)

	if not os.path.isfile(INDEX_PATH+'/SILVA_132_SSURef_Nr99_tax_silva.clean.1.bt2'):
		to_download.append('https://bitbucket.org/CibioCM/viromeqc/downloads/SILVA_132_LSURef_tax_silva.clean.zip')

	if not os.path.isfile(INDEX_PATH+'/SILVA_132_LSURef_tax_silva.clean.1.bt2'):
		to_download.append('https://bitbucket.org/CibioCM/viromeqc/downloads/SILVA_132_SSURef_Nr99_tax_silva.clean_1.zip')
		to_download.append('https://bitbucket.org/CibioCM/viromeqc/downloads/SILVA_132_SSURef_Nr99_tax_silva.clean_2.zip')
		to_download.append('https://bitbucket.org/CibioCM/viromeqc/downloads/SILVA_132_SSURef_Nr99_tax_silva.clean_3.zip')
 
	if not os.path.isfile(INDEX_PATH+'/amphora_bacteria.dmnd'):
		to_download.append('https://bitbucket.org/CibioCM/viromeqc/downloads/amphora_bacteria.dmnd.zip')


	fancy_print("Checking Database Files",'OK',bcolors.OKGREEN,reline=True,newLine=True)

	if(to_download):
		fancy_print("Need to download "+str(len(to_download))+' files','...',bcolors.OKBLUE,reline=True)
		for downloadable in to_download:
			download(downloadable, INDEX_PATH+'/tmp.zip')
			zipDB = zipfile.ZipFile(INDEX_PATH+'/tmp.zip', 'r')
			zipDB.extractall(INDEX_PATH)
			zipDB.close()
			os.remove(INDEX_PATH+'/tmp.zip')
		fancy_print("Need to download "+str(len(to_download))+' files','DONE',bcolors.OKGREEN,reline=True,newLine=True)
	 

except IOError: 
	print("Failed to retrieve DB")
	fancy_print("Failed to retrieve DB",'FAIL',bcolors.FAIL)
	sys.exit(1)


tmpdir = tempfile.TemporaryDirectory()
tmpdirname = tmpdir.name

if len(args.input) > 1:
	fancy_print('Merging '+str(len(args.input))+' files','...',bcolors.OKBLUE,reline=True)
	with open(tmpdirname+'/combined.fastq','a') as combinedFastq:
		for infile in args.input:
			#print(['cat',infile,'>>',tmpdirname+'/combined.fastq'])
			if infile.endswith('.gz'):
				uncompression_cmd = 'zcat'
			elif infile.endswith('.bz2'):
				uncompression_cmd = 'bzcat'
			else:
				uncompression_cmd = 'cat'

			subprocess.check_call([uncompression_cmd,infile],stdout=combinedFastq)

	inputFile = tmpdirname+'/combined.fastq'

	fancy_print('Merging '+str(len(args.input))+' files','DONE',bcolors.OKGREEN,reline=True,newLine=True)
else:
	inputFile = args.input[0]

fileName = os.path.basename(inputFile).replace('.fastq','').replace('.fq','')

fastq_len_cmd = [CHECKER_PATH+'/fastq_len_filter.py', '--min_len',args.minlen, '--min_qual',args.minqual,'--count',tmpdirname+'/'+fileName+'.nreads','-i',inputFile,'-o',tmpdirname+'/'+fileName+'.filter.fastq']
try: 
	fancy_print('[fastq_len_filter] | filtering HQ reads','...',bcolors.OKBLUE)

	subprocess.check_call(fastq_len_cmd)


	with open(tmpdirname+'/'+fileName+'.nreads') as readCounts:
		HQReads, totalReads = [line.strip().split('\t')[0:2] for line in readCounts][0]

	fancy_print('[fastq_len_filter] | '+HQReads+' / '+totalReads+' ('+str(round(float(HQReads)/float(totalReads),2)*100)+'%) reads selected','DONE',bcolors.OKGREEN)

except Exception as e: 
	fancy_print('Fatal error running fastq_len_filter. Error message: '+str(e),'FAIL',bcolors.FAIL)
	sys.exit(1)


try: 
	
	filteredFile = tmpdirname+'/'+fileName+'.filter.fastq' 
	

	fancy_print('[SILVA_LSU] | Bowtie2 Aligning','...',bcolors.OKBLUE,reline=True)
	
	bt2_command = ['bowtie2','--quiet','-p',args.bowtie2_threads,'--very-sensitive-local','-x',INDEX_PATH+'/SILVA_132_LSURef_tax_silva.clean','--no-unal','-U',filteredFile,'-S','-']
	#print("AAA")
	
	p4 = subprocess.Popen(['wc','-l'], stdin=subprocess.PIPE,stdout=subprocess.PIPE)
	p3 = subprocess.Popen(['samtools','view','-'], stdin=subprocess.PIPE,stdout=p4.stdin)
	p2 = subprocess.Popen([CHECKER_PATH+'/cmseq/filter.py','--minlen','50','--minqual','20','--maxsnps','0.075'],stdin=subprocess.PIPE,stdout=p3.stdin)
	p1 = subprocess.Popen(bt2_command, stdout=p2.stdin)

	p1.wait() 
	p2.communicate()
	p3.communicate()
	
	LSU_READS = int(p4.communicate()[0])
	LSU_READS_RATE = float(LSU_READS)/float(HQReads)*100
	enrichment_LSU = min(100,float(medians.loc[medians['parameter']=='rRNA_LSU',args.enrichment_preset]) / float(LSU_READS_RATE))


	fancy_print('[SILVA_LSU] | Bowtie2 Alignment rate: '+str(round(LSU_READS_RATE,4))+'% (~'+str(round(enrichment_LSU,1))+'x)','DONE',bcolors.OKGREEN,reline=True,newLine=True)

except Exception as e: 
	fancy_print('Fatal error running Bowtie2 on LSU rRNA. Error message: '+str(e),'FAIL',bcolors.FAIL)
	sys.exit(1)


try: 
	
	filteredFile = tmpdirname+'/'+fileName+'.filter.fastq' 
	

	fancy_print('[SILVA_SSU] | Bowtie2 Aligning','...',bcolors.OKBLUE,reline=True)
	
	bt2_command = ['bowtie2','--quiet','-p',args.bowtie2_threads,'--very-sensitive-local','-x',INDEX_PATH+'/SILVA_132_SSURef_Nr99_tax_silva.clean','--no-unal','-U',filteredFile,'-S','-']
	#print("AAA")
	
	p4 = subprocess.Popen(['wc','-l'], stdin=subprocess.PIPE,stdout=subprocess.PIPE)
	p3 = subprocess.Popen(['samtools','view','-'], stdin=subprocess.PIPE,stdout=p4.stdin)
	p2 = subprocess.Popen([CHECKER_PATH+'/cmseq/filter.py','--minlen','50','--minqual','20','--maxsnps','0.075'],stdin=subprocess.PIPE,stdout=p3.stdin)
	p1 = subprocess.Popen(bt2_command, stdout=p2.stdin)

	p1.wait() 
	p2.communicate()
	p3.communicate()
	
	SSU_READS = int(p4.communicate()[0])
	SSU_READS_RATE = float(SSU_READS)/float(HQReads)*100
	enrichment_SSU = min(100,float(medians.loc[medians['parameter']=='rRNA_SSU',args.enrichment_preset]) / float(SSU_READS_RATE))

	fancy_print('[SILVA_SSU] | Bowtie2 Alignment rate: '+str(round(SSU_READS_RATE,4))+'% (~'+str(round(enrichment_SSU,1))+'x)','DONE',bcolors.OKGREEN,reline=True,newLine=True)

except Exception as e: 
	fancy_print('Fatal error running Bowtie2 on SSU rRNA. Error message: '+str(e),'FAIL',bcolors.FAIL)
	sys.exit(1)


try:
	
	filteredFile = tmpdirname+'/'+fileName+'.filter.fastq' 
	
	fancy_print('[AMPHORA]   | Diamond Aligning','...',bcolors.OKBLUE,reline=True)
	
	diamond_command = ['diamond','blastx','-q',filteredFile,'--threads',args.diamond_threads,'--outfmt','6','--db',INDEX_PATH+'/amphora_bacteria.dmnd','--id','50','--max-hsps','35','-k','0']
	#print("AAA")
	 
	cut_command='wc -l'.split(' ')
	p2 = subprocess.Popen('cut -f1 | sort | uniq | wc -l',shell=True, stdin=subprocess.PIPE,stdout=subprocess.PIPE)
	p1 = subprocess.Popen(diamond_command, stdout=p2.stdin)

	AMPHORA_READS = int(p2.communicate()[0])
	AMPHORA_READS_RATE = float(AMPHORA_READS)/float(HQReads)*100
	enrichment_AMPHORA = min(100,float(medians.loc[medians['parameter']=='AMPHORA2',args.enrichment_preset]) / float(AMPHORA_READS_RATE))

	fancy_print('[AMPHORA]   | Diamond Alignment rate: '+str(round(AMPHORA_READS_RATE,4))+'% (~'+str(round(enrichment_AMPHORA,1))+'x)','DONE',bcolors.OKGREEN,reline=True,newLine=True)

except Exception as e: 
	fancy_print('Fatal error running Diamond on Single-Copy-Proteins. Error message: '+str(e),'FAIL',bcolors.FAIL)
	sys.exit(1)


 
OVERALL_ENRICHM_RATE = min(enrichment_SSU,enrichment_LSU,enrichment_AMPHORA)
to_out=[fileName,totalReads,HQReads,SSU_READS_RATE,LSU_READS_RATE,AMPHORA_READS_RATE,OVERALL_ENRICHM_RATE]


fancy_print('Finished: overall enrichment score is ~'+str(round(OVERALL_ENRICHM_RATE,1))+'x','END',bcolors.OKGREEN,reline=True,newLine=True)


outFile = open(args.output,'w')
outFile.write("Sample\tReads\tReads_HQ\tSSU rRNA alignment rate\tLSU rRNA alignment rate\tAmphora2_bacterial_markers alignment rate\ttotal enrichmnet score\n")
outFile.write('\t'.join([str(x) for x in to_out])+'\n')
outFile.close()
tmpdir.cleanup()