# ViromeQC #
 
## Description ##
 
* Provides an enrichment score for VLP viromes with respect to metagenomes
* Useful benchmark for the quality of enrichment of a virome
* Tested on Linux Ubuntu Server 16.04 LTS and on Linux Mint 19

**Requires:**

* [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) >= v. 2.3.4
* [Samtools](http://samtools.sourceforge.net/) >= 1.3.1
* [Biopython](https://github.com/biopython/biopython) >= 1.69
* [Pysam](http://pysam.readthedocs.io/en/latest/) >= 0.14
* [Diamond](http://github.com/bbuchfink/diamond) (tested on v.0.9.9 and 0.9.29)
* Python3 (tested on 3.6)
* [pandas](https://pandas.pydata.org) >= 0.20

**Update:** _ViromeQC_ now works with newer versions of diamond (e.g. v0.9.29) 
Thanks to Ryan Cook ([@RyanCookAMR](https://twitter.com/RyanCookAMR)) for the new diamond db

## Usage ##

### Step 1: clone or download the repository ###

`git clone --recurse-submodules https://github.com/SegataLab/viromeqc.git`

or download the repository from the **[releases](https://github.com/SegataLab/viromeqc/releases)** page

### Step 2: install the database: ###

This steps downloads the database file. This needs to be done only the first time you run ViromeQC. This may require a few minutes, depending on your internet connection.

`viromeQC.py --install`

Alternatively, you can also download the database files from [Zenodo](https://zenodo.org/record/4020594#.X1jxgGMzZDM). Once downloaded the files, create a folder named `index/` in the ViromeQC installation folder and unzip all the files in this folder.

### Step 3: Run on your sample ###

`viromeQC.py -i <input_virome_file(s)> -o <report_file.txt>`

*Please Note:* 
You can pass more than one file as input (e.g. for multiple runs or paired end reads). However, you can process only one sample at a time with this command. If you want to parallelize the execution, this can be easily done with [Parallel](https://www.gnu.org/software/parallel/) or equivalent tools.

You can try the test example (`test/test.sh`) which analyzes 10'000 reads from the sample `SRR829034`. This should take approximately 1 or 2 minutes.

Parameters:

```
usage: viromeQC.py -i <input_virome_file> -o <report_file.txt>

optional arguments:
  -h, --help            show this help message and exit
  -i [INPUT [INPUT ...]], --input [INPUT [INPUT ...]]
                        Raw Reads in FASTQ format. Supports multiple inputs
                        (plain, gz o bz2) (default: None)
  -o OUTPUT, --output OUTPUT
                        output file (default: None)
  --minlen MINLEN       Minimum Read Length allowed (default: 75)
  --minqual MINQUAL     Minimum Read Average Phred quality (default: 20)
  --bowtie2_threads BOWTIE2_THREADS
                        Number of Threads to use with Bowtie2 (default: 4)
  --diamond_threads DIAMOND_THREADS
                        Number of Threads to use with Diamond (default: 4)
  -w {human,environmental}, --enrichment_preset {human,environmental}
                        Calculate the enrichment basing on human or
                        environmental metagenomes. Defualt: human-microbiome
                        (default: human)
  --bowtie2_path BOWTIE2_PATH
                        Full path to the bowtie2 command to use, deafult
                        assumes that bowtie2 is present in the system path
                        (default: bowtie2)
  --diamond_path DIAMOND_PATH
                        Full path to the diamond command to use, deafult
                        assumes that diamond is present in the system path
                        (default: diamond)
  --version             Prints version informations (default: False)
  --install             Downloads database files (default: False)
  --sample_name SAMPLE_NAME
                        Optional label for the sample to be included in the
                        output file (default: None)
  --tempdir TEMPDIR     Temporary Directory override (default is the system
                        temp. directory) (default: None)
```

### Pipeline structure ###

ViromeQC starts from FASTQ files (compressed files are supported), and will:

1. Elimitate short and low quality reads
    - *adjust the `minqual` and `minlen` parameters if you want to change the thresholds*
2. Map the reads against a curated collection of rRNAs and single-copy bacteral markers
3. Filter the reads to remove short and dlsivergent alignments
4. Compute the enrichment value of the sample, compared to the median observed in human metagenomes
    - use `-w environmental` for envronmental reads
    - reference medians for un-enriched metagenomes are taken from `medians.csv`, you can provide your own data to ViromeQC by changing this file accordingly
5. Produce a report file with the alignment rates and the final enrichment score (which is the minimum enrichment observed across SSU-rRNA, LSU-rRNA and single-copy markers)


### Output ###

Output is given as a TSV file with the following structure:


|    Sample    |    Reads    |    Reads_HQ    |    SSU rRNA alignment (%)    |    LSU rRNA alignment (%)   |    Bacterial_Markers alignment (%)   |    total enrichmnet score
|---|---|---|---|---|---|---|
|    your_sample.fq | 40000 | 39479 | 0.00759898  | 0.0227969 | 0.01266496  | 5.795329


- An alignment score of 5.8 means that the virome is 5.8 times more enriched than a comparable metagenome
- High score (e.g. 10-50) reflect high VLP enrichment 


## Citation ##

If you find this tool useful, please cite:

*Zolfo, M., Pinto, F., Asnicar, F., Manghi, P., Tett A., Segata N.* **[Detecting contamination in viromes using ViromeQC](https://www.nature.com/articles/s41587-019-0334-5)**, *Nature Biotechnology* 37, 1408–1412 (2019)

