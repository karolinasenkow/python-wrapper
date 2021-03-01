import os
from Bio import Entrez
from Bio import SeqIO
file = open('input.txt').readlines()

# Question 1
# retrieve HCMV transcriptomes 2- and 6-days post-infection (dpi)
def fastq(file):
	file = open(file).readlines()
	for i in file:
		# retrieve  2dpi (Donors 1, 3) & 6dpi (Donors 1, 3)
		os.system('wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/' + i[:-3] + '/' + i)
		# convert to paired end fastq files
		os.system('fastq-dump -I --split-files ' + i)
		print('hellow world')
# Question 2
# build a transcriptome index for HCMV (NCBI accession EF999921)
# use biopython to retrieve and generate the appropriate input and then build the index with kallisto
# you will need to extract the CDS features from the GenBank format

