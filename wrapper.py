import os
from Bio import Entrez

# Question 1
# retrieve HCMV transcriptomes 2- and 6-days post-infection (dpi)

# retrieve  2dpi (Donor 1)
os.system('wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660030/SRR5660030.1')
# convert to paired end fastq files
os.system('fastq-dump -I --split-files SRR5660030.1')

# retrieve 6dpi (Donor 1)
os.system('wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660033/SRR5660033.1')
# convert to paired end fastq files
os.system('fastq-dump -I --split-files SRR5660033.1')

# retrieve 2dpi (Donor 3)
os.system('wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660044/SRR5660044.1')
# convert to paired end fastq files
os.system('fastq-dump -I --split-files SRR5660044.1')

# retrieve 6dpi (Donor 3)
os.system('wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660045/SRR5660045.1')
# convert to paired end fastq files
os.system('fastq-dump -I --split-files SRR5660045.1')

# Question 2
# build a transcriptome index for HCMV (NCBI accession EF999921)
# use biopython to retrieve and generate the appropriate input and then build the index with kallisto
# you will need to extract the CDS features from the GenBank format

# write output to outfile
outfile = open('miniProject.log','w')

# 
Entrez.email = 'ksenkow@luc.edu'
handle = Entrez.esearch(db="nucleotide", term=)

