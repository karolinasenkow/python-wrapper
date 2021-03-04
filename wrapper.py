import os
import csv
import argparse
from Bio import Entrez
from Bio import SeqIO

# declare args
args = argparse.Namespace

file = open('input.txt').readlines()

# create outfiles
output = open('miniProject.log','w')
EF999921_fasta = open('EF999921_fasta.fasta', 'w')
EF999921_CDS = open('EF999921_CDS.fasta','w')
longest_contig = dict()
blast_input = open('blast_input.fasta', 'w')
blast_db = open('blast_db.fasta', 'w')

# set current path
path = os.getcwd()

# Question 1
# retrieve HCMV transcriptomes 2- and 6-days post-infection (dpi)

def fastq(file):
	file = open(file).readlines()
	for i in file:
		# retrieve  2dpi (Donors 1, 3) & 6dpi (Donors 1, 3)
		os.system('wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/' + i[:-3] + '/' + i)

# convert to paired end fastq files
def split_lines(file):
	file = open(file).read().splitlines()
	for i in file:
		os.system('fastq-dump -I --split-files ' + i)


''' build a transcriptome index for HCMV (NCBI accession EF999921)
 use biopython to retrieve and generate the appropriate input and then build the index with kallisto'''
def CDS_record():
	Entrez.email = 'ksenkow@luc.edu'
	count = 0
	# write out CDS of EF999921
	handle = Entrez.efetch(db='nucleotide', id='EF999921',rettype='gb', retmode='text')
	record = SeqIO.read(handle, 'genbank') # retrieve genbank record for EF999921
	for i in record.features:
		if i.type == 'CDS': # find CDS in record
			count += 1 # count each CDS and write to fasta file
			coordinates = i.location 
			EF999921_CDS.write('>' + ' '.join(map(str, i.qualifiers['protein_id'])))
			EF999921_CDS.write('\n')
			EF999921_CDS.write(str(coordinates.extract(record.seq)))
			EF999921_CDS.write('\n')
	output.write('The HCMV genome (EF999921) has ' + str(count) + ' CDS.')
	output.write('\n')

def kallisto(file): # run kallisto
	os.system('time kallisto index -i index.idx EF999921_CDS.fasta')
	file = open(file).read().splitlines() # loop through all accessions in input file
	for i in file:
		os.system('time kallisto quant -i index.idx -o ' + path + '/results_' + i + ' -b 30 -t 2 ' + path + '/' + i + '_1.fastq ' + path + '/' + i + '_2.fastq')

def sleuth():
	os.system('python3 create_table.py') # run code to create sample data
	os.system('Rscript sleuth.R') # run sleuth
	sleuth_results = open('sleuth_results.txt').readlines()
	for i in sleuth_results: # write out results of sleuth.R to miniProject.log
		output.write(i)
	output.write('\n')

''' Using Bowtie2, create an index for HCMV (NCBI accession EF999921)
    Next, save the reads that map to the HCMV index for use in assembly.
    Write to your log file the number of reads in each transcriptome before and after the Bowtie2 mapping '''

def bowtie2_index():
        handle = Entrez.efetch(db='nucleotide', id='EF999921', rettype='fasta', retmode='text')
        record = list(SeqIO.parse(handle, 'fasta')) # retrieve fasta record for EF999921
        for i in record:
                SeqIO.write(i, 'EF999921_fasta.fasta', "fasta")
        #record = SeqIO.read(handle, 'fasta')
        #EF999921_fasta.write('>' + str(record.description)) + EF999921_fasta.write('\n' + str(record.seq)

def bowtie2_map(file): # run bowtie2
	os.system('bowtie2-build EF999921_fasta.fasta EF999921') 
	file = open(file).read().splitlines() # loop through accession IDs in input file
	for i in file:
		os.system('bowtie2 --quiet -x EF999921 -1 ' + path + '/' + i + '_1.fastq -2 ' + path + '/' + i + '_2.fastq -S ' + path + '/' + i + '_EF999921_map.sam --al-conc ' + path + '/bt2_' + i + '.fastq')

def transcriptome_reads(file): # count reads before and after bowtie2 filtering
	file = open(file).read().splitlines()
	before = []
	after = []
	count = 0
	for i in file:
		donor = path + '/' + i + '_1.fastq'
		f = open(donor).readlines()
		for x in f:
			if x.startswith('@'): # identifies new fastq read
				count +=2 # 1 for _1.fastq and 1 for _2.fastq
		before.append(count) # append to list before
		count = 0

	count = 0
	for i in file:
		donor = path + '/bt2_' + i + '.1.fastq'
		f = open(donor).readlines()
		for x in f:
                        if x.startswith('@'): # identifies new fastq read
                                count +=2 # 1 for _1.fastq and 1 for _2.fastq
		after.append(count) # append to list after
		count = 0

	output.write('Donor 1 (2dpi) had ' + str(before[0]) + ' read pairs before Bowtie2 filtering and ' + str(after[0]) + ' read pairs after. \n')
	output.write('Donor 1 (6dpi) had ' + str(before[1]) + ' read pairs before Bowtie2 filtering and ' + str(after[1]) + ' read pairs after. \n')
	output.write('Donor 3 (2dpi) had ' + str(before[2]) + ' read pairs before Bowtie2 filtering and ' + str(after[2]) + ' read pairs after. \n')
	output.write('Donor 3 (6dpi) had ' + str(before[3]) + ' read pairs before Bowtie2 filtering and ' + str(after[3]) + ' read pairs after. \n')

def spades(file):
	file = open(file).read().splitlines()
	first = 'bt2_' + file[0] # rename output reads
	second = 'bt2_' + file[1]
	third = 'bt2_' + file[2]
	fourth = 'bt2_' + file[3]
	command = 'spades -k 55,77,99,127 -t 2 --only-assembler --pe1-1 ' + path + '/' + first + '.1.fastq --pe1-2 ' + path + '/' + first + '.2.fastq --pe2-1 ' + path + '/' + second + '.1.fastq --pe2-2 ' +  path + '/' + second + '.2.fastq --pe3-1 ' + path + '/' + third + '.1.fastq --pe3-2 ' + path + '/' + third + '.2.fastq --pe4-1 ' + path + '/' + fourth + '.1.fastq --pe4-2 ' + path + '/' + fourth + '.2.fastq -o ' + path + '/spades_assembly/'
	os.system(command) # run spades
	output.write(command)
	output.write('\n')


def contig_calc():
	os.chdir(path + '/spades_assembly') # navigate to spades_assembly folder
	record = SeqIO.parse('contigs.fasta','fasta')
	count = 0
	seq_len = 0
	for i in record:
		if len(i.seq) > 1000: # count all sequences > 1000 bp
			count += 1
			seq_len += len(i.seq) # count len of each seq > 10000 bp
			longest_contig[len(i.seq)] = (i.description, i.seq) # key = seq len, value = (description, seq)
	output.write('There are ' + str(count) + ' contigs > 1000 bp in the assembly.') # write to output file
	output.write('\n')
	output.write('There are ' + str(seq_len) + ' bp in the assembly.')
	output.write('\n')

def blast_inputs():
	# create fna file of longest contig (query)
	key_max = max(longest_contig, key=int) # retrieve seq of longest contig
	blast_input.write('>' + longest_contig[key_max][0]) # write >description of longest contig
	blast_input.write('\n')
	blast_input.write(str(longest_contig[key_max][1]) + '\n') #write seq of longest contig
	blast_input.close()

	# create db file - pull from NCBI
	Entrez.email = 'ksenkow@luc.edu'
	handle = Entrez.esearch(db='nucleotide', term='Betaherpesvirinae [Organism] AND refseq[filter]')
	record = Entrez.read(handle)
	handle = Entrez.esearch(db='nucleotide', term='Betaherpesvirinae [Organism] AND refseq[filter]', retmax =record['Count']) # retmax displays all hits
	record = Entrez.read(handle) # read all records that have betaherpesvirinae refseq
	for i in record['IdList']:
		handle = Entrez.efetch(db='nucleotide', id=i, rettype='fasta') # write them out to fasta file
		record = SeqIO.read(handle, 'fasta')
		blast_db.write('>' + str(record.description))
		blast_db.write('\n')
		blast_db.write(str(record.seq))
		blast_db.write('\n')
	blast_db.close()

def blast():
	# create db in blast+
	os.system('makeblastdb -in ' + path + '/blast_db.fasta -out ' + path + '/betaherpesvirinae -title betaherpesvirinae -dbtype nucl')
	os.system('blastn -query ' + path + '/blast_input.fasta -db betaherpesvirinae -max_target_seqs 10 -out ' + path + '/blast_results.txt -outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitle"')
	output.write('sacc' + '\t' + 'pident' + '\t' + 'length' + '\t' + 'qstart' + '\t' + 'qend' + '\t' + 'sstart' + '\t' + 'send' + '\t' + 'bitscore' + '\t' + 'eval' + '\t' + 'stitle')
	output.write('\n') # write desired headers to output file
	read_blast_results = open('blast_results.txt').read().splitlines() # read output of blast+ output file and write to miniProject log
	for i in read_blast_results:
		output.write(str(i))
		output.write('\n')

parser = argparse.ArgumentParser(description='Run pipeline.') # use argparse to have option to run test data or download data using wget
parser.add_argument('--download','-d', type=int, required=False, dest="download",help='0 for download 1 for test data')

args = parser.parse_args()

if __name__ == '__main__':
	choice = args.download # if user chooses to download, run all methods
	if choice == 0 :
		fastq('input.txt')
		split_lines('input.txt')
		CDS_record()
		kallisto('input.txt')
		sleuth()
		bowtie2_index()
		bowtie2_map('input.txt')
		transcriptome_reads('input.txt')
		spades('input.txt')
		contig_calc()
		blast_inputs()
		blast()

	elif choice == 1: # if user uses to use test data, run methods after split_lines()
		CDS_record()
		kallisto('input.txt')
		sleuth()
		bowtie2_index()
		bowtie2_map('input.txt')
		transcriptome_reads('input.txt')
		spades('input.txt')
		contig_calc()
		blast_inputs()
		blast()
	else:
		print("Please choose 0 for download and 1 for test data") 

