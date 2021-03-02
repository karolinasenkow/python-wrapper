import os
from Bio import Entrez
from Bio import SeqIO
file = open('input.txt').readlines()

# create outfile
output = open('miniProject.log','w')
EF999921_fasta = open('EF999921_fasta.fasta', 'w')
EF999921_CDS = open('EF999921_CDS.fasta','w')
longest_contig = dict()
blast_input = open('blast_input.fna', 'w')
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
		os.system('fastq-dump -I --split-files ' + i)

# Question 2
# build a transcriptome index for HCMV (NCBI accession EF999921)
# use biopython to retrieve and generate the appropriate input and then build the index with kallisto
# you will need to extract the CDS features from the GenBank format
def CDS_record():
	Entrez.email = 'ksenkow@luc.edu'
	count = 0
	# write out CDS of EF999921
	handle = Entrez.efetch(db='nucleotide', id='EF999921',rettype='gb', retmode='text')
	record = SeqIO.read(handle, 'genbank')
	for i in record.features:
		if i.type == 'CDS':
			count += 1
			coordinates = i.location
			EF999921_CDS.write('>' + ' '.join(map(str, i.qualifiers['protein_id'])))
			EF999921_CDS.write('\n')
			EF999921_CDS.write(str(coordinates.extract(record.seq)))
			EF999921_CDS.write('\n')
	output.write('The HCMV genome (EF999921) has ' + str(count) + ' CDS.')
	output.write('\n')

def kallisto(file):
	os.system('time kallisto index -i index.idx EF999921_CDS.fasta')
	file = open(file).read().splitlines()
	for i in file:
		os.system('time kallisto quant -i index.idx -o ' + path + '/results_' + i + ' -b 30 -t 2 ' + path + '/' + i + '_1.fastq ' + path + '/' + i + '_2.fastq')

def sleuth():
	os.system('python3 create_table.py')
	os.system('Rscript sleuth.R')
	sleuth_results = open('sleuth_results.txt').readlines()
	for i in sleuth_results: # write out results of sleuth.R to miniProject.log
		output.write(i)
	output.write('\n')

# Using Bowtie2, create an index for HCMV (NCBI accession EF999921)
# Next, save the reads that map to the HCMV index for use in assembly.
# Write to your log file the number of reads in each transcriptome before and after the Bowtie2 mapping

def bowtie2_index():
        handle = Entrez.efetch(db='nucleotide', id='EF999921', rettype='fasta')
        record = SeqIO.read(handle, 'fasta')
        EF999921_fasta.write('>' + str(record.description)) + EF999921_fasta.write('\n' + str(record.seq))

def bowtie2_map(file):
	os.system('bowtie2-build EF999921_fasta.fasta EF999921')
	file = open(file).read().splitlines()
	for i in file:
		os.system('bowtie2 --quiet -x EF999921 -1 ' + i + '_1.fastq -2 ' + i + '_2.fastq -S ' + i + '_EF999921_map.sam')

def file_version(): # convert sam to fastq
        pass

def spades(file):
	file = open(file).read().splitlines()
	first = file[0]
	second = file[1]
	third = file[2]
	fourth = file[3]
	command = 'spades -k 55,77,99,127 -t 2 --only-assembler --pe1-1 ' + path + '/' + first + '_1.fastq --pe1-2 ' + path + '/' + first + '_2.fastq --pe2-1 ' + path + '/' + second + '_1.fastq --pe2-2 ' +  path + '/' + second + '_2.fastq --pe3-1 ' + path + '/' + third + '_1.fastq --pe3-2 ' + path + '/' + third + '_2.fastq --pe4-1 ' + path + '/' + fourth + '_1.fastq --pe4-2 ' + path + '/' + fourth + '_2.fastq -o ' + path + '/spades_assembly/'
	os.system(command)
	output.write(command)
	output.write('\n')


def contig_calc():
	os.chdir(path + '/spades_assembly') # navigate to spades_assembly folder
	record = SeqIO.parse('contigs.fasta','fasta')
	count = 0
	seq_len = 0
	for i in record:
		if len(i.seq) > 1000:
			count += 1
			seq_len += len(i.seq)
			longest_contig[len(i.seq)] = (i.description, i.seq)
	output.write('There are ' + str(count) + ' contigs > 1000 bp in the assembly.')
	output.write('\n')
	output.write('There are ' + str(seq_len) + ' bp in the assembly.')
	output.write('\n')

def blast_inputs():
	# create fna file of longest contig (query)
	key_max = max(longest_contig, key=int) # retrieve seq of longest contig
	blast_input.write('>' + longest_contig[key_max][0])
	blast_input.write('\n')
	blast_input.write(str(longest_contig[key_max][1]) + '\n')

	# create db file
	Entrez.email = 'ksenkow@luc.edu'
	handle = Entrez.esearch(db='nucleotide', term='Betaherpesvirinae [Organism] AND refseq[filter]')
	record = Entrez.read(handle)
	handle = Entrez.esearch(db='nucleotide', term='Betaherpesvirinae [Organism] AND refseq[filter]', retmax =record['Count'])
	record = Entrez.read(handle)
	for i in record['IdList']:
		handle = Entrez.efetch(db='nucleotide', id=i, rettype='fasta')
		record = SeqIO.read(handle, 'fasta')
		blast_db.write('>' + str(record.description))
		blast_db.write('\n')
		blast_db.write(str(record.seq))
		blast_db.write('\n')

def blast():
	# create db in blast+
	os.system('makeblastdb -in ' + path + '/blast_db.fasta -out betaherpesvirinae -title betaherpresvirinae -dbtype nucl')

if __name__ == '__main__':
	#fastq('input.txt')
	#CDS_record()
	#kallisto('input.txt')
	#sleuth()
	#bowtie2_index()
	#bowtie2_map('input.txt')
	#spades('input.txt')
	contig_calc()
	blast_inputs()
	blast()
