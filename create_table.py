import csv

# create table describing samples and kallisto output

def create_table(file):
	file = open(file).read().splitlines()

	# append header and all samples to list
	sample = ['sample']
	for i in file:
		sample.append(i)

	# append header and all conditions to list
	condition = ['condition', '2dpi', '6dpi', '2dpi', '6dpi']

	# append header and all paths to list
	path = ['path']
	for i in file:
		path.append('results_' + i)

	# write out data to tab-delimited txt file
	with open('sample_table.txt', 'w') as f:
		writer = csv.writer(f, delimiter='\t')
		writer.writerows(zip(sample, condition, path))

create_table('input.txt')
