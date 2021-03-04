library(sleuth)

#read in sample_table.txt which describes samples and kallisto output
stab <- read.table("sample_table.txt",header=TRUE,stringsAsFactors=FALSE)

# initialize object
so <- sleuth_prep(stab)

#fit model comparing conditions
so <- sleuth_fit(so, ~condition, 'full')

#fit reduced model
so <- sleuth_fit(so, ~1, 'reduced')

#perform likelihood ratio test for differential expression between conditions
so <- sleuth_lrt(so, 'reduced', 'full')

#load dplyr package
library(dplyr)

#pull the test results
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)

#filter most significant results (FDR/qval < 0.05)
#sort by pval
test_stat <- dplyr::filter(sleuth_table, qval <= 0.05) %>% dplyr::arrange(pval)

#select variables to display 
select = dplyr::select(test_stat, target_id, test_stat, pval, qval)

#write FDR < 0.05 transcripts to file
write.table(select, file="sleuth_results.txt",quote = FALSE,row.names = FALSE, sep = '\t')

