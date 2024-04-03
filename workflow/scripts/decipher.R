suppressMessages(library(DECIPHER))

sample <- snakemake@wildcards[['sample']]
multifasta <- snakemake@input[['fa']]
html_out <- snakemake@output[['html_out']]
aln_out <- snakemake@output[['aln_out']]

seqs <- readDNAStringSet(multifasta)
aligned <- AlignSeqs(seqs, iterations = 2, refinements = 1, verbose=FALSE)

BrowseSeqs(aligned,htmlFile = html_out, openURL = FALSE,colWidth = 100)
writeXStringSet(aligned,file=aln_out)
