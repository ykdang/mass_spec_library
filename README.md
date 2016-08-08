# mass_spec_library
create  a reference library to search for any abnormal sequence caused by frame-shift

### The code is written in perl. requires two reference files:
1, the mRNA sequences (i.e. ORF sequence with defined position of start codon, here I added 100nt upstream of start codon)
2, the sorted gene list, ranking the emPAI value from high to low. 

### The result is a fasta file, containing 
1, the number of in-frame products translated from ORF
2, the non-redundant frame-shift fragment, defined by trypsin digestion. (i.e. either end with R/K or stop)

### Use, need to define the number of genes and type of frame-shift 
 i.e. "+1, -1 or -2" in the code, at first a few lines. 
