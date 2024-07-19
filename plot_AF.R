# Read in the allele frequency data
allele_freqs <- read.table("tmp.txt", header = TRUE, sep="\t", row.names=NULL)

# Extract the allele frequency column for plotting
frequencies <- allele_freqs$col5

png(filename = "allele_frequency_histogram.png")

# Create a histogram of allele frequencies
hist(frequencies, 
     main = "Histogram of Allele Frequencies", 
     xlab = "Allele Frequency", 
     ylab = "Frequency", 
     col = "pink", 
     border = "black")

dev.off()




