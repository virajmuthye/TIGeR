# Read the FisherStrand Phred-scaled p-values
FStest <- read.table("tmpFS.INFO", header = TRUE, sep="\t", row.names=NULL)

# Extract the column for plotting
data <- FStest$FS

png(filename = "allele_FS_histogram.png")

# Create a histogram of FS p-values
hist(data, 
     main = "Histogram of Fisher Strand Phred-scaled p-value", 
     xlab = "Phred-scaled p-value", 
     ylab = "Frequency", 
     col = "green", 
     border = "black")

dev.off()



