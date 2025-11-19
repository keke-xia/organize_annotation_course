# Load ggplot2
# If ggplot2 is not installed in this R environment, run install.packages("ggplot2") once.
library(ggplot2)

# 1) Read AED cumulative distribution table
aed <- read.table("assembly.all.maker.renamed.gff.AED.txt",
                  header = TRUE,
                  sep = "\t",
                  stringsAsFactors = FALSE)

# The second column name is long: "/data/users/..."
# Rename it to a simpler name: "cum_frac"
colnames(aed)[2] <- "cum_frac"

# 2) Open a PDF device to save the plot (no interactive window on HPC)
pdf("AED_CDF.pdf", width = 6, height = 5)

# 3) Plot cumulative AED distribution
p <- ggplot(aed, aes(x = AED, y = cum_frac)) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = c(0.25, 0.5), linetype = "dashed") +
  theme_bw() +
  labs(title = "Cumulative distribution of Annotation Edit Distance (AED)",
       x = "AED",
       y = "Cumulative fraction of gene models")

print(p)

# 4) Close the PDF device
dev.off()
