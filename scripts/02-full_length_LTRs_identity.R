# In this script, we analyze the identity of full-length LTR retrotransposons (LTR-RTs)
# annotated in a genome. We read a GFF file with LTR-RT annotations,
# extract the LTR identity values, merge with classification data from TEsorter,
# and generate plots showing the distribution of LTR identities per clade within
# the Copia and Gypsy superfamilies.

library(tidyverse)
library(data.table)
library(cowplot)

#-------------------------------------------------
# Input files (edit paths if needed)
#-------------------------------------------------
gff_file <- "genomic.fna.mod.LTR.intact.raw.gff3"
cls_file <- "genomic.fna.mod.LTR.raw.fa.rexdb-plant.cls.tsv"
# cls_file is the output from TEsorter on the raw LTR-RT fasta file
#-------------------------------------------------
# Read and preprocess input data
#-------------------------------------------------
message("Reading GFF: ", gff_file)
anno <- read.table(gff_file, sep = "\t", header = FALSE)

# Remove subfeatures (Terminal repeats, TSDs) so we keep top-level TE annotations
exclude_feats <- c("long_terminal_repeat", "repeat_region", "target_site_duplication")
anno <- anno %>% filter(!V3 %in% exclude_feats)

# Extract Name and ltr_identity from the ninth column (attributes). This uses regex.

anno <- anno %>%
  rowwise() %>%
  mutate(
    # extract Name=... from attributes (V9)
    Name = str_extract(V9, "(?<=Name=)[^;]+"),
    # extract ltr_identity=... 
    Identity = as.numeric(str_extract(V9, "(?<=ltr_identity=)[^;]+")),
    # compute length as end - start
    length = as.numeric(V5) - as.numeric(V4)
  ) %>%
  # keep only the columns used downstream
  select(V1, V4, V5, V3, Name, Identity, length)

message("Reading classification: ", cls_file)
# Read classification table (TE name in first column). If your file doesn't have a header, set header=FALSE.
cls <- fread(cls_file, sep = "\t", header = TRUE)
setnames(cls, 1, "TE")

# TEsorter outputs encode the internal domain classification as TEName_INT#Classification. We split on '#',
# then keep only rows that correspond to internal-domain matches (Name ends with _INT), and strip _INT.
cls <- cls %>%
  separate(TE, into = c("Name", "Classification"), sep = "#", fill = "right") %>%
  filter(str_detect(Name, "_INT")) %>%
  mutate(Name = str_remove(Name, "_INT$"))

## Merge annotation with classification table
# Use a left join so all annotated TEs are kept even if they have no classification match
anno_cls <- merge(anno, cls, by = "Name", all.x = TRUE)

# Quick checks: how many per Superfamily/Clade (may be NA if classification missing)
message("Counts per Superfamily")
print(table(anno_cls$Superfamily, useNA = "ifany"))
message("Counts per Clade")
print(table(anno_cls$Clade, useNA = "ifany"))

#-------------------------------------------------
# Plot setup
#-------------------------------------------------
# binwidth controls histogram resolution around identity values 
binwidth <- 0.005
# x axis limits, these can be adjusted if needed, minimum identity in your data may differ
xlims <- c(0.80, 1.00)

# Compute a single y-max across ALL Copia and Gypsy clades. This ensures consistent y-axis scaling.
global_ymax <- anno_cls %>%
  filter(Superfamily %in% c("Copia", "Gypsy"), !is.na(Identity)) %>%
  # bin Identity into consistent breaks and count occurrences
  count(Superfamily, Clade, Identity = cut(Identity, seq(xlims[1], xlims[2], by = binwidth))) %>%
  pull(n) %>%
  max(na.rm = TRUE)

message("Global y-limit (shared for overview plots): ", global_ymax)

#-------------------------------------------------
# Plot function for one superfamily
#-------------------------------------------------
plot_by_clade <- function(df, sf, ymax) {
  df %>%
    filter(Superfamily == sf, !is.na(Identity), !is.na(Clade)) %>%
    ggplot(aes(x = Identity)) +
    # histogram with color coding per superfamily
    geom_histogram(binwidth = binwidth,
                   color = "black",
                   fill = ifelse(sf == "Copia", "#1b9e77", "#d95f02")) +
    # vertical stacking: one facet per Clade (ncol = 1) and fixed y scale so bars are comparable
    facet_wrap(~Clade, ncol = 1, scales = "fixed") +
    # x axis focused around xlims values
    scale_x_continuous(limits = xlims, breaks = seq(xlims[1], xlims[2], 0.05)) +
    # set y limit to the provided ymax (useful for consistent overview plots)
    scale_y_continuous(limits = c(0, ymax), expand = c(0, 0)) +  
    theme_cowplot() +
    theme(strip.background = element_rect(fill = "#f0f0f0"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold", hjust = 0.5)) +
    labs(title = sf, x = "Identity", y = "Count")
}

#-------------------------------------------------
# Generate Copia and Gypsy plots
#-------------------------------------------------
p_copia <- plot_by_clade(anno_cls, "Copia", global_ymax)
p_gypsy <- plot_by_clade(anno_cls, "Gypsy", global_ymax)

# Combine with cowplot side-by-side 
combined <- plot_grid(p_copia, p_gypsy, ncol = 2, rel_widths = c(1, 1))

ggsave("plots/01_LTR_Copia_Gypsy_cladelevel.png", combined, width = 12, height = 10, dpi = 300)
ggsave("plots/01_LTR_Copia_Gypsy_cladelevel.pdf", combined, width = 12, height = 10)