#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(readr)
})

opt_list <- list(
  make_option("--workdir", type="character",
              default="/data/users/yliu2/Organization_and_annotation"),
  make_option("--outdir",  type="character",
              default=NULL),
  make_option("--prefix",  type="character",
              default="assembly_p_ctg"),
  make_option("--rate",    type="double",
              default=8.22e-9, help="Substitution rate per site per year")
)
opt <- parse_args(OptionParser(option_list = opt_list))

workdir <- opt$workdir
edta_dir <- file.path(workdir, "results/EDTA_annotation")
land    <- if (is.null(opt$outdir)) file.path(edta_dir, "landscape") else opt$outdir
dir.create(land, showWarnings = FALSE, recursive = TRUE)

divsum_file <- file.path(land, paste0(opt$prefix, ".parseRM.divsum.tsv"))
summary_file <- file.path(land, paste0(opt$prefix, ".parseRM.summary.tsv"))

read_divsum <- function(f){
  dt <- suppressWarnings(read_tsv(f, show_col_types = FALSE))
  # Try to normalize columns: expect something like columns for class/family and divergence bins/counts
  # Typical parseRM .divsum has columns like: "class","family","bin","count" or similar.
  # We reshape to have: superfamily, bin_div, count
  names(dt) <- tolower(names(dt))
  # Heuristics
  if (!("bin" %in% names(dt)) && ("divergence" %in% names(dt))) {
    dt$bin <- dt$divergence
  }
  if (!("count" %in% names(dt)) && ("copies" %in% names(dt))) {
    dt$count <- dt$copies
  }
  # Build superfamily from class/family fields if present
  sup <- if ("class" %in% names(dt)) dt$class else if ("order" %in% names(dt)) dt$order else NA
  fam <- if ("family" %in% names(dt)) dt$family else if ("superfamily" %in% names(dt)) dt$superfamily else NA
  dt$superfamily <- ifelse(!is.na(fam), fam, sup)
  dt$superfamily <- toupper(dt$superfamily)
  # Keep reasonable
  keep <- c("superfamily","bin","count")
  dt <- dt[, intersect(keep, names(dt))]
  dt <- dt %>% filter(!is.na(bin), !is.na(count))
  # Coerce bin to numeric (percentage divergence 0..100)
  dt$bin <- as.numeric(dt$bin)
  dt$count <- as.numeric(dt$count)
  dt
}

read_summary <- function(f){
  dt <- suppressWarnings(read_tsv(f, show_col_types = FALSE))
  names(dt) <- tolower(names(dt))
  # Fallback: if we only have per-copy rows with %div, aggregate
  cand_div <- intersect(c("percdiv","div","k2p","kimura","pdiv"), names(dt))
  cand_sup <- intersect(c("class","order","superfamily","family"), names(dt))
  if (length(cand_div) == 0 || length(cand_sup) == 0) {
    stop("Cannot infer columns from parseRM summary.")
  }
  divcol <- cand_div[1]
  supcol <- cand_sup[1]
  dt$bin <- floor(as.numeric(dt[[divcol]])) # 0..100 bins
  dt$superfamily <- toupper(as.character(dt[[supcol]]))
  dt <- dt %>% filter(!is.na(bin))
  dt <- dt %>% count(superfamily, bin, name="count")
  dt
}

if (file.exists(divsum_file)) {
  ds <- read_divsum(divsum_file)
} else if (file.exists(summary_file)) {
  ds <- read_summary(summary_file)
} else {
  stop("No parseRM divsum or summary found in: ", land)
}

# Clean superfamily labels into broad groups
map_sup <- function(x){
  x <- toupper(x)
  x <- gsub("_INT$", "", x)
  x <- ifelse(grepl("GYPSY", x), "LTR/GYPSY",
       ifelse(grepl("COPIA", x), "LTR/COPIA",
       ifelse(grepl("LINE", x), "LINE",
       ifelse(grepl("SINE", x), "SINE",
       ifelse(grepl("DNA|TIR|CACTA|MULE|HELITRON", x), "DNA/TIR",
       ifelse(grepl("RC/HELITRON|HELITRON", x), "RC/HELITRON",
              "OTHERS"))))))
  x
}
ds$group <- map_sup(ds$superfamily)

# Compute normalized abundance per group (per-bin fraction)
ds_sum <- ds %>% group_by(group, bin) %>% summarize(n = sum(count), .groups="drop")
totals <- ds_sum %>% group_by(group) %>% summarize(total = sum(n), .groups="drop")
ds_norm <- ds_sum %>% left_join(totals, by="group") %>% mutate(frac = ifelse(total>0, n/total, 0))

# Convert divergence (%) to age (Myr)
r <- opt$rate
ds_norm <- ds_norm %>%
  mutate(K = pmax(bin, 0)/100,
         age_years = K/(2*r),
         age_Myr = age_years/1e6)

# Save binned table
write_csv(ds_norm, file.path(land, "TE_landscape_bins.tsv"))

groups_order <- c("LTR/GYPSY","LTR/COPIA","DNA/TIR","LINE","SINE","RC/HELITRON","OTHERS")

p1 <- ggplot(ds_norm %>% filter(group %in% groups_order),
             aes(x = bin, y = frac, color = group)) +
  geom_line(size=0.8) +
  facet_wrap(~ factor(group, levels=groups_order), scales="free_y", ncol=2) +
  labs(x = "Divergence (%)",
       y = "Abundance (fraction per group)",
       title = "TE landscape by superfamily (divergence)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

# Secondary axis for age (approximate mapping via linear transform on K=bin/100)
# age = (bin/100)/(2r) -> Myr
scale_factor <- (1.0/100.0)/(2*r)/1e6
p1 <- p1 + scale_x_continuous(sec.axis = sec_axis(~ . * scale_factor, name = "Approx. age (Myr)"))

ggsave(file.path(land, "TE_landscape_by_superfamily.pdf"), p1, width=10, height=8)
ggsave(file.path(land, "TE_landscape_by_superfamily.png"), p1, width=10, height=8, dpi=200)

p2 <- ggplot(ds_norm, aes(x = bin, y = frac)) +
  geom_area(fill="grey60") +
  labs(x = "Divergence (%)",
       y = "Abundance (fraction, all TEs)",
       title = "TE landscape (overall)") +
  theme_bw(base_size = 12) +
  scale_x_continuous(sec.axis = sec_axis(~ . * scale_factor, name = "Approx. age (Myr)"))

ggsave(file.path(land, "TE_landscape_overall.pdf"), p2, width=8, height=5)
ggsave(file.path(land, "TE_landscape_overall.png"), p2, width=8, height=5, dpi=200)

message("Saved plots and tables into: ", land)