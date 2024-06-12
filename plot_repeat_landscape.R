#!/usr/bin/env Rscript

#
# Plot Repeat Landscape
#

# Set environment
library(ggplot2)
library(dplyr)

# Set working directory 
setwd("~/path/to/working/directory")

# Load the data containing all species and clusters
targets_f <- './cluster_sizes.txt'
target_seqs <- read.delim(targets_f)

# For troubleshooting purposes
# target_seqs <- target_seqs %>%
#   filter(hb == 'MN') %>%
#   filter(spp == 'treBer')

# Default vars
min.length <- 10

# Loop over the targets
for (i in 1:nrow(target_seqs)){
  row <- target_seqs[i,]
  # The target species-cluster pair
  spp <- paste(row$spp, '_', row$hb, sep='')
  # length of the sequence
  assm.len <- row$bp
  # Report
  message(paste('Working on', spp, '...'))
  
  # Directory variables
  work.dir <- paste('./species_database/', spp, sep='')
  input_f <- paste(work.dir, '/', spp, '.divsum.tsv', sep='')
  
  # Read input file
  rep.annot <- read.delim(input_f)
  message('    Read input file with ', nrow(rep.annot), ' elements')
  
  # Define classes for filtering
  target.classes <- c('DNA', 'LINE', 'LTR', 'SINE', 'Unknown')
  rna.classes <- c('tRNA', 'rRNA', 'scRNA', 'snRNA')
  remove.classes <- c('Satellite')
  final.classes <- c('DNA', 'LTR', 'LINE', 'SINE', 'small RNA',
                     'Other', 'Unclass.')
  
  # Filter the file for:
  message('    Filtering input file')
  rep.annot <- rep.annot %>%
    # Create a repeat family column
    mutate(Family = Class) %>%
    # Modify the class to remove family info (everything before the "/")
    mutate(Class = sapply(strsplit(Class, '/'), function(x){x[1]})) %>%
    # Remove satellite
    filter(!Class %in% remove.classes) %>%
    # Convert to a general RNA class
    mutate(Class = replace(Class, Class %in% rna.classes, 'small RNA')) %>%
    # Convert non-target classes to Other
    mutate(Class = replace(Class, !Class %in% target.classes, 'Other')) %>%
    # Rename Unknown to Unclassified
    mutate(Class = replace(Class, Class=='Unknown', 'Unclass.')) %>%
    # Larger than a min length
    filter(wellCharLen >= min.length) %>%
    # Positive Kimura divergence only
    filter(Kimura >= 0)
  
  # Create binned repeats classes proportional to genome size
  message('    Binning repeats')
  binned.reps <- rep.annot %>%
    # Select only some columns
    select(Class, wellCharLen, Kimura) %>%
    # Remove repeated rows
    distinct() %>%
    # Convert Kimura to percentages
    mutate(Kimura = Kimura*100) %>%
    # Create bins of the Kimura percentages, bins=1
    mutate(Kimura.bin = cut(Kimura, breaks=seq(-1,max(Kimura)+1,1),
                            labels=FALSE)) %>%
    # Add some rows manually to ensure all elements are present in the legend
    add_row(Kimura = rep(0,length(final.classes)),
            Class = final.classes,
            wellCharLen = rep(0,length(final.classes))) %>%
    # Group the bins by each Repeat Class
    group_by(Kimura.bin, Class) %>%
    # And summarize the lengths of the elements
    summarize(Len = sum(wellCharLen), .groups='drop_last') %>%
    # Calculate a proportion based on the size of the genome
    mutate(Prop.Genome = (Len/assm.len)*100) %>%
    # Reorder some factors for plotting
    mutate(Class = factor(Class, levels=final.classes))
  
  
  # Make the plot
  message('    Making plot')
  fig <- ggplot(data=binned.reps, aes(x=Kimura.bin, y=Prop.Genome,
                                      fill=Class, color=Class)) +
    geom_bar(position="stack", stat="identity", alpha=0.9, linewidth=0.1) +
    xlim(0,50) +
    scale_y_continuous(limits=c(0,50), breaks=seq(0,50,10)) +
    scale_fill_brewer(palette='Dark2', drop=FALSE) +
    scale_color_brewer(palette='Dark2', drop=FALSE) +
    labs(title=spp, x='Kimura Divergence %', y='Locus Coverage %') + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  print(fig)
  
  # Export the PDF
  pdf_f <- paste(work.dir, '/', spp, '.repeat_landscape.pdf', sep='')
  ggsave(pdf_f, plot=fig, device='pdf', width=4, height=3)
}
