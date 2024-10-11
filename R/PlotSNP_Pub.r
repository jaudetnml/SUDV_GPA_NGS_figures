library(tidyverse)
library(cowplot)
library(patchwork)
library(bioseq) # from Bioconductor
library(ggrepel)
library(writexl)

# Identify which protein a position belongs to
# Information output in Geneious VCF is not helpful for this
Label_protein <- function(position){
    output <- character(length(position))
    for(i in seq_along(position)){
        output[i] <- case_when(
            position[i] >= 458 & position[i] <= 2674 ~ "NP",
            position[i] >= 3138 & position[i] <= 4127 ~ "VP35",
            position[i] >= 4454 & position[i] <= 5434 ~ "VP40",
            position[i] >= 5998 & position[i] <= 8027 ~ "GP",
            position[i] >= 8441 & position[i] <= 9307 ~ "VP30",
            position[i] >= 10299 & position[i] <= 11054 ~ "VP24",
            position[i] > 11535 & position[i] <= 18167 ~ "L",
            .default = "Non-coding"
        )
    }
    # Set up a factor before returning the values to ensure the order is correct.
    output <- factor(output, levels = c("NP", "VP35", "VP40", "GP", "VP30", "VP24", "L", "Non-coding"))
    return(output)
}

# Label the location of non-coding mutations

Label_NC <- function(position){
  output <- character(length(position))
  for(i in seq_along(position)){
    output[i] <- case_when(
      position[i] < 458 ~ "3' Leader",
      position[i] > 2674 & position[i] <= 3138 ~ "NP-VP35",
      position[i] > 4127 & position[i] <= 4454 ~ "VP35-VP40",
      position[i] > 5434 & position[i] <= 5998 ~ "VP40-GP",
      position[i] > 8027 & position[i] <= 8441 ~ "GP-VP30",
      position[i] > 9307 & position[i] <= 10299 ~ "VP30-VP24",
      position[i] > 11054 & position[i] <= 11535 ~ "VP24-L",
      position[i] > 18167 ~ "5' Trailer",
      .default = "coding"
    )
  }
  # Set up a factor before returning the values to ensure the order is correct.
  output <- factor(output, levels = c("3' Leader", "NP-VP35", "VP35-VP40", "VP40-GP",
                                      "GP-VP30", "VP30-VP24", "VP24-L", "5' Trailer", "coding"))
  return(output)
}

# Retrieve a value from the info section of the VCF
# annots: vector of annotation strings
# fieldname: name of the field to retrieve
# numeric: convert to numeric?

Extract_info_field <- function(annots, fieldname, numeric = FALSE){
  strings <- str_extract(annots, paste0(";",fieldname, "=[a-zA-Z0-9->\\.]+;"))
  vals <- str_remove_all(str_split_fixed(strings, "=", 2)[, 2], ";")
  if (numeric) {
    vals <- as.numeric(vals)
  }
  return(vals)
}


# Set diverging colors for each gene (and non-coding)
# Based on the order of the factor returned by Label_protein
gene_colors <- c("forestgreen","navy","gold", "orange",
                 "steelblue", "purple","skyblue","grey92")

# Find and read all the depth files

depth_files <- list.files("../Data/Depth/",
                          pattern = " Coverage.csv", recursive = T, full.names = T)

depth <- read_csv(depth_files, id = "FileName", comment = "#")

depth_clean <- depth %>%
  mutate(Passage = str_extract(FileName, "[0-9]+"),
         Passage = if_else(is.na(Passage), "Stock", Passage),   # The file for the stock is just "Coverage.csv"
         Passage = factor(Passage, levels = c(20:0, "Stock")),  # Set up for vertical display with 0 at the top and Stock at the bottom
         Type = if_else(Passage == "Stock", "Stock", "Passage"))

# Find and read all the vcf files

variant_files <- list.files("../Data/VCF/",
                             pattern = "vcf", recursive = T, full.names = T)

variants <- read_tsv(variant_files, col_names =  c("CHR", "Position", "ID",
                                                   "Reference", "Alternate", "Quality", "Filter",
                                                   "Annotations", "Names", "Depth_nums"),
                     id = "FileName", comment = "#") %>%
  distinct() # Some mutations in GP are listed twice due to overlapping CDF, but not labelled differently

# Clean up the variant information
# Geneious outputs variable fields depending on the nature of the mutation (snp, ins, del)

variants_pre <- variants %>%
  mutate(Type = if_else(str_detect(FileName, "-p"), "Passage", "Stock")) %>%  # Define whether a VCF file was for a passage or a stock
  separate_wider_delim(FileName, "-p", names = c(NA, "PassageDot"), too_few = "align_start") %>%   # Get the passage number (along with next line)
  separate_wider_delim(PassageDot, ".", names = c("Passage", NA), too_few = "align_start") %>%
  mutate(Passage = if_else(is.na(Passage), "Stock", Passage),               # Set the passage number (or "Stock")
         Passage = factor(Passage, levels = c(as.character(20:0), "Stock")),
         Codon_Chg = Extract_info_field(Annotations, "CDNCHG")) %>%         # Get the codon change (VCF does not contain single letter aa change)
  separate_wider_delim(Codon_Chg, names = c("From", "To"), delim = "->", too_few = "align_start") %>% 
  mutate(Alt_Freq = Extract_info_field(Annotations, "VF", T),     # Get the frequency of the mutation
         ProtPos = Extract_info_field(Annotations, "CDSCN"),      # Get the codon position (aa position)
         Effect = Extract_info_field(Annotations, "PE"),          # Get the effect (Insertion, Deletion, etc.)
         Effect = if_else(is.na(Effect), "Non-coding", Effect),
         From = as.vector(seq_translate(rna(str_replace_all(From, "T", "U")))),    # Translate codons to aa, needs to be RNA so switch T to U
         To = as.vector(seq_translate(rna(str_replace_all(To, "T", "U")))),
         Prot = if_else(Effect %in% c("Non-coding", "Insertion", "Deletion", "") | is.na(Effect), NA, paste0(From, ProtPos, To)))  # For missense mutations, create the E64K format

# Get the list of position with at least 1 reported mutation
unique_pos <- sort(unique(variants$Position)) 

# Create a tibble where we have information on the depth of all samples for
# every position where a mutation was detected. if a mutation was detected in
# a passage, that info is joined. If the mutation was not detected, the mutation information
# is going to be NA.
variants_clean <- depth_clean %>%
            filter(Position %in% unique_pos) %>% 
            left_join(variants_pre, by = c("Position", "Type", "Passage")) %>%
            rename(Depth = Coverage) %>% 
            mutate(#Cases: 1- Depth >= 1000 &   Freq >= 0.05 :: Freq
                   #       2- Depth >= 1000 &   Freq < 0.05  :: **Missing** despite sufficient data (considered not there)
                   #       3- Depth >= 1000 &   Freq == NA   :: **Missing** despite sufficient data (not detected at all)
                   #       4- Depth < 1000                   :: NA  (cannot tell due to insufficient data, greyed out)
                   #Missing despite sufficient data is set to 0 for the moment and will be filtered out before plotting the heatmap
                   Depth = as.numeric(Depth),
                   Alt_Freq = if_else(Depth >= 1000 & (is.na(Alt_Freq) | Alt_Freq <= 0.05), 0, 
                                      if_else(Depth < 1000 , NA, as.numeric(Alt_Freq)))) %>%
            # First forward primer ends at position 205, so no reliable sequencing before that
            # Except for the stock, which was shotgun sequenced
            filter(Position > 205 | (Position < 205 & Type == "Stock")) %>%        
            select(Type, Passage, Depth, Position, Reference, Alternate, Alt_Freq, Prot, From, To, ProtPos)

# Spread out the variant information, specifically for locations where more than one variant exists,
# the information for passages where one or more variants is missing needs to be assigned to all the
# correct variants

variants_spread <- tibble()

for(pos in unique_pos){
  temp <- filter(variants_clean, Position == pos)
  # Frequencies < 0.05 were set to 0, so the sum should be at least that if one passage
  # was detected at that level. If the mutation(s) never got above that, ignore them
  if (sum(temp$Alt_Freq, na.rm = T) <= 0.05) next
  # Count the number of different mutations at that site
  # Need to count ref as well due to indels
  Alts <- unique(filter(temp, !is.na(Alternate), !is.na(Alt_Freq), Alt_Freq > 0) %>%
                   .$Alternate)
  ref <- unique(filter(temp, !is.na(Reference), !is.na(Alt_Freq), Alt_Freq > 0) %>%
                  .$Reference)
  
  # if there is a single mutation at the site, add the data to the output.
  # if there are more, we need to have information for all the passages for each mutation
  # Also, recheck individual mutations to meet the frequency threshold.
  if(length(Alts) == 1 & length(ref) == 1){
    temp$Nuc <- paste0(ref, pos, Alts)
    variants_spread <- bind_rows(variants_spread, filter(temp, Alternate == Alts | is.na(Alternate)))
  } else {
    for(alt in Alts){
      for(subref in ref){
        temp3 <- temp %>%
          filter(Alternate == alt | is.na(Alternate), Reference == subref | is.na(Reference)) %>%
          mutate(Alt_Freq = if_else(Alternate == alt & Reference == subref, Alt_Freq, 0, 0))
        if (sum(temp3$Alt_Freq, na.rm = T) <= 0.05) next
        temp3$Nuc <- paste0(subref, pos, alt)
        variants_spread <- bind_rows(variants_spread, temp3)
      }
    }
  }
}

variants_spread <- distinct(variants_spread)

# Output a list of the DNA mutations and their aa effect.
variants_clean %>%
  filter(!is.na(Reference)) %>% 
  transmute(Position = Position,
            Nuc = paste0(Reference, Position, Alternate),
            Prot = Prot) %>%
  distinct() %>%
  arrange(Position, Nuc) %>%
  write_csv("../Results/SNP_effects.csv")

#### Multiples only ######

# The heatmap will only report mutations seen in at least 2 passages
# The code below is based on the code for variants_spread above,
# but also checks that the mutations are present in at least 2 passages or
# in the final stock
variants_spread_Multiples <- tibble()

for(pos in unique_pos){
    temp <- filter(variants_clean, Position == pos)
    rows <- nrow(temp)
    # All frequencies below 0.05 were set to 0, so if the sum is 0.05 or more, there is at least one passage with the mutation
    if (sum(temp$Alt_Freq, na.rm = T) <= 0.05) next
    #       if there is only one passage with the mutation                               and the mutation is not found in the final stock                         then skip to the next mutation
    if ((sum(temp$Alt_Freq > 0 & !is.na(temp$Alt_Freq)) < 2) & (sum(is.na(temp$Alt_Freq[temp$Passage == "Stock"])) | sum(temp$Alt_Freq[temp$Type == "Stock"] == 0))) next
    
    # Get all unique mutations at this location
    Alts <- unique(filter(temp, !is.na(Alternate), !is.na(Alt_Freq), Alt_Freq > 0) %>%
                       .$Alternate)
    ref <- unique(filter(temp, !is.na(Reference), !is.na(Alt_Freq), Alt_Freq > 0) %>%
                      .$Reference)
    if(length(Alts) == 1 & length(ref) == 1){
        # Only one mutation and it already passed the minimum thresholds, so we add the data to the output
        temp$Nuc <- paste0(ref, pos, Alts)
        variants_spread_Multiples <- bind_rows(variants_spread_Multiples, filter(temp, Alternate == Alts | is.na(Alternate)))
    } else {
      # Multiple mutations at this position, redo the checks for each unique mutation
        for(alt in Alts){
            for(subref in ref){
                temp3 <- temp %>%
                    filter(Alternate == alt | is.na(Alternate), Reference == subref | is.na(Reference)) %>%
                    mutate(Alt_Freq = if_else(Alternate == alt & Reference == subref, Alt_Freq, 0, 0))
                if (sum(temp3$Alt_Freq, na.rm = T) <= 0.05) next
                if ((sum(temp3$Alt_Freq > 0 & !is.na(temp3$Alt_Freq)) < 2)){
                  if (sum(str_detect("Stock", temp3$Type))){
                    if (is.na(temp3$Alt_Freq[temp3$Type == "Stock"]) | temp3$Alt_Freq[temp3$Type == "Stock"] == 0) {
                      next
                    }
                  } else {
                    next
                  }
                }
                temp3$Nuc <- paste0(subref, pos, alt)
                variants_spread_Multiples <- bind_rows(variants_spread_Multiples, temp3)
            }
        }
    }
}

# Get each unique mutation that passed the filters above and label them with the appropriate protein
# This will be used to make the "gene bar" below the heatmap

nuc_order_Multiples <- select(variants_spread_Multiples, Position, Nuc) %>%
    arrange(Position, Nuc) %>%
    distinct() %>%
    mutate(Protein = Label_protein(Position),
           AnimalID = "1")

# Create the "gene bar" plot
gene_legend_Multiples <- ggplot(mutate(nuc_order_Multiples, Nuc = factor(Nuc, levels = nuc_order_Multiples$Nuc)), aes(x = Nuc, y = AnimalID, fill = Protein)) +
    theme_minimal() +
    geom_raster() +
    scale_fill_manual(values = gene_colors) +
    guides(fill = "none") +
    theme(legend.position = "bottom",
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank(),
          plot.margin = margin(t = -10))

# Create the main heatmap
Multiples_main <- variants_spread_Multiples %>%
  # Remove rows where Alt_Freq is 0 (The mutation is not there)
  # This will make those heatmap entries transparent
  # NA (insufficient depth) will show up in grey
    filter(Alt_Freq > 0 | is.na(Alt_Freq)) %>% 
    mutate(Nuc = factor(Nuc, levels = nuc_order_Multiples$Nuc)) %>%
    ggplot(aes(Nuc, Passage, fill = Alt_Freq)) +
    theme_minimal() +
    geom_tile() +
    # NA values are set to a very light gray, blocking the grid and showing where
    # data was insufficient.
    scale_fill_viridis_c(na.value = "grey98", begin = 0, end = 1, limits = c(0, 1)) +
    facet_grid(Type~., scales = "free_y", space = "free_y", switch = "y") +
    guides(fill = guide_colorbar(direction = "horizontal")) +
    theme(axis.text.x = element_text(angle = 090, vjust = 0.5, hjust = 1,
                                     size = unit(4, "point")),
          axis.text.y = element_text(size = unit(6, "point")),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          legend.title = element_text(size = unit(7, "point")),
          legend.text = element_text(size = unit(5, "point")),
          legend.key.height = unit(0.2, "cm"),
          panel.grid = element_line(linewidth = unit(0.18, "point")),
          legend.position = "top",
          plot.margin = margin(b = -1, l = -1, unit = "pt"),
          strip.text.y.left = element_text(size = unit(5, "point"), angle = 90,
                                           hjust = 0.5),
          strip.placement = "outside")

# Assemble the plots with patchwork

Multiples_figure <- Multiples_main +
    gene_legend_Multiples +
    plot_layout(nrow = 2, ncol = 1, heights = c(24, 1))

ggsave("../Results/Multiples_1000.pdf", Multiples_figure, width = 8, height = 3, device = cairo_pdf())

######## Line Graph non-coding #######

# Create the line graph showing the frequency of non-coding mutations over passages

noncoding <- variants_spread %>% 
  filter(Label_protein(Position) == "Non-coding", !is.na(Alt_Freq)) %>% 
  mutate(Nuc = tolower(Nuc),
         Passage = as.numeric(Passage),
         Passage = if_else(Passage == 22, 21, 21 - Passage),
         Region = droplevels(Label_NC(Position)), # Drops the "coding" level
         Region = factor(Region, levels = c("3' Leader", "NP-VP35", "VP35-VP40", "VP40-GP", "GP-VP30", "VP30-VP24", "VP24-L", "5' Trailer"))) %>% 
  arrange(Passage) %>% 
  group_by(Nuc) %>% 
  # Does the mutation exist in the final passage?
  mutate(LastPassage = which.max(cumsum(Alt_Freq > 0)),
         LastPassage = Passage[LastPassage] == 21) %>% 
  ungroup


# Set the mutations in increasing order, so that the colour scale runs in the leader -> trailer direction
noncoding_order <- select(noncoding, Nuc, Position) %>% 
  distinct() %>% 
  arrange(Position)

noncoding <- mutate(noncoding, Nuc = factor(Nuc, levels = noncoding_order$Nuc))

regions <-  c("3' Leader", "NP-VP35", "VP35-VP40", "VP40-GP", "GP-VP30", "VP30-VP24", "VP24-L", "5' Trailer")

# Get the frequency of each mutation in the last passage
noncod_annot <- noncoding %>% 
  filter(LastPassage) %>% 
  group_by(Nuc) %>% 
  top_n(1, Passage)

# Create a list of plots (one per region)
noncoding_indiv <- noncoding %>% 
  group_by(Region) %>% 
  group_map(~ ggplot(data = .x, aes(x = Passage, y = Alt_Freq, colour = Nuc)) +
        theme_classic() +
        geom_point(size = 0.4) +
        geom_line(aes(linetype = LastPassage), linewidth = 0.3) +
        geom_text_repel(data = filter(noncod_annot, Region == .y$Region),
                        aes(label = Nuc), size = 2) +
          labs(title = regions[.y$Region], y = "Mutation Frequency") +
        scale_linetype_manual(values = c(3, 1)) +
        scale_x_continuous(breaks = seq(0, 21, 3), labels = c("SUDV\nstock", seq(3, 18, 3), "SUDV\nAdapted")) +
        scale_y_continuous(limits = c(0, 1)) +
        guides(linetype = "none", colour = guide_legend(ncol = 3, override.aes = list(linewidth = 0.1, size = 0.2))) +
        theme(legend.position = "bottom",
              legend.direction = "horizontal",
              legend.title = element_blank(),
              legend.text = element_text(size = unit(5, "pt")),
              legend.key.size = unit(2, "mm"),
              legend.key.spacing = unit(0.5, "mm"),
              axis.text = element_text(size = unit(5, "pt")),
              axis.title = element_text(size = unit(7, "pt")),
              plot.title = element_text(size = unit(10, "pt")),
              axis.line = element_line(linewidth = unit(0.25, "mm")),
              axis.ticks = element_line(linewidth = unit(0.2, "mm")),
              plot.margin = unit(c(0, 0, 0, 0), "mm")))

# Set up the plots in a grid
noncoding_comb <- wrap_plots(noncoding_indiv, nrow = 2, byrow = T) +
  plot_annotation(tag_levels = "A") +
  plot_layout(axis_titles = "collect")
ggsave("../Results/noncoding_split.pdf", noncoding_comb, width = 6, height = 6)


######## Line Graph coding nonsynonymous #######

# Line graph of missense mutations, same logic as for non-coding

coding_mis <- variants_spread %>% 
  filter(Label_protein(Position) != "Non-coding", !is.na(Alt_Freq), !is.na(Reference), From != To) %>% 
  mutate(Passage = as.numeric(Passage),
         Passage = if_else(Passage == 22, 21, 21 - Passage),
         Region = droplevels(Label_protein(Position)), # Drops the "Non-coding" level
         Region = factor(Region, levels = c("NP", "VP35", "VP40", "GP", "VP30", "VP24", "L")),
         Prot = if_else(is.na(Prot), paste0(From, ProtPos, "del"), Prot)) %>% 
  arrange(Passage) %>% 
  group_by(Nuc) %>% 
  mutate(LastPassage = which.max(cumsum(Alt_Freq > 0)),
         LastPassage = Passage[LastPassage],
         LastPassage = LastPassage == 21) %>% 
  ungroup


coding_mis_order <- select(coding_mis, Prot, Position) %>% 
  distinct() %>% 
  arrange(Position)

coding_mis <- mutate(coding_mis, Prot = factor(Prot, levels = coding_mis_order$Prot))

regions <-  c("NP", "VP35", "VP40", "GP", "VP30", "VP24", "L")

coding_mis_annot <- coding_mis %>% 
  filter(LastPassage) %>% 
  group_by(Nuc) %>% 
  top_n(1, Passage)

coding_all <- ggplot(coding_mis, aes(Passage, Alt_Freq, colour = Prot)) +
  theme_classic() +
  geom_point() +
  geom_line(aes(linetype = LastPassage)) +
  scale_linetype_manual(values = c(3, 1)) +
  scale_x_continuous(breaks = 0:21, labels = c("SUDV\nstock", 1:20, "SUDV\nAdapted")) +
  guides(linetype = "none", colour = "none") +
  theme(legend.position = "bottom",
        legend.direction = "horizontal")

coding_mis_indiv <- coding_mis %>% 
  group_by(Region) %>% 
  group_map(~ ggplot(data = .x, aes(x = Passage, y = Alt_Freq, colour = Prot)) +
              theme_classic() +
              geom_point(size = 0.4) +
              geom_line(aes(linetype = LastPassage), linewidth = 0.3) +
              geom_text_repel(data = filter(coding_mis_annot, Region == .y$Region),
                              aes(label = Prot), size = 2) +
              labs(title = regions[.y$Region], y = "Mutation Frequency") +
              scale_linetype_manual(values = c(3, 1)) +
              scale_x_continuous(limits = c(0, 21), breaks = seq(0, 21, 3), labels = c("SUDV\nstock", seq(3, 18, 3), "SUDV\nAdapted")) +
              scale_y_continuous(limits = c(0, 1)) +
              guides(linetype = "none", colour = guide_legend(ncol = 7, override.aes = list(linewidth = 0.1, size = 0.2))) +
              theme(legend.position = "bottom",
                    legend.direction = "horizontal",
                    legend.title = element_blank(),
                    legend.text = element_text(size = unit(5, "pt")),
                    legend.key.size = unit(2, "mm"),
                    legend.key.spacing = unit(0.5, "mm"),
                    axis.text = element_text(size = unit(5, "pt")),
                    axis.title = element_text(size = unit(7, "pt")),
                    plot.title = element_text(size = unit(10, "pt")),
                    axis.line = element_line(linewidth = unit(0.25, "mm")),
                    axis.ticks = element_line(linewidth = unit(0.2, "mm")),
                    plot.margin = unit(c(0, 0, 0, 0), "mm"),
                    legend.margin =margin(0, 0, 0, 0, "mm")))
coding_misc_comb <- wrap_plots(coding_mis_indiv, ncol = 2, byrow = T) +
  plot_annotation(tag_levels = "A") +
  plot_layout(axis_titles = "collect")
ggsave("../Results/coding_missense_split.pdf", coding_misc_comb, width = 6, height = 9)

######## Line Graph coding synonymous #######

# Line graphs of synonymous mutations, same logic as non-coding

coding_syn <- variants_spread %>% 
  filter(Label_protein(Position) != "Non-coding", !is.na(Alt_Freq), !is.na(Reference), From == To) %>% 
  mutate(Passage = as.numeric(Passage),
         Passage = if_else(Passage == 22, 21, 21 - Passage),
         Region = droplevels(Label_protein(Position)),
         Region = factor(Region, levels = c("NP", "VP35", "VP40", "GP", "VP30", "VP24", "L"))) %>% 
  arrange(Passage) %>% 
  group_by(Nuc) %>% 
  mutate(LastPassage = which.max(cumsum(Alt_Freq > 0)),
         LastPassage = Passage[LastPassage],
         LastPassage = LastPassage == 21) %>% 
  ungroup


coding_syn_order <- select(coding_syn, Nuc, Position) %>% 
  distinct() %>% 
  arrange(Position)

coding_syn <- mutate(coding_syn, Nuc = factor(tolower(Nuc), levels = tolower(coding_syn_order$Nuc)))

regions <-  c("NP", "VP35", "VP40", "GP", "VP30", "VP24", "L")

coding_syn_annot <- coding_syn %>% 
  filter(LastPassage) %>% 
  group_by(Nuc) %>% 
  top_n(1, Passage)

coding_all <- ggplot(coding_syn, aes(Passage, Alt_Freq, colour = Nuc)) +
  theme_classic() +
  geom_point() +
  geom_line(aes(linetype = LastPassage)) +
  scale_linetype_manual(values = c(3, 1)) +
  scale_x_continuous(breaks = 0:21, labels = c("SUDV\nstock", 1:20, "SUDV\nAdapted")) +
  guides(linetype = "none", colour = "none") +
  theme(legend.position = "bottom",
        legend.direction = "horizontal")

coding_syn_indiv <- coding_syn %>% 
  group_by(Region) %>% 
  group_map(~ ggplot(data = .x, aes(x = Passage, y = Alt_Freq, colour = Nuc)) +
              theme_classic() +
              geom_point(size = 0.4) +
              geom_line(aes(linetype = LastPassage), linewidth = 0.3) +
              geom_text_repel(data = filter(coding_syn_annot, Region == .y$Region),
                              aes(label = Nuc), size = 2) +
              labs(title = regions[.y$Region], y = "Mutation Frequency") +
              scale_linetype_manual(values = c(3, 1)) +
              scale_x_continuous(limits = c(0, 21), breaks = seq(0, 21, 3), labels = c("SUDV\nstock", seq(3, 18, 3), "SUDV\nAdapted")) +
              scale_y_continuous(limits = c(0, 1)) +
              guides(linetype = "none", colour = guide_legend("Mutation", ncol = 7, override.aes = list(linewidth = 0.1, size = 0.2))) +
              theme(legend.position = "bottom",
                    legend.direction = "horizontal",
                    legend.title = element_blank(),
                    legend.text = element_text(size = unit(5, "pt")),
                    legend.key.size = unit(2, "mm"),
                    legend.key.spacing = unit(0.5, "mm"),
                    axis.text = element_text(size = unit(5, "pt")),
                    axis.title = element_text(size = unit(7, "pt")),
                    plot.title = element_text(size = unit(10, "pt")),
                    axis.line = element_line(linewidth = unit(0.25, "mm")),
                    axis.ticks = element_line(linewidth = unit(0.2, "mm")),
                    plot.margin = unit(c(0, 0, 0, 0), "mm"),
                    legend.margin =margin(0, 0, 0, 0, "mm")))
coding_syn_comb <- wrap_plots(coding_syn_indiv, ncol = 2, byrow = T) +
  plot_annotation(tag_levels = "A") +
  plot_layout(axis_titles = "collect")
ggsave("../Results/coding_synonymous_split.pdf", coding_syn_comb, width = 6, height = 6)


##### Write the data for the graphs ######

all_dat <- list(Missense = coding_mis,
                Synonymous = coding_syn,
                "Non-coding" = noncoding)
write_xlsx(all_dat, "../Results/LineGraph_data.xlsx")
