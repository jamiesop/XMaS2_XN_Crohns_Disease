

# Project: Xanthohumol Microbiome and Signature (XMaS) 2 Clinical Trial
# Description: Figures for XMAS2 Manuscript
# Author: Paige Jamieson

# Load Libraries --------------

library(phyloseq)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(vegan)


# Functions -------------------

# Helper function to impute missing values as NA
NA_impute_helper <- function(x) {
        res <- rep(NA_real_,length(x))
        for (i in seq_along(x)) 
                res[i] = ifelse(is.na(x[i]), mean(c(x[i-1], x[i+1])), x[i])
        return(res)
}

# Helper function to impute missing values as non-real zero 
zero_impute_helper <- function(x) {
        res <- rep(NA_real_,length(x))
        for (i in seq_along(x)) 
                res[i] = ifelse((x[i]) == 0, mean(c(x[i-1], x[i+1])), x[i])
        return(res)
}

# Funciton to calculate standard error
ser <- function(x){
        sd(x)/sqrt(length(x))
}



# Datasets -------------------

# METADATA
meta <- read.csv("~/Documents/Xmas 2/Data Analysis/meta_full.csv") %>% 
        mutate_at('treatment', ~ ifelse(treatment == 'X', 'B', 'A')) %>% 
        mutate_at('treatment', ~ factor(.x, levels = c('A', 'B'))) %>% 
        mutate_at('time', ~ factor(.x, levels = c('v2', 'v3', 'v4', 'v5', 'v6'))) %>% 
        # Remove participant that did not complete trial
        filter(part.id != 'p218')

# XN METABOLITES
xn <- read.csv("~/Documents/Xmas 2/Data Analysis/xmas2_xn_metabs_ng_per_g.csv") %>%  
        rename_with(.cols = 2:dim(.)[2], ~ paste0("XN_", .x))

# CONJUGATED XN METABOLITES
xn_conj <- read.csv("~/Documents/Xmas 2/Data Analysis/xn_conjugates_response.csv")

# CDAI SCORE
cdai <- read.csv("~/Documents/Xmas 2/Data Analysis/XMaSCrohnsDisease_CDAI clean.csv") %>% 
        pivot_longer(cols = 2:6, names_to = 'time', values_to = 'cdai_score') %>% 
        rowwise() %>% 
        mutate(time = str_split(time, "_")[[1]][[2]]) %>% 
        mutate(sample.id = paste0('s', str_split(part.id, 'p')[[1]][[2]], time))

# FATTY ACIDS
# umol/g dry stool
scfa <- read.csv("~/Documents/Xmas 2/Data Analysis/xmas2_scfa.csv") %>% 
        rename_with(.cols = 2:dim(.)[2], ~ paste0('FA_', .x))

# BILE ACIDS
# ug/g dry stool
biles <- read.csv("~/Documents/Xmas 2/Data Analysis/xmas2_bileacids_ug_per_g.csv") %>% 
        rename_with(.cols = 2:dim(.)[2], ~ paste0('BA_', .x)) %>% 
        mutate(BA_secondary = BA_DCA + BA_LCA)

# INFLAMMATORY BIOMARKERS
elisas <- read.csv("~/Documents/Xmas 2/Data Analysis/xmas2_ELISAs.csv") %>% 
        rename_with(.cols = 2:dim(.)[2], ~ paste0('EL_', .x)) %>% 
        # Change negative values from platereader to zero
        mutate(across(.cols = 2:dim(.)[2], ~ ifelse(.x < 0, 0, .x))) %>% 
        # Impute NAs as mean of participant values at time point before and after missing value
        mutate(across(.cols = 2:dim(.)[2], ~ NA_impute_helper(.x)))



# 16S rRNA Sequencing Data --------

# Load phyloseq object
ps_raw <- readRDS("~/Documents/16S/xmas2/xmas2_ps.rds") 

# Change treatment variable for consistency with other datasets
mdata <- sample_data(ps_raw) %>% 
        data.frame() %>% 
        rownames_to_column("sample.id") %>% 
        mutate(treatment = ifelse(treatment == "X", "B", "A")) %>% 
        mutate(part.id = paste0('p', part.id)) %>% 
        rename(sex = sex..0..female..1...male.) %>% 
        mutate(across(c(part.id, treatment, time), ~ factor(.x))) %>% 
        column_to_rownames('sample.id')
asv <- otu_table(ps_raw) 
tax <- tax_table(ps_raw)

# New Phyloseq object with adjusted metadata
psnew <- phyloseq(otu_table(asv),
                  tax_table(tax),
                  sample_data(mdata))



# Create placeholder names for un-annotated genera using family
renames <- rownames(tax_table(psnew)[is.na(tax_table(psnew)[, 'Genus'])])
taxdf <- tax_table(psnew)[renames,]
renamed_genus <- unname(sapply(taxa_names(taxdf), function(x) paste0('f_', taxdf[x, 'Family'], '_', x)))
# Add Family placeholder names
tax_table(psnew)[renames, 'Genus'] <- renamed_genus
# Create placeholder names for un-annotated genera using order
rename_order <- tax_table(psnew) %>% data.frame() %>% filter(., grepl("f_NA", Genus)) %>% rownames()
taxadf <- tax_table(psnew)[rename_order,]
renamed_from_order <- unname(sapply(taxa_names(taxadf), function(x) paste0('o_', taxadf[x, 'Order'], '_', x)))
# Add Order placeholder names 
tax_table(psnew)[rename_order, 'Genus'] <- renamed_from_order




# Filter low abundance taxa - keep at ASV level

# Filter out taxa seen fewer than 5 times in less than 15% of samples
ps_counts <- psnew %>% filter_taxa(function(x) sum(x > 3) > (0.20*length(x)), TRUE)
# retained 218 ASVs

ps_counts <- ps_counts %>% filter_taxa(function(x) mean(x) > 1e-5, TRUE) # ASVs = 110

# Transform to relative abundance
psrelab <- ps_counts %>% transform_sample_counts(function(x) x / sum(x) )

# Taxa table
tax <- ps_counts %>% 
        tax_table() %>% 
        data.frame() %>% 
        rownames_to_column('ASV')


## ---------------------------------------------------------------------------##
## FIGURE 1: Crohn's Disease Activity Index (CDAI) score between intervention arms.

# CDAI scores
cdaidata <- left_join(meta, cdai) %>% 
        dplyr::filter(part.id != 'p218') %>% 
        # Impute missing values with mean of measurement before and after (p213v4)
        mutate_at('cdai_score', ~ zero_impute_helper(.x)) %>% 
        filter(!c(part.id == 'p220' & time == 'v5'))

# Mean and Std error statistics
cdai_stats <- cdaidata %>% 
        group_by(treatment, time) %>% 
        summarize(mScore = mean(cdai_score),
                  SER = ser(cdai_score))

# Create palette for XN vx placebo
xanpal <- c('#F78D5C', '#55C7BE')
names(xanpal) <- c('A', 'B')


# CDAI score plot
cdaiplot <- ggplot(cdai_stats, aes(x = time, y = mScore, color = treatment, group = treatment)) +
        geom_point(size = 2.2, position = position_dodge(width=0.2)) +
        geom_path(size = 2.2, show.legend = F, position = position_dodge(width=0.2), alpha = 1) +
        geom_errorbar(aes(ymin = mScore - SER, ymax = mScore + SER), position = position_dodge(width=0.2), 
                      show.legend = F, width = 0.2, size = 2.2) +
        geom_point(data = cdaidata, aes(time, cdai_score, color = treatment, group = part.id), show.legend = F, alpha = 0.3, size = 0.2) +
        geom_path(data = cdaidata, aes(time, cdai_score, color = treatment, group = part.id), show.legend = F, alpha = 0.5, size = 1.5, linetype = 8) +
        theme_minimal() +
        scale_color_manual(values = xanpal, labels = c('Placebo', 'XN'),
                           name = '') +
        scale_x_discrete(labels = c(0, 2, 4, 6, 8)) +
        labs(x = "week", y = "CDAI score", title = '') +
        theme(axis.line = element_line(size = 1),
              axis.title = element_text(size = 30),
              legend.text = element_text(size = 30),
              axis.text = element_text(size = 30),
              legend.title = element_text(size = 30),
              strip.text = element_text(size = 30),
              plot.title = element_text(hjust = 0.5, size = 30),
              legend.position = 'none',
              strip.placement = 'outside') +
        guides(colour = guide_legend(override.aes = list(size=8))) +
        geom_hline(yintercept = 150, linetype = 5, size = 1)


## Create Legend for Placebo vs XN-treated  ---------------------------
saveleg <- ggplot(cdai_stats, aes(time, mScore, fill = treatment, color = treatment)) +
        geom_point(aes(group = treatment), size = 2.2, position = position_dodge(width=0.2), show.legend = F) +
        geom_path(aes(group = treatment), size = 12, position = position_dodge(width=0.2), alpha = 1) +
        #geom_path(aes(group = treatment), size = 1) + 
        geom_errorbar(aes(ymin = mScore - SER, ymax = mScore + SER), position = position_dodge(width=0.2), 
                      show.legend = F, width = 0.3, size = 3.2) +
        #geom_point(size = 3) +
        #geom_point(shape = 21, fill = 'white', size = 3, show.legend = F) +
        scale_fill_manual(values = xanpal) +
        scale_color_manual(values = xanpal, name = 'Group:', labels = c('Placebo', 'Treatment')) +
        scale_x_discrete(labels = c('v2' = 0, 'v3' = 2, 'v4' = 4, 'v5' = 6, 'v6' = 8)) +
        theme_minimal() +
        labs(x = 'Week', y = "") +
        theme(panel.grid.minor = element_blank(), 
              panel.grid.major.x = element_blank(),
              strip.placement = "outside",
              panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
              legend.position = "top",
              axis.title = element_text(size = 30),
              axis.text = element_text(size = 30),
              legend.text = element_text(size = 24),
              legend.title = element_text(size = 24),
              strip.text = element_text(size = 30))

legend <- cowplot::get_plot_component(saveleg, 'guide-box-top', return_all = T)
legend <- cowplot::ggdraw(legend)


# SAVE 1200w x 900h
plot_grid(legend, cdaiplot, ncol = 1, rel_heights = c(1, 10))



## ---------------------------------------------------------------------------##
## FIGURE 2: Biomarkers of inflammation and gut barrier function between intervention arms.

# ELISA Biomarkers
el_long <- left_join(meta, elisas) %>% 
        left_join(cdai) %>% 
        pivot_longer(starts_with('EL'), names_to = 'EL', values_to = 'value') %>% 
        # change pg/mL to ng/mL
        mutate(across(value, ~ .x/1000)) %>% 
        mutate_at('EL', ~ factor(.x, levels = c('EL_IL1B', 'EL_IL6', 'EL_IL10',
                                                'EL_TNF', 'EL_CAL', 'EL_CD14', 
                                                'EL_FABP2', 'EL_LBP', 'EL_PCSK9'))) %>% 
        arrange(EL) 

# Mean and std error statistics
el.stats <- el_long %>% 
        group_by(treatment, time, EL) %>% 
        summarise(mMeas = mean(value),
                  SER = ser(value)) %>% 
        arrange(EL)

# Facet titles
EL_labels <- c('EL_CAL' = 'Calprotectin',
               'EL_CD14' = 'sCD14',
               'EL_FABP2' = 'I-FABP',
               'EL_LBP' = 'LBP',
               'EL_IL10' = 'IL-10',
               'EL_IL6' = 'IL-6',
               'EL_IL1B' = 'IL-1β',
               'EL_TNF' = 'TNF-α',
               'EL_PCSK9' = 'PCSK9')

# Significance annotation
ann_text <- data.frame(time = 'v3', mMeas = 2.5, EL = factor("EL_CD14" ), treatment = c('X'))

# ELISA plots
elplot <- ggplot(el.stats, aes(x = time, y = mMeas, color = treatment, group = treatment)) +
        geom_point(size = 3, position = position_dodge(width=0.2), show.legend = F) +
        geom_path(size = 3, show.legend = F, position = position_dodge(width=0.2), alpha = 1) +
        geom_errorbar(aes(ymin = mMeas - SER, ymax = mMeas + SER), position = position_dodge(width=0.2), 
                      show.legend = F, width = 0.4, size = 4) +
        geom_point(data = el_long, aes(time, value, color = treatment, group = part.id), alpha = 0.3, size = 0.2) +
        geom_path(data = el_long, aes(time, value, color = treatment, group = part.id), show.legend = F, alpha = 0.5, size = 2, linetype = 8) +
        theme_minimal()+
        #cowplot::theme_cowplot() +
        scale_color_manual(values = xanpal, labels = c('Placebo', 'XN'),
                           name = '') +
        scale_x_discrete(labels = c(0, 2, 4, 6, 8)) +
        labs(x = "week", y = "ng/mL") +
        theme(axis.line = element_line(size = 1.4),
              axis.title = element_text(size = 40),
              axis.text = element_text(size = 40),
              strip.text = element_text(size = 40),
              strip.placement = 'outside', 
              legend.position = 'none'
              ) +
        guides(colour = guide_legend(override.aes = list(size=12))) +
        facet_wrap(~EL, scales = 'free',
                   labeller = labeller(EL = EL_labels),
                   strip.position = "top") +
        geom_point(data = ann_text, aes(time, mMeas), shape = "#", size = 13, color = 'grey20')

# Create new legend
saveleg <- ggplot(cdai_stats, aes(time, mScore, fill = treatment, color = treatment)) +
        geom_point(aes(group = treatment), size = 2.2, position = position_dodge(width=0.2), show.legend = F) +
        geom_path(aes(group = treatment), size = 12, position = position_dodge(width=0.2), alpha = 1) +
        #geom_path(aes(group = treatment), size = 1) + 
        geom_errorbar(aes(ymin = mScore - SER, ymax = mScore + SER), position = position_dodge(width=0.2), 
                      show.legend = F, width = 0.3, size = 3.2) +
        #geom_point(size = 3) +
        #geom_point(shape = 21, fill = 'white', size = 3, show.legend = F) +
        scale_fill_manual(values = xanpal) +
        scale_color_manual(values = xanpal, name = 'Group:', labels = c('Placebo', 'Treatment')) +
        scale_x_discrete(labels = c('v2' = 0, 'v3' = 2, 'v4' = 4, 'v5' = 6, 'v6' = 8)) +
        theme_minimal() +
        labs(x = 'Week', y = "") +
        theme(panel.grid.minor = element_blank(), 
              panel.grid.major.x = element_blank(),
              strip.placement = "outside",
              panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
              legend.position = "top",
              axis.title = element_text(size = 30),
              axis.text = element_text(size = 30),
              legend.text = element_text(size = 38),
              legend.title = element_text(size = 38),
              strip.text = element_text(size = 30))

legend <- cowplot::get_plot_component(saveleg, 'guide-box-top', return_all = T)
legend <- cowplot::ggdraw(legend)

# SAVE [3000w x 2200h]
cowplot::plot_grid(legend, elplot, ncol = 1, rel_heights = c(1, 6))




##----------------------------------------------------------------------------##
## FIGURE 3: Signatures of gut microbiota composition and metabolism between intervention arms. 

# Alpha Diversity (Figure 3A) ---------------------

# alpha diversity measures
adiv <- psnew %>% 
        estimate_richness(measures = c("Observed", "Shannon", "Simpson")) %>% 
        rownames_to_column("sample.id") %>% 
        left_join(meta) %>% 
        # Remove participant lost to follow up
        filter(part.id != 'p218') 

# Mean and Std error statistics
adiv4plot <- adiv %>% 
        pivot_longer(cols = c(Observed, Shannon, Simpson), names_to = 'measure', values_to = 'values') %>% 
        group_by(treatment, time, measure) %>% 
        summarize(meanADIV = mean(values), 
                  serADIV = ser(values))


# alpha div plots
adiv_treat <- ggplot(adiv4plot, aes(time, meanADIV, fill = treatment, color = treatment)) +
        geom_point(aes(group = treatment), size = 2.2, position = position_dodge(width=0.2)) +
        geom_path(aes(group = treatment), size = 2, show.legend = F, position = position_dodge(width=0.2), alpha = 1) +
        #geom_path(aes(group = treatment), size = 1) + 
        geom_errorbar(aes(ymin = meanADIV - serADIV, ymax = meanADIV + serADIV), size =2, width = 0.3,
                      position = position_dodge(width=0.2)) +
        #geom_point(size = 3) +
        #geom_point(shape = 21, fill = 'white', size = 3, show.legend = F) +
        facet_wrap(~measure, scales = 'free', strip.position = 'top') +
        scale_fill_manual(values = xanpal) +
        scale_color_manual(values = xanpal, name = 'Group:', labels = c('nonresponders', 'responders')) +
        scale_x_discrete(labels = c('v2' = 0, 'v3' = 2, 'v4' = 4, 'v5' = 6, 'v6' = 8)) +
        theme_minimal() +
        labs(x = 'Week', y = "") +
        theme(panel.grid.minor = element_blank(), 
              panel.grid.major.x = element_blank(),
              strip.placement = "outside",
              panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
              legend.position = "none",
              axis.title = element_text(size = 30),
              axis.text = element_text(size = 30),
              legend.text = element_text(size = 30),
              legend.title = element_text(size = 30),
              strip.text = element_text(size = 30))



## Beta Diversity (Figure 3B) -----------------------

# Subset counts by each time point
counts_by_time <- ps_counts %>% 
        subset_samples(time == 'v6')

clr <- counts_by_time %>% 
        rarefy_even_depth(rngseed = 11) %>% 
        otu_table() %>% 
        as.data.frame() %>% 
        decostand(method = 'clr', pseudocount = 1) 

euclid_dist <- clr %>% 
        stats::prcomp()

summary(euclid_dist)$importance[2,] 

mdata <- sample_data(counts_by_time) %>% 
        data.frame() %>% 
        rownames_to_column("sample.id")

micro.comps <- euclid_dist$x %>% 
        as.data.frame() %>% 
        rownames_to_column('sample.id') %>% 
        left_join(., mdata) %>% 
        left_join(clr %>% rownames_to_column('sample.id'))



# BASELINE -----

# PERMANOVA
distmat <- stats::dist(clr, diag=TRUE)
bdisp <- betadisper(distmat, mdata$treatment)
# Evaluate using a permutation test
permutest(bdisp) # p = 0.57
set.seed(123)
adonis2(distmat~treatment, data = mdata) # p = 0.70


# SAVE 600w x 800h
beta.v2 <- ggplot(micro.comps, aes(PC1, PC2, color = treatment, group = treatment)) +
        geom_point(size = 4) + 
        stat_ellipse(linewidth = 2) +
        labs(title = 'Baseline', x = '', y = '') +
        annotate("text", x = 5, y = -20, label = 'adonis p = 0.70', size = 8) +
        theme_minimal() +
        scale_color_manual(values = xanpal) +
        theme(panel.grid.minor = element_blank(), 
              panel.grid.major.x = element_blank(),
              strip.placement = "outside",
              panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
              legend.position = "none",
              axis.title = element_text(size = 24),
              axis.text = element_text(size = 26),
              legend.text = element_text(size = 24),
              legend.title = element_text(size = 24),
              strip.text = element_text(size = 24),
              plot.title = element_text(hjust = 0.5, size = 30),
              plot.margin = unit(c(0.5,0,0,0), "cm"))


# WEEK 2 -----

# PERMANOVA
distmat <- stats::dist(clr, diag=TRUE)

bdisp <- betadisper(distmat, mdata$treatment)
# Evaluate using a permutation test
permutest(bdisp) # p = 0.44
set.seed(123)
adonis2(distmat~treatment, data = mdata) # p = 0.92


# SAVE 600w x 800h
beta.v3 <- ggplot(micro.comps, aes(PC1, PC2, color = treatment, group = treatment)) +
        geom_point(size = 4) + 
        stat_ellipse(linewidth = 2) +
        labs(title = 'Week 2', x = '', y = '') +
        annotate("text", x = 5, y = -20, label = 'adonis p = 0.92', size = 8) +
        theme_minimal() +
        scale_color_manual(values = xanpal) +
        theme(panel.grid.minor = element_blank(), 
              panel.grid.major.x = element_blank(),
              strip.placement = "outside",
              panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
              legend.position = "none",
              axis.title = element_text(size = 24),
              axis.text = element_text(size = 26),
              legend.text = element_text(size = 24),
              legend.title = element_text(size = 24),
              strip.text = element_text(size = 24),
              plot.title = element_text(hjust = 0.5, size = 30),
              plot.margin = unit(c(0.5,0,0,0), "cm"))


# WEEK 4 -----

# PERMANOVA
distmat <- stats::dist(clr, diag=TRUE)

bdisp <- betadisper(distmat, mdata$treatment)
# Evaluate using a permutation test
permutest(bdisp) # p = 0.73
set.seed(123)
adonis2(distmat~treatment, data = mdata) # p = 0.86


# SAVE 600w x 800h
beta.v4 <- ggplot(micro.comps, aes(PC1, PC2, color = treatment, group = treatment)) +
        geom_point(size = 4) + 
        stat_ellipse(linewidth = 2) +
        labs(title = 'Week 4', x = '', y = '') +
        annotate("text", x = 5, y = -20, label = 'adonis p = 0.86', size = 8) +
        theme_minimal() +
        scale_color_manual(values = xanpal) +
        theme(panel.grid.minor = element_blank(), 
              panel.grid.major.x = element_blank(),
              strip.placement = "outside",
              panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
              legend.position = "none",
              axis.title = element_text(size = 24),
              axis.text = element_text(size = 26),
              legend.text = element_text(size = 24),
              legend.title = element_text(size = 24),
              strip.text = element_text(size = 24),
              plot.title = element_text(hjust = 0.5, size = 30),
              plot.margin = unit(c(0.5,0,0,0), "cm"))


# WEEK 6 -----

# PERMANOVA
distmat <- stats::dist(clr, diag=TRUE)

bdisp <- betadisper(distmat, mdata$treatment)
# Evaluate using a permutation test
permutest(bdisp) # p = 0.44
set.seed(123)
adonis2(distmat~treatment, data = mdata) # p = 0.31


# SAVE 600w x 800h
beta.v5 <- ggplot(micro.comps, aes(PC1, PC2, color = treatment, group = treatment)) +
        geom_point(size = 4) + 
        stat_ellipse(linewidth = 2) +
        labs(title = 'Week 6', x = '', y = '') +
        annotate("text", x = 5, y = -20, label = 'adonis p = 0.31', size = 8) +
        theme_minimal() +
        scale_color_manual(values = xanpal) +
        theme(panel.grid.minor = element_blank(), 
              panel.grid.major.x = element_blank(),
              strip.placement = "outside",
              panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
              legend.position = "none",
              axis.title = element_text(size = 24),
              axis.text = element_text(size = 26),
              legend.text = element_text(size = 24),
              legend.title = element_text(size = 24),
              strip.text = element_text(size = 24),
              plot.title = element_text(hjust = 0.5, size = 30),
              plot.margin = unit(c(0.5,0,0,0), "cm"))

# WEEK 8 -----

# PERMANOVA
distmat <- stats::dist(clr, diag=TRUE)

bdisp <- betadisper(distmat, mdata$treatment)
# Evaluate using a permutation test
permutest(bdisp) # p = 0.67
set.seed(123)
adonis2(distmat~treatment, data = mdata) # p = 0.12


# SAVE 600w x 800h
beta.v6 <- ggplot(micro.comps, aes(PC1, PC2, color = treatment, group = treatment)) +
        geom_point(size = 4) + 
        stat_ellipse(linewidth = 2) +
        labs(title = 'Week 8', x = '', y = '') +
        annotate("text", x = 12, y = -22, label = 'adonis p = 0.12', size = 8) +
        theme_minimal() +
        scale_color_manual(values = xanpal) +
        theme(panel.grid.minor = element_blank(), 
              panel.grid.major.x = element_blank(),
              strip.placement = "outside",
              panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
              legend.position = "none",
              axis.title = element_text(size = 24),
              axis.text = element_text(size = 26),
              legend.text = element_text(size = 24),
              legend.title = element_text(size = 24),
              strip.text = element_text(size = 24),
              plot.title = element_text(hjust = 0.5, size = 30),
              plot.margin = unit(c(0.5,1,0,0), "cm"))

# Join PCA plots
betadiv <- cowplot::plot_grid(beta.v2, beta.v3, beta.v4, beta.v5, beta.v6, ncol = 5)




# Acetic Acid (Figure 3C) -------------------------

# Full dataset
fa_tidy <- left_join(meta, scfa) 

# Mean and std error
fa_stats <- fa_tidy %>% 
        group_by(treatment, time) %>% 
        summarise(mConc = mean(FA_acetic_acid),
                  mErr = ser(FA_acetic_acid))

# Annotate significant feature
ann_text <- data.frame(time = 'v6', mConc = 500, FA = factor("FA_acetic_acid" ), treatment = c('X'))


# longitudinal plot of acetic acid by intervention group
ace.plot <- ggplot(fa_stats, aes(x = time, y = mConc, color = treatment, group = treatment)) +
        geom_point(size = 2.2, position = position_dodge(width=0.2)) +
        geom_path(size = 2.2, show.legend = F, position = position_dodge(width=0.2), alpha = 1) +
        geom_errorbar(aes(ymin = mConc - mErr, ymax = mConc + mErr), position = position_dodge(width=0.2), 
                      show.legend = F, width = 0.3, size = 2.2) +
        geom_point(data = fa_tidy, aes(time, FA_acetic_acid, color = treatment, group = part.id), alpha = 0.3, size = 0.2) +
        geom_path(data = fa_tidy, aes(time, FA_acetic_acid, color = treatment, group = part.id), show.legend = F, alpha = 0.5, size = 1.5, linetype = 8) +
        theme_minimal() +
        scale_color_manual(values = xanpal, labels = c('Placebo', 'XN'),
                           name = '') +
        scale_x_discrete(labels = c(0, 2, 4, 6, 8)) +
        labs(x = "week", y = "µmol/g dry wt", title = 'Acetic Acid') +
        guides(colour = guide_legend(override.aes = list(size=10))) +
        theme(axis.line = element_line(size = 1.4),
              axis.title = element_text(size = 30),
              legend.text = element_text(size = 34),
              axis.text = element_text(size = 30),
              legend.title = element_text(size = 30),
              strip.text = element_text(size = 30),
              plot.title = element_text(hjust = 0.5, size = 30),
              legend.position = 'none') +
        guides(colour = guide_legend(override.aes = list(size=14))) +
        geom_point(data = ann_text, aes(time, mConc), shape = "*", size = 20, color = 'grey20')



# Secondary Bile Acids (Figure 3D) ---------------------

# Join bile acid data
badata <- left_join(meta, biles) %>% 
        mutate(SEC.UNC = rowSums(across(c(BA_DCA, BA_LCA)))) %>% 
        # change to mg/g dry wt
        mutate(across(SEC.UNC, ~ .x/1000))

# Mean and Std error statistics
sec_stats <- badata %>% 
        group_by(treatment, time) %>% 
        summarize(mConc = mean(SEC.UNC),
                  SER = ser(SEC.UNC))

# Significance annotation
ann_text <- data.frame(time = 'v6', mConc = 5, BA = factor("SEC.UNC" ), treatment = c('X'))

# SAVE 1400 x 700
sec.plot <- ggplot(sec_stats, aes(x = time, y = mConc, color = treatment, group = treatment)) +
        geom_point(size = 2.2, position = position_dodge(width=0.2)) +
        geom_path(size = 2.2, show.legend = F, position = position_dodge(width=0.2), alpha = 1) +
        geom_errorbar(aes(ymin = mConc - SER, ymax = mConc + SER), position = position_dodge(width=0.2), 
                      show.legend = F, width = 0.3, size = 2.2) +
        geom_point(data = badata, aes(time, SEC.UNC, color = treatment, group = part.id), alpha = 0.3, size = 0.2) +
        geom_path(data = badata, aes(time, SEC.UNC, color = treatment, group = part.id), show.legend = F, alpha = 0.5, size = 1.5, linetype = 8) +
        #cowplot::theme_cowplot() +
        theme_minimal()+
        scale_color_manual(values = xanpal, labels = c('Placebo', 'XN'),
                           name = '') +
        scale_x_discrete(labels = c(0, 2, 4, 6, 8)) +
        labs(x = "week", y = "mg/g dry wt", title = 'Secondary Bile Acids') +
        guides(colour = guide_legend(override.aes = list(size=10))) +
        theme(axis.line = element_line(size = 1.4),
              axis.title = element_text(size = 30),
              legend.text = element_text(size = 34),
              axis.text = element_text(size = 30),
              legend.title = element_text(size = 30),
              strip.text = element_text(size = 30),
              plot.title = element_text(hjust = 0.5, size = 30),
              legend.position = 'none') +
        guides(colour = guide_legend(override.aes = list(size=14))) +
        geom_point(data = ann_text, aes(time, mConc), shape = "#", size = 10, color = 'grey20')


# Join, legend, alpha adn beta div plots
ABdivplot <- cowplot::plot_grid(legend, adiv_treat, betadiv, cols = 1, align = 'v', 
                                rel_heights = c(1,4, 6), labels = c("", 'A', 'B'), label_size = 32)
# Join acetic acid and secondary bile acid plots
ace.sec.plot <- cowplot::plot_grid(NULL, ace.plot, NULL, sec.plot, NULL, cols = 5, rel_widths = c(1,3,1,3,1), 
                                   labels = c("","C", "", "D", ""), label_size = 32)
# JOIN all GM plots [SAVE 2800w x 2100h]
gm.plot <- cowplot::plot_grid(ABdivplot, ace.sec.plot, ncol = 1, align = 'v', rel_heights = c(2,1))




##----------------------------------------------------------------------------##
## Figure 5: CDAI score between responders and nonresponders


# Determine severe and no/low-grade inflammatory categories:
inflam.summ <-  left_join(meta, elisas) %>% 
        left_join(cdai) %>%  
        pivot_longer(starts_with('EL_'), names_to = 'marker', values_to = 'conc') %>% 
        group_by(part.id, marker, treatment) %>% 
        summarize(mInflam = mean(conc)) %>% 
        pivot_wider(names_from = 'marker', values_from = 'mInflam') %>% 
        dplyr::select(part.id, treatment, EL_IL6, EL_IL10, EL_TNF) %>% 
        # "high" inflammation indicates severe vs "low" indicates low-grade or no inflammaiton
        mutate(inflam.status = ifelse(c(EL_IL6 <= 30 & EL_TNF <= 30 & EL_IL10 <= 30), "low", "high"))

# Tidy data for plot CDAI by individual participant
inf.status.part <- left_join(meta, elisas) %>% 
        left_join(cdai) %>% 
        left_join(inflam.summ %>% dplyr::select(part.id, inflam.status)) %>%
        # Filter out for noncompliance
        filter(!c(part.id == 'p220' & time == 'v5')) %>%
        mutate(inf.group = paste0(treatment, "_", inflam.status))


# Tidy data for plot mean CDAI by inflammation status*treatment
inf.status.mean <- left_join(meta, elisas) %>% 
        left_join(cdai) %>% 
        left_join(inflam.summ %>% dplyr::select(part.id, inflam.status)) %>%
        # Filter out for noncompliance
        filter(!c(part.id == 'p220' & time == 'v5')) %>%
        mutate(inf.group = paste0(treatment, "_", inflam.status)) %>% 
        #pivot_longer(cols = c('EL_IL6', 'EL_IL10', 'EL_TNF'), names_to = 'marker', values_to = 'conc') %>% 
        group_by(treatment, inflam.status, inf.group, time) %>% 
        summarize(mCDAI = mean(cdai_score),
                  serCDAI = ser(cdai_score))

# color palette for response groups in XN-treated and placebo
response_vals <- c('cyan4', '#D1283F', '#B4BA45', '#D6773F')
#names(response_vals) <- c('responder', 'nonresponder', 'placebo-responder', 'placebo-nonresponder')
names(response_vals) <- c('B_high', 'B_low', 'A_high', 'A_low')

# Plot CDAI by inflammatory status*intervention group
resplot <- ggplot(inf.status.mean, aes(time, mCDAI, group = inf.group, color = inf.group)) +
        geom_point(size = 2.2, position = position_dodge(width=0.2)) +
        geom_path(size = 2.2, show.legend = F, position = position_dodge(width=0.2), alpha = 1) +
        geom_errorbar(aes(ymin = mCDAI - serCDAI, ymax = mCDAI + serCDAI), position = position_dodge(width=0.2),
                      show.legend = F, width = 0.2, size = 2.2) +
        geom_point(data = inf.status.part, aes(time, cdai_score, color = inf.group, group = part.id), 
                   show.legend = F, alpha = 0.3, size = 0.2) +
        geom_path(data = inf.status.part, aes(time, cdai_score, color = inf.group, group = part.id),
                  show.legend = F, alpha = 0.5, size = 1.5, linetype = 8) +
        theme_minimal() +
        scale_color_manual(values = response_vals, 
                           labels = c('B_high' = 'XN-treated (severe inflammation)', 
                                      'B_low' = 'XN-treated (low inflammation)', 
                                      'A_high' = 'Placebo (severe inflammation)', 
                                      'A_low' = 'Placebo (low inflammation)'),
                           name = 'Response Group:') +
        scale_x_discrete(labels = c(0, 2, 4, 6, 8)) +
        labs(x = "week", y = "CDAI score", title = '') +
        theme(axis.line = element_line(size = 1),
              axis.title = element_text(size = 30),
              legend.text = element_text(size = 30),
              axis.text = element_text(size = 30),
              legend.title = element_text(size = 30),
              strip.text = element_text(size = 30),
              plot.title = element_text(hjust = 0.5, size = 30),
              legend.position = 'none',
              strip.placement = 'outside') +
        guides(colour = guide_legend(override.aes = list(size=8))) +
        geom_hline(yintercept = 150, linetype = 5, size = 1)

# Create legend ------------------------------
saveleg <- ggplot(inf.status.mean, aes(time, mCDAI, group = inf.group, color = inf.group)) +
        geom_point(size = 2.2, position = position_dodge(width=0.2), show.legend = F) +
        geom_path(size = 8, position = position_dodge(width=0.2), alpha = 1) +
        geom_errorbar(aes(ymin = mCDAI - serCDAI, ymax = mCDAI + serCDAI), position = position_dodge(width=0.2),
                      show.legend = F, width = 0.2, size = 2.2) +
        #geom_point(data = inf.status.part, aes(time, cdai_score, color = inf.group, group = part.id), 
        #           show.legend = F, alpha = 0.3, size = 0.2) +
        #geom_path(data = inf.status.part, aes(time, cdai_score, color = inf.group, group = part.id),
        #          show.legend = F, alpha = 0.5, size = 1.5, linetype = 8) +
        theme_minimal() +
        scale_color_manual(values = response_vals, 
                           labels = c('B_high' = 'XN-treated (severe inflammation)', 
                                      'B_low' = 'XN-treated (low inflammation)', 
                                      'A_high' = 'Placebo (severe inflammation)', 
                                      'A_low' = 'Placebo (low inflammation)'),
                           name = 'Group:') +
        scale_x_discrete(labels = c(0, 2, 4, 6, 8)) +
        labs(x = "week", y = "CDAI score", title = '') +
        theme(axis.line = element_line(size = 1),
              axis.title = element_text(size = 30),
              legend.text = element_text(size = 20),
              axis.text = element_text(size = 30),
              legend.title = element_text(size = 20),
              strip.text = element_text(size = 30),
              plot.title = element_text(hjust = 0.5, size = 30),
              legend.position = 'top',
              strip.placement = 'outside') +
        guides(color = guide_legend(nrow = 2))

legend <- cowplot::get_plot_component(saveleg, 'guide-box-top', return_all = T)
legend <- cowplot::ggdraw(legend)

# SAVE 1200w x 900h
cowplot::plot_grid(legend, resplot, ncol = 1, rel_heights = c(1, 10))



## Create Legend ---------------------------
saveleg <- ggplot(cdai_stats, aes(time, mScore, fill = treatment, color = treatment)) +
        geom_point(aes(group = treatment), size = 2.2, position = position_dodge(width=0.2), show.legend = F) +
        geom_path(aes(group = treatment), size = 12, position = position_dodge(width=0.2), alpha = 1) +
        #geom_path(aes(group = treatment), size = 1) + 
        geom_errorbar(aes(ymin = mScore - SER, ymax = mScore + SER), position = position_dodge(width=0.2), 
                      show.legend = F, width = 0.3, size = 3.2) +
        #geom_point(size = 3) +
        #geom_point(shape = 21, fill = 'white', size = 3, show.legend = F) +
        scale_fill_manual(values = xanpal) +
        scale_color_manual(values = xanpal, name = 'Group:', labels = c('Placebo', 'Treatment')) +
        scale_x_discrete(labels = c('v2' = 0, 'v3' = 2, 'v4' = 4, 'v5' = 6, 'v6' = 8)) +
        theme_minimal() +
        labs(x = 'Week', y = "") +
        theme(panel.grid.minor = element_blank(), 
              panel.grid.major.x = element_blank(),
              strip.placement = "outside",
              panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
              legend.position = "top",
              axis.title = element_text(size = 30),
              axis.text = element_text(size = 30),
              legend.text = element_text(size = 24),
              legend.title = element_text(size = 24),
              strip.text = element_text(size = 30))

legend <- cowplot::get_plot_component(saveleg, 'guide-box-top', return_all = T)
legend <- cowplot::ggdraw(legend)

# SAVE 1200w x 900h
plot_grid(legend, cdaiplot, ncol = 1, rel_heights = c(1, 10))




##----------------------------------------------------------------------------##
## Figure 6. Significantly different features between responders that achieved 
## remission and nonresponders in the XN-treated group. 


# Join data sets:
metab.data <- left_join(meta, xn) %>% 
        left_join(cdai) %>% 
        left_join(scfa) %>% 
        left_join(biles) %>% 
        left_join(elisas)

# color palette for XN-treated responders and nonresponders
response_vals <- c('cyan4', '#D1283F')
names(response_vals) <- c('responder', 'nonresponder')


# Figure 6A: Inflammatory Markers -------------------
el.long <- metab.data %>% 
        mutate(response = ifelse(part.id %in% c('p207', 'p209', 'p212', 'p220'), 
                                 'responder', 'nonresponder')) %>% 
        filter(treatment == 'B') %>% 
        pivot_longer(cols = c('EL_IL10', 'EL_IL6', 'EL_TNF'), names_to = 'marker', values_to = 'conc') %>% 
        mutate(across(conc, ~ .x / 1000))



# mean and std error statistics
el.summ <- el.long %>% 
        group_by(marker, time, response) %>% 
        summarize(mMarker = mean(conc),
                  serMarker = ser(conc))


# Significance annotation
el_ann_text <- data.frame(time = c('v2', 'v3', 'v4', 'v5', 'v6'), 
                          mMarker = c(rep(2.7, 5), rep(0.5, 5), rep(0.5, 5)),
                          marker = c("EL_IL10", "EL_IL10", "EL_IL10", "EL_IL10", "EL_IL10",
                                     "EL_IL6", "EL_IL6", "EL_IL6", "EL_IL6", "EL_IL6", 
                                     "EL_TNF", "EL_TNF", "EL_TNF", "EL_TNF", "EL_TNF"),
                          response = rep('responder', 15)) %>% 
        mutate(across(marker, ~ factor(.x)))

# inflammatory markers plot
inf.plot <- ggplot(el.summ, aes(time, mMarker, group = response, color = response)) +
        geom_point(size = 2.2, position = position_dodge(width=0.2)) +
        geom_path(size = 2.2, show.legend = F, position = position_dodge(width=0.2), alpha = 1) +
        geom_errorbar(aes(ymin = mMarker - serMarker, ymax = mMarker + serMarker), size = 2.2, width = 0.3,
                      show.legend = F, position = position_dodge(width=0.2)) +
        facet_wrap(~marker, strip.position = "top", scales = 'free',
                   labeller = as_labeller(c('EL_IL10' = 'IL10',
                                            'EL_IL6' = 'IL6',
                                            'EL_TNF' = 'TNF-α')),
                   ncol = 3) +
        scale_color_manual(values = response_vals, name = 'Group:', labels = c('Nonresponders', 'Responders')) +
        scale_x_discrete(labels = c('v2' = 0, 'v3' = 2, 'v4' = 4, 'v5' = 6, 'v6' = 8)) +
        theme_minimal() +
        labs(x = 'week', y = "ng/mL") +
        guides(colour = guide_legend(override.aes = list(size=10))) +
        theme(axis.line = element_line(size = 1.4),
              axis.title = element_text(size = 30),
              legend.text = element_text(size = 34),
              axis.text = element_text(size = 30),
              legend.title = element_text(size = 30),
              strip.text = element_text(size = 30),
              plot.title = element_text(hjust = 0.5, size = 30),
              legend.position = 'none') +
        guides(colour = guide_legend(override.aes = list(size=14))) +
        geom_point(data = el_ann_text, aes(time, mMarker), color = 'grey25', shape = "*", size = 20)




# Figure 6B: Gut Microbiota -----------------------
metadata <- psrelab %>% 
        subset_samples(treatment == 'B') %>% 
        sample_data() %>% 
        data.frame() %>% 
        rownames_to_column('sample.id') %>% 
        mutate(response = ifelse(part.id %in% c('p207', 'p209', 'p212', 'p220'), 'responder', 'nonresponder')) 

# Extract data from phyloseq object
relabs <- psrelab %>% 
        subset_samples(treatment == 'B') %>% 
        otu_table() %>% 
        data.frame() %>% 
        rownames_to_column('sample.id') %>% 
        right_join(metadata)

# Pivot longer
relab.long <- relabs %>% 
        dplyr::select(sample.id, part.id, treatment, time, response, ASV77, ASV57) %>% 
        pivot_longer(starts_with('ASV'), names_to = 'ASV', values_to = 'relab') %>% 
        left_join(tax) %>% 
        mutate(taxa.name = ifelse(is.na(Species), 
                                  paste0('g_', Genus, '_', ASV), 
                                  paste0('g_', Genus, '_', Species, '_', ASV))) 

# Significance Annotations
relab_ann_asv77 <- data.frame(time = c('v2', 'v5', 'v6'), 
                              relab = rep(0.01, 3), 
                              taxa.name = factor(rep("g_Coprococcus_comes_ASV77", 3)),
                              response = rep('responder', 3)) %>% 
        mutate(across(response, ~ factor(.x)))

relab_ann_asv57 <- data.frame(time = 'v6',
                              relab = 0.03,
                              taxa.name = factor("g_[Ruminococcus] torques group_ASV57"),
                              response = 'responder') %>% 
        mutate(across(response, ~ factor(.x, levels = c('nonresponder', 'responder'))))

# Facet titles
asv_names <- list(
        "g_Coprococcus_comes_ASV77" = "Coprococcus comes",
        "g_[Ruminococcus] torques group_ASV57" = "Ruminococcus torques group"
)

# Relative abundance of sig ASVs
relab.plot <- ggplot(relab.long, aes(time, relab, color = response, fill = response)) +
        geom_boxplot(alpha = 0.5, , linewidth = 2.2) +
        geom_point(shape = 21, fill = 'white', size = 3, show.legend = F, position = position_dodge(width=0.75)) +
        facet_wrap(~ taxa.name, scales = 'free', ncol = 2, strip.position = 'top',
                   labeller = labeller(taxa.name = c("g_Coprococcus_comes_ASV77" = "Coprococcus comes",
                                                     "g_[Ruminococcus] torques group_ASV57" = "Ruminococcus torques group"))) +
        scale_fill_manual(values = response_vals) +
        scale_color_manual(values = response_vals, name = 'Group:', labels = c('Nonresponders', 'Responders')) +
        scale_x_discrete(labels = c('v2' = 0, 'v3' = 2, 'v4' = 4, 'v5' = 6, 'v6' = 8)) +
        theme_minimal() +
        labs(x = 'week', y = "Relative Abundance") +
        theme(axis.line = element_line(size = 1),
              axis.title = element_text(size = 30),
              legend.text = element_text(size = 34),
              axis.text = element_text(size = 30),
              legend.title = element_text(size = 30),
              strip.text = element_text(size = 30),
              plot.title = element_text(hjust = 0.5, size = 30),
              legend.position = 'none') +
        geom_point(data = relab_ann_asv77, aes(time, relab), shape = '*', 
                   color = 'grey25', size = 20) +
        geom_text(data = relab_ann_asv57, aes(time, relab, label = "#", group = response), 
                  nudge_x = 0.2, size = 10, color = 'grey25')



# Figure 6C: Secondary Bile Acids -------------------
ba.long <- metab.data %>% 
        mutate(response = ifelse(part.id %in% c('p207', 'p209', 'p212', 'p220'), 
                                 'responder', 'nonresponder')) %>% 
        mutate(SEC.UNC = rowSums(across(c(BA_DCA, BA_LCA)))) %>% 
        filter(treatment == 'B') %>% 
        # change to mg/g dry wt
        mutate(across(SEC.UNC, ~ .x/1000))

# mean and std error statistics
ba.summ <- ba.long %>% 
        group_by(time, response) %>% 
        summarize(mBA = mean(SEC.UNC),
                  serBA = ser(SEC.UNC))

# Significance Annotation
ba_ann <- data.frame(time = 'v6', mBA = 1.5, response = 'responder') 

# secondary bile acids plot
ba.plot <- ggplot(ba.summ, aes(time, mBA, group = response, color = response)) +
        geom_errorbar(aes(ymin = mBA - serBA, ymax = mBA + serBA), size = 2.2, 
                      width = 0.3, position = position_dodge(width = 0.2)) +
        geom_point(size = 2.2, show.legend = F, position = position_dodge(width = 0.2)) +
        geom_path(size = 2.2, position = position_dodge(width = 0.2)) +
        scale_color_manual(values = response_vals, name = 'Group:', labels = c('Nonresponders', 'Responders')) +
        scale_x_discrete(labels = c('v2' = 0, 'v3' = 2, 'v4' = 4, 'v5' = 6, 'v6' = 8)) +
        theme_minimal() +
        labs(x = 'week', y = "mg/g dry wt", title = 'Secondary Bile Acids') +
        theme(axis.line = element_line(size = 1),
              axis.title = element_text(size = 30),
              legend.text = element_text(size = 34),
              axis.text = element_text(size = 30),
              legend.title = element_text(size = 30),
              strip.text = element_text(size = 30),
              plot.title = element_text(hjust = 0.5, size = 30),
              legend.position = 'none') +
        geom_text(data = ba_ann, aes(time, mBA, label = "#", group = response), 
                  nudge_x = 0.2, size = 10, color = 'grey25')
        


# Create Legend ------------
leg <- ggplot(ba.summ, aes(time, mBA, group = response, color = response)) +
        geom_path(size = 12) +
        geom_point(shape = 21, fill = 'white', size = 3, show.legend = F) +
        scale_color_manual(values = response_vals, name = 'XN-treated Group:', labels = c('Nonresponder', 'Responder')) +
        scale_x_discrete(labels = c('v2' = 0, 'v3' = 2, 'v4' = 4, 'v5' = 6, 'v6' = 8)) +
        theme_minimal() +
        labs(x = 'Week', y = "µg/g") +
        theme(panel.grid.minor = element_blank(), 
              panel.grid.major.x = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
              legend.position = "bottom",
              axis.title = element_text(size = 16),
              axis.text = element_text(size = 16),
              legend.text = element_text(size = 35),
              legend.title = element_text(size = 35),
              strip.text = element_text(size = 16), 
              plot.title = element_text(size = 18, hjust = 0.5))

leg <- cowplot::get_plot_component(leg, 'guide-box-bottom', return_all = T)
leg <- cowplot::ggdraw(leg)


## Join plots
plot1 <- cowplot::plot_grid(relab.plot, ba.plot, rel_widths = c(2,1), ncol = 2, labels = c('B', 'C'), label_size = 34)

# SAVE [2800w x 1600h]
response.plot <- cowplot::plot_grid(leg, inf.plot, plot1, align = 'v', ncol = 1, rel_heights = c(1,4,4), labels = c("", 'A', ''), label_size = 34)





## ---------------------------------------------------------------------------##
## FIGURE 6B: Signature of XN metabolism in XN-treated participants. 

# XN Metabolites
xn_tidy <- left_join(meta, xn) %>% 
        filter(treatment == 'B') %>% 
        pivot_longer(cols = c(starts_with('XN_')), names_to = 'XN', values_to = 'conc') %>% 
        # change ng/g to ug/g
        mutate(across(conc, ~ .x /1000)) %>% 
        mutate(across(XN, ~ factor(.x, 
                                   levels = c('XN_XN', 'XN_IXN', 'XN_DXN', 'XN_x8PN', 'XN_DDXN', 'XN_x6PN'))))

# Mean and Std error statistics
xn_stats <- xn_tidy %>% 
        group_by(XN, time) %>% 
        summarize(mConc = mean(conc),
                  mErr = ser(conc))

# Create facet titles
xn_labels <- c('XN_XN' = 'XN',
               'XN_IXN' = 'IXN', 
               'XN_DXN' = 'DXN',
               'XN_x8PN' = '8PN',
               'XN_DDXN' = 'DDXN',
               'XN_x6PN' = '6PN')

# Create palette for each subject id
part_pal <- c('#960472', '#D60969', '#E0DD3C', '#2E294E', '#95D0EE', '#00A6A6', '#D0E562', 
              '#87E18E', '#E85548', '#9A94BC')
names(part_pal) <- c("p202", "p205", "p207", "p209", "p210", "p211", "p212", "p219", "p220")


# XN Metabolite plots
plot1 <- ggplot(xn_tidy, aes(time, conc, group = part.id, color = part.id)) +
        geom_point(size = 6) +
        geom_point(data = xn_stats, aes(time, mConc, group = XN, color = XN), 
                   show.legend = F, size = 2, alpha = 0.3, inherit.aes = F) +
        geom_path(data = xn_stats, aes(time, mConc, group = XN, color = XN), size = 4,
                  show.legend = FALSE) +
        facet_wrap(~ XN, scales = 'free', labeller = labeller(XN = xn_labels),
                   strip.position = "top", ncol = 2) + 
        labs(x = 'week', y = 'µg/g dry wt') +
        theme_minimal() +
        scale_x_discrete(labels = c(0, 2, 4, 6, 8)) +
        theme(strip.placement = 'outside') +
        scale_color_manual(values = part_pal, name = 'Subject:', guide = 'none') +
        theme(axis.line = element_line(size = 2.2),
              axis.title = element_text(size = 40),
              axis.text = element_text(size = 40),
              strip.text = element_text(size = 40),
              plot.title = element_text(hjust = 0.5, size = 40),
              legend.position = 'top',
              legend.text = element_text(size = 30)) 

# Legend for participants
leg1 <- ggplot(xn_tidy, aes(time, conc, group = part.id, color = part.id)) +
        geom_point(size = 9) +
        facet_wrap(~ XN, scales = 'free', labeller = labeller(XN = xn_labels),
                   strip.position = "top", ncol = 2) + 
        labs(x = 'week', y = 'µg/g dry wt') +
        theme_minimal() +
        scale_x_discrete(labels = c(0, 2, 4, 6, 8)) +
        theme(strip.placement = 'outside') +
        scale_color_manual(values = part_pal, name = '') +
        theme(axis.line = element_line(size = 2.2),
              axis.title = element_text(size = 40),
              axis.text = element_text(size = 40),
              strip.text = element_text(size = 40),
              plot.title = element_text(hjust = 0.5, size = 45),
              legend.position = 'top',
              legend.text = element_text(size = 25)) +
        guides(color = guide_legend(nrow =2))

# Extract legend
legend = cowplot::get_plot_component(leg1, 'guide-box-top', return_all = TRUE)
save.leg <- cowplot::ggdraw(legend)

# SAVE 1250w x 1800h
cowplot::plot_grid(save.leg, plot1, 
                   ncol = 1, rel_widths = c(1, 4),
                   rel_heights = c(1, 8))






