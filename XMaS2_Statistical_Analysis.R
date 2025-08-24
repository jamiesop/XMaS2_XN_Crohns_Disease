
# Project: Xanthohumol Microbiome and Signature (XMaS) 2 Clinical Trial
# Description: Placebo-controlled RCT with XN in Crohn's Disease
# Author: Paige Jamieson


# Environment ---------------------

library(tidyverse)
library(phyloseq)
library(vegan)
library(microbiome)
library(nortest)
library(lmerTest)
library(mixOmics)
library(ggplot2)
library(stats)
library(GLMMadaptive)
library(multcomp)


# Functions ---------------------

# Helper function to impute missing values (imputes mean of timepoints before and after missing value)
NA_impute_helper <- function(x) {
        res <- rep(NA_real_,length(x))
        for (i in seq_along(x)) 
                res[i] = ifelse(is.na(x[i]), mean(c(x[i-1], x[i+1])), x[i])
        return(res)
}

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

#Help functions:
log_helper <- function(x, min.val){
        log((x + sqrt(x^2 + min.val^2))/2)
}
#Pareto Scaling:
PS_helper <- function(x){
        (x - mean(x))/sqrt(sd(x, na.rm = T))
}	

#Transformation Functions:
#Log Scaling:
log_transform <- function(mtb){
        mtb_nz <- mtb[ ,which(apply(mtb, 2, sum) != 0)]
        min.val <- min(abs(mtb_nz[mtb_nz!=0]))/10
        mtb_log_trans <- apply(mtb_nz, 2, log_helper, min.val)
        return(mtb_log_trans)
}

#Pareto Scaling:
pareto_scale <- function(mtb){
        mtb_scaled <- apply(mtb, 2, PS_helper) 
        return(mtb_scaled)
}

#Function to find min value not zero
nzmin <- function(x){
        min(x[x > 0 ])
}

#Log2FC function to use on a dataframe
l2fcmin <- function(df, counts, var, condA, condB) {
        #Replace 0 values with 0.1% of min val
        df[[counts]] <- replace(df[[counts]], df[[counts]] <= 0, nzmin(df[[counts]]*0.001))
        c1 <- df[[counts]][sapply(df[[var]], function(x) x == condA)]
        c2 <- df[[counts]][sapply(df[[var]], function(x) x == condB)]
        l2fc <- mean(log2(c2[sapply(c2, is.finite)])) - mean(log2(c1[sapply(c1, is.finite)]))
        return(l2fc)
}

# Helper function to pull participant beta diversity values from distance matrix
bdiv_helper <- function(df, sample, part.var, time.var) {
        # Convert to symbols for tidy evaluation
        sample <- ensym(sample)
        part.var <- ensym(part.var)
        time.var <- ensym(time.var)
        
        # Initialize dataframe output
        out_df <- df %>% dplyr::select(sample, part.var, time.var) %>% mutate(bdiv = NA)
        
        # Create empty vector for bdiv values
        bdiv.vec <- c()
        
        # Loop through each participant and append bdiv values to a vector
        for (i in unique(df %>% pull(part.var))) {
                
                col.var <- paste0(i, 'v2')
                bdivvals <- df[[col.var]][sapply(df[[part.var]], function(x) x == i)]
                bdiv.vec <- c(bdiv.vec, bdivvals)
                
        }
        
        out_df$bdiv <- bdiv.vec
        return(out_df)
}



# Load Datasets --------------------

# DEMOGRAPHICS for Table 1
demo <- read.csv("~/Documents/Xmas 2/XMaS2_DEMO.csv")

# METADATA
meta <- read.csv("~/Documents/Xmas 2/Data Analysis/meta_full.csv") %>% 
        mutate_at('treatment', ~ ifelse(treatment == 'X', 'B', 'A')) %>% 
        mutate_at('treatment', ~ factor(.x, levels = c('A', 'B'))) %>% 
        mutate_at('time', ~ factor(.x, levels = c('v2', 'v3', 'v4', 'v5', 'v6')))

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



# GUT MICROBIOME -------------------------

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

# New Phyloseq object
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

# Taxa Table
tax <- ps_counts %>% 
        tax_table() %>% 
        data.frame() %>% 
        rownames_to_column('ASV')




## Participant Characteristics & Demographics ----------------------------------

# Create Demographics Table
TB1 <- demo %>% 
        # remove participant lost to attrition
        filter(part.id != 'p218') %>% 
        rename(sex = "Biological.sex..0...Female..1...Male.") %>% 
        group_by(group) %>% 
        summarize(mAge = mean(Age),
                  sdAge = sd(Age),
                  pSex = sum(sex == 0)/ sum(sex == 0 | sex == 1),
                  mCDAI = mean(cdai_base), 
                  sdCDAI = sd(cdai_base),
                  mBMI = mean(BMI),
                  sdBMI = sd(BMI),
                  mHt = mean(Height_meters),
                  sdHt = sd(Height_meters),
                  mWt = mean(Weight_kg),
                  sdWt = sd(Weight_kg),
                  pEth.His = sum(Ethnicity == 'Hispanic/Latino'),
                  pEth.Non = sum(Ethnicity == 'Non-Hispanic/Latino'),
                  pRac.White = sum(Race == 'White/Caucasian'),
                  pRac.Asian = sum(Race == 'Asian'),
                  pRac.Black = sum(Race == 'Black/African American'),
                  pFem = sum(sex == 0)
        )

# One-Way ANOVA test for continuous demographic variables
anov.summ <- demo %>% 
        filter(part.id != 'p218') %>% 
        rename(sex = "Biological.sex..0...Female..1...Male.") %>% 
        pivot_longer(cols = c('Age', 'cdai_base', 'BMI', 'Height_meters', 'Weight_kg'),
                     names_to = 'demos', values_to = 'measurement') %>% 
        group_by(demos) %>% 
        nest() %>% 
        mutate(aov.mod = purrr::map(data, function(x) aov(x$measurement ~ x$group))) %>% 
        mutate(pval = purrr::map_dbl(aov.mod, function(x) summary(x)[[1]][[1,5]]))

# Chi-square test for proportion Ethnicity (Non-hispanic/Hispanic)
eth.tab <- TB1 %>% 
        column_to_rownames('group') %>% 
        dplyr::select(starts_with('pEth'))
chisq.test(eth.tab) # pval = 0.17

# Chi-square test for proportion Race (Asian/Black/White)
rac.tab <- TB1 %>% 
        column_to_rownames('group') %>% 
        dplyr::select(starts_with('pRac'))
chisq.test(rac.tab) # pval = 0.37

# Chi-square test for proportion female
fem.tab <- TB1 %>% 
        column_to_rownames('group') %>% 
        dplyr::select(pFem)
chisq.test(fem.tab) # pval = 0.37



# CDAI scores ------------------------------------------------------------------

cdaidata <- left_join(meta, cdai) %>% 
        # Remove participant lost to attrition
        dplyr::filter(part.id != 'p218') %>% 
        # Impute missing values with mean of measurement before and after (p213v4) 
        mutate_at('cdai_score', ~ zero_impute_helper(.x))

# Normality test (null = data follow a normal distribution)
cdai_norm <- cdaidata %>% 
        group_by(treatment, time) %>% 
        summarize(normtest = nortest::ad.test(cdai_score)$p.value)
# passes normality test

# Assess QQ Plots
qqnorm(cdaidata$cdai_score)

# LMER Model
cdaiLMM <- lmerTest::lmer(cdai_score ~ treatment*time + (1|part.id), data = cdaidata)

# Assess Contrasts
cdai_contrasts <- data.frame(
        v2_BvA = lmerTest::contest(cdaiLMM, L = c(0,1,0,0,0,0,0,0,0,0))[[6]],
        v3_BvA = lmerTest::contest(cdaiLMM, L = c(0,1,0,0,0,0,1,0,0,0))[[6]],
        v4_BvA = lmerTest::contest(cdaiLMM, L = c(0,1,0,0,0,0,0,1,0,0))[[6]],
        v5_BvA = lmerTest::contest(cdaiLMM, L = c(0,1,0,0,0,0,0,0,1,0))[[6]],
        v6_BvA = lmerTest::contest(cdaiLMM, L = c(0,1,0,0,0,0,0,0,0,1))[[6]],
        B_v2v3 = lmerTest::contest(cdaiLMM, L = c(0,0,1,0,0,0,1,0,0,0))[[6]],
        B_v2v4 = lmerTest::contest(cdaiLMM, L = c(0,0,0,1,0,0,0,1,0,0))[[6]],
        B_v2v5 = lmerTest::contest(cdaiLMM, L = c(0,0,0,0,1,0,0,0,1,0))[[6]],
        B_v2v6 = lmerTest::contest(cdaiLMM, L = c(0,0,0,0,0,1,0,0,0,1))[[6]]
)

# Biomarkers for Inflammation and Gut Barrier ----------------------------------

# Join Data
eldata <- left_join(meta, elisas) %>% 
        # Remove participant that did not complete study
        filter(part.id != 'p218')

# Check normality
el.normtest <- eldata %>% 
        pivot_longer(starts_with('EL'), names_to = 'EL') %>% 
        group_by(treatment, time, EL) %>% 
        nest() %>% 
        mutate(norm = purrr::map(data, function(x) ad.test(x$value ))) %>%
        mutate(pval = purrr::map_dbl(norm, function(x) x$p.value))

# LMER Model
els.lmm <- eldata %>% 
        mutate(across(starts_with('EL'), ~ log(.x + 1))) %>% 
        pivot_longer(starts_with('EL'), names_to = 'EL') %>% 
        group_by(EL) %>% 
        nest() %>% 
        mutate(model = purrr::map(data, function(x) 
                lmerTest::lmer(value ~ treatment*time + (1|part.id), data = x)))

# Assess Contrasts
els.contrasts <- els.lmm %>% 
        mutate(v2_BvA = map_dbl(model, function(x) 
                as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,0,0)))$`Pr(>F)`)) %>% 
        mutate(v3_BvA = map_dbl(model, function(x) 
                as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,1,0,0,0)))$`Pr(>F)`)) %>%
        mutate(v4_BvA = map_dbl(model, function(x) 
                as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,1,0,0)))$`Pr(>F)`)) %>%
        mutate(v5_BvA = map_dbl(model, function(x) 
                as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,1,0)))$`Pr(>F)`)) %>%
        mutate(v6_BvA = map_dbl(model, function(x) 
                as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,0,1)))$`Pr(>F)`)) %>%
        mutate(B_v2v3 = map_dbl(model, function(x) 
                as.data.frame(lmerTest::contest(x, L = c(0,0,1,0,0,0,1,0,0,0)))$`Pr(>F)`)) %>%
        mutate(B_v2v4 = map_dbl(model, function(x) 
                as.data.frame(lmerTest::contest(x, L = c(0,0,0,1,0,0,0,1,0,0)))$`Pr(>F)`)) %>%
        mutate(B_v2v5 = map_dbl(model, function(x) 
                as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,1,0,0,0,1,0)))$`Pr(>F)`)) %>% 
        mutate(B_v2v6 = map_dbl(model, function(x) 
                as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,0,1,0,0,0,1)))$`Pr(>F)`)) %>% 
        ungroup() %>% 
        mutate(across(c(everything(), -EL, -data, -model), ~ p.adjust(.x, method = 'BH')))



# Gut Microbiome & Metabolites -------------------------------------------------


# Alpha Diversity Analysis -----------------------------------------------------

# Determine alpha diversity metrics
adiv <- psnew %>% 
        estimate_richness(measures = c("Observed", "Shannon", "Simpson")) %>% 
        rownames_to_column("sample.id") %>% 
        left_join(meta) %>% 
        # Remove participant lost to follow up
        filter(part.id != 'p218')

# Test Normality
adiv_sum <- adiv %>% 
        pivot_longer(cols = c(Observed, Shannon, Simpson), 
                     names_to = 'measures', values_to = 'values') %>%
        group_by(treatment, measures, time) %>% 
        nest() %>% 
        mutate(norm.test = purrr::map(data, function(x) (ad.test(x$values)))) %>% 
        mutate(pval = purrr::map_dbl(norm.test, function(x) x$p.value))
# Normality failed at a couple time points for Simpson (3) and shannon (1)

# LMER Model & contrasts
LMMadiv <- adiv %>% 
        mutate(Simpson = log(Simpson),
               Shannon = log(Shannon)) %>% 
        pivot_longer(cols = c(Observed, Shannon, Simpson), 
                     names_to = 'measures', values_to = 'values') %>% 
        group_by(measures) %>% 
        nest() %>% 
        mutate(mod = purrr::map(data, function(x) 
                lmerTest::lmer(values ~ treatment*time + (1|part.id), data = x))) %>% 
        mutate(T1_BvA = map_dbl(mod, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,0,0)))$`Pr(>F)`)) %>% 
        mutate(T2_BvA = map_dbl(mod, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,1,0,0,0)))$`Pr(>F)`)) %>%
        mutate(T3_BvA = map_dbl(mod, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,1,0,0)))$`Pr(>F)`)) %>%
        mutate(T4_BvA = map_dbl(mod, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,1,0)))$`Pr(>F)`)) %>%
        mutate(T5_BvA = map_dbl(mod, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,0,1)))$`Pr(>F)`)) %>%
        mutate(T2_BvB = map_dbl(mod, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,1,0,0,0,1,0,0,0)))$`Pr(>F)`)) %>%
        mutate(T3_BvB = map_dbl(mod, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,1,0,0,0,1,0,0)))$`Pr(>F)`)) %>%
        mutate(T4_BvB = map_dbl(mod, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,1,0,0,0,1,0)))$`Pr(>F)`)) %>% 
        mutate(T5_BvB = map_dbl(mod, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,0,1,0,0,0,1)))$`Pr(>F)`)) %>% 
        ungroup() %>% 
        mutate(across(c(everything(), -data, -measures, -mod), ~ p.adjust(.x, method = 'BH')))
# No significant differences



# Beta Diversity (Aitchison Distance) ----------------------------


# Pull meta data
mdata <- sample_data(ps_counts) %>% 
        data.frame() %>% 
        rownames_to_column("sample.id")

# Treatment VS Control
ps_counts <- ps_counts %>% 
        # Remove participant lost to follow up
        subset_samples(part.id != 'p218')


# Time: baseline

ct_v2 <- ps_counts %>% 
        subset_samples(time == "v2")

clr <- ct_v2 %>% 
        rarefy_even_depth(rngseed = 11) %>% 
        otu_table() %>% 
        as.data.frame() %>% 
        decostand(method = 'clr', pseudocount = 1)

euclid_dist <- clr %>% 
        stats::prcomp()

summary(euclid_dist)$importance[2,] 

comps <- euclid_dist$x %>% 
        as.data.frame() %>% 
        rownames_to_column('sample.id') %>% 
        left_join(., mdata)


# PERMANOVA
distmat <- stats::dist(clr, diag=TRUE)
meta <- data.frame(sample_data(ct_v2))
bdisp <- betadisper(distmat, meta$treatment)
# Evaluate using a permutation test
permutest(bdisp) # p = 0.57 (treatment)
set.seed(123)
adonis2(distmat~treatment, data = meta) 
# p = 0.70


# Time: week 2

ct_v3 <- ps_counts %>% 
        subset_samples(time == "v3")

clr <- ct_v3 %>% 
        rarefy_even_depth(rngseed = 11) %>% 
        otu_table() %>% 
        as.data.frame() %>% 
        decostand(method = 'clr', pseudocount = 1)

euclid_dist <- clr %>% 
        stats::prcomp()

summary(euclid_dist)$importance[2,] # 12.9% and 11.1%

comps <- euclid_dist$x %>% 
        as.data.frame() %>% 
        rownames_to_column('sample.id') %>% 
        left_join(., mdata)


# PERMANOVA
distmat <- stats::dist(clr, diag=TRUE)
meta <- data.frame(sample_data(ct_v3))
bdisp <- betadisper(distmat, meta$treatment)
# Evaluate using a permutation test
permutest(bdisp) # p = 0.44 (treatment)
set.seed(123)
adonis2(distmat~treatment, data = meta) 
# p = 0.92


# Time: week 4

ct_v4 <- ps_counts %>% 
        subset_samples(time == "v4")

clr <- ct_v4 %>% 
        rarefy_even_depth(rngseed = 11) %>% 
        otu_table() %>% 
        as.data.frame() %>% 
        decostand(method = 'clr', pseudocount = 1)

euclid_dist <- clr %>% 
        stats::prcomp()

summary(euclid_dist)$importance[2,] 

comps <- euclid_dist$x %>% 
        as.data.frame() %>% 
        rownames_to_column('sample.id') %>% 
        left_join(., mdata)


# PERMANOVA
distmat <- stats::dist(clr, diag=TRUE)
meta <- data.frame(sample_data(ct_v4))
bdisp <- betadisper(distmat, meta$treatment)
# Evaluate using a permutation test
permutest(bdisp) # p = 0.73 (treatment)
set.seed(123)
adonis2(distmat~treatment, data = meta) 
# p = 0.86


# Time: week 6

ct_v5 <- ps_counts %>% 
        subset_samples(time == "v5")

clr <- ct_v5 %>% 
        rarefy_even_depth(rngseed = 11) %>% 
        otu_table() %>% 
        as.data.frame() %>% 
        decostand(method = 'clr', pseudocount = 1)

euclid_dist <- clr %>% 
        stats::prcomp()

summary(euclid_dist)$importance[2,] 

comps <- euclid_dist$x %>% 
        as.data.frame() %>% 
        rownames_to_column('sample.id') %>% 
        left_join(., mdata)


# PERMANOVA
distmat <- stats::dist(clr, diag=TRUE)
meta <- data.frame(sample_data(ct_v5))
bdisp <- betadisper(distmat, meta$treatment)
# Evaluate using a permutation test
permutest(bdisp) # p = 0.44 (treatment)
set.seed(123)
adonis2(distmat~treatment, data = meta) 
# p = 0.31


# Time: week 8

ct_v6 <- ps_counts %>% 
        subset_samples(time == "v6")

clr <- ct_v6 %>% 
        rarefy_even_depth(rngseed = 11) %>% 
        otu_table() %>% 
        as.data.frame() %>% 
        decostand(method = 'clr', pseudocount = 1)

euclid_dist <- clr %>% 
        stats::prcomp()

summary(euclid_dist)$importance[2,] 

comps <- euclid_dist$x %>% 
        as.data.frame() %>% 
        rownames_to_column('sample.id') %>% 
        left_join(., mdata)


# PERMANOVA
distmat <- stats::dist(clr, diag=TRUE)
meta <- data.frame(sample_data(ct_v6))
bdisp <- betadisper(distmat, meta$treatment)
# Evaluate using a permutation test
permutest(bdisp) # p = 0.67(treatment)
set.seed(123)
adonis2(distmat~treatment, data = meta) 
# p = 0.118



# LM Model Beta Diversity by Time ---------------

# Aitchison Distance
ait.dist <- ps_counts %>% 
        rarefy_even_depth(rngseed = 11) %>% 
        otu_table() %>% 
        as.data.frame() %>% 
        vegdist(method = 'aitchison', pseudocount = 1) %>% 
        as.matrix() %>% 
        data.frame() %>% 
        rownames_to_column('sample.id') %>% 
        mutate(part.id = str_split(sample.id, 'v', simplify = T)[,1]) %>% 
        mutate(time = str_split(sample.id, pattern = 's2[0-9]{2}', simplify = T)[,2]) 


# Helper function to pull participant beta diversity values from distance matrix
bdiv_helper <- function(df, sample, part.var, time.var) {
        # Convert to symbols for tidy evaluation
        sample <- ensym(sample)
        part.var <- ensym(part.var)
        time.var <- ensym(time.var)
        # Initialize dataframe output
        out_df <- df %>% dplyr::select(sample, part.var, time.var) %>% mutate(bdiv = NA)
        # Create empty vector for bdiv values
        bdiv.vec <- c()
        # Loop through each participant and append bdiv values to a vector
        for (i in unique(df %>% pull(part.var))) {
                col.var <- paste0(i, 'v2')
                bdivvals <- df[[col.var]][sapply(df[[part.var]], function(x) x == i)]
                bdiv.vec <- c(bdiv.vec, bdivvals)
        }
        out_df$bdiv <- bdiv.vec
        return(out_df)
}

# Reorder
ait.dist.reorder <- bdiv_helper(ait.dist, sample.id, part.id, time) %>% 
        mutate(part.id = str_replace(part.id, 's', 'p')) %>% 
        left_join(meta %>% dplyr::select(sample.id, treatment))


# Fit LMM to Beta Diversity ~ Treatment * Time
ait.dist.mod <- lmerTest::lmer(bdiv ~ treatment*time + (1|part.id), data = ait.dist.reorder, REML = F)

# Assess Contrasts
AITcontsTvP <- data.frame(
        v2_TvP = as.data.frame(lmerTest::contest(ait.dist.mod, L = c(0,1,0,0,0,0,0,0,0,0)))$`Pr(>F)`,
        v3_TvP = as.data.frame(lmerTest::contest(ait.dist.mod, L = c(0,1,0,0,0,0,1,0,0,0)))$`Pr(>F)`,
        v4_TvP = as.data.frame(lmerTest::contest(ait.dist.mod, L = c(0,1,0,0,0,0,0,1,0,0)))$`Pr(>F)`,
        v5_TvP = as.data.frame(lmerTest::contest(ait.dist.mod, L = c(0,1,0,0,0,0,0,0,1,0)))$`Pr(>F)`,
        v6_TvP = as.data.frame(lmerTest::contest(ait.dist.mod, L = c(0,1,0,0,0,0,0,0,0,1)))$`Pr(>F)`
        
) %>% 
        pivot_longer(everything(), names_to = 'comparisons', values_to = 'pvals') %>% 
        mutate(across(pvals, ~ p.adjust(.x, method = 'BH')))
# No Significant differences



# Differential Abundance Analysis ----------------

# clr-transform counts
clr <- ps_counts %>% 
        rarefy_even_depth(rngseed = 11) %>% 
        otu_table() %>% 
        as.data.frame() %>% 
        decostand(method = 'clr', pseudocount = 1) %>% 
        rownames_to_column('sample.id') %>% 
        left_join(meta) %>% 
        # Remove participant lost to follow up
        filter(part.id != 'p218')

# Pivot Long
clr_long <- clr %>% 
        pivot_longer(cols = starts_with('ASV'), names_to = 'ASV', values_to = 'clr')

# LMER Model 
clr_lm <- clr_long %>% 
        group_by(ASV) %>% 
        nest() %>% 
        mutate(model = purrr::map(data, function(x) lmerTest::lmer(clr ~ treatment*time + (1|part.id), data = x) )) %>% 
        mutate(T1_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,0,0)))$`Pr(>F)`)) %>% 
        mutate(T2_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,1,0,0,0)))$`Pr(>F)`)) %>%
        mutate(T3_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,1,0,0)))$`Pr(>F)`)) %>%
        mutate(T4_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,1,0)))$`Pr(>F)`)) %>%
        mutate(T5_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,0,1)))$`Pr(>F)`)) %>%
        mutate(T2_BvB = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,1,0,0,0,1,0,0,0)))$`Pr(>F)`)) %>%
        mutate(T3_BvB = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,1,0,0,0,1,0,0)))$`Pr(>F)`)) %>%
        mutate(T4_BvB = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,1,0,0,0,1,0)))$`Pr(>F)`)) %>% 
        mutate(T5_BvB = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,0,1,0,0,0,1)))$`Pr(>F)`)) %>% 
        ungroup() %>% 
        mutate(across(c(everything(), -ASV, -data, -model), ~ p.adjust(.x, method = 'BH'))) %>% 
        left_join(., tax)
# No Significant differences



# Short-chain Fatty Acids ------------------

# Join Data
fa_tidy <- left_join(meta, scfa) %>% 
        # Remove participant lost to follow up
        filter(part.id != 'p218')

# Test Normality
fa_norm <- fa_tidy %>% 
        pivot_longer(starts_with('FA'), names_to = 'FA', values_to = 'conc') %>% 
        group_by(FA, treatment, time) %>% 
        nest() %>% 
        mutate(norm_test = purrr::map(data, function(x) nortest::ad.test(x$conc)),
               pval = purrr::map_dbl(norm_test, function(x) x$p.value))

# LMER Model
fa_LMM <- fa_tidy %>% 
        mutate(across(c(starts_with('FA'), -FA_acetic_acid), ~ log(.x))) %>% 
        pivot_longer(cols = starts_with('FA'), names_to = 'FA', values_to = 'conc') %>% 
        group_by(FA) %>% 
        nest() %>% 
        dplyr::mutate(model = purrr::map(data, function(x) 
                lmerTest::lmer(conc ~ treatment*time + (1|part.id), data = x)))

# Assess Contrasts
fa.contrasts <- fa_LMM %>% 
        mutate(v2_BvA = map_dbl(model, function(x) 
                as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,0,0)))$`Pr(>F)`)) %>% 
        mutate(v3_BvA = map_dbl(model, function(x) 
                as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,1,0,0,0)))$`Pr(>F)`)) %>%
        mutate(v4_BvA = map_dbl(model, function(x) 
                as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,1,0,0)))$`Pr(>F)`)) %>%
        mutate(v5_BvA = map_dbl(model, function(x) 
                as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,1,0)))$`Pr(>F)`)) %>%
        mutate(v6_BvA = map_dbl(model, function(x) 
                as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,0,1)))$`Pr(>F)`)) %>%
        mutate(v3_BvB = map_dbl(model, function(x) 
                as.data.frame(lmerTest::contest(x, L = c(0,0,1,0,0,0,1,0,0,0)))$`Pr(>F)`)) %>%
        mutate(v4_BvB = map_dbl(model, function(x) 
                as.data.frame(lmerTest::contest(x, L = c(0,0,0,1,0,0,0,1,0,0)))$`Pr(>F)`)) %>%
        mutate(v5_BvB = map_dbl(model, function(x) 
                as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,1,0,0,0,1,0)))$`Pr(>F)`)) %>% 
        mutate(v6_BvB = map_dbl(model, function(x) 
                as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,0,1,0,0,0,1)))$`Pr(>F)`)) %>% 
        ungroup() %>% 
        mutate(across(c(everything(), -FA, -data, -model), ~ p.adjust(.x, method = 'BH')))

# acetic acid between groups q = 0.039 at close-out



# BILE ACIDS ----------------


# Join Data (individual bile acids)
badata <- left_join(meta, biles) %>% 
        # Remove participant that did not complete the study
        dplyr::filter(part.id != 'p218')

# Test Normality
normcheck <- badata %>% 
        #remove DHCA because conc = 0 in many individuals and 
        #throwing error when trying to perform norm stat test
        dplyr::select(-BA_DHCA) %>% 
        pivot_longer(cols = starts_with("BA_"), names_to = "BA", values_to = "conc") %>% 
        group_by(BA, treatment, time) %>% 
        nest() %>% 
        mutate(norm_test = purrr::map(data, function(x) nortest::ad.test(x$conc)),
               pval = purrr::map_dbl(norm_test, function(x) x$p.value))

# LMER Model - individuals bile acids
ba.LMM <- badata %>% 
        mutate(across(starts_with('BA'), ~ log( .x + 1 ))) %>% 
        pivot_longer(cols = starts_with("BA_"), names_to = "BA", values_to = "conc") %>% 
        group_by(BA) %>% 
        nest() %>% 
        dplyr::mutate(model = purrr::map(data, function(x) 
                lmerTest::lmer(conc ~ treatment*time + (1|part.id), data = x)))

# Assess Contrasts
ba.contrasts <- ba.LMM %>% 
        mutate(v2_BvA = map_dbl(model, function(x) 
                as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,0,0)))$`Pr(>F)`)) %>% 
        mutate(v3_BvA = map_dbl(model, function(x) 
                as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,1,0,0,0)))$`Pr(>F)`)) %>%
        mutate(v4_BvA = map_dbl(model, function(x) 
                as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,1,0,0)))$`Pr(>F)`)) %>%
        mutate(v5_BvA = map_dbl(model, function(x) 
                as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,1,0)))$`Pr(>F)`)) %>%
        mutate(v6_BvA = map_dbl(model, function(x) 
                as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,0,1)))$`Pr(>F)`)) %>%
        mutate(v3_BvB = map_dbl(model, function(x) 
                as.data.frame(lmerTest::contest(x, L = c(0,0,1,0,0,0,1,0,0,0)))$`Pr(>F)`)) %>%
        mutate(v4_BvB = map_dbl(model, function(x) 
                as.data.frame(lmerTest::contest(x, L = c(0,0,0,1,0,0,0,1,0,0)))$`Pr(>F)`)) %>%
        mutate(v5_BvB = map_dbl(model, function(x) 
                as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,1,0,0,0,1,0)))$`Pr(>F)`)) %>% 
        mutate(v6_BvB = map_dbl(model, function(x) 
                as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,0,1,0,0,0,1)))$`Pr(>F)`)) %>% 
        ungroup() %>% 
        mutate(across(c(everything(), -BA, -data, -model), ~ p.adjust(.x, method = 'BH')))

# Borderline within-group non-significant difference in DCA at close-out (q = 0.06)





# BILE ACID GROUPS --------------------

# Create BA groups based on chemical similarity 
BA_groups <- badata %>% 
        # Total
        mutate(TOTAL = rowSums(across(c('BA_CA':'BA_MCA')))) %>% 
        # Primary unconjugated
        mutate(PR.UNC = rowSums(across(c(BA_CA, BA_CDCA)))) %>% 
        # Primary conjugated0
        mutate(PR.CON = rowSums(across(c(BA_GCA, BA_TCA, BA_GCDCA, BA_TCDCA)))) %>% 
        # Secondary unconjugated
        mutate(SEC.UNC = rowSums(across(c(BA_DCA, BA_LCA)))) %>% 
        # Secondary conjugated
        mutate(SEC.CON = rowSums(across(c(BA_GDCA, BA_TDCA, BA_GLCA, BA_GLCA)))) %>% 
        # Glycine conjugates
        mutate(GLY = rowSums(across(c(BA_GCA,, BA_GCDCA, BA_GUDCA, BA_GLCA, BA_GDCA)))) %>% 
        # Taurine conjugates
        mutate(TAU = rowSums(across(c(BA_TCA, BA_TCDCA, BA_TUDCA)))) %>% 
        dplyr::select(sample.id, part.id, treatment, time, 'TOTAL':'TAU')

# Test Normality
BA.G.norm <- BA_groups %>% 
        pivot_longer(cols = 'PR.UNC':'TAU.FL', names_to = "BA", values_to = "conc") %>% 
        group_by(BA, treatment, time) %>% 
        nest() %>% 
        mutate(norm_test = purrr::map(data, function(x) nortest::ad.test(x$conc)),
               pval = purrr::map_dbl(norm_test, function(x) x$p.value))

# LMER Model
ba_groups_LMM <- BA_groups %>% 
        mutate(across(c( 'TOTAL':'TAU' ), ~ log(.x + 1))) %>% 
        pivot_longer(cols = c('TOTAL':'TAU'), names_to = 'BA', values_to = 'conc') %>% 
        group_by(BA) %>% 
        nest() %>% 
        dplyr::mutate(model = purrr::map(data, function(x) 
                lmerTest::lmer(conc ~ treatment*time + (1|part.id), data = x))) 

# Assess Contrasts
ba.g.contrasts <- ba_groups_LMM %>% 
        mutate(v2_BvA = map_dbl(model, function(x) 
                as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,0,0)))$`Pr(>F)`)) %>% 
        mutate(v3_BvA = map_dbl(model, function(x) 
                as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,1,0,0,0)))$`Pr(>F)`)) %>%
        mutate(v4_BvA = map_dbl(model, function(x) 
                as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,1,0,0)))$`Pr(>F)`)) %>%
        mutate(v5_BvA = map_dbl(model, function(x) 
                as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,1,0)))$`Pr(>F)`)) %>%
        mutate(v6_BvA = map_dbl(model, function(x) 
                as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,0,1)))$`Pr(>F)`)) %>%
        mutate(v3_BvB = map_dbl(model, function(x) 
                as.data.frame(lmerTest::contest(x, L = c(0,0,1,0,0,0,1,0,0,0)))$`Pr(>F)`)) %>%
        mutate(v4_BvB = map_dbl(model, function(x) 
                as.data.frame(lmerTest::contest(x, L = c(0,0,0,1,0,0,0,1,0,0)))$`Pr(>F)`)) %>%
        mutate(v5_BvB = map_dbl(model, function(x) 
                as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,1,0,0,0,1,0)))$`Pr(>F)`)) %>% 
        mutate(v6_BvB = map_dbl(model, function(x) 
                as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,0,1,0,0,0,1)))$`Pr(>F)`)) %>% 
        ungroup() %>% 
        mutate(across(c(everything(), -BA, -data, -model), ~ p.adjust(.x, method = 'BH'))) 

# Significant within-group difference in secondary bile acids at close-out (q = 0.032)




# Differences between Responders & Non-responders ------------------------------


# CDAI Score stratified by Inflammatory Status ---------------------------------

# Join data sets:
metab.data <- left_join(meta, xn) %>% 
        left_join(cdai) %>% 
        left_join(scfa) %>% 
        left_join(biles) %>% 
        left_join(elisas)

# Create variable for Severe Inflammation and No/Low-grade inflammation
inflam.summ <- metab.data %>% 
        pivot_longer(starts_with('EL_'), names_to = 'marker', values_to = 'conc') %>% 
        group_by(part.id, marker, treatment) %>% 
        summarize(mInflam = mean(conc)) %>% 
        pivot_wider(names_from = 'marker', values_from = 'mInflam') %>% 
        dplyr::select(part.id, treatment, EL_IL6, EL_IL10, EL_TNF) %>% 
        # "High" inflammation indicates severe vs low-grade or no inflammaiton
        mutate(inflam.status = ifelse(c(EL_IL6 <= 30 & EL_TNF <= 30 & EL_IL10 <= 30), "low", "severe"))

# Summary statistics of inflammatory markers at baseline
inf_stats <-  metab.data %>% 
        left_join(inflam.summ %>% dplyr::select(part.id, inflam.status)) %>% 
        filter(time == 'v2') %>% 
        pivot_longer(starts_with('EL_'), names_to = 'marker', values_to = 'conc') %>% 
        group_by(inflam.status, marker) %>% 
        nest() %>% 
        mutate(mod = purrr::map(data, function(x) t.test(data = x, conc ~ treatment))) %>% 
        mutate(pval = purrr::map_dbl(mod, function(x) x$p.value))
# By inflammatory category, no significant difference in inflammatory cytokines 
# between XN and placebo at baseline
# -- suggests inflammatory status is driving difference between XN and placebo at week 4

# Model CDAI score by placebo vs treatment in severe inflammation group
resp_lm <- metab.data.scaled %>% 
        left_join(inflam.summ %>% dplyr::select(part.id, inflam.status)) %>%
        mutate(res.group = paste0(treatment, '_', inflam.status)) %>%
        dplyr::select(sample.id, part.id, treatment, time, inflam.status, res.group, cdai_score, EL_IL6, EL_IL10, EL_TNF) %>% 
        filter(inflam.status == 'severe') %>% 
        filter(part.id != 'p218') %>% 
        filter(!c(part.id == 'p220' & time == 'v5')) %>%
        pivot_longer(starts_with('EL'), names_to = 'marker', 
                     values_to = 'conc') %>% 
        group_by(marker) %>% 
        nest() %>% 
        mutate(mod = purrr::map(data, function(x) 
                lmerTest::lmer(cdai_score ~ treatment*time + (1|part.id), data = x))) %>% 
        mutate(v2_BvA = map_dbl(mod, function(x) 
                as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,0,0)))$`Pr(>F)`)) %>% 
        mutate(v3_BvA = map_dbl(mod, function(x) 
                as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,1,0,0,0)))$`Pr(>F)`)) %>%
        mutate(v4_BvA = map_dbl(mod, function(x) 
                as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,1,0,0)))$`Pr(>F)`)) %>%
        mutate(v5_BvA = map_dbl(mod, function(x) 
                as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,1,0)))$`Pr(>F)`)) %>%
        mutate(v6_BvA = map_dbl(mod, function(x) 
                as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,0,1)))$`Pr(>F)`)) %>%
        mutate(B_v2v3 = map_dbl(mod, function(x) 
                as.data.frame(lmerTest::contest(x, L = c(0,0,1,0,0,0,1,0,0,0)))$`Pr(>F)`)) %>%
        mutate(B_v2v4 = map_dbl(mod, function(x) 
                as.data.frame(lmerTest::contest(x, L = c(0,0,0,1,0,0,0,1,0,0)))$`Pr(>F)`)) %>%
        mutate(B_v2v5 = map_dbl(mod, function(x) 
                as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,1,0,0,0,1,0)))$`Pr(>F)`)) %>% 
        mutate(B_v2v6 = map_dbl(mod, function(x) 
                as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,0,1,0,0,0,1)))$`Pr(>F)`)) %>% 
        ungroup() %>% 
        mutate(across(c(everything(), -mod, -data, -marker), ~p.adjust(.x, method = 'BH')))




# Phenotypic differences btw Responders & Non-responders ------------------------------


# Biomarkers of Inflammation and Gut Barrier --------------

# LMER Model
el_lm <- eldata %>% 
        filter(treatment == 'B') %>% 
        mutate(response = ifelse(part.id %in% c('p207', 'p209', 'p212', 'p220'), 'responder', 'nonresponder')) %>% 
        # Removed due to treatment noncompliance
        filter(!c(part.id == 'p220' & time == 'v5')) %>% 
        mutate(across(starts_with('EL'), ~ log(.x + 1))) %>% 
        pivot_longer(cols = starts_with('EL'), names_to = 'EL', values_to = 'value') %>% 
        group_by(EL) %>% 
        nest() %>% 
        mutate(model = purrr::map(data, function(x) lmerTest::lmer(value ~ response*time + (1|part.id), data = x) )) %>% 
        mutate(T1_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,0,0)))$`Pr(>F)`)) %>% 
        mutate(T2_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,1,0,0,0)))$`Pr(>F)`)) %>%
        mutate(T3_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,1,0,0)))$`Pr(>F)`)) %>%
        mutate(T4_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,1,0)))$`Pr(>F)`)) %>%
        mutate(T5_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,0,1)))$`Pr(>F)`)) %>%
        mutate(T2_BvB = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,1,0,0,0,1,0,0,0)))$`Pr(>F)`)) %>%
        mutate(T3_BvB = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,1,0,0,0,1,0,0)))$`Pr(>F)`)) %>%
        mutate(T4_BvB = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,1,0,0,0,1,0)))$`Pr(>F)`)) %>% 
        mutate(T5_BvB = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,0,1,0,0,0,1)))$`Pr(>F)`)) %>% 
        ungroup() %>% 
        mutate(across(c(everything(), -EL, -data, -model), ~ p.adjust(.x, method = 'BH')))



# Alpha Diversity -----------------------------

# LMER Model
response_adiv_lm <- adiv %>% 
        filter(treatment == 'B') %>% 
        mutate(response = ifelse(part.id %in% c('p207', 'p209', 'p212', 'p220'), 'responder', 'nonresponder')) %>% 
        # Removed due to treatment noncompliance
        filter(!c(part.id == 'p220' & time == 'v5')) %>% 
        mutate(Simpson = log(Simpson), 
               Shannon = log(Shannon)) %>% 
        pivot_longer(cols = c(Observed, Shannon, Simpson), names_to = 'measures', values_to = 'values') %>% 
        group_by(measures) %>% 
        nest() %>% 
        mutate(model = purrr::map(data, function(x) lmerTest::lmer(values ~ response*time + (1|part.id), data = x))) %>% 
        mutate(T1_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,0,0)))$`Pr(>F)`)) %>% 
        mutate(T2_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,1,0,0,0)))$`Pr(>F)`)) %>%
        mutate(T3_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,1,0,0)))$`Pr(>F)`)) %>%
        mutate(T4_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,1,0)))$`Pr(>F)`)) %>%
        mutate(T5_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,0,1)))$`Pr(>F)`)) %>%
        mutate(T2_BvB = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,1,0,0,0,1,0,0,0)))$`Pr(>F)`)) %>%
        mutate(T3_BvB = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,1,0,0,0,1,0,0)))$`Pr(>F)`)) %>%
        mutate(T4_BvB = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,1,0,0,0,1,0)))$`Pr(>F)`)) %>% 
        mutate(T5_BvB = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,0,1,0,0,0,1)))$`Pr(>F)`)) %>% 
        ungroup() %>% 
        mutate(across(c(everything(), -data, -measures, -model), ~ p.adjust(.x, method = 'BH')))
# No significant differences



# Differential Abundance -------------------

ct_treat <- ps_counts %>% 
        subset_samples(treatment == 'B')

mdata <- sample_data(ct_treat) %>% 
        data.frame() %>% 
        rownames_to_column("sample.id") %>% 
        mutate(response = ifelse(part.id %in% c('p207', 'p209', 'p212', 'p220'), 'responder', 'nonresponder'))

clr <- ct_treat %>% 
        rarefy_even_depth(rngseed = 11) %>% 
        otu_table() %>% 
        as.data.frame() %>% 
        decostand(method = 'clr', pseudocount = 1) %>% 
        rownames_to_column('sample.id')

# LMER Model
clr_lm <- left_join(mdata, clr) %>% 
        pivot_longer(cols = starts_with('ASV'), names_to = 'ASV', values_to = 'clr') %>% 
        # Removed due to treatment noncompliance
        filter(!c(part.id == 'p220' & time == 'v5')) %>% 
        group_by(ASV) %>% 
        nest() %>% 
        mutate(model = purrr::map(data, function(x) lmerTest::lmer(clr ~ response*time + (1|part.id), data = x) )) %>% 
        mutate(T1_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,0,0)))$`Pr(>F)`)) %>% 
        mutate(T2_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,1,0,0,0)))$`Pr(>F)`)) %>%
        mutate(T3_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,1,0,0)))$`Pr(>F)`)) %>%
        mutate(T4_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,1,0)))$`Pr(>F)`)) %>%
        mutate(T5_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,0,1)))$`Pr(>F)`)) %>%
        mutate(T2_BvB = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,1,0,0,0,1,0,0,0)))$`Pr(>F)`)) %>%
        mutate(T3_BvB = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,1,0,0,0,1,0,0)))$`Pr(>F)`)) %>%
        mutate(T4_BvB = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,1,0,0,0,1,0)))$`Pr(>F)`)) %>% 
        mutate(T5_BvB = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,0,1,0,0,0,1)))$`Pr(>F)`)) %>% 
        ungroup() %>% 
        mutate(across(c(everything(), -ASV, -data, -model), ~ p.adjust(.x, method = 'BH'))) %>% 
        left_join(., tax)

# Significant differences:
# ASV77 - coprococcus comes -- between group
# ASV29 - Lachnoclostridum -- within-group



# Short-chain fatty acids --------------------------------

# LMER Model
fa_lm <- fa_tidy %>% 
        filter(treatment == 'B') %>% 
        mutate(response = ifelse(part.id %in% c('p207', 'p209', 'p212', 'p220'), 'responder', 'nonresponder')) %>% 
        # Removed due to treatment noncompliance
        filter(!c(part.id == 'p220' & time == 'v5')) %>% 
        mutate(across(starts_with('FA'), ~ log(.x + 1))) %>%
        pivot_longer(cols = starts_with('FA'), names_to = 'FA', values_to = 'conc') %>% 
        group_by(FA) %>% 
        nest() %>% 
        mutate(model = purrr::map(data, function(x) lmerTest::lmer(conc ~ response*time + (1|part.id), data = x) )) %>% 
        mutate(T1_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,0,0)))$`Pr(>F)`)) %>% 
        mutate(T2_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,1,0,0,0)))$`Pr(>F)`)) %>%
        mutate(T3_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,1,0,0)))$`Pr(>F)`)) %>%
        mutate(T4_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,1,0)))$`Pr(>F)`)) %>%
        mutate(T5_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,0,1)))$`Pr(>F)`)) %>%
        mutate(T2_BvB = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,1,0,0,0,1,0,0,0)))$`Pr(>F)`)) %>%
        mutate(T3_BvB = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,1,0,0,0,1,0,0)))$`Pr(>F)`)) %>%
        mutate(T4_BvB = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,1,0,0,0,1,0)))$`Pr(>F)`)) %>% 
        mutate(T5_BvB = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,0,1,0,0,0,1)))$`Pr(>F)`)) %>% 
        ungroup() %>% 
        mutate(across(c(everything(), -FA, -data, -model), ~ p.adjust(.x, method = 'BH')))
# No significant differences


# Bile Acids ------------------------

# Create BA groups based on chemical similarity 
BA_groups <- badata %>% 
        mutate(BA_TOTAL = rowSums(across(c('BA_CA':'BA_MCA')))) %>% 
        mutate(BA_PR.UNC = rowSums(across(c(BA_CA, BA_CDCA)))) %>% 
        mutate(BA_PR.CON = rowSums(across(c(BA_GCA, BA_TCA, BA_GCDCA, BA_TCDCA)))) %>% 
        mutate(BA_SEC.UNC = rowSums(across(c(BA_DCA, BA_LCA)))) %>% 
        mutate(BA_SEC.CON = rowSums(across(c(BA_GDCA, BA_TDCA, BA_GLCA, BA_GLCA)))) %>% 
        mutate(BA_GLY = rowSums(across(c(BA_GCA,, BA_GCDCA, BA_GUDCA, BA_GLCA, BA_GDCA)))) %>% 
        mutate(BA_TAU = rowSums(across(c(BA_TCA, BA_TCDCA, BA_TUDCA)))) %>% 
        dplyr::select(sample.id, part.id, treatment, time, 'BA_TOTAL':'BA_TAU')

# LMER Model
ba_lm <- BA_groups %>% 
        mutate(across(starts_with('BA'), ~ log( .x + 1 ))) %>% 
        mutate(response = ifelse(part.id %in% c('p207', 'p209', 'p212', 'p220'), 'responder', 'nonresponder')) %>% 
        dplyr::filter(treatment == 'B') %>% 
        # Removed due to treatment noncompliance
        filter(!c(part.id == 'p220' & time == 'v5')) %>% 
        pivot_longer(cols = starts_with('BA'), names_to = 'BA', values_to = 'value') %>% 
        group_by(BA) %>% 
        nest() %>% 
        mutate(model = purrr::map(data, function(x) lmerTest::lmer(value ~ response*time + (1|part.id), data = x) )) %>% 
        mutate(T1_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,0,0)))$`Pr(>F)`)) %>% 
        mutate(T2_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,1,0,0,0)))$`Pr(>F)`)) %>%
        mutate(T3_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,1,0,0)))$`Pr(>F)`)) %>%
        mutate(T4_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,1,0)))$`Pr(>F)`)) %>%
        mutate(T5_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,0,1)))$`Pr(>F)`)) %>%
        mutate(T2_BvB = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,1,0,0,0,1,0,0,0)))$`Pr(>F)`)) %>%
        mutate(T3_BvB = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,1,0,0,0,1,0,0)))$`Pr(>F)`)) %>%
        mutate(T4_BvB = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,1,0,0,0,1,0)))$`Pr(>F)`)) %>% 
        mutate(T5_BvB = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,0,1,0,0,0,1)))$`Pr(>F)`)) %>% 
        ungroup() %>%
        mutate(across(c(everything(), -BA, -data, -model), ~ p.adjust(.x, method = 'BH')))


# XN Metabolites ---------------------------

# LMER Model
xn_lm <- left_join(meta, xn) %>% 
        filter(treatment == 'B') %>% 
        mutate(response = ifelse(part.id %in% c('p207', 'p209', 'p212', 'p220'), 'responder', 'nonresponder')) %>% 
        # Removed due to treatment noncompliance
        filter(!c(part.id == 'p220' & time == 'v5')) %>% 
        mutate(across(starts_with('XN'), ~ log(.x + 1))) %>%
        pivot_longer(cols = starts_with('XN'), names_to = 'XN', values_to = 'value') %>% 
        group_by(XN) %>% 
        nest() %>% 
        mutate(model = purrr::map(data, function(x) lmerTest::lmer(value ~ response*time + (1|part.id), data = x) )) %>% 
        mutate(T1_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,0,0)))$`Pr(>F)`)) %>% 
        mutate(T2_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,1,0,0,0)))$`Pr(>F)`)) %>%
        mutate(T3_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,1,0,0)))$`Pr(>F)`)) %>%
        mutate(T4_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,1,0)))$`Pr(>F)`)) %>%
        mutate(T5_BvA = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,1,0,0,0,0,0,0,0,1)))$`Pr(>F)`)) %>%
        mutate(T2_BvB = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,1,0,0,0,1,0,0,0)))$`Pr(>F)`)) %>%
        mutate(T3_BvB = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,1,0,0,0,1,0,0)))$`Pr(>F)`)) %>%
        mutate(T4_BvB = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,1,0,0,0,1,0)))$`Pr(>F)`)) %>% 
        mutate(T5_BvB = map_dbl(model, function(x) as.data.frame(lmerTest::contest(x, L = c(0,0,0,0,0,1,0,0,0,1)))$`Pr(>F)`)) %>% 
        ungroup() %>% 
        mutate(across(c(everything(), -XN, -data, -model), ~ p.adjust(.x, method = 'BH')))





# Stepwise Backward Regression: CDAI ~ phenotype features ----------------------

# Log Transform & Join data sets:
metab.data.scaled <- left_join(meta, xn %>% mutate(across(starts_with('XN'), ~ log(.x + 1)))) %>% 
        left_join(cdai) %>% 
        left_join(scfa %>% mutate(across(starts_with('FA'), ~log(.x + 1)))) %>% 
        left_join(biles %>% mutate(across(starts_with('BA'), ~ log(.x + 1)))) %>% 
        left_join(elisas %>% mutate(across(starts_with('EL'), ~ log(.x + 1))))

# Microbiome Data
ct_treat <- ps_counts %>% 
        subset_samples(treatment == 'B')

# metadata
mdata <- sample_data(ct_treat) %>% 
        data.frame() %>% 
        rownames_to_column("sample.id") %>% 
        mutate(response = ifelse(part.id %in% c('p207', 'p209', 'p212', 'p220'), 'responder', 'nonresponder'))
# clr transform
clr <- ct_treat %>% 
        rarefy_even_depth(rngseed = 11) %>% 
        otu_table() %>% 
        as.data.frame() %>% 
        decostand(method = 'clr', pseudocount = 1) %>% 
        rownames_to_column('sample.id')
# JOIN
full.data.treat <- left_join(mdata, clr) %>% 
        left_join(metab.data.scaled) %>% 
        filter(treatment == 'B') %>% 
        filter(part.id != 'p218') %>% 
        filter(!c(part.id == 'p220' & time == 'v5')) 

# Full model
fullMod <- lmerTest::lmer(cdai_score ~ BA_secondary + ASV77 + ASV57 + EL_TNF + EL_IL10 + EL_IL6 + (1|part.id),
                          data = full.data.treat, REML = F)

# backward stepwise for fixed effects
stepMod <- lmerTest::step(fullMod, reduce.random=F)

# Final model
finalMod <- lmerTest::lmer(cdai_score ~ BA_secondary + EL_IL10 + (1|part.id),
                           data = full.data.treat, REML = F)

summary(finalMod)
MuMIn::r.squaredGLMM(finalMod) # R2 = 0.34




# Stepwise Backward Regression: DELTA CDAI ~ phenotypes features -------------

# Create new variable (delta CDAI)
delta.cdai <- full.data.treat %>% 
        filter(time %in% c('v2', 'v6')) %>% 
        dplyr::select(part.id, time, cdai_score, BA_secondary, 
                      EL_IL10, EL_IL6, EL_TNF, ASV77, ASV57) %>% 
        pivot_wider(names_from = 'time', values_from = c('cdai_score', 'ASV77', 'ASV57', 
                                                         'BA_secondary', 'EL_IL10', 'EL_IL6',
                                                         'EL_TNF')) %>% 
        mutate(delt.cdai = cdai_score_v6 - cdai_score_v2)

# Full model
fullMod <- lm(delt.cdai ~ EL_IL10_v2 + EL_IL6_v2 + EL_TNF_v2 + 
                      ASV77_v2 + ASV57_v2 + BA_secondary_v2, data = delta.cdai)

# Backward stepwise
stepMod <- stats::step(fullMod, direction = 'backward')

# Final model
finalMod <- lm(delt.cdai ~ EL_TNF_v2 + ASV77_v2 + ASV57_v2, data = delta.cdai)
summary(finalMod)
# adjusted R-squared = 0.33

# BASELINE levels of plasma TNF and abundance of ASV57 and ASV77 best predictor change in CDAI score




# sPLS Analysis to integrate Microbiota & XN Metabolites ----------------------

# Join XN + metadata
xndata <- left_join(meta, xn) %>% 
        mutate(across(starts_with('XN'), ~ log(.x + 1)))

# Microbiome data (CLR transformed)
microdata <- ps_counts %>% 
        #set.seed(711)
        rarefy_even_depth(rngseed = 711) %>% 
        otu_table() %>% 
        as.data.frame() %>% 
        decostand(method = 'clr', pseudocount = 1) %>% 
        rownames_to_column('sample.id')

# Join all by sample ids (remove baseline without XN exposure)
fulldata <- left_join(xndata, microdata) %>% 
        filter(treatment == 'B' & time != 'v2') %>% 
        column_to_rownames('sample.id')


# Extract Blocks so everything aligns correctly
micro.block <- fulldata %>% dplyr::select(starts_with('ASV'))
metab.block <- fulldata %>% dplyr::select(starts_with('XN'))


# run pca on each block
pca.micro <- pca(micro.block, multilevel = as.factor(fulldata$part.id),
                 ncomp = 10, center = TRUE, scale = TRUE)
pca.xns <- pca(metab.block, multilevel = as.factor(fulldata$part.id),
               ncomp = 5, center = TRUE, scale = TRUE)

# plot the eigenvalues (explained variance per component)
plot(pca.micro)
plot(pca.xns)

# plot samples projected on pca subspace
plotIndiv(pca.micro, comp = c(1,2))
plotIndiv(pca.xns, comp = c(1,2))


# Tune sPLS Model

# starting model
spls.xns <- spls(X = micro.block, Y = metab.block, multilevel = as.factor(fulldata$part.id),
                 ncomp = 6, mode = 'regression') # regression mode allows to decomposition by repeated measures


#Tune component parameter
perf.spls.xns <- perf(spls.xns, validation = 'Mfold', folds = 10, nrepeat = 5)
plot(perf.spls.xns, criterion = 'Q2.total') # 1 component suitable - keep 2 in final model for ploting purposes

# tune keepX parameter
list.keepX <- c(seq(10, 50, 5))
list.keepY <- c(6) # include all XNs variables since there are only 6

# Model tune
tune.spls.xns <- tune.spls(X = micro.block, Y = metab.block, ncomp = 1,
                           multilevel = as.factor(fulldata$part.id),
                           test.keepX = list.keepX,
                           test.keepY = list.keepY,
                           nrepeat = 1, folds = 10, 
                           mode = 'regression', measure = 'cor'
)
plot(tune.spls.xns)

# Extract optimal number of features to keep on component 1
optimal.keepX <- tune.spls.xns$choice.keepX # keep 45
optimal.keepY <- tune.spls.xns$choice.keepY # keep 6


#Final model with tuned parameters
final.spls.xns <- spls(micro.block, metab.block,
                       multilevel = fulldata$part.id,
                       ncomp = 2, # keep 2 for graphical purposes
                       keepX = c(45, 1),
                       keepY = c(6, 1),
                       mode = 'regression')

plotIndiv(final.spls.xns,
          rep.space = 'X-variate',
          group = fulldata$time,
          pch = as.factor(fulldata$time),
          legend = TRUE)
plotIndiv(final.spls.xns,
          rep.space = 'Y-variate',
          group = fulldata$time,
          pch = as.factor(fulldata$time),
          legend = TRUE)
plotIndiv(final.spls.xns,
          rep.space = 'XY-variate',
          group = fulldata$time,
          pch = as.factor(fulldata$time),
          legend = TRUE)

plotArrow(final.spls.xns,
          group = fulldata$time)

# Assess performance
perf.spls.xns <- perf(final.spls.xns,
                      folds = 5, nrepeat = 10,
                      validation = 'Mfold',
                      dist = 'max.dist',
                      progressBar = FALSE)


# plot the stability of each feature for the first two componets ('h' refers to histogram)
plot(perf.spls.xns$features$stability.X$comp1, type = 'h',
     ylab = 'Stability',
     xlab = 'Features',
     main = '(a) Comp 1', las = 2,
     xlim = c(0, 10))
plot(perf.spls.xns$features$stability.X$comp2, type = 'h',
     ylab = 'Stability',
     xlab = 'Features',
     main = '(b) Comp 2', las = 2,
     xlim = c(0, 10))

#correlation circle
plotVar(final.spls.xns, cex = c(3,4), var.names = c(T,T))

perf.spls.xns$measures$MSEP$summary

xn.names <- c('XN', 'IXN', '8PN', '6PN', 'DXN', 'DDXN')

asv.names <- data.frame(ASV = colnames(final.spls.xns$X)) %>% 
        left_join(tax) %>% 
        mutate(taxa = ifelse(is.na(Species), 
                             paste0('g_', Genus, '_', ASV),
                             paste0('g_', Genus, '_', Species, '_', ASV)))

#Heatmap
cim(final.spls.xns, comp = 1:2, 
    xlab = "Metabolites", ylab = 'ASVs', 
    margins = c(5,21),
    scale = T,
    transpose=F,
    zoom = T,
    cutoff = 0.25,
    col.names = xn.names,
    row.names = asv.names$taxa,
    keysize = c(1,1))
















