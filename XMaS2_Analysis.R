
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



# Functions -----

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



# Load Datasets -----

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


# Join data sets:
metab.data <- left_join(meta, xn) %>% 
        left_join(cdai) %>% 
        left_join(scfa) %>% 
        left_join(biles) %>% 
        left_join(elisas)

# Log Transform & Join data sets:
metab.data.scaled <- left_join(meta, xn %>% mutate(across(starts_with('XN'), ~ log(.x + 1)))) %>% 
        left_join(cdai) %>% 
        left_join(scfa %>% mutate(across(starts_with('FA'), ~log(.x + 1)))) %>% 
        left_join(biles %>% mutate(across(starts_with('BA'), ~ log(.x + 1)))) %>% 
        left_join(elisas %>% mutate(across(starts_with('EL'), ~ log(.x + 1))))




# GUT MICROBIOME

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


tax <- ps_counts %>% 
        tax_table() %>% 
        data.frame() %>% 
        rownames_to_column('ASV')




## TABLE 1: DEMOGRAPHICS

# Create Demographics Table
TB1 <- demo %>% 
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





# Analyses by Treatment vs Placebo --------------------


# Beta Diversity by Time -----------------------------------

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


ait.dist.reorder <- bdiv_helper(ait.dist, sample.id, part.id, time) %>% 
        mutate(part.id = str_replace(part.id, 's', 'p')) %>% 
        left_join(meta %>% dplyr::select(sample.id, treatment))
        

# Fit LMM to Beta Diversity by Response category * Time
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



bdivTvP <- ait.dist.reorder %>% 
        group_by(treatment, time) %>% 
        summarize(meanBDIV = mean(bdiv), 
                  serBDIV = ser(bdiv))


# Set color palette
xanpal <- c('#F78D5C', '#55C7BE')
names(xanpal) <- c('A', 'B')


# SAVE 
bdivxtime <- ggplot(bdivTvP, aes(time, meanBDIV, fill = treatment, color = treatment)) +
        geom_path(aes(group = treatment), size = 1) + 
        geom_errorbar(aes(ymin = meanBDIV - serBDIV, ymax = meanBDIV + serBDIV), size =1, width = 0.2) +
        geom_point(shape = 21, fill = 'white', size = 3, show.legend = F) +
        scale_fill_manual(values = xanpal) +
        scale_color_manual(values = xanpal, name = 'Group:', labels = c('Placebo', 'XN-treated')) +
        scale_x_discrete(labels = c('v2' = 0, 'v3' = 2, 'v4' = 4, 'v5' = 6, 'v6' = 8)) +
        theme_minimal() +
        labs(x = 'Week', y = "Aitchison Distance") +
        theme(panel.grid.minor = element_blank(), 
              panel.grid.major.x = element_blank(),
              strip.placement = "outside",
              panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
              legend.position = "none",
              axis.title = element_text(size = 16),
              axis.text = element_text(size = 16),
              legend.text = element_text(size = 16),
              legend.title = element_text(size = 16),
              strip.text = element_text(size = 16))







# Analyses by Responders vs Nonresponders ----------------------


# sPLS Analysis to integrate Microbiota & XN Metabolites ------------------

xndata <- left_join(meta, xn) %>% 
        mutate(across(starts_with('XN'), ~ log(.x + 1)))


microdata <- ps_counts %>% 
        #set.seed(711)
        rarefy_even_depth(rngseed = 711) %>% 
        otu_table() %>% 
        as.data.frame() %>% 
        decostand(method = 'clr', pseudocount = 1) %>% 
        rownames_to_column('sample.id')

# Join all by sample ids
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
                       ncomp = 2,
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
#feature comp      mean         sd
#1  XN_DDXN    1 0.6795890 0.04680926
#2  XN_DDXN    2 0.7509224 0.09061140
#3   XN_DXN    1 0.7144525 0.06668228
#4   XN_DXN    2 0.8125420 0.18807042
#5   XN_IXN    1 0.9928588 0.05363215
#6   XN_IXN    2 1.1594508 0.20767292
#7    XN_XN    1 0.9125406 0.04361350
#8    XN_XN    2 1.0411001 0.14708978
#9  XN_x6PN    1 0.8053840 0.04101209
#10 XN_x6PN    2 0.8517705 0.08110533
#11 XN_x8PN    1 0.5636066 0.04791971
#12 XN_x8PN    2 0.6755338 0.09880040


dev.off()

xn.names <- c('XN', 'IXN', '8PN', '6PN', 'DXN', 'DDXN')

asv.names <- data.frame(ASV = colnames(final.spls.xns$X)) %>% 
        left_join(tax) %>% 
        mutate(taxa = ifelse(is.na(Species), 
                             paste0('g_', Genus, '_', ASV),
                             paste0('g_', Genus, '_', Species, '_', ASV)))



#Heatmap
cim.ob <- cim(final.spls.xns, comp = 1:2, 
    xlab = "Metabolites", ylab = 'ASVs', 
    margins = c(5,20),
    scale = T,
    transpose=F,
    zoom = T,
    cutoff = 0.25,
    col.names = xn.names,
    row.names = asv.names$taxa,
    keysize = c(1,1))

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

#Extract matrix of correlation coeffiecients
cor.mat <- mat$mat

pheatmap::pheatmap(cor.mat)





# Responder Analysis ------





# GUT MICROBIOME -----------------------------------------------


## FIGURE ## -- BETA DIVERSITY OF GUT MICROBIOME (RESPONDERS VS NONRESPONDERS)
## CAPTION ## -- 


ct_treat <- ps_counts %>% 
        subset_samples(treatment == 'B') 

clr <- ct_treat %>% 
        rarefy_even_depth(rngseed = 11) %>% 
        otu_table() %>% 
        as.data.frame() %>% 
        decostand(method = 'clr', pseudocount = 1) 

euclid_dist <- clr %>% 
        stats::prcomp()

summary(euclid_dist)$importance[2,] 

mdata <- sample_data(ct_treat) %>% 
        data.frame() %>% 
        rownames_to_column("sample.id") %>% 
        mutate(response = ifelse(part.id %in% c('p207', 'p209', 'p212', 'p220'), 'responder', 'nonresponder'))

micro.comps <- euclid_dist$x %>% 
        as.data.frame() %>% 
        rownames_to_column('sample.id') %>% 
        left_join(., mdata) %>% 
        left_join(clr %>% rownames_to_column('sample.id'))



# BASELINE -----

# PERMANOVA
distmat <- stats::dist(clr, diag=TRUE)

bdisp <- betadisper(distmat, mdata$response)
# Evaluate using a permutation test
permutest(bdisp) # p = 0.82 
set.seed(123)
adonis2(distmat~response, data = mdata) # p = 0.17


# PC1 22.4% PC2 15.4%

# SAVE 600w x 800h
beta.v2 <- ggplot(micro.comps, aes(PC1, PC2, color = response, group = response)) +
        geom_point(size = 4) + 
        stat_ellipse(linewidth = 2) +
        labs(title = 'Baseline', x = '', y = '') +
        annotate("text", x = 15, y = -38, label = 'adonis p = 0.17', size = 10) +
        theme_minimal() +
        scale_color_manual(values = response_vals) +
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
              plot.title = element_text(hjust = 0.5, size = 30))


# WEEK 2

# PERMANOVA
distmat <- stats::dist(clr, diag=TRUE)

bdisp <- betadisper(distmat, mdata$response)
# Evaluate using a permutation test
permutest(bdisp) # p = 0.73 
set.seed(123)
adonis2(distmat~response, data = mdata) # p = 0.17


# PC1 22.1% PC2 17.0%

# SAVE 600w x 800h
beta.v3 <- ggplot(micro.comps, aes(PC1, PC2, color = response, group = response)) +
        geom_point(size = 4) + 
        stat_ellipse(linewidth = 2) +
        labs(title = 'Week 2', x = '', y = '') +
        annotate("text", x = 15, y = -38, label = 'adonis p = 0.17', size = 10) +
        theme_minimal() +
        scale_color_manual(values = response_vals) +
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
              plot.title = element_text(hjust = 0.5, size = 30))




# WEEK 4 ------------

# PERMANOVA
distmat <- stats::dist(clr, diag=TRUE)

bdisp <- betadisper(distmat, mdata$response)
# Evaluate using a permutation test
permutest(bdisp) # p = 0.796
set.seed(123)
adonis2(distmat~response, data = mdata) # p = 0.26


# PC1 22.1% PC2 17.0%

# SAVE 600w x 800h
beta.v4 <- ggplot(micro.comps, aes(PC1, PC2, color = response, group = response)) +
        geom_point(size = 4) + 
        stat_ellipse(linewidth = 2) +
        labs(title = 'Week 4', x = '', y = '') +
        annotate("text", x =10, y = -38, label = 'adonis p = 0.26', size = 10) +
        theme_minimal() +
        scale_color_manual(values = response_vals) +
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
              plot.title = element_text(hjust = 0.5, size = 30))





# WEEK 6 -----------

# PERMANOVA
distmat <- stats::dist(clr, diag=TRUE)

bdisp <- betadisper(distmat, mdata$response)
# Evaluate using a permutation test
permutest(bdisp) # p = 0.25
set.seed(123)
adonis2(distmat~response, data = mdata) # p = 0.12


# PC1 26.5% PC2 16.6%

# SAVE 600w x 800h
beta.v5 <- ggplot(micro.comps, aes(PC1, PC2, color = response, group = response)) +
        geom_point(size = 4) + 
        stat_ellipse(linewidth = 2) +
        labs(title = 'Week 6', x = '', y = '') +
        annotate("text", x = 12, y = -35, label = 'adonis p = 0.12', size = 10) +
        theme_minimal() +
        scale_color_manual(values = response_vals) +
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
              plot.title = element_text(hjust = 0.5, size = 30))




# WEEK 8 --------

# PERMANOVA
distmat <- stats::dist(clr, diag=TRUE)

bdisp <- betadisper(distmat, mdata$response)
# Evaluate using a permutation test
permutest(bdisp) # p = 0.44
set.seed(123)
adonis2(distmat~response, data = mdata) # p = 0.12


# PC1 22.1% PC2 16.8%

# SAVE 600w x 800h
beta.v6 <- ggplot(micro.comps, aes(PC1, PC2, color = response, group = response)) +
        geom_point(size = 4) + 
        stat_ellipse(linewidth = 2) +
        labs(title = 'Week 8', x = '', y = '') +
        annotate("text", x = 12, y = -35, label = 'adonis p = 0.12', size = 10) +
        theme_minimal() +
        scale_color_manual(values = response_vals) +
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
              plot.title = element_text(hjust = 0.5, size = 30))

betadiv <- cowplot::plot_grid(beta.v2, beta.v3, beta.v4, beta.v5, beta.v6, ncol = 5)


# BETA DIV by TIME -----------------------------------

ct_treat <- ps_counts %>% 
        subset_samples(treatment == 'B') 

# Bray curtis distance -------
rlb_treat <- psrelab %>% 
        subset_samples(treatment == 'B')

BCdist <- distance(rlb_treat, method = 'bray') %>% as.matrix()

AITdist <- ct_treat %>% 
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



AITdist_reorder <- bdiv_helper(AITdist, sample.id, part.id, time) %>% 
        mutate(part.id = str_replace(part.id, 's', 'p')) %>% 
        mutate(response = ifelse(part.id %in% c('p207', 'p209', 'p212', 'p220'), 'responder', 'nonresponder'))

# Fit LMM - Beta Diversity as a function of Response Category * Time
AITmod <- lmerTest::lmer(bdiv ~ response*time + (1|part.id), data = AITdist_reorder, REML = F)

# Assess Contrasts
AITconts <- data.frame(
        v2_RvN = as.data.frame(lmerTest::contest(AITmod, L = c(0,1,0,0,0,0,0,0,0,0)))$`Pr(>F)`,
        v3_RvN = as.data.frame(lmerTest::contest(AITmod, L = c(0,1,0,0,0,0,1,0,0,0)))$`Pr(>F)`,
        v4_RvN = as.data.frame(lmerTest::contest(AITmod, L = c(0,1,0,0,0,0,0,1,0,0)))$`Pr(>F)`,
        v5_RvN = as.data.frame(lmerTest::contest(AITmod, L = c(0,1,0,0,0,0,0,0,1,0)))$`Pr(>F)`,
        v6_RvN = as.data.frame(lmerTest::contest(AITmod, L = c(0,1,0,0,0,0,0,0,0,1)))$`Pr(>F)`


) %>% 
        pivot_longer(everything(), names_to = 'comparisons', values_to = 'pvals') %>% 
        mutate(across(pvals, ~ p.adjust(.x, method = 'BH')))



bdiv4plot <- AITdist_reorder %>% 
        group_by(response, time) %>% 
        summarize(meanBDIV = mean(bdiv), 
                  serBDIV = ser(bdiv))


response_vals <- c('orange', 'cyan4')
names(response_vals) <- c('responder', 'nonresponder')


# SAVE 
bdivplot <- ggplot(bdiv4plot, aes(time, meanBDIV, fill = response, color = response)) +
        geom_path(aes(group = response), size = 1) + 
        geom_errorbar(aes(ymin = meanBDIV - serBDIV, ymax = meanBDIV + serBDIV), size =1, width = 0.2) +
        geom_point(shape = 21, fill = 'white', size = 3, show.legend = F) +
        scale_fill_manual(values = response_vals) +
        scale_color_manual(values = response_vals, name = 'Group:', labels = c('nonresponders', 'responders')) +
        scale_x_discrete(labels = c('v2' = 0, 'v3' = 2, 'v4' = 4, 'v5' = 6, 'v6' = 8)) +
        theme_minimal() +
        labs(x = 'Week', y = "Aitchison Distance") +
        theme(panel.grid.minor = element_blank(), 
              panel.grid.major.x = element_blank(),
              strip.placement = "outside",
              panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
              legend.position = "none",
              axis.title = element_text(size = 16),
              axis.text = element_text(size = 16),
              legend.text = element_text(size = 16),
              legend.title = element_text(size = 16),
              strip.text = element_text(size = 16))












## FIGURE -- ALPHA DIVERSITY OF GUT MICROBIOME
## CAPTION -- 

adiv <- psnew %>% 
        estimate_richness(measures = c("Observed", "Shannon", "Simpson")) %>% 
        rownames_to_column("sample.id") %>% 
        left_join(meta) %>% 
        # Remove participant lost to follow up
        filter(part.id != 'p218') %>% 
        filter(!c(part.id == 'p220' & time == 'v5')) %>% 
        mutate(response = ifelse(part.id %in% c('p207', 'p209', 'p212', 'p220'), 
                                 'responder', 'nonresponder'))

response_adiv_lm <- adiv %>% 
        filter(treatment == 'B') %>% 
        filter(part.id != 'p220') %>% 
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


# significant or borderline significant changes in all adiv measures in week 6 compared to baseline 
# -- interestingly this corresponds to a dip in XN metabolism. 
# -- or does the dip in XN metabolism related to the one person who likely stopped taking XN
# -- YES SIG CHANGE DUE TO p220 WHO STOPPED TAKING XN



adiv4plot <- adiv %>% 
        pivot_longer(cols = c(Observed, Shannon, Simpson), names_to = 'measure', values_to = 'values') %>% 
        group_by(response, time, measure) %>% 
        summarize(meanADIV = mean(values), 
                  serADIV = ser(values))

response_vals <- c('orange', 'cyan4')
names(response_vals) <- c('responder', 'nonresponder')


# SAVE 1200w x 400h
adivplot <- ggplot(adiv4plot, aes(time, meanADIV, fill = response, color = response)) +
        geom_path(aes(group = response), size = 1) + 
        geom_errorbar(aes(ymin = meanADIV - serADIV, ymax = meanADIV + serADIV), size =1, width = 0.2) +
        geom_point(shape = 21, fill = 'white', size = 3, show.legend = F) +
        facet_wrap(~measure, scales = 'free', strip.position = 'top') +
        scale_fill_manual(values = response_vals) +
        scale_color_manual(values = response_vals, name = 'Group:', labels = c('nonresponders', 'responders')) +
        scale_x_discrete(labels = c('v2' = 0, 'v3' = 2, 'v4' = 4, 'v5' = 6, 'v6' = 8)) +
        theme_minimal() +
        labs(x = 'Week', y = "Relative Abundance") +
        theme(panel.grid.minor = element_blank(), 
              panel.grid.major.x = element_blank(),
              strip.placement = "outside",
              panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
              legend.position = "none",
              axis.title = element_text(size = 16),
              axis.text = element_text(size = 16),
              legend.text = element_text(size = 16),
              legend.title = element_text(size = 16),
              strip.text = element_text(size = 16))


cowplot::plot_grid(leg, adivplot, cols = 1, align = 'v', rel_heights = c(1,2.5))







## Linear Mixed Effect Model of CLR-transformed Data




## FIGURE --
## CAPTION -- 

metadata <- psrelab %>% 
        subset_samples(treatment == 'B') %>% 
        sample_data() %>% 
        data.frame() %>% 
        rownames_to_column('sample.id') %>% 
        mutate(response = ifelse(part.id %in% c('p207', 'p209', 'p212', 'p220'), 'responder', 'nonresponder'))

relabs <- psrelab %>% 
        subset_samples(treatment == 'B') %>% 
        otu_table() %>% 
        data.frame() %>% 
        rownames_to_column('sample.id') %>% 
        right_join(metadata)

relab.long <- relabs %>% 
        dplyr::select(sample.id, part.id, treatment, time, response, ASV77, ASV57) %>% 
        pivot_longer(cols = starts_with('ASV'), names_to = 'ASV', values_to = 'relab') %>% 
        left_join(tax) %>% 
        mutate(taxa.name = ifelse(is.na(Species), 
                                  paste0('g_', Genus, '_', ASV), 
                                  paste0('g_', Genus, '_', Species, '_', ASV))) 

response_vals <- c('orange', 'cyan4')
names(response_vals) <- c('responder', 'nonresponder')


relab_ann_asv77 <- data.frame(time = c('v2', 'v5', 'v6'), 
                              relab = rep(0.016, 3), 
                              taxa.name = factor(rep("g_Coprococcus_comes_ASV77", 3)),
                              response = rep('responder', 3))
relab_ann_asv57 <- data.frame(time = 'v6',
                              relab = 0.03,
                              taxa.name = factor("g_[Ruminococcus] torques group_ASV57"),
                              response = 'responder')

# SAVE 900w x 450h
relab.plot <- ggplot(relab.long, aes(time, relab, color = response, fill = response)) +
        geom_path(aes(group = part.id), size = 1) +
        #geom_boxplot(alpha = 0.5, , linewidth = 1) +
        geom_point(shape = 21, fill = 'white', size = 3, show.legend = F) +
        facet_wrap(~ taxa.name, scales = 'free', ncol = 1, strip.position = 'top') +
        scale_fill_manual(values = response_vals) +
        scale_color_manual(values = response_vals, name = 'Group:', labels = c('nonresponders', 'responders')) +
        scale_x_discrete(labels = c('v2' = 0, 'v3' = 2, 'v4' = 4, 'v5' = 6, 'v6' = 8)) +
        theme_minimal() +
        labs(x = 'Week', y = "Relative Abundance") +
        theme(panel.grid.minor = element_blank(), 
              panel.grid.major.x = element_blank(),
              strip.placement = "outside",
              panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
              legend.position = "none",
              axis.title = element_text(size = 16),
              axis.text = element_text(size = 16),
              legend.text = element_text(size = 16),
              legend.title = element_text(size = 16),
              strip.text = element_text(size = 16)) + 
        geom_point(data = relab_ann_asv57, aes(time, relab), shape = '#', color = 'grey25', size = 4) +
        geom_point(data = relab_ann_asv77, aes(time, relab), shape = '*', color = 'grey25', size = 8)



## Presence/Absence of Microbiota between Responder // Nonresponder Groups ------

ct.treat.melt <- ct_treat %>% 
        subset_samples(part.id != 'p218') %>% 
        rarefy_even_depth(rngseed = 11) %>% 
        psmelt() %>% 
        mutate(response = ifelse(part.id %in% c('p207', 'p209', 'p212', 'p220'), 'responder', 'nonresponder')) %>% 
        mutate(sample.id = Sample) %>% 
        mutate(presence = ifelse(Abundance > 0, 1, 0))

NBmod.presence <- ct.treat.melt %>% 
        group_by(OTU) %>% 
        nest() %>% 
        mutate(model = purrr::map(data, function(x) try(glmmTMB(presence ~ response + (1|part.id), family=nbinom2, data = x)))) 


NBmod.presence.clean <- NBmod.presence[sapply(NBmod.presence$model, function(x) class(x[[1]]) !=  'character'),] 

NBmod.presence.contrast <- NBmod.presence.clean %>% 
        filter(!c(OTU %in% c('ASV215', 'ASV32'))) %>% 
        mutate(pres_resp_pval = purrr::map_dbl(model, function(x) summary(x)$coefficients$cond[[2,4]])) %>%
        ungroup() %>% 
        mutate(padj = p.adjust(pres_resp_pval, method = 'BH')) %>% 
        rename(ASV = OTU) %>% 
        left_join(tax, by = 'ASV')
# NO SIGNIFICANT ASVs when considering Presence and Absence


        
## Negative Binomial Differential Abundance Analysis ------------

library(glmmTMB)
library(emmeans)
library(multcomp)

ct.treat.melt <- ct_treat %>% 
        subset_samples(part.id != 'p218') %>% 
        rarefy_even_depth(rngseed = 11) %>% 
        psmelt() %>% 
        mutate(response = ifelse(part.id %in% c('p207', 'p209', 'p212', 'p220'), 'responder', 'nonresponder')) %>% 
        mutate(sample.id = Sample) #%>% 
        #mutate(Abundance = Abundance + 1)

NBmod.res <- ct.treat.melt %>% 
        group_by(OTU) %>% 
        nest() %>% 
        mutate(model = purrr::map(data, function(x) try(glmmTMB(Abundance ~ response*time + (1|part.id), family=nbinom2, data = x)))) 


NBmod.clean <- NBmod.res[sapply(NBmod.res$model, function(x) class(x[[1]]) !=  'character'),] 

NBmod.contrast <- NBmod.clean %>% 
        filter(!c(OTU %in% c('ASV105', 'ASV113', 'ASV147', 'ASV150', 'ASV178', 'ASV191', 'ASV200', 'ASV66', 'ASV92', 'ASV7'))) %>% 
        mutate(T1_BvA = map_dbl(model, function(x) broom::tidy(multcomp::glht(x, linfct = matrix(c(0,1,0,0,0,0,0,0,0,0), nrow = 1)))[[1,6]][1])) %>% 
        mutate(T2_BvA = map_dbl(model, function(x) broom::tidy(multcomp::glht(x, linfct = matrix(c(0,1,0,0,0,0,1,0,0,0), nrow = 1)))[[1,6]][1])) %>% 
        mutate(T3_BvA = map_dbl(model, function(x) broom::tidy(multcomp::glht(x, linfct = matrix(c(0,1,0,0,0,0,0,1,0,0), nrow = 1)))[[1,6]][1])) %>%
        mutate(T4_BvA = map_dbl(model, function(x) broom::tidy(multcomp::glht(x, linfct = matrix(c(0,1,0,0,0,0,0,0,1,0), nrow = 1)))[[1,6]][1])) %>%
        mutate(T5_BvA = map_dbl(model, function(x) broom::tidy(multcomp::glht(x, linfct = matrix(c(0,1,0,0,0,0,0,0,0,1), nrow = 1)))[[1,6]][1])) %>%
        mutate(B_T1vT2 = map_dbl(model, function(x) broom::tidy(multcomp::glht(x, linfct = matrix(c(0,0,1,0,0,0,1,0,0,0), nrow = 1)))[[1,6]][1])) %>%
        mutate(B_T1vT3 = map_dbl(model, function(x) broom::tidy(multcomp::glht(x, linfct = matrix(c(0,0,0,1,0,0,0,1,0,0), nrow = 1)))[[1,6]][1])) %>%
        mutate(B_T1vT4 = map_dbl(model, function(x) broom::tidy(multcomp::glht(x, linfct = matrix(c(0,0,0,0,1,0,0,0,1,0), nrow = 1)))[[1,6]][1])) %>%
        mutate(B_T1vT5 = map_dbl(model, function(x) broom::tidy(multcomp::glht(x, linfct = matrix(c(0,0,0,0,0,1,0,0,0,1), nrow = 1)))[[1,6]][1])) %>%
        ungroup() %>% 
        mutate(T1_BvA_adj = p.adjust(T1_BvA, method = 'BH')) %>% 
        mutate(T2_BvA_adj = p.adjust(T2_BvA, method = 'BH')) %>%
        mutate(T3_BvA_adj = p.adjust(T3_BvA, method = 'BH')) %>%
        mutate(T4_BvA_adj = p.adjust(T4_BvA, method = 'BH')) %>%
        mutate(T5_BvA_adj = p.adjust(T5_BvA, method = 'BH')) %>%
        mutate(B_T1vT2_adj = p.adjust(B_T1vT2, method = 'BH')) %>%
        mutate(B_T1vT3_adj = p.adjust(B_T1vT3, method = 'BH')) %>%
        mutate(B_T1vT4_adj = p.adjust(B_T1vT4, method = 'BH')) %>%
        mutate(B_T1vT5_adj = p.adjust(B_T1vT5, method = 'BH')) %>% 
        rename(ASV = OTU) %>% 
        left_join(tax, by = 'ASV')

sig.asvs <- NBmod.contrast %>% 
        dplyr::select(-c(T1_BvA:T5_BvA_adj), -data) %>% 
        filter(if_any(B_T1vT2_adj:B_T1vT5_adj, ~ .x <= 0.0505)) %>% 
        filter(!c(ASV %in% c('ASV32', 'ASV40', 'ASV112'))) %>% 
        mutate(taxa.name = ifelse(is.na(Species), paste0('g_', Genus, '_', ASV), paste0('g_', Genus, "_", Species, '_', ASV)))

# T2: ASV32, ASV40 - T3: ASV32, ASV40, ASV108 ----- accept 108
# T4: ASV112, ASV139, ASV32, ASV68 ----- accept 139, 68
# T5: ASV32, ASV161, ASV112, ASV57, ASV49, ASV85, ASV139 ---- accept 161, 57, 49, 85, 139


asv77 <- ct.treat.melt %>% filter(OTU == 'ASV77')

model77 <- glmmTMB(Abundance ~ response*time + (1|part.id), family = nbinom2, data = asv77)

model77.2 <- glmer.nb(Abundance ~ response*time + (1|part.id), data = asv77)



## FIGURE ## --------------------
## CAPTION ## -------------------


relab.long <- relabs %>% 
        dplyr::select(sample.id, part.id, treatment, time, response, starts_with('ASV')) %>% 
        pivot_longer(starts_with('ASV'), names_to = 'ASV', values_to = 'relab') %>% 
        filter(ASV %in% sig.asvs$ASV) %>% 
        left_join(tax) %>% 
        mutate(taxa.name = ifelse(is.na(Species), 
                                  paste0('g_', Genus, '_', ASV), 
                                  paste0('g_', Genus, '_', Species, '_', ASV))) %>% 
        filter(response == 'responder')

response_vals <- c('orange', 'cyan4')
names(response_vals) <- c('responder', 'nonresponder')


# SAVE 900w x 450h
ggplot(relab.long, aes(time, relab, group = part.id, color = response)) +
        geom_path(size = 1) +
        geom_boxplot() +
        geom_point(shape = 21, fill = 'white', size = 3, show.legend = F) +
        scale_color_manual(values = response_vals, name = 'Group:', labels = c('nonresponders', 'responders')) +
        scale_x_discrete(labels = c('v2' = 0, 'v3' = 2, 'v4' = 4, 'v5' = 6, 'v6' = 8)) +
        theme_minimal() +
        facet_wrap(~ taxa.name, strip.position = "top", scales = 'free', ncol = 4) +
        labs(x = 'Week', y = "Relative Abundance") +
        theme(panel.grid.minor = element_blank(), 
              panel.grid.major.x = element_blank(),
              strip.placement = "outside",
              panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
              legend.position = "top",
              axis.title = element_text(size = 16),
              axis.text = element_text(size = 16),
              legend.text = element_text(size = 16),
              legend.title = element_text(size = 16),
              strip.text = element_text(size = 12))



# With Boxplots overlaid [800w x 1400h] or [1600w x 800h]
ggplot(relab.long, aes(time, relab, color = response, fill = response)) +
        geom_path(aes(group = part.id), size = 1) +
        geom_boxplot(aes(alpha = 0.5), , linewidth = 1) +
        geom_point(shape = 21, fill = 'white', size = 3, show.legend = F) +
        facet_wrap(~ taxa.name, scales = 'free', ncol = 4, strip.position = 'top') +
        theme_minimal() +
        scale_fill_manual(values = response_vals) +
        scale_color_manual(values = response_vals, name = 'Group:', labels = c('nonresponders', 'responders')) +
        scale_x_discrete(labels = c('v2' = 0, 'v3' = 2, 'v4' = 4, 'v5' = 6, 'v6' = 8)) +
        labs(x = 'Week', y = "Relative Abundance") +
        theme(panel.grid.minor = element_blank(), 
              panel.grid.major.x = element_blank(),
              strip.placement = "outside",
              panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
              legend.position = "none",
              axis.title = element_text(size = 16),
              axis.text = element_text(size = 12),
              legend.text = element_text(size = 16),
              legend.title = element_text(size = 16),
              strip.text = element_text(size = 12))



# -----------------------------------------------------------------------------#





counts <- ct_treat %>% 
        rarefy_even_depth(rngseed = 11) %>% 
        otu_table() %>% 
        as.data.frame() %>% 
        rownames_to_column('sample.id') %>% 
        left_join(mdata) 

# All Sig ASVs -- boxplots
counts.long <- counts %>% 
        pivot_longer(starts_with('ASV'), names_to = 'ASV', values_to = 'counts') %>% 
        filter(ASV %in% sig.asvs$ASV) %>% 
        left_join(tax) %>% 
        mutate(taxa.name = ifelse(is.na(Species), paste0('g_', Genus, '_', ASV), paste0('g_', Genus, '_', Species, '_', ASV)))

ggplot(counts.long, aes(time, counts, fill = response)) + 
        geom_boxplot() +
        facet_wrap(~taxa.name, scales = 'free')

# BTW Sig ASVs -- boxplots
counts.btw.long <- counts %>% 
        pivot_longer(starts_with('ASV'), names_to = 'ASV', values_to = 'counts') %>% 
        filter(ASV %in% sig.asvs.btw$ASV) %>% 
        left_join(tax) %>% 
        mutate(taxa.name = ifelse(is.na(Species), paste0('g_', Genus, '_', ASV), paste0('g_', Genus, '_', Species, '_', ASV)))

ggplot(counts.btw.long, aes(time, counts, fill = response)) + 
        geom_boxplot() +
        facet_wrap(~taxa.name, scales = 'free')

# BTW -- Longitudinal 
ggplot(counts.btw.long, aes(time, counts, color = response, group = part.id)) +
        geom_point() +
        geom_path() +
        facet_wrap(~taxa.name, scales = 'free')



# WITHIN Sig ASVs -- boxplots
counts.within.long <- counts %>% 
        pivot_longer(starts_with('ASV'), names_to = 'ASV', values_to = 'counts') %>% 
        filter(ASV %in% sig.asvs.within$ASV) %>% 
        filter(response == 'responder') %>% 
        left_join(tax) %>% 
        mutate(taxa.name = ifelse(is.na(Species), 
                                  paste0('g_', Genus, '_', ASV), 
                                  paste0('g_', Genus, '_', Species, '_', ASV))) 


############## can i filter by ASVs with > 50% zeros as baseline ###############
counts.within.long %>% 
        filter(time == 'v2') %>% 
        group_by(ASV, part.id) %>% 
        summarize(sum(counts))
################################################################################


ggplot(counts.within.long, aes(time, counts)) + 
        geom_boxplot() +
        facet_wrap(~taxa.name, scales = 'free')

# WITHIN -- Longitudinal 
ggplot(counts.within.long, aes(time, counts, group = part.id, color = part.id)) +
        geom_point() +
        geom_path() +
        facet_wrap(~taxa.name, scales = 'free')



# Mean Longitudinal Plots
counts.summ <- counts.long %>% 
        group_by(taxa.name, time, response) %>% 
        summarize(meanCounts = mean(counts),
                  sdCounts = sd(counts))

ggplot(counts.summ, aes(time, meanCounts, group = response, color = response)) +
        geom_point() +
        geom_path() + 
        geom_errorbar(aes(ymin = meanCounts - sdCounts, ymax = meanCounts + sdCounts), width = 0.2) +
        facet_wrap(~taxa.name, scales = 'free')




# INFLAMMATORY BIOMARKERS --------------------------


el.transformed <- metab.data.scaled %>% 
        dplyr::select(starts_with('EL'), sample.id, part.id, treatment, time) %>% 
        mutate(response = ifelse(part.id %in% c('p207', 'p209', 'p212', 'p220'), 
                                 'responder', 'nonresponder')) %>% 
        filter(!c(part.id == 'p220' & time == 'v5')) %>% 
        filter(treatment == 'B')

el_tidy <- el.transformed %>% 
        pivot_longer(cols = starts_with('EL'), names_to = 'EL', values_to = 'value')

el_lm <- el_tidy %>% 
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
# between the responders and nonresponders -- IL6, IL10, and TNF are significantly different


## FIGURE -- Inflammatory Markers Significantly different between Responders and Nonresponders
## CAPTION -- 

el.long <- metab.data %>% 
        mutate(response = ifelse(part.id %in% c('p207', 'p209', 'p212', 'p220'), 
                                 'responder', 'nonresponder')) %>% 
        filter(treatment == 'B') %>% 
        pivot_longer(cols = c('EL_IL10', 'EL_IL6', 'EL_TNF'), names_to = 'marker', values_to = 'conc')

el.summ <- el.long %>% 
        group_by(marker, time, response) %>% 
        summarize(mMarker = mean(conc),
                  serMarker = ser(conc))


response_vals <- c('orange', 'cyan4')
names(response_vals) <- c('responder', 'nonresponder')


# Significance annotation
el_ann_text <- data.frame(time = c('v2', 'v3', 'v4', 'v5', 'v6'), 
                          mMarker = c(rep(2700, 5), rep(460, 5), rep(420, 5)),
                          marker = c("EL_IL10", "EL_IL10", "EL_IL10", "EL_IL10", "EL_IL10",
                                     "EL_IL6", "EL_IL6", "EL_IL6", "EL_IL6", "EL_IL6", 
                                     "EL_TNF", "EL_TNF", "EL_TNF", "EL_TNF", "EL_TNF"),
                          response = rep('responder', 15)) %>% 
        mutate(across(marker, ~ factor(.x)))


inf.plot <- ggplot(el.summ, aes(time, mMarker, group = response, color = response)) +
        geom_errorbar(aes(ymin = mMarker - serMarker, ymax = mMarker + serMarker), size = 1, width = 0.2) +
        geom_path(size = 1) +
        geom_point(shape = 21, fill = 'white', size = 3, show.legend = F) +
        facet_wrap(~marker, strip.position = "top", scales = 'free',
                   labeller = as_labeller(c('EL_IL10' = 'IL10',
                                            'EL_IL6' = 'IL6',
                                            'EL_TNF' = 'TNF')),
                   ncol = 1) +
        scale_color_manual(values = response_vals, name = 'Group:', labels = c('nonresponders', 'responders')) +
        scale_x_discrete(labels = c('v2' = 0, 'v3' = 2, 'v4' = 4, 'v5' = 6, 'v6' = 8)) +
        theme_minimal() +
        labs(x = 'Week', y = "pg/mL") +
        theme(panel.grid.minor = element_blank(), 
              panel.grid.major.x = element_blank(),
              strip.placement = "outside",
              panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
              legend.position = "none",
              axis.title = element_text(size = 16),
              axis.text = element_text(size = 16),
              legend.text = element_text(size = 16),
              legend.title = element_text(size = 16),
              strip.text = element_text(size = 16)) +
        geom_point(data = el_ann_text, aes(time, mMarker), color = 'grey25', shape = "*", size = 8)







# BILE ACIDS ----------------------------------------------------

ba.transformed <- metab.data.scaled %>% 
        dplyr::select(sample.id, part.id, treatment, time, starts_with('BA')) %>% 
        mutate(response = ifelse(part.id %in% c('p207', 'p209', 'p212', 'p220'), 
                                 'responder', 'nonresponder')) %>% 
        filter(!c(part.id == 'p220' & time == 'v5')) %>% 
        filter(treatment == 'B')


ba.transformed_tidy <- ba.transformed %>% 
        pivot_longer(cols = starts_with('BA'), names_to = 'BA', values_to = 'value')


ba_lm <- ba.transformed_tidy %>% 
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
        #mutate(across(T5_BvB, ~ p.adjust(.x, method = 'BH'))) #%>% 
        mutate(across(c(everything(), -BA, -data, -model), ~ p.adjust(.x, method = 'BH')))

# borderline significant differentce in DCA and secondary BA within responders groups
# -- no sig changes in nonresponders compared to baseline

# between groups -- DHCA (v4) and GUDCA (v6)



## FIGURE -- Bile Acids Significantly different between Responders and Nonresponders
## CAPTION -- 

ba.long <- metab.data %>% 
        mutate(response = ifelse(part.id %in% c('p207', 'p209', 'p212', 'p220'), 
                                 'responder', 'nonresponder')) %>% 
        filter(treatment == 'B') %>% 
        pivot_longer(cols = c('BA_DCA', 'BA_DHCA', 'BA_GUDCA', 'BA_secondary', 'BA_LCA'), names_to = 'BA', values_to = 'conc')

ba.summ <- ba.long %>% 
        group_by(BA, time, response) %>% 
        summarize(mBA = mean(conc),
                  serBA = ser(conc))


response_vals <- c('orange', 'cyan4')
names(response_vals) <- c('responder', 'nonresponder')

ggplot(ba.summ, aes(time, mBA, group = response, color = response)) +
        geom_errorbar(aes(ymin = mBA - serBA, ymax = mBA + serBA), size = 1, width = 0.2) +
        geom_path(size = 1) +
        geom_point(shape = 21, fill = 'white', size = 3) +
        facet_wrap(~BA, strip.position = "top", scales = 'free',
                   labeller = as_labeller(c('BA_DCA' = 'DCA',
                                            'BA_DHCA' = 'DHCA',
                                            'BA_GUDCA' = 'GUDCA', 
                                            'BA_secondary' = 'Secondary Bile Acids',
                                            'BA_LCA' = 'LCA')),
                   ncol = 2) +
        scale_color_manual(values = response_vals, name = 'Group:', labels = c('nonresponders', 'responders')) +
        scale_x_discrete(labels = c('v2' = 0, 'v3' = 2, 'v4' = 4, 'v5' = 6, 'v6' = 8)) +
        theme_minimal() +
        labs(x = 'Week', y = "µg/g") +
        theme(panel.grid.minor = element_blank(), 
              panel.grid.major.x = element_blank(),
              strip.placement = "outside",
              panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
              legend.position = "top",
              axis.title = element_text(size = 16),
              axis.text = element_text(size = 16),
              legend.text = element_text(size = 16),
              legend.title = element_text(size = 16),
              strip.text = element_text(size = 16))


## FIGURE --
## CAPTION --

ba.long <- metab.data %>% 
        mutate(response = ifelse(part.id %in% c('p207', 'p209', 'p212', 'p220'), 
                                 'responder', 'nonresponder')) %>% 
        filter(treatment == 'B') 

ba.summ <- ba.long %>% 
        group_by(time, response) %>% 
        summarize(mBA = mean(BA_secondary),
                  serBA = ser(BA_secondary))


response_vals <- c('orange', 'cyan4')
names(response_vals) <- c('responder', 'nonresponder')

ba_ann_text <- data.frame(time = 'v6', response = 'responder', mBA = 1800)

ba.plot <- ggplot(ba.summ, aes(time, mBA, group = response, color = response)) +
        geom_errorbar(aes(ymin = mBA - serBA, ymax = mBA + serBA), size = 1, width = 0.2) +
        geom_path(size = 1) +
        geom_point(shape = 21, fill = 'white', size = 3) +
        scale_color_manual(values = response_vals, name = 'Group:', labels = c('nonresponders', 'responders')) +
        scale_x_discrete(labels = c('v2' = 0, 'v3' = 2, 'v4' = 4, 'v5' = 6, 'v6' = 8)) +
        theme_minimal() +
        labs(x = 'Week', y = "µg/g", title = 'Secondary Bile Acids') +
        theme(panel.grid.minor = element_blank(), 
              panel.grid.major.x = element_blank(),
              strip.placement = "outside",
              panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
              legend.position = "none",
              axis.title = element_text(size = 16),
              axis.text = element_text(size = 16),
              legend.text = element_text(size = 16),
              legend.title = element_text(size = 16),
              strip.text = element_text(size = 16), 
              plot.title = element_text(size = 18, hjust = 0.5)) +
        geom_point(data = ba_ann_text, aes(time, mBA), shape = '#', size = 5, color = 'grey25')




## Extract LEGEND --------------------------------------------------------


legend1 <- ggplot(ba.summ, aes(time, mBA, group = response, color = response)) +
        #geom_errorbar(aes(ymin = mBA - serBA, ymax = mBA + serBA), size = 1, width = 0.2) +
        geom_path(size = 10) +
        geom_point(shape = 21, fill = 'white', size = 3, show.legend = F) +
        scale_color_manual(values = response_vals, name = 'Group:', labels = c('nonresponders', 'responders')) +
        scale_x_discrete(labels = c('v2' = 0, 'v3' = 2, 'v4' = 4, 'v5' = 6, 'v6' = 8)) +
        theme_minimal() +
        labs(x = 'Week', y = "µg/g") +
        theme(panel.grid.minor = element_blank(), 
              panel.grid.major.x = element_blank(),
              #strip.placement = "outside",
              panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
              legend.position = "bottom",
              axis.title = element_text(size = 16),
              axis.text = element_text(size = 16),
              legend.text = element_text(size = 14),
              legend.title = element_text(size = 14),
              strip.text = element_text(size = 16), 
              plot.title = element_text(size = 18, hjust = 0.5))

leg <- cowplot::get_plot_component(legend1, 'guide-box-bottom', return_all = T)
leg <- cowplot::ggdraw(leg)



## Combine plots into one grid -------------------------------------------

layer1 <- cowplot::plot_grid(ba.plot, relab.plot, ncol = 1, 
                             rel_heights = c(1,2), labels = c("B", "C"),
                             label_size = 16)

layer2 <- cowplot::plot_grid(inf.plot, layer1, ncol = 2, labels = c('A', ''),
                             label_size = 16)
# SAVE 1200 x 1200 
final.plot <- cowplot::plot_grid(leg, layer2, ncol = 1, rel_heights = c(1, 8))






# REGRESSION ANALYSIS -- METABOLITES AND CDAI SCORE -----------------------


ct_treat <- ps_counts %>% 
        subset_samples(treatment == 'B') 

clr <- ct_treat %>% 
        rarefy_even_depth(rngseed = 11) %>% 
        otu_table() %>% 
        as.data.frame() %>% 
        decostand(method = 'clr', pseudocount = 1) 

euclid_dist <- clr %>% 
        stats::prcomp()

summary(euclid_dist)$importance[2,] 

mdata <- sample_data(ct_treat) %>% 
        data.frame() %>% 
        rownames_to_column("sample.id") %>% 
        mutate(response = ifelse(part.id %in% c('p207', 'p209', 'p212', 'p220'), 'responder', 'nonresponder'))

micro.comps <- euclid_dist$x %>% 
        as.data.frame() %>% 
        rownames_to_column('sample.id') %>% 
        left_join(., mdata) %>% 
        left_join(clr %>% rownames_to_column('sample.id'))



full.data.treat <- metab.data.scaled %>% 
        left_join(micro.comps) %>% 
        filter(treatment == 'B') %>% 
        filter(part.id != 'p218') %>% 
        filter(!c(part.id == 'p220' & time == 'v5')) 


fullMod <- lmerTest::lmer(cdai_score ~ BA_secondary + ASV77 + ASV57 + EL_TNF + EL_IL10 + EL_IL6 + (1|part.id),
                         data = full.data.treat, REML = F)

stepMod <- lmerTest::step(fullMod, reduce.random=F)

finalMod <- lmerTest::lmer(cdai_score ~ BA_secondary + EL_IL10 + ASV77 + (1|part.id),
                           data = full.data.treat, REML = F)


summary(finalMod)
MuMIn::r.squaredGLMM(finalMod)


# Model found:
# cdai_score ~ BA_secondary + EL_IL10 + (1 | part.id)




#nullMod <- lme4::lmer(cdai_score ~ (1|part.id), data = full.data.treat, REML = F)
# Stepwise Backward Regression
#stepFit <- MASS::stepAIC(fullMod, direction = 'backward', trace = T, scope = list(upper = fullMod, lower = nullMod))






# DELTA CDAI as a function of BASELINE BIOMARKERS ------------------------------

delta.cdai <- full.data.treat %>% 
        filter(time %in% c('v2', 'v6')) %>% 
        dplyr::select(part.id, time, cdai_score, BA_secondary, 
                      EL_IL10, EL_IL6, EL_TNF, ASV77, ASV57) %>% 
        pivot_wider(names_from = 'time', values_from = c('cdai_score', 'ASV77', 'ASV57', 
                                                         'BA_secondary', 'EL_IL10', 'EL_IL6',
                                                         'EL_TNF')) %>% 
        mutate(delt.cdai = cdai_score_v6 - cdai_score_v2)



fullMod <- lm(delt.cdai ~ EL_IL10_v2 + EL_IL6_v2 + EL_TNF_v2 + 
                      ASV77_v2 + ASV57_v2 + BA_secondary_v2, data = delta.cdai)


stepMod <- stats::step(fullMod, direction = 'backward')


finalMod <- lm(delt.cdai ~ EL_TNF_v2 + ASV77_v2 + ASV57_v2, data = delta.cdai)
# adjusted R-squared = 0.33

# BASELINE levels of plasma TNF and abundance of ASV57 and ASV77 best predictor change in CDAI score




# revision analysis - categorize individuals by chronic inflammation - plot cdai_score

inflam.summ <- metab.data %>% 
        pivot_longer(starts_with('EL_'), names_to = 'marker', values_to = 'conc') %>% 
        group_by(part.id, marker, treatment) %>% 
        summarize(mInflam = mean(conc)) %>% 
        pivot_wider(names_from = 'marker', values_from = 'mInflam') %>% 
        dplyr::select(part.id, treatment, EL_IL6, EL_IL10, EL_TNF) %>% 
        # "High" inflammation indicates severe vs low-grade or no inflammaiton
        mutate(inflam.status = ifelse(c(EL_IL6 <= 30 & EL_TNF <= 30 & EL_IL10 <= 30), "A_low", "B_high"))

# Normal Range:
# IL6: 0-5 pg/mL
# TNF: 0-8.1 pg/mL
# IL10: 0-10 pg/mL


# maybe borderline: 
# B: 205, 210
# A: 203


plot <- metab.data %>% 
        left_join(inflam.summ %>% dplyr::select(part.id, inflam.status)) %>% 
        filter(part.id != 'p218') %>% 
        filter(!c(part.id == 'p220' & time == 'v5'))



ggplot(plot, aes(time, cdai_score, group = part.id, color = inflam.status, shape = treatment)) +
        geom_point() +
        geom_path()


plot2 <- metab.data %>% 
        left_join(inflam.summ %>% dplyr::select(part.id, inflam.status)) %>% 
        filter(part.id != 'p218') %>% 
        filter(!c(part.id == 'p220' & time == 'v5')) %>% 
        group_by(treatment, inflam.status, time) %>% 
        summarize(mCDAI = mean(cdai_score))

ggplot(plot2, aes(time, mCDAI, group = interaction(treatment, inflam.status), 
                  color = inflam.status,
                  shape = treatment)) +
        geom_point(size = 3) +
        geom_path()

# Does XN work in synch with immune system in severe inflammatory participants. 


inflam.summ.summ <- metab.data %>% 
        left_join(inflam.summ %>% dplyr::select(part.id, inflam.status)) %>% 
        filter(time == 'v2') %>% 
        pivot_longer(starts_with('EL_'), names_to = 'marker', values_to = 'conc') %>% 
        group_by(treatment, inflam.status, marker) %>% 
        summarize(mInflam = mean(conc),
                  sdInflam = sd(conc)) %>% 
        pivot_wider(names_from = 'marker', values_from = c('mInflam', 'sdInflam')) %>% 
        dplyr::select(treatment, ends_with(c('EL_IL6', 'EL_IL10', 'EL_TNF')))

stats <-  metab.data %>% 
        left_join(inflam.summ %>% dplyr::select(part.id, inflam.status)) %>% 
        filter(time == 'v2') %>% 
        pivot_longer(starts_with('EL_'), names_to = 'marker', values_to = 'conc') %>% 
        group_by(inflam.status, marker) %>% 
        nest() %>% 
        mutate(mod = purrr::map(data, function(x) t.test(data = x, conc ~ treatment))) %>% 
        mutate(pval = purrr::map_dbl(mod, function(x) x$p.value))
        

# Model CDAI score

lmm <- metab.data.scaled %>% 
        left_join(inflam.summ %>% dplyr::select(part.id, inflam.status)) %>%
        #filter(treatment == 'B') %>% 
        filter(inflam.status == 'B_high') %>% 
        #mutate(treatment = ifelse(treatment == 'A', 'PL', 'A')) %>% 
        #mutate(group = paste0(treatment, '_', inflam.status)) %>% 
        filter(part.id != 'p218') %>% 
        filter(!c(part.id == 'p220' & time == 'v5')) %>%
        pivot_longer(c('EL_IL6', 'EL_IL10', 'EL_TNF'), names_to = 'marker', 
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
        ungroup() #%>% 
        mutate(across(starts_with(c('v' | 'B_'), ~p.adjust(.x, method = 'BH'))))
        
        
        
        
        
        
        
        
        
        
        
        


