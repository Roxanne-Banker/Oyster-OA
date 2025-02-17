
# July 2023
# Document to analyze oyster larvae/juvenile data from 
# Oyster multi-omics project run in Aug-Oct 2019
# Following some of this tutorial: https://benjjneb.github.io/dada2/tutorial.html

# loading in required packages
library(phyloseq)
library(dada2)
library(decontam)
library(phangorn)
library(DECIPHER)
library(reshape)
library(coin)
library(FSA)
library(Biostrings)
library(DESeq2)
library(seqinr)
library(ALDEx2)
library(ape)
library(biomformat)

library(edgeR)
library(KEGGREST)
library(SummarizedExperiment)
library(metagenomeSeq)
library(Maaslin2)
library(limma)
library(lefser)

library(readr)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)

library(ggpicrust2)

library(ggplot2)
library(RColorBrewer)

library(viridis)
library(gtable)
library(grid)
library(gridExtra)

library(vegan)
library(rmarkdown)
library(dplyr)
library(tidyr)
library(lme4)
library(brms)
library(lmerTest)
library(bayesplot)


set.seed(5311)





# I. Import and prep data ----



load("~/Dropbox/Oyster_microbes/dna_data/oyster_oa_ms-v2.RData")
#save.image("~/Dropbox/Oyster_microbes/dna_data/oyster_oa_ms-v2.RData")
#save.image("~/Dropbox/Oyster_microbes/dna_data/oyster_oa_ms-v3.RData")
#save.image("~/Dropbox/Oyster_microbes/dna_data/oyster_oa_ms-v2-copy.RData")


## A. Create ASV table ----
# Using Dada2 to amplicon sequence variant (ASV) table

# Setting the path for the fastq files
path <- "~/Dropbox/Oyster_microbes/dna_data/raw_reads/" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

# Sort and get sample names
fnFs <- sort(list.files(path, pattern="_R1_001.fastq"))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq"))
sample.names <- sapply(strsplit(fnFs, "_S"), `[`, 1)

# specify full paths to the data
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)

# Quality Inspection
# Looking at our data to try to see at what point the read quality drops below Q20

#plotQualityProfile(fnFs[1:4])
# For the forward reads, the controls look a bit bad, but the other samples
# look good. How about the reverse reads?

#plotQualityProfile(fnRs[1:4])
# These are maybe worse, but tutorials indicate that is okay for reverse reads

# Note: might have to be careful in subsequent trimming steps
# The tutorial this script is following indicates that the tutorial 
# was written for 2x250 V4 sequence data, which means that the forward and
# reverse reads almost completely overlap, and can use quality scores to 
# guide filtering. I am not sure how well our primer set, V4-V5, overlaps. 
# I am going to try what the tutorial uses, then adjust
# the truncLen from there.

# tutorial "your truncLen must be large enough to maintain 20 + 
# biological.length.variation nucleotides of overlap between them."

## B. Quality Filtering ----

# Assign the filenames for the filtered fastq.gz files.
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(300,240),
                     maxN=0, maxEE=c(3,6), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE,
                     trimLeft = c(19, 20)) # On Windows set multithread=FALSE
head(out)

# this is the first attempt with the default truncLen values in the tutorial
# truncLen=c(240,160)
# reads.in reads.out
# b-A_S179_L001_R1_001.fastq       358       109
# b-B_S191_L001_R1_001.fastq       329        94
# C1-A_S332_L001_R1_001.fastq    61807     56665
# C1-B_S344_L001_R1_001.fastq    82331     75507
# C1-C_S108_L001_R1_001.fastq   103194     95612
# C2-A_S120_L001_R1_001.fastq   111945    103285

# truncLen=c(200,250)
# reads.in reads.out
# b-A_S179_L001_R1_001.fastq       358        72
# b-B_S191_L001_R1_001.fastq       329        66
# C1-A_S332_L001_R1_001.fastq    61807     50486
# C1-B_S344_L001_R1_001.fastq    82331     67597
# C1-C_S108_L001_R1_001.fastq   103194     88307
# C2-A_S120_L001_R1_001.fastq   111945     95590

# truncLen=c(300,300)
# reads.in reads.out
# b-A_S179_L001_R1_001.fastq       358        41
# b-B_S191_L001_R1_001.fastq       329        42
# C1-A_S332_L001_R1_001.fastq    61807     37846
# C1-B_S344_L001_R1_001.fastq    82331     52180
# C1-C_S108_L001_R1_001.fastq   103194     70509
# C2-A_S120_L001_R1_001.fastq   111945     74158

# Settled on these parameters for now:
# truncLen=c(300,240)
# maxEE=c(3,6)
# reads.in reads.out
# b-A_S179_L001_R1_001.fastq       358       129
# b-B_S191_L001_R1_001.fastq       329       111
# C1-A_S332_L001_R1_001.fastq    61807     57342
# C1-B_S344_L001_R1_001.fastq    82331     76414
# C1-C_S108_L001_R1_001.fastq   103194     96472
# C2-A_S120_L001_R1_001.fastq   111945    104226



## C. Learn the Error Rates ----
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# It is always worthwhile, as a sanity check if nothing else, to visualize the estimated error rates:
plotErrors(errF, nominalQ=TRUE)


# Sample Inference

# We are now ready to apply the core sample inference algorithm 
# to the filtered and trimmed sequence data.
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

# Inspecting the returned dada-class object:
dadaFs[[5]]

# Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)


# Inspect the merger data.frame from the first sample
head(mergers[[4]])

## D. Construct sequence table ----


# We can now construct an amplicon sequence variant table (ASV) table, 
# a higher-resolution version of the OTU table produced 
# by traditional methods.
seqtab <- makeSequenceTable(mergers)

dim(seqtab)
# [1]   38 7916


# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# 287  312  326  341  351  352  356  359  362  364  365  366  367  368  369  370 
# 1    1    1    3    1    1    1    1    4    2    6   12    7  107   63  623 
# 371  372  373  374  375  376  377  378  379  380  381  382  383  389  390  392 
# 145 2445 1639 2592  142   75    8    2    1    4    2    1    4    1    1    1 
# 394  396  398  399  400  401  402  403  406  408  440  449  459 
# 2    1    1    3    1    2    1    1    2    2    1    1    1 

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

dim(seqtab.nochim)
# [1]   38 2043

sum(seqtab.nochim)/sum(seqtab)
# [1] 0.8725215


# Track reads through the pipeline. As a final check of our progress, we’ll look at the number of 
# reads that made it through each step in the pipeline:
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
print(track)

## E. Assign taxonomy ----

taxa <- assignTaxonomy(seqtab.nochim, "/home/rbanker/Dropbox/Oyster_microbes/dna_data/raw_reads/silva_nr_v132_train_set.fa.gz", multithread=TRUE)

# Assign species identity
taxa <- addSpecies(taxa, "/home/rbanker/Dropbox/Oyster_microbes/dna_data/raw_reads/silva_species_assignment_v132.fa.gz")

# Let’s inspect the taxonomic assignments:
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


## F. Build Phylogeny ----

# All from Ettinger Harbor Notebook
# March 26 2020: This is computationally expensive and I am thus going to skip right now.
# So far, it seems like the need for the phylogeny comes from calculating Unifrac distances
# I plan to use Bray-Curtis, so perhaps not needed right now

#get DNA sequences from chimera-removed ASV table
seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs # This propagates to the tip labels of the tree

#Align DNA sequences
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)

#Turn alignment into a matrix
phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
#calculate maximum liklihood values for alignment
dm <- dist.ml(phangAlign)


#build neighbor joining tree
#a GTR+G+I (Generalized time-reversible with Gamma rate variation)
#maximum likelihood tree using the neighbor-joining tree
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE,
                    optGamma=TRUE, rearrangement = "stochastic", control = pml.control(trace = 0))
# need to rerun, took too long ^

#saveRDS(fitGTR, "Harbor_seqtab2_tree.rds")
#Now to root the tree
#fitGTR <- readRDS("Harbor_seqtab2_tree.rds")
#using most abundant archaea as outgroup

# find most abundant archea to use as outgroup
# taxa.df <- data.frame(taxa)
# archaea <- filter(taxa.df, Kingdom == "Archaea")
# archaea.names <- row.names(archaea)
# 
# archaea.count <- data.frame(seqtab.nochim) %>% select(archaea.names)
# colSums(archaea.count)

outgroupseq <- "TACCGGCAGCTCAAGTGGTCGTCGCTTTTATTGGGCCTAAAACGTCCGTAGCCTGTTTGGTAAATCTGTGGGTAAATCAACCAGCTTAACTGGTTGAATTCTGCAGAGACTGCCAGACTAGGGACCGGGAGAGGTGTGGGGTACTCGAGGGGTAGGGGTAAAATCCTGTCATCCTTCGAGGACCACCAGTTGCGAAGGCGCCACACTGGAACGGATCCGACGGTCAGGGACGAAGCCTAGGGGCACGAACCGGATTAGATACCCGGGTAGTCCTAGGTGTAAACGCTGTGAACTTGGTGTCGGGGGTCCGCAAGGGGTCCCCGGTGCCGGAGTGAAGATGTTAAGTTCACTGCCTGGGGAGTACGGTCGCAAGGCTG"

phyloseq.tree <- root(fitGTR$tree, outgroup = outgroupseq, resolve.root = TRUE)


## G. Create phyloseq object ----

mapping <- read.csv("mapping_file.csv" , row.names = NULL) 
names(mapping)
# transform sample_names column from mapping to character class
mapping$sample_names <- as.character(mapping$sample_names)
row.names(mapping) <- mapping$sample_names
mapping$group <- gsub('LP','OA',mapping$group)
otu_table <- otu_table (seqtab.nochim, taxa_are_rows= FALSE ) 
mapping_file <- sample_data(mapping) 
taxa_table <- tax_table(taxa) 

ps <- phyloseq(otu_table, mapping_file, taxa_table)




## H. Identify Contaminants ----



sample_names(ps) # bead blanks: 1,2; tube blanks: 36,37,38
vector_for_decontam <- c(rep(TRUE, 2), rep(FALSE, 33), rep(TRUE, 3))

contam_df <- isContaminant(ps, neg=vector_for_decontam)

table(contam_df$contaminant) # identified 1 as contaminant

#which ASV is contaminant
head(which(contam_df$contaminant))
#seq 1
contam_df[1,]
# TACGTAGGTGGCAAGCGTTGTCCGGAATTATTGGGCGTAAAGGGCTCGCAGGCGGTTCCTTAAGTCTGATGTGAAAGCCCCCGGCTCAACCGGGGAGGGTCATTGGAAACTGGGGAACTTGAGTGCAGAAGAGGAGAGTGGAATTCCACGTGTAGCGGTGAAATGCGTAGAGATGTGGAGGAACACCAGTGGCGAAGGCGACTCTCTGGTCTGTAACTGACGCTGAGGAGCGAAAGCGTGGGGAGCGAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAGTGCTAAGTGTTAGGGGGTTTCCGCCCCTTAGTGCTGCAGCTAACGCATTAAGCACTCCGCCTGGGGAGTACGGTCGCAAGACTG 

# make a list of the contaminant
contaminant <- ("TACGTAGGTGGCAAGCGTTGTCCGGAATTATTGGGCGTAAAGGGCTCGCAGGCGGTTCCTTAAGTCTGATGTGAAAGCCCCCGGCTCAACCGGGGAGGGTCATTGGAAACTGGGGAACTTGAGTGCAGAAGAGGAGAGTGGAATTCCACGTGTAGCGGTGAAATGCGTAGAGATGTGGAGGAACACCAGTGGCGAAGGCGACTCTCTGGTCTGTAACTGACGCTGAGGAGCGAAAGCGTGGGGAGCGAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAGTGCTAAGTGTTAGGGGGTTTCCGCCCCTTAGTGCTGCAGCTAACGCATTAAGCACTCCGCCTGGGGAGTACGGTCGCAAGACTG")


# Remove contaminant
#get all taxa names
allTaxa <- taxa_names(ps)
#returns list of all taxa, except contaminants
allTaxa <- allTaxa[!(allTaxa %in% contaminant)]
#keeps taxa in allTaxa and gets rid of anthing else
ps_filt <- prune_taxa(allTaxa, ps)


## I. Remove Chloroplasts, Mitochondria, and Animal sequences ----
# removing chloroplasts, mitochondria FOR SILVA 132 - silva 128 has chloro as a Class
ps_noChloro <- subset_taxa(ps_filt, Order !="o__Chloroplast")
ps_noChloro_noMito <- subset_taxa(ps_noChloro, Family != "f__Mitochondria")
ps_noChloro_noMito_noMollusc <- subset_taxa(ps_noChloro_noMito, Phylum != "p__Mollusca")
ps_noChloro_noMito_noMollusc_noAni <- subset_taxa(ps_noChloro_noMito_noMollusc, Kingdom != "p__Animalia")
# checking to see how many sequences lost during these filtering steps
filter_seq_loss <- (1-(sample_sums(ps_noChloro_noMito_noMollusc_noAni)/sample_sums(ps)))*100
print(filter_seq_loss)
# need to remove samples that have a zero count, I think it is disrupting deseq2 down the line
data <- prune_samples( sample_sums(ps_noChloro_noMito_noMollusc_noAni) > 0, ps_noChloro_noMito_noMollusc_noAni )



# creating a list to subset phyloseq object to remove negative and positive controls
names <- sample_names(data)
sample_names(data)
# list of the samples we want to keep
names <- names[-c(1,2,12:14,24:36)]


# II. Taxonomic diversity analyses ----

## A. Alpha diversity analysis ----

# prune samples to just keep those from low ph and control groups
# use to perform alpha diversity statistics on raw reads
ps.alpha <- prune_samples(names,data)

# apply estimate richness function from phyloseq as observed asvs and using the shannon index
alpha.measure <- estimate_richness(ps.alpha, measures=c("Observed", "Shannon", "Simpson", "Chao1"))

# add new key column with sample names
alpha.measure$sample_name <- names 

# join these two dataframes so now we have observed and shannon index along with appropriate key columns for easy graphing
alpha.measure.final <- full_join(alpha.measure, mapping, by = c("sample_name" = "sample_names"))
# removes samples from the final dataframe that we are not dealing with here
alpha.measure.final <- na.omit(alpha.measure.final)

### i. Alpha statistics ----

kruskal.test(Observed ~ bucket.id, data=alpha.measure.final)
# Kruskal-Wallis chi-squared = 5.4678, df = 5, p-value = 0.3615

kruskal.test(Shannon ~ bucket.id, data=alpha.measure.final)
# Kruskal-Wallis chi-squared = 10.029, df = 5, p-value = 0.07441

kruskal.test(Chao1 ~ bucket.id, data=alpha.measure.final)
# Kruskal-Wallis chi-squared = 5.4678, df = 5, p-value = 0.3615


### ii. Alpha figure ----

alpha.cols <- c("#721F81FF", "#B63679FF", "#F1605DFF", "#FEAF77FF", "#FCFDBFFF")
#alpha.breaks <- as.character(observed.bucket$bucket.id)
alpha.ids <- c("C38-1",  "C38-2",  "C51", 
               "LP38-1", "LP38-2", "LP51")

# Creating vectors of names that will be used as labels in the ggplot commands below
supp.labs <- c("38-Days", "51-Days")
names(supp.labs) <- c("38", "51")

ab <-ggplot(alpha.measure.final, aes(x=bucket.id, y=Observed, fill=treatment))+
  geom_boxplot(outlier.color="white", width = .9) + geom_point() +
  theme(legend.title=element_blank()) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 12)) +
  xlab("Treatment") + ylab("Observed ASVs") +
  facet_grid(.~age, scales="free", labeller = labeller(age = supp.labs), space = "free") +
  scale_fill_manual(values = magma(5)[3:5] ) +
  scale_x_discrete(breaks=c("C1","C2", "C3", "OA1", "OA2", "OA3"),
                   labels=c("C-38D-1","C-38D-2", "C-51D-3", "LP-38D-1", "LP-38D-2","LP-51D-3")) +
  ylim(150,400)
ab


bc <-ggplot(alpha.measure.final, aes(x=bucket.id, y=Shannon, fill=treatment))+
  geom_boxplot(outlier.color="white", width = .9) + geom_point() +
  theme(legend.title=element_blank()) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 12)) +
  xlab("Treatment") + ylab("Shannon Index") +
  facet_grid(.~age, scale="free", labeller = labeller(age = supp.labs), space = "free") +
  scale_fill_manual(values = magma(5)[3:5] ) +
  scale_x_discrete(breaks=c("C1","C2", "C3", "OA1", "OA2", "OA3"),
                   labels=c("C-38D-1","C-38D-2", "C-51D-3", "LP-38D-1", "LP-38D-2","LP-51D-3")) +
  ylim(2.7,4.5)
bc

# cd <-ggplot(alpha.measure.final, aes(x=bucket.id, y=Chao1, fill=treatment))+
#   geom_boxplot(outlier.color="white", width = .9) + geom_point() +
#   theme(legend.title=element_blank()) + 
#   theme_bw() + 
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         legend.position="none") +
#   xlab("Treatment") + ylab("Chao1") +
#   facet_grid(.~age, scales="free", labeller = labeller(age = supp.labs), space = "free") +
#   scale_fill_manual(values = magma(5)[3:5] ) +
#   scale_x_discrete(breaks=c("C1","C2", "C3", "OA1", "OA2", "OA3"),
#                    labels=c("C-38D-1","C-38D-2", "C-51D-3", "LP-38D-1", "LP-38D-2","LP-51D-3")) 
# cd

# get onto one image, both plots same size
xx <- ggplotGrob(ab)
yy <- ggplotGrob(bc)
#zz <- ggplotGrob(cd)


yy$heights[10] <- xx$heights[10] 

xyz <- rbind(xx, yy)

grid.newpage()
grid.draw(xyz, recording=TRUE)





#### A. end ####


## B. Normalizing sample depth ----
# Using variance stabalizing transformation 

# Will be using Deseq2 to variance stabilizing transformation, instead of the classic rarefaction or
# turning counts into proportions. This is recommended by the McMurdie and Holmes 2014 Plos Computational Biology paper, Waste not Want not.
# using this tutorial: https://astrobiomike.github.io/amplicon/dada2_workflow_ex#analysis-in-r

# first we need to make a DESeq2 object
# the design I set, bucket.id, will create groups based on both treatment and time (point sampled)
data_deseq <- phyloseq_to_deseq2(data, ~ bucket.id)

gm_mean <- function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans <- apply(counts(data_deseq), 1, gm_mean)
data_deseq <- estimateSizeFactors(data_deseq, geoMeans = geoMeans, type="poscount")

data_deseq_vst <- varianceStabilizingTransformation(data_deseq, blind = FALSE, fitType = "parametric")

# and here is pulling out our transformed table
vst_transform <- assay(data_deseq_vst)

## C. Hierarchical clustering ----

euc_dist <- dist(t(vst_transform))
euc_clust <- hclust(euc_dist, method="ward.D2")

# hclust objects like this can be plotted with the generic plot() function
plot(euc_clust) 



# okay. the VST log-like transformation produces negative values
# for counts that are <1
# therefore I am going to set these to zero so that distance metrics can be applied
# https://www.bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#negative-numbers-in-my-transformed-data-table

vst_transform[vst_transform < 0.0] <- 0.0

# making our phyloseq object with transformed table
vst_transform_phy <- otu_table(vst_transform, taxa_are_rows=T)

# creating phyloseq object
ps.vst <- phyloseq(vst_transform_phy, mapping_file, taxa_table)

# Calculate RA per sample
ps.vst.RA <- transform_sample_counts(ps.vst, function(x) 100 * x/sum(x))







##### C. end ####



## D. Taxonomic Ordination ----



ps.control <- subset_samples(ps.vst.RA, sample_data(ps.vst.RA)$treatment == "control")
ps.lp <- subset_samples(ps.vst.RA, sample_data(ps.vst.RA)$treatment == "low-ph")

ps.final <- merge_phyloseq(ps.control,ps.lp)

# creating a factor so I can use a discrete scale later for figure
sample_data(ps.final)$age <- as.factor(sample_data(ps.final)$age)

# Ordination between control and molybdate samples
ord.nmds.bray <- ordinate(ps.final, method="NMDS", distance="bray")
ord.nmds.bray

bucket.labs <- c("C-38D-1","C-38D-2", "C-51D-3", "LP-38D-1", "LP-38D-2","LP-51D-3")

bucket.shapes <-c(17, 17, 17,
                  19, 19, 19)

bucket.cols <- c("#fa9fb5", "#c51b8a", "#7a0177",
                 "#fed976", "#fd8d3c", "#f03b20")

plot_ordination(ps.final, ord.nmds.bray, shape="bucket.id", color="bucket.id") +
  geom_point(size = 4) +
  theme(legend.title=element_blank()) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_colour_manual(name = "Bucket ID",
                      labels = bucket.labs,
                      values = bucket.cols ) +   
  scale_shape_manual(name = "Bucket ID",
                     labels = bucket.labs,
                     values = bucket.shapes) +
  annotate(geom="text", x=.5, y=0.4, label="Stress = 0.064", color="black") 

#### D. end ####



## E. Ordination Statistics ----


### i. Overall - Treatment ----

# distance call
Dist.cntrl.lp = phyloseq::distance(ps.final, method = "bray", type="samples")

sample.data <- data.frame(sample_data(ps.final))

# betadisper: dispersion test
cntrl.lp.beta <- betadisper(Dist.cntrl.lp, sample.data$treatment)
permutest(cntrl.lp.beta)
#           Df   Sum Sq   Mean Sq     F N.Perm Pr(>F)
# Groups     1 0.001179 0.0011790 0.419    999  0.538
# Residuals 16 0.045018 0.0028136       

# Adonis test
adonis2(Dist.cntrl.lp ~ treatment, data = sample.data)
#           Df SumOfSqs      R2      F Pr(>F)    
# treatment  1  0.72363 0.37236 9.4924  0.001 ***
# Residual  16  1.21972 0.62764                  
# Total     17  1.94335 1.00000   



### ii. Overall - Treatment/Time ----

# betadisper: dispersion test
cntrl.lp.beta.group <- betadisper(Dist.cntrl.lp, sample.data$group)
permutest(cntrl.lp.beta.group)
#           Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)  
# Groups     3 0.024918 0.0083060 4.4403    999  0.024 *
# Residuals 14 0.026188 0.0018706  

# Adonis test
adonis2(Dist.cntrl.lp ~ group, data = sample.data)
#           Df SumOfSqs      R2      F Pr(>F)    
# group     3  1.26449 0.65068 8.6926  0.001 ***
# Residual 14  0.67885 0.34932                  
# Total    17  1.94335 1.00000     


### iii. Among treatment-time groups ----
# Pairwise tests between treatment-time points (groups column in sample.data

ps.c38 <- subset_samples(ps.final, group=="C38")
ps.c51 <- subset_samples(ps.final, group=="C51")
ps.lp38 <- subset_samples(ps.final, group=="OA38")
ps.lp51 <- subset_samples(ps.final, group=="OA51")


# C38:C51
ps.c38.c51 <- merge_phyloseq(ps.c38,ps.c51)
# make a data frame from the sample_data
data.c38.c51 <- data.frame(sample_data(ps.c38.c51))
Dist.c38.c51 = phyloseq::distance(ps.c38.c51, method = "bray", type="samples")

# betadisper: dispersion test
c38.c51.beta <- betadisper(Dist.c38.c51, data.c38.c51$group)
permutest(c38.c51.beta)
#           Df    Sum Sq   Mean Sq     F N.Perm Pr(>F)
# Groups     1 0.0012086 0.0012086 0.492    999  0.515
# Residuals  7 0.0171944 0.0024563             

# Adonis test
adonis2(Dist.c38.c51 ~ group, data = data.c38.c51)
#           Df SumOfSqs      R2      F Pr(>F)  
# group     1  0.27857 0.48002 6.4619  0.016 *
# Residual  7  0.30177 0.51998                
# Total     8  0.58034 1.00000   

# C38:lp38
ps.c38.lp38 <- merge_phyloseq(ps.c38,ps.lp38)
# make a data frame from the sample_data
data.c38.lp38 <- data.frame(sample_data(ps.c38.lp38))
Dist.c38.lp38 = phyloseq::distance(ps.c38.lp38, method = "bray", type="samples")

# betadisper: dispersion test
c38.lp38.beta <- betadisper(Dist.c38.lp38, data.c38.lp38$group)
permutest(c38.lp38.beta)
#           Df    Sum Sq   Mean Sq      F N.Perm Pr(>F)  
# Groups     1 0.0057013 0.0057013 5.0745    999  0.047 *
# Residuals 10 0.0112352 0.0011235  

# Adonis test
adonis2(Dist.c38.lp38 ~ group, data = data.c38.lp38)
#         Df SumOfSqs      R2      F Pr(>F)   
# group     1  0.61271 0.52975 11.265  0.003 **
# Residual 10  0.54389 0.47025                 
# Total    11  1.15660 1.00000  


# C51:lp51
ps.c51.lp51 <- merge_phyloseq(ps.c51,ps.lp51)
# make a data frame from the sample_data
data.c51.lp51 <- data.frame(sample_data(ps.c51.lp51))
Dist.c51.lp51 = phyloseq::distance(ps.c51.lp51, method = "bray", type="samples")

# betadisper: dispersion test
c51.lp51.beta <- betadisper(Dist.c51.lp51, data.c51.lp51$group)
permutest(c51.lp51.beta)
#           Df    Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups     1 0.0021518 0.0021518 0.5757    719 0.6014
# Residuals  4 0.0149509 0.0037377   

# Adonis test
adonis2(Dist.c51.lp51 ~ group, data = data.c51.lp51)
#           Df SumOfSqs      R2      F Pr(>F)
# group     1  0.29558 0.68653 8.7603    0.1
# Residual  4  0.13496 0.31347              
# Total     5  0.43055 1.00000     


# lp38:lp51
ps.lp38.lp51 <- merge_phyloseq(ps.lp38,ps.lp51)
# make a data frame from the sample_data
data.lp38.lp51 <- data.frame(sample_data(ps.lp38.lp51))
Dist.lp38.lp51 = phyloseq::distance(ps.lp38.lp51, method = "bray", type="samples")

# betadisper: dispersion test
lp38.lp51.beta <- betadisper(Dist.lp38.lp51, data.lp38.lp51$group)
permutest(lp38.lp51.beta)
#           Df    Sum Sq   Mean Sq      F N.Perm Pr(>F)    
# Groups     1 0.0224936 0.0224936 17.507    999  0.001 ***
# Residuals  7 0.0089938 0.0012848     

# Adonis test
adonis2(Dist.lp38.lp51 ~ group, data = data.lp38.lp51)
#           Df SumOfSqs      R2      F Pr(>F)  
# group     1  0.26230 0.41024 4.8692  0.011 *
# Residual  7  0.37708 0.58976                
# Total     8  0.63938 1.00000  



# lp38:c51
ps.lp38.c51 <- merge_phyloseq(ps.lp38,ps.c51)
# make a data frame from the sample_data
data.lp38.c51 <- data.frame(sample_data(ps.lp38.c51))
Dist.lp38.c51 = phyloseq::distance(ps.lp38.c51, method = "bray", type="samples")

# betadisper: dispersion test
lp38.c51.beta <- betadisper(Dist.lp38.c51, data.lp38.c51$group)
permutest(lp38.c51.beta)
#            Df    Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups     1 0.0091655 0.0091655 2.6529    999   0.12
# Residuals  7 0.0241841 0.0034549     

# Adonis test
adonis2(Dist.lp38.c51 ~ group, data = data.lp38.c51)
#           Df SumOfSqs    R2      F Pr(>F)  
# group     1  0.40270 0.489 6.6985  0.013 *
# Residual  7  0.42083 0.511                
# Total     8  0.82353 1.000  




#### E. end ####




## F. Taxonomic Differential Expression with DESeq2 ----

# using no Relative Abundance data
ps.control2 <- subset_samples(ps.vst, sample_data(ps.vst)$treatment == "control")
ps.lp2 <- subset_samples(ps.vst, sample_data(ps.vst)$treatment == "low-ph")

ps.final_des <- merge_phyloseq(ps.control2,ps.lp2)

group_des <- phyloseq_to_deseq2(ps.final_des, ~ group)
group_des <- DESeq(group_des)

### i. C38:C51 ----
# pulling out our results table, we specify the object, the p-value we are going to use to filter our results, and what contrast we want to consider by first naming the column, then the two groups we care about
deseq_res_C38_vs_C51 <- results(group_des, contrast = c("group", "C38", "C51"),  alpha=0.01)

sigtab <- deseq_res_C38_vs_C51

sigtab <- cbind(as(sigtab, "data.frame"), as(tax_table(ps.final_des)[rownames(sigtab), ], "matrix"))
sigtab <- sigtab %>% filter(padj <= 0.05)
head(sigtab)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))
ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))


### ii. C38:OA38 ----

# pulling out our results table, we specify the object, the p-value we are going to use to filter our results, and what contrast we want to consider by first naming the column, then the two groups we care about
deseq_res_C38_vs_OA38 <- results(group_des, contrast = c("group", "C38", "OA38"),  alpha=0.01)

sigtab2 <- deseq_res_C38_vs_OA38

sigtab2 <- cbind(as(sigtab2, "data.frame"), as(tax_table(ps.final_des)[rownames(sigtab2), ], "matrix"))
sigtab2 <- sigtab2 %>% filter(padj <= 0.05)


theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x2 = tapply(sigtab2$log2FoldChange, sigtab2$Phylum, function(x2) max(x2))
x2 = sort(x2, TRUE)
sigtab2$Phylum = factor(as.character(sigtab2$Phylum), levels=names(x2))
# Genus order
x2 = tapply(sigtab2$log2FoldChange, sigtab2$Family, function(x2) max(x2))
x2 = sort(x2, TRUE)
sigtab2$Family = factor(as.character(sigtab2$Family), levels=names(x2))
ggplot(sigtab2, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))


### iii. OA38:OA51 ----

# pulling out our results table, we specify the object, the p-value we are going to use to filter our results, and what contrast we want to consider by first naming the column, then the two groups we care about
deseq_res_OA51_vs_OA38 <- results(group_des, contrast = c("group", "OA51", "OA38"),  alpha=0.01)

sigtab3 <- deseq_res_OA51_vs_OA38

sigtab3 <- cbind(as(sigtab3, "data.frame"), as(tax_table(ps.final_des)[rownames(sigtab3), ], "matrix"))
sigtab3 <- sigtab3 %>% filter(padj <= 0.05)


theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x3 = tapply(sigtab3$log2FoldChange, sigtab3$Phylum, function(x3) max(x3))
x3 = sort(x3, TRUE)
sigtab3$Phylum = factor(as.character(sigtab3$Phylum), levels=names(x3))
# Genus order
x3 = tapply(sigtab3$log2FoldChange, sigtab3$Family, function(x3) max(x3))
x3 = sort(x3, TRUE)
sigtab3$Family = factor(as.character(sigtab3$Family), levels=names(x3))
ggplot(sigtab3, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

## G. Taxonomy tests ----


### i. Treatments/Age ----
# Comparison on means, KW, Dunn 

#collapse ASVs at Order level
AvgRA_o = tax_glom(ps.final, taxrank="Order", NArm = FALSE)

#filter out taxa that don't vary > 0.2 %
AvgRA99_O = filter_taxa(AvgRA_o, function(x) var(x) > .002, TRUE)
df_o <- psmelt(AvgRA99_O)
#group and calculate mean, sd and se for different taxonomic levels
grouped_o <- group_by(df_o, group, Phylum, Class, Order)
avgs_o <- dplyr::summarise(grouped_o, mean=mean(Abundance), sd=sd(Abundance), se=se(Abundance))

avgs_filt <- filter(avgs_o, mean > 0)
# I want to read out this table for sup data
write.csv(avgs_filt, "OA_treat-time_fam_avg_sd_se.csv")


avgs_filt$ST <- as.factor(avgs_filt$group)
avgs_filt_st <- melt(data.frame(avgs_filt), id.vars=c("Order", "group"),
                     measure.vars=c("mean"))

#start
chisq = NULL
pvals = NULL
listofcats_sig = NULL
DataSet = df_o
taxa = as.factor(unique(df_o$Order))


for (cat in taxa) {
  new_df <- subset(DataSet, DataSet$Order == cat)
  new_df$ST <- as.factor(new_df$group)
  kw = kruskal.test(Abundance ~ ST, data=new_df)
  chisq = c(kw$statistic, chisq)
  pvals = c(kw$p.value, pvals)
  if (kw$p.value <= 0.05) {
    listofcats_sig = c(cat, listofcats_sig)
  }
}

pvals.bonf = p.adjust(pvals, method="bonferroni")
df_taxa = data.frame(rev(taxa), chisq, pvals, pvals.bonf)
write.csv(df_taxa, 'OA_treat-time_VST_Mean_FAM_KW.csv')

# Order summary
pvals.dunn = NULL
pvals.dunn.bonf = NULL
Zsc.dunn = NULL
comparison.dunn = NULL
cats.dunn = NULL

for (cat in listofcats_sig){
  new_df <- subset(DataSet, DataSet$Order == cat)
  new_df$ST <- as.factor(new_df$group)
  dT = dunnTest(Abundance ~ ST, data = new_df, method = "bonferroni")
  for (i in 1:length(dT$res$Comparison)) {
    print(cat)
    pvals.dunn = c(dT$res$P.unadj[i], pvals.dunn)
    pvals.dunn.bonf = c(dT$res$P.adj[i], pvals.dunn.bonf)
    Zsc.dunn = c(dT$res$Z[i], Zsc.dunn)
    comparison.dunn = c(dT$res$Comparison[i], comparison.dunn)
    cats.dunn = c(cat, cats.dunn)
  }
}

df_taxa = data.frame(cats.dunn, comparison.dunn, Zsc.dunn, pvals.dunn, pvals.dunn.bonf)

tax_oi <- filter(df_taxa, pvals.dunn.bonf < 0.05)

#write.csv(df_taxa, 'OA_treat-time_VST_Mean_FAM_KW_Dunn_OI.csv')

write.csv(df_taxa, file = paste0("OA_treat-time_VST_Mean_FAM_KW_Dunn_OI", format(Sys.time(), "%d-%b-%Y %H.%M"), ".csv"), row.names = F)



## H. Taxonomic bar chart figure ----

#group and calculate mean, sd and se for different taxonomic levels
groups <- group_by(df_o, group, Phylum, Class, Order)
buckets <- group_by(df_o, bucket.id, Phylum, Class, Order)
head(groups)

groups_avgs <- dplyr::summarise(groups, mean=mean(Abundance), sd=sd(Abundance), se=se(Abundance))

buckets_avgs <- dplyr::summarise(buckets, mean=mean(Abundance), sd=sd(Abundance), se=se(Abundance))


tax_oi_filt <- tax_oi %>% filter(comparison.dunn != "C38 - OA51") %>% filter(comparison.dunn != "C51 - OA38") 

orders <- tax_oi_filt$cats.dunn


# groups_avgs was created above, groups melted by order
sig_orders <- filter(buckets_avgs, Order %in% orders)

sig_orders <- filter(buckets_avgs, Order %in% orders)

ggplot(sig_orders, aes(x=bucket.id, y=mean, fill=Order))+
  geom_bar(stat="identity", position = "dodge", color="black") +
  geom_errorbar(aes(ymin=(mean-se), ymax=(mean+se)),
                width=.4, position=position_dodge(.9)) +
  theme(legend.title=element_blank()) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(face="bold"),
        legend.position="bottom") +
  xlab("Treatment-Age") + ylab("Mean Percent Abundance") +
  scale_fill_manual(values = viridis(19)) +
  scale_x_discrete(breaks=c("C1","C2", "C3", "OA1", "OA2", "OA3"),
                   labels=c("C-38D-1", "C-38D-2", "C-51D-3", "LP-38D-1", "LP-38D-2", "LP-51D-3")) 



 


# III. Functional diversity analyses ----

## A. Prepare export files ----

# Need to run picrust2 in terminal, prep files for export

# ps.final is the main data set with vst transformed asv counts
# re extract information to feed into picrust in terminal using python

# OTU/ASV (BIOM) Table 
# rows: ASV1, ASV2, etc
# columns: samples


# Export feature/OTU table
#otu <- t(as(otu_table(ps.final),"matrix")) # 't' to transform if taxa_are_rows=FALSE
#if taxa_are_rows=TRUE
otu <- as(otu_table(ps.alpha),"matrix")
otu_t <- t(otu) 
otu_biom <- make_biom(data=otu_t)
write_biom(otu_biom,"asv_abundance_alpha-2.biom")
# this seems to work




# Taxonomy table (FASTA/FNA)

# character list
# names are ASV labels

# filter out sequences that are not present in dataset
seqs2 <- data.frame(seqs)

seqs3 <- filter(seqs2, seqs %in% asv.key$sequence)
names(seqs3) <- "sequence"
# join id column by sequence column
seqs3_key <- left_join(seqs3, asv.key, by = "sequence")
row.names(seqs3_key) <- seqs3_key$id
seqs4 <- seqs3_key %>% select(-c(id))


# Prepare the header lines for the FASTA file
header_lines <- paste(">", seqs4$sequence, "\n", seqs4$sequence, sep = "")

# Write the sequences to a FASTA file
writeLines(header_lines, "asv_sequences.fna")
asv_sequences_out <- file("asv_sequences.fna", "a")
#writeLines(seqs4$sequence, asv_sequences_out)
close(asv_sequences_out)
# I think this actually works


## B. Import and analyses ----

### iv. Kegg redo using manual module numbers ----

# January 30 2025

# import stratified and unstratified ko output
ko_strat_all <- read.delim("picrust2_out_pipeline_4/KO_metagenome_out/pred_metagenome_contrib.tsv")
ko_unstrat <- read.delim("picrust2_out_pipeline_4/KO_metagenome_out/pred_metagenome_unstrat.tsv")


# prep tables
link_table <- read.csv(file = "picrust2_out_pipeline_4/KO_metagenome_out/link_table.csv")
module_genes <- read.csv("picrust2_out_pipeline_4/KO_metagenome_out/module_genes.csv")
module_genes_full <- left_join(module_genes, link_table, by = "module")


# # Perform pathway differential abundance analysis (DAA) using ALDEx2 method
colnames(ko_unstrat) <- paste0(c("function", metadata$sample_names))
rownames(ko_unstrat) <- ko_unstrat$`function`
ko_final <- ko_unstrat[,c(2:19)]

# run daa
daa_results_deseq_df_2 <- pathway_daa(abundance = ko_final, metadata = metadata, group = "group", daa_method = "DESeq2", select = NULL, reference = NULL, p.adjust = "bonferroni")


# filter for significance
daa_results_deseq_sig <- daa_results_deseq_df_2 %>% filter(p_adjust < 0.05)
# keep only results from modules of interest
daa_results_deseq_sig_mods <- daa_results_deseq_sig %>% filter(feature %in% module_genes_full$knumber)
names(daa_results_deseq_sig_mods)[1] <- "knumber" # making names match for another join

# joining by k number so now we associate sig tests with modules that have useful bioligical info
daa_results_deseq_sig_mods_combined <- left_join(daa_results_deseq_sig_mods, module_genes_full, by = "knumber")


# joining strat output and module metadata
names(ko_strat_all)[2] <- "knumber"

# creat new storage object for main taxa.df so we do not corrupt it
taxa.ko <- taxa.df
# pull sequences as row names into new col
taxa.ko$sequence <- row.names(taxa.ko)

# merge with metacyc data by seq to add taxonomic information to metacyc
ko_strat_all_tax <- ko_strat_all %>% left_join(taxa.ko, by = c("taxon" = "sequence"))
# join the rest of the metadata for group info
ko_strat_all_tax_meta <- ko_strat_all_tax %>% left_join(metadata, by = c("sample" = "sample_names"))

# Going to keep sequences that have a significant DAA between C38-LP38 and C51-LP51
# create list of genes of interest
goi.38 <- daa_results_deseq_sig_mods %>% filter(group1 == "C38" & group2 == "OA38")
goi.51 <- daa_results_deseq_sig_mods %>% filter(group1 == "C51" & group2 == "OA51")

# now filter large dataset usin these genes
ko_strat_meta_38 <- ko_strat_all_tax_meta %>% filter(knumber %in% goi.38$knumber) %>% filter(age == 38)
ko_strat_meta_51 <- ko_strat_all_tax_meta %>% filter(knumber %in% goi.51$knumber) %>% filter(age == 51)

# Summarize by group (e.g., C38), Phylum, and Gene (KO number)
ko_strat_meta_38_sum <- ko_strat_meta_38 %>%
  group_by(group, Phylum, knumber) %>%
  dplyr::summarise(
    taxon_abun = sum(taxon_abun),
    taxon_function_abun = sum(taxon_function_abun),
    knumber = first(knumber)
    #.groups = "drop" # Optional: prevents grouping in the result
  )

# Summarize by group (e.g., C51) 
ko_strat_meta_51_sum <- ko_strat_meta_51 %>%
  group_by(group, Phylum, knumber) %>%
  dplyr::summarise(
    taxon_abun = sum(taxon_abun),
    taxon_function_abun = sum(taxon_function_abun),
    knumber = first(knumber)
    #.groups = "drop" # Optional: prevents grouping in the result
  )


# finally need to add M numbers on to associate genes with specific pathways
module_genes_full_38 <- module_genes_full %>% filter(knumber %in% ko_strat_meta_38_sum$knumber)
module_genes_full_51 <- module_genes_full %>% filter(knumber %in% ko_strat_meta_51_sum$knumber)


# now need to merge M numbers onto filtered strat outputs
ko_strat_module_genes_38 <- left_join(ko_strat_meta_38_sum, module_genes_full_38, by = "knumber")
ko_strat_module_genes_51 <- left_join(ko_strat_meta_38_sum, module_genes_full_51, by = "knumber")






ko_strat_module_genes_38$process <- factor(ko_strat_module_genes_38$process, levels = c("denitrification", "nitrate reduction", "sulfate reduction",              
                                                                    "acetoclastic methanogenesis", "hydrogenotrophic methanogenesis", "photosynthesis",                 
                                                                    "photoautotrophy", "sulfide oxidation", "methylotrophic methanogenesis",
                                                                    "ammonium oxidation"))

# plot 
z_fxn_abun <- ggplot(ko_strat_module_genes_38, aes(x = process, y = taxon_function_abun, fill = group)) +
  geom_col(position = "dodge") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "bottom") +
  coord_flip() +
  xlab("Biochemical Reactions (Kegg Modules)") +
  ylab("Function (KO) Abundance") +
  scale_x_discrete(breaks=c("denitrification", "nitrate reduction", "sulfate reduction",              
                            "acetoclastic methanogenesis", "hydrogenotrophic methanogenesis", "photosynthesis",                 
                            "photoautotrophy", "sulfide oxidation", "methylotrophic methanogenesis",
                            "ammonium oxidation"),
                   labels=c("Denitrification; \nM00529", "Nitrate reduction; \nM00530, M00531", "Sulfate reduction; \nM00176, M00569",              
                            "Acetoclastic methanogenesis; \nM00357", "Hydrogenotrophic methanogenesis; \nM00567", "Photosynthesis; \nM00161, M00163",                 
                            "Photoautotrophy; \nM00597, M00598", "Sulfide oxidation;	\nM00595", "Methylotrophic methanogenesis; \nM00356",
                            "Ammonium oxidation; \nM00528, M00973")) +
  scale_fill_discrete(name = "Treatment",
                      breaks = c("C38", "OA38"),
                      labels = c("Control", "Low-pH"))
z_fxn_abun





# first panel, going to summarize again 
ko_38_orth_sum <- ko_strat_module_genes_38 %>%
  group_by(group, process) %>%
  dplyr::summarise(unique_orthologs = n_distinct(knumber),
            orthologs_count = n()
            )

# make columns numeric instead of integer
ko_38_orth_sum$unique_orthologs <- as.numeric(ko_38_orth_sum$unique_orthologs)
ko_38_orth_sum$orthologs_count <- as.numeric(ko_38_orth_sum$orthologs_count)

ko_38_orth_sum$process <- factor(ko_38_orth_sum$process, levels = c("denitrification", "nitrate reduction", "sulfate reduction",              
                                                         "acetoclastic methanogenesis", "hydrogenotrophic methanogenesis", "photosynthesis",                 
                                                         "photoautotrophy", "sulfide oxidation", "methylotrophic methanogenesis",
                                                         "ammonium oxidation"))


# plot for number of orthologs
a_orthologs <- ggplot(ko_38_orth_sum, aes(x = process, y = orthologs_count, fill = group)) +
  geom_col(position = "dodge") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "bottom") +
  coord_flip() +
  xlab("Biochemical Reactions (Kegg Modules)") +
  ylab("Orthologs Represented") +
  scale_x_discrete(breaks=c("denitrification", "nitrate reduction", "sulfate reduction",              
                             "acetoclastic methanogenesis", "hydrogenotrophic methanogenesis", "photosynthesis",                 
                             "photoautotrophy", "sulfide oxidation", "methylotrophic methanogenesis",
                             "ammonium oxidation"),
                   labels=c("Denitrification; \nM00529", "Nitrate reduction; \nM00530, M00531", "Sulfate reduction; \nM00176, M00569",              
                            "Acetoclastic methanogenesis; \nM00357", "Hydrogenotrophic methanogenesis; \nM00567", "Photosynthesis; \nM00161, M00163",                 
                            "Photoautotrophy; \nM00597, M00598", "Sulfide oxidation;	\nM00595", "Methylotrophic methanogenesis; \nM00356",
                            "Ammonium oxidation; \nM00528, M00973")) +
  scale_fill_discrete(name = "Treatment",
                      breaks = c("C38", "OA38"),
                      labels = c("Control", "Low-pH"))
a_orthologs

# plot for number of unique orthologs
# ggplot(ko_38_orth_sum, aes(x = process, y = unique_orthologs, fill = group)) +
#   geom_col(position = "dodge")




# Summarize and calculate relative proportions
ko_38_tax_sum_relative <- ko_strat_module_genes_38 %>%
  group_by(process, Phylum) %>%  # Group by both Group1 and Group2
  summarise(TotalAbundance = sum(taxon_abun), .groups = "drop") %>%  # Get total per Group2 within Group1
  group_by(process) %>%  # Group by Group1 to calculate relative proportions
  mutate(RelativeProportion = TotalAbundance / sum(TotalAbundance)) %>%
  ungroup()  # Remove grouping for cleaner output

palette_12 <- c("#a50026",
                "#d73027",
                "#f46d43",
                "#fdae61",
                "#fee090",
                "#ffffbf",
                "#ffffff",
                "#e0f3f8",
                "#abd9e9",
                "#74add1",
                "#4575b4",
                "#313695")


ko_38_tax_sum_relative$process <- factor(ko_38_tax_sum_relative$process, levels = c("denitrification", "nitrate reduction", "sulfate reduction",              
                                                                    "acetoclastic methanogenesis", "hydrogenotrophic methanogenesis", "photosynthesis",                 
                                                                    "photoautotrophy", "sulfide oxidation", "methylotrophic methanogenesis",
                                                                    "ammonium oxidation"))


b_taxa <- ggplot(ko_38_tax_sum_relative, aes(x = process, y = RelativeProportion, fill = Phylum)) +
  geom_col() +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "bottom",
        axis.text.y=element_blank()) +
  scale_fill_manual( values = palette_12) + 
  coord_flip() +
  xlab("") +
  ylab("Relative Abundance") +
  scale_x_discrete(breaks=c("denitrification", "nitrate reduction", "sulfate reduction",              
                            "acetoclastic methanogenesis", "hydrogenotrophic methanogenesis", "photosynthesis",                 
                            "photoautotrophy", "sulfide oxidation", "methylotrophic methanogenesis",
                            "ammonium oxidation"),
                   labels=c("Denitrification; \nM00529", "Nitrate reduction; \nM00530, M00531", "Sulfate reduction; \nM00176, M00569",              
                            "Acetoclastic methanogenesis; \nM00357", "Hydrogenotrophic methanogenesis; \nM00567", "Photosynthesis; \nM00161, M00163",                 
                            "Photoautotrophy; \nM00597, M00598", "Sulfide oxidation;	\nM00595", "Methylotrophic methanogenesis; \nM00356",
                            "Ammonium oxidation; \nM00528, M00973")) 
b_taxa


# adjusting link table for final panel of calcification for this fig
# Add a new column with 0 or 1 based on Group1 identity
link_table_graph <- link_table %>%
  mutate(binary = if_else(effect == "increase", 1, 0))

process.red <- unique(ko_38_tax_sum_relative$process)
link_table_graph_red <- link_table_graph %>% filter(process %in% process.red)

link_table_graph_red$process <- factor(link_table_graph_red$process, levels = c("denitrification", "nitrate reduction", "sulfate reduction",              
                                                                                    "acetoclastic methanogenesis", "hydrogenotrophic methanogenesis", "photosynthesis",                 
                                                                                    "photoautotrophy", "sulfide oxidation", "methylotrophic methanogenesis",
                                                                                    "ammonium oxidation"))
link_table_graph_red$new <- rep(1, length(link_table_graph_red))

c_carb <- ggplot(link_table_graph_red, aes(x = process, y = new, fill = effect)) +
  geom_tile(position = , color = "black") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "bottom",
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        plot.margin = margin(10, 100, 10, 100)) +
  coord_flip() +
  ylab("Effect on \nCalcium Carbonate") +
  scale_x_discrete(breaks=c("denitrification", "nitrate reduction", "sulfate reduction",              
                            "acetoclastic methanogenesis", "hydrogenotrophic methanogenesis", "photosynthesis",                 
                            "photoautotrophy", "sulfide oxidation", "methylotrophic methanogenesis",
                            "ammonium oxidation"),
                   labels=c("Denitrification; \nM00529", "Nitrate reduction; \nM00530, M00531", "Sulfate reduction; \nM00176, M00569",              
                            "Acetoclastic methanogenesis; \nM00357", "Hydrogenotrophic methanogenesis; \nM00567", "Photosynthesis; \nM00161, M00163",                 
                            "Photoautotrophy; \nM00597, M00598", "Sulfide oxidation;	\nM00595", "Methylotrophic methanogenesis; \nM00356",
                            "Ammonium oxidation; \nM00528, M00973")) +
  scale_fill_manual(name = "Calcium Carbonate",
                    values = c("lightgray","darkgray"),
                    breaks = c("decrease", "increase"),
                    labels = c("Decrease", "Increase"))
c_carb






# get onto one image, both plots same size
grob_orth <- ggplotGrob(a_orthologs)
grob_tax <- ggplotGrob(b_taxa)
grob_carb <- ggplotGrob(c_carb)


fig1 <- cbind(grob_orth,
              grob_tax,
              grob_carb)

grid.newpage()
grid.draw(fig1, recording=TRUE)


















# IV. Shell analysis ----


shell.area.data <- read.csv("shell-data/oyster_seed_shell_data.csv") %>% filter(treatment != "molybdate") %>% filter(stage != "larval")
shell.mass.data <- read.csv("shell-data/shell_mass.csv") %>% filter(treatment != "molybdate")

shell.data <- cbind(shell.area.data, shell.mass.data)
# remove replicate columns
shell.data <- shell.data[,c(1:9)]
# convert area um to area mm

## A. Shell figures ----

### i. Shell area ----

supp.labs <- c("38-Days", "51-Days")
names(supp.labs) <- c("38", "52")

fig1a <- ggplot(shell.data, aes(bucket, area.um2, fill=treatment)) +
  geom_boxplot() + geom_jitter(width = 0.2) +
  theme(legend.title=element_blank()) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="none") +
  xlab("Treatment") + 
  ylab( expression(paste( "area", " (", mu,"m"^2, ")") ) ) +
  facet_grid(.~age, scales = "free", space = "free", labeller = labeller(age = supp.labs)) +
  scale_fill_manual(values = magma(5)[3:5] ) +
  scale_x_discrete(breaks=c("C1",  "C2",  "C3",  "OA1", "OA2", "OA3"),
                   labels=c("C-38D-1",  "C-38D-2",  "C-51D-3",  "LP-38D-1", "LP-38D-2", "LP-51D-3"))


### ii. Shell mass ----


fig1b <- ggplot(shell.data, aes(bucket, mass.mg, fill=treatment)) +
  geom_boxplot() + geom_jitter(width = 0.2) +
  theme(legend.title=element_blank()) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="none") +
  xlab("Treatment") + 
  ylab(bquote('mass (mg)')) +
  facet_grid(.~age, scales = "free", space = "free", labeller = labeller(age = supp.labs)) +
  scale_fill_manual(values = magma(5)[3:5] ) +
  scale_x_discrete(breaks=c("C1",  "C2",  "C3",  "OA1", "OA2", "OA3"),
                   labels=c("C-38D-1",  "C-38D-2",  "C-51D-3",  "LP-38D-1", "LP-38D-2", "LP-51D-3"))



# get onto one image, both plots same size
grob1a <- ggplotGrob(fig1a)
grob1b <- ggplotGrob(fig1b)
#zz <- ggplotGrob(cd)

grob1a$heights[10] <- grob1b$heights[10]

fig1 <- rbind(grob1a, grob1b)

grid.newpage()
grid.draw(fig1, recording=TRUE)



## B. Shell statistical comparison ----

### i. Shell area ----

# need to keep only 38 day old oysters because of our replicates issue
names(shell.data)
shell.data.38 <- shell.data %>% filter(age < 40)

#### a. raw data ----
area_mixed_model <- lmer(area.um2 ~ treatment + (1 | bucket/treatment), data = shell.data.38, REML = T)
summary(area_mixed_model)
anova(area_mixed_model)

# Extract fitted values and residuals
fitted_values <- fitted(area_mixed_model)
residuals <- resid(area_mixed_model)

# Plot fitted values vs residuals
plot(fitted_values, residuals, xlab = "Fitted Values", ylab = "Residuals", main = "Fitted Values vs Residuals")
abline(h = 0, col = "red")  # Add a horizontal line at y = 0 for reference
# this shows that variance changes among our fitted values -- log tansform below



#### b. log transform ----

area_mixed_model_log <- lmer(log(area.um2) ~ treatment + (1 | bucket/treatment), data = shell.data.38, REML = T)
summary(area_mixed_model_log)
anova(area_mixed_model_log)
qqnorm(residuals(area_mixed_model_log))
qqline(residuals(area_mixed_model_log))

plot(residuals(area_mixed_model_log)) 
abline(h = 0.6, col = "blue", lwd = 4, lty = 4)  
abline(h = -0.6, col = "blue", lwd = 4, lty = 4)  


which(residuals(area_mixed_model_log) > 0.6)
which(residuals(area_mixed_model_log) < -0.6)
shell.data.38[c(37:38),]


# Extract fitted values and residuals
fitted_values <- fitted(area_mixed_model_log)
residuals <- resid(area_mixed_model_log)

# Plot fitted values vs residuals
plot(fitted_values, residuals, xlab = "Fitted Values", ylab = "Residuals", main = "Fitted Values vs Residuals")
abline(h = 0, col = "red")  # Add a horizontal line at y = 0 for reference
# this shows that variance changes among our fitted values -- log tansform below

data_plot <- data.frame(residuals = residuals, fitted_values = fitted_values, treatment = shell.data.38$treatment)

# Plot residuals vs. fitted values by treatment
ggplot(data_plot, aes(x = fitted_values, y = residuals, color = treatment)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE) +
  facet_wrap(~ treatment) +
  labs(x = "Fitted Values", y = "Residuals") +
  theme_minimal()

treatments <- shell.data.38$treatment  # Replace 'your_data' with the actual name of your dataset

# Combine residuals and treatment information into a data frame
residuals_data <- data.frame(residuals = residuals, treatment = treatments)

# Create a boxplot
boxplot(residuals ~ treatment, data = residuals_data, 
        xlab = "Treatment", ylab = "Residuals", main = "Residuals by Treatment")




##### Levene's test ----
library(car)
library(rcompanion)
# Assuming 'group1' and 'group2' are vectors containing your data
# Run Levene's test
levene_area <- car::leveneTest(area.um2 ~ bucket, data = shell.data.38)
print(levene_area)

levene_area <- car::leveneTest(area.um2 ~ bucket, data = shell.data.38)

# pairwise 

lgroup1.1 <- ((shell.data.38 %>% filter(treatment == "control")) %>% filter(bucket == "C1"))
lgroup1.2 <- ((shell.data.38 %>% filter(treatment == "control")) %>% filter(bucket == "C2"))

lgroup2.1 <- ((shell.data.38 %>% filter(treatment == "low-ph")) %>% filter(bucket == "OA1"))
lgroup2.2 <- ((shell.data.38 %>% filter(treatment == "low-ph")) %>% filter(bucket == "OA2"))

# C1 v C2
levene.df.1 <- rbind(lgroup1.1, lgroup1.2)
levene_area1 <- car::leveneTest(area.um2 ~ bucket, data = levene.df.1)
levene_area1

# C1 v LP1
levene.df.2 <- rbind(lgroup1.1, lgroup2.1)
levene_area2 <- car::leveneTest(area.um2 ~ bucket, data = levene.df.2)
levene_area2

# C1 v LP2
levene.df.3 <- rbind(lgroup1.1, lgroup2.2)
levene_area3 <- car::leveneTest(area.um2 ~ bucket, data = levene.df.3)
levene_area3



# C2 v LP1
levene.df.4 <- rbind(lgroup1.2, lgroup2.1)
levene_area4 <- car::leveneTest(area.um2 ~ bucket, data = levene.df.4)
levene_area4

# C2 v LP2
levene.df.5 <- rbind(lgroup1.2, lgroup2.2)
levene_area5 <- car::leveneTest(area.um2 ~ bucket, data = levene.df.5)
levene_area5

# LP1 v LP2
levene.df.6 <- rbind(lgroup2.1, lgroup2.2)
levene_area6 <- car::leveneTest(area.um2 ~ bucket, data = levene.df.6)
levene_area6$`Pr(>F)`[1]


l.pvals <- c(levene_area1$`Pr(>F)`[1], levene_area2$`Pr(>F)`[1],
             levene_area3$`Pr(>F)`[1], levene_area4$`Pr(>F)`[1],
             levene_area5$`Pr(>F)`[1], levene_area6$`Pr(>F)`[1])
l.pvals.bonf <- l.pvals * 6
l.pvals.bonf 



##### F-test ----
names(shell.data.38)
group1 <- (shell.data.38 %>% filter(treatment == "control"))$area.um2
group2 <- (shell.data.38 %>% filter(treatment == "low-ph"))$area.um2

F_test_result <- var.test(group1, group2)
print(F_test_result)

# pairwise F test

group1.1 <- ((shell.data.38 %>% filter(treatment == "control")) %>% filter(bucket == "C1"))$area.um2
group1.2 <- ((shell.data.38 %>% filter(treatment == "control")) %>% filter(bucket == "C2"))$area.um2

group2.1 <- ((shell.data.38 %>% filter(treatment == "low-ph")) %>% filter(bucket == "OA1"))$area.um2
group2.2 <- ((shell.data.38 %>% filter(treatment == "low-ph")) %>% filter(bucket == "OA2"))$area.um2


pairwise_F_test_result_1 <- var.test(group1.1, group1.2)
print(pairwise_F_test_result_1)

pairwise_F_test_result_2 <- var.test(group1.1, group2.1)
print(pairwise_F_test_result_2)

pairwise_F_test_result_3 <- var.test(group1.1, group2.2)
print(pairwise_F_test_result_3)


pairwise_F_test_result_4 <- var.test(group1.2, group2.1)
print(pairwise_F_test_result_4)

pairwise_F_test_result_5 <- var.test(group1.2, group2.2)
print(pairwise_F_test_result_5)

pairwise_F_test_result_6 <- var.test(group2.1, group2.2)
print(pairwise_F_test_result_6)



f.pvals <- c(pairwise_F_test_result_1$p.value, pairwise_F_test_result_2$p.value, 
             pairwise_F_test_result_3$p.value, pairwise_F_test_result_4$p.value, 
             pairwise_F_test_result_5$p.value, pairwise_F_test_result_6$p.value)
f.pvals.bonf <- f.pvals * 6
f.pvals
f.pvals.bonf 



# Perform pairwise F-tests
pairwise_F_test_results <- list()
pairwise_F_test_results[[1]] <- pairwise_F_test_result_1
pairwise_F_test_results[[2]] <- pairwise_F_test_result_2
pairwise_F_test_results[[3]] <- pairwise_F_test_result_3
pairwise_F_test_results[[4]] <- pairwise_F_test_result_4
pairwise_F_test_results[[5]] <- pairwise_F_test_result_5
pairwise_F_test_results[[6]] <- pairwise_F_test_result_6

# Calculate variances
variances <- c(var(group1.1), var(group1.2), var(group2.1), var(group2.2))

critical_values <- qf(c(0.025, 0.975), 2, 3)

# Calculate confidence intervals for variances (95% confidence)
conf_intervals <- sapply(pairwise_F_test_results, function(x) {
  qf(c(0.025, 0.975), df1 = x$parameter[1], df2 = x$parameter[2]) * x$statistic / variances
})

conf_intervals <- sapply(pairwise_F_test_results, function(x) {
  qf(c(0.025, 0.975), df1 = x$parameter[1], df2 = x$parameter[2]) / x$statistic * variances
})

# Create dataframe for confidence intervals
ci_df <- data.frame(Group = c("Group 1", "Group 2", "Group 3", "Group 4"),
                    Variance = variances,
                    CI_Lower = variances / conf_intervals[1, ],
                    CI_Upper = variances / conf_intervals[2, ])

ci_df <- data.frame(Group = c("Group 1", "Group 2", "Group 3", "Group 4"),
                    Variance = variances,
                    CI_Lower = variances - conf_intervals[2 ],
                    CI_Upper = variances + conf_intervals[1 ])


# Plot
library(ggplot2)

ggplot(ci_df, aes(x = Group, y = Variance, fill = Group)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper),
                width = 0.2, position = position_dodge(width = 0.9)) +
  labs(x = "Groups", y = "Variance", title = "Group Variances with 95% Confidence Intervals") +
  theme_minimal() +
  theme(legend.position = "none")



# pairwise F test

group1.1 <- ((shell.data.38 %>% filter(treatment == "control")) %>% filter(bucket == "C1"))$area.um2
group1.2 <- ((shell.data.38 %>% filter(treatment == "control")) %>% filter(bucket == "C2"))$area.um2

group2.1 <- ((shell.data.38 %>% filter(treatment == "low-ph")) %>% filter(bucket == "OA1"))$area.um2
group2.2 <- ((shell.data.38 %>% filter(treatment == "low-ph")) %>% filter(bucket == "OA2"))$area.um2


pairwise_F_test_result_1 <- var.test(group1.1, group1.2)
print(pairwise_F_test_result_1)

pairwise_F_test_result_2 <- var.test(group1.1, group2.1)
print(pairwise_F_test_result_2)

pairwise_F_test_result_3 <- var.test(group1.1, group2.2)
print(pairwise_F_test_result_3)


pairwise_F_test_result_4 <- var.test(group1.2, group2.1)
print(pairwise_F_test_result_4)

pairwise_F_test_result_5 <- var.test(group1.2, group2.2)
print(pairwise_F_test_result_5)

pairwise_F_test_result_6 <- var.test(group2.1, group2.2)
print(pairwise_F_test_result_6)



f.pvals <- c(pairwise_F_test_result_1$p.value, pairwise_F_test_result_2$p.value, 
             pairwise_F_test_result_3$p.value, pairwise_F_test_result_4$p.value, 
             pairwise_F_test_result_5$p.value, pairwise_F_test_result_6$p.value)
f.pvals.bonf <- f.pvals * 6
f.pvals.bonf 










### ii. Shell mass ----

#### a. raw data ----
# need to keep only 38 day old oysters because of our replicates issue

mass_mixed_model <- lmer(mass.mg ~ treatment + (1 | bucket/treatment), data = shell.data.38, REML = T)
summary(mass_mixed_model)
anova(mass_mixed_model)

# Extract fitted values and residuals
fitted_values <- fitted(mass_mixed_model)
residuals <- resid(mass_mixed_model)

# Plot fitted values vs residuals
plot(fitted_values, residuals, xlab = "Fitted Values", ylab = "Residuals", main = "Fitted Values vs Residuals")

abline(h = 0, col = "red")  # Add a horizontal line at y = 0 for reference
# this shows that variance changes among our fitted values -- log tansform below

data_plot <- data.frame(residuals = residuals, fitted_values = fitted_values, treatment = shell.data.38$treatment)

# Plot residuals vs. fitted values by treatment
ggplot(data_plot, aes(x = fitted_values, y = residuals, color = treatment)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE) +
  facet_wrap(~ treatment) +
  labs(x = "Fitted Values", y = "Residuals") +
  theme_minimal()

treatments <- shell.data.38$treatment  # Replace 'your_data' with the actual name of your dataset

# Combine residuals and treatment information into a data frame
residuals_data <- data.frame(residuals = residuals, treatment = treatments)

# Create a boxplot
boxplot(residuals ~ treatment, data = residuals_data, 
        xlab = "Treatment", ylab = "Residuals", main = "Residuals by Treatment")


# is bucket important?
mass_mixed_model_g <- lm(mass.mg ~ treatment, data = shell.data.38, REML = T)
anova(mass_mixed_model, mass_mixed_model_g)


# now lets try the bayesian equivalent as implemented in brms
# mass_mixed_bayesian_model <- brm(
#   mass.mg ~ treatment + (1 | bucket/treatment),
#   data = shell.data.38,
#   cores = 4,  # Number of CPU cores to use for parallelization
#   chains = 4,  # Number of MCMC chains
#   iter = 4000  # Number of iterations
# )
# 
# summary(mass_mixed_bayesian_model)




# Obtain marginal effects
# mass_conditional_effects_plot <- conditional_effects(mass_mixed_bayesian_model)
# 
# # Plot the marginal effects
# plot(mass_conditional_effects_plot)
# # spruce up later, but shows that the 95% credible intervals overlap - size is not significantly different at 38 days
# 
# 
# # How much variance does bucket contribute to area?
# 
# icc(mass_mixed_bayesian_model, by_group = T)

# Group               |   ICC
# bucket.no           | 0.928
# bucket.no:treatment | 0.066


#### b. log transform ----


mass_mixed_model_log <- lmer(log10(mass.mg*1000) ~ treatment + (1 | bucket/treatment), data = shell.data.38, REML = T)
summary(mass_mixed_model_log)
qqnorm(residuals(mass_mixed_model_log))
qqline(residuals(mass_mixed_model_log))

plot(residuals(mass_mixed_model_log)) 
abline(h = 0.5, col = "blue", lwd = 4, lty = 4)  
abline(h = -0.5, col = "blue", lwd = 4, lty = 4)  

which(residuals(mass_mixed_model_log) > 0.5)
which(residuals(mass_mixed_model_log) < -0.5)
shell.data.38[c(7,28,37,25,39),]


# Extract fitted values and residuals
fitted_values <- fitted(mass_mixed_model_log)
residuals <- resid(mass_mixed_model_log)

# Plot fitted values vs residuals
plot(fitted_values, residuals, xlab = "Fitted Values", ylab = "Residuals", main = "Fitted Values vs Residuals")

abline(h = 0, col = "red")  # Add a horizontal line at y = 0 for reference
# this shows that variance changes among our fitted values -- log tansform below





##### Levene's test ----
library(car)
library(rcompanion)
# Assuming 'group1' and 'group2' are vectors containing your data
# Run Levene's test
levene_mass <- car::leveneTest(mass.mg ~ bucket, data = shell.data.38)
print(levene_mass)

# pairwise 

lgroup1.1 <- ((shell.data.38 %>% filter(treatment == "control")) %>% filter(bucket == "C1"))
lgroup1.2 <- ((shell.data.38 %>% filter(treatment == "control")) %>% filter(bucket == "C2"))

lgroup2.1 <- ((shell.data.38 %>% filter(treatment == "low-ph")) %>% filter(bucket == "OA1"))
lgroup2.2 <- ((shell.data.38 %>% filter(treatment == "low-ph")) %>% filter(bucket == "OA2"))

# C1 v C2
levene.df.1 <- rbind(lgroup1.1, lgroup1.2)
levene_area1 <- car::leveneTest(mass.mg ~ bucket, data = levene.df.1)
levene_area1

# C1 v LP1
levene.df.2 <- rbind(lgroup1.1, lgroup2.1)
levene_area2 <- car::leveneTest(mass.mg ~ bucket, data = levene.df.2)
levene_area2

# C1 v LP2
levene.df.3 <- rbind(lgroup1.1, lgroup2.2)
levene_area3 <- car::leveneTest(mass.mg ~ bucket, data = levene.df.3)
levene_area3



# C2 v LP1
levene.df.4 <- rbind(lgroup1.2, lgroup2.1)
levene_area4 <- car::leveneTest(mass.mg ~ bucket, data = levene.df.4)
levene_area4

# C2 v LP2
levene.df.5 <- rbind(lgroup1.2, lgroup2.2)
levene_area5 <- car::leveneTest(mass.mg ~ bucket, data = levene.df.5)
levene_area5

# LP1 v LP2
levene.df.6 <- rbind(lgroup2.1, lgroup2.2)
levene_area6 <- car::leveneTest(mass.mg ~ bucket, data = levene.df.6)
levene_area6$`Pr(>F)`[1]


l.pvals <- c(levene_area1$`Pr(>F)`[1], levene_area2$`Pr(>F)`[1],
             levene_area3$`Pr(>F)`[1], levene_area4$`Pr(>F)`[1],
             levene_area5$`Pr(>F)`[1], levene_area6$`Pr(>F)`[1])
l.pvals.bonf <- l.pvals * 6
l.pvals.bonf 



##### F-test ----
#names(shell.data.38)
mgroup1 <- (shell.data.38 %>% filter(treatment == "control"))$mass.mg
mgroup2 <- (shell.data.38 %>% filter(treatment == "low-ph"))$mass.mg

mF_test_result <- var.test(group1, group2)
print(mF_test_result)

# pairwise F test

mgroup1.1 <- ((shell.data.38 %>% filter(treatment == "control")) %>% filter(bucket == "C1"))$mass.mg
mgroup1.2 <- ((shell.data.38 %>% filter(treatment == "control")) %>% filter(bucket == "C2"))$mass.mg

mgroup2.1 <- ((shell.data.38 %>% filter(treatment == "low-ph")) %>% filter(bucket == "OA1"))$mass.mg
mgroup2.2 <- ((shell.data.38 %>% filter(treatment == "low-ph")) %>% filter(bucket == "OA2"))$mass.mg


mpairwise_F_test_result_1 <- var.test(mgroup1.1, mgroup1.2)
print(mpairwise_F_test_result_1)

mpairwise_F_test_result_2 <- var.test(mgroup1.1, mgroup2.1)
print(mpairwise_F_test_result_2)

mpairwise_F_test_result_3 <- var.test(mgroup1.1, mgroup2.2)
print(mpairwise_F_test_result_3)


mpairwise_F_test_result_4 <- var.test(mgroup1.2, mgroup2.1)
print(mpairwise_F_test_result_4)

mpairwise_F_test_result_5 <- var.test(mgroup1.2, mgroup2.2)
print(mpairwise_F_test_result_5)

mpairwise_F_test_result_6 <- var.test(mgroup2.1, mgroup2.2)
print(mpairwise_F_test_result_6)



mf.pvals <- c(mpairwise_F_test_result_1$p.value, mpairwise_F_test_result_2$p.value, 
             mpairwise_F_test_result_3$p.value, mpairwise_F_test_result_4$p.value, 
             mpairwise_F_test_result_5$p.value, mpairwise_F_test_result_6$p.value)
mf.pvals.bonf <- mf.pvals * 6
mf.pvals.bonf 





## C. Shell percent growth ----


# calculate percent area change between 38 control and low ph shells 

#low ph at 38 days
((mean((shell.data.38 %>% filter(treatment == "low-ph"))$area.um2) -
  mean((shell.data.38 %>% filter(treatment == "control"))$area.um2)) / mean((shell.data.38 %>% filter(treatment == "control"))$area.um2)) * 100
# [1] 136.3381
# This means that the low ph shells were 136 percent larger than the control at 38 days


# now at 51 days
shell.data.51 <- shell.data %>% filter(age > 40)

((mean((shell.data.51 %>% filter(treatment == "low-ph"))$area.um2) -
    mean((shell.data.51 %>% filter(treatment == "control"))$area.um2)) / mean((shell.data.51 %>% filter(treatment == "control"))$area.um2)) * 100
# [1] 44.92739
# This means that the low ph shells were 44 percent larger than the control at 51 days



# mass

#low ph at 38 days
((mean((shell.data.38 %>% filter(treatment == "low-ph"))$mass.mg) -
    mean((shell.data.38 %>% filter(treatment == "control"))$mass.mg)) / mean((shell.data.38 %>% filter(treatment == "control"))$mass.mg)) * 100
# [1] 38.32853
# This means that the low ph shells were 38 percent heavier than the control at 38 days


# now at 51 days

((mean((shell.data.51 %>% filter(treatment == "low-ph"))$mass.mg) -
    mean((shell.data.51 %>% filter(treatment == "control"))$mass.mg)) / mean((shell.data.51 %>% filter(treatment == "control"))$mass.mg)) * 100
# [1] 186.351
# This means that the low ph shells were 186 percent larger than the control at 51 days





## D. Main Shell growth statistics ----

### i. NMHP ----

# area.38.ttest <- t.test( (shell.data.38 %>% filter(treatment == "low-ph"))$area.um2 , (shell.data.38 %>% filter(treatment == "control"))$area.um2 )
# area.38.ttest
# t = 4.6314, df = 20.695, p-value = 0.0001487

log10((shell.data.51 %>% filter(treatment == "low-ph"))$area.um2)

area.51.ttest <- t.test( log10((shell.data.51 %>% filter(treatment == "low-ph"))$area.um2) , log10((shell.data.51 %>% filter(treatment == "control"))$area.um2) )
#area.51.ttest <- t.test( (shell.data.51 %>% filter(treatment == "low-ph"))$area.um2 , (shell.data.51 %>% filter(treatment == "control"))$area.um2 )
area.51.ttest
# t = 2.2103, df = 10.925, p-value = 0.04936

# So statistically different, but cannot account for bucket variance, which may be important


# mass.38.ttest <- t.test( (shell.data.38 %>% filter(treatment == "low-ph"))$mass.mg , (shell.data.38 %>% filter(treatment == "control"))$mass.mg )
# mass.38.ttest
# t = 1.3093, df = 36.015, p-value = 0.1987

mass.51.ttest <- t.test( log10((shell.data.51 %>% filter(treatment == "low-ph"))$mass.mg) , log10((shell.data.51 %>% filter(treatment == "control"))$mass.mg) )
#mass.51.ttest <- t.test( (shell.data.51 %>% filter(treatment == "low-ph"))$mass.mg , (shell.data.51 %>% filter(treatment == "control"))$mass.mg )
mass.51.ttest
# t = 3.2996, df = 9.5882, p-value = 0.00848
# significant, and bucket does not appear to influence shell weight based on the mixed model test





# V. Water analysis ----

water.dat.raw <- read.csv("water-data/avg_alk_results-v2.csv")
unique(water.dat.raw$treatment)

keep <- c("control", "oa")

water.dat.filt <- filter(water.dat.raw, treatment %in% keep)
names(water.dat.filt)


t.test(pH.CO2SyS.correct.avg.alk ~ treatment, water.dat.filt)

# t = 10.526, df = 38.621, p-value = 6.611e-13

# Oct 31 2024
# Rather than t test, doing mixed model

# remove rows with NA in column of interest
water.dat.filt.na <- water.dat.filt %>%
  filter(!is.na(pH.CO2SyS.correct.avg.alk))

names(water.dat.filt.na)
ph_mixed_model <- lmer(ph.tris.correct ~ treatment + (1 | sample.ID/treatment), data = water.dat.filt.na, REML = T)

summary(ph_mixed_model)

(8.679e-05 / (8.679e-05 + 4.098e-07 + 1.512e-02))*100

(4.098e-07 / (8.679e-05 + 4.098e-07 + 1.512e-02))*100

(1.512e-02 / (8.679e-05 + 4.098e-07 + 1.512e-02))*100





