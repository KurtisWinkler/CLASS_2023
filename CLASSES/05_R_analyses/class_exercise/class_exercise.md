Class\_exercise
================
Kurt Winkler
3/16/2023

# Load the libraries you need

# Load functions you need “my\_class\_functions”

# load in your peak files for each replicate of each protein

# Here I am starting to analyze my data for my proteins of interest:

# proteinX, Y, Z …..

# First I will read in each replicate file

``` r
# filepath to import peaks
basepath <- "/scratch/Shares/rinnclass/CLASS_2023/kurt"
peak_path <- "CLASS_2023/group_chip/kurt/results/bwa/mergedLibrary/macs/broadPeak"
broadpeakfilepath <- file.path(basepath, peak_path)

# import peaks
peak_list <- import_peaks(consensus_file_path = broadpeakfilepath)

# printing out a table of the number of peaks in each file:
peak_num <- sapply(peak_list, length) %>% as.data.frame()
names(peak_num) <- c("num_peaks")
peak_num
```

    ##            num_peaks
    ## ARID3A_R1      47962
    ## ARID3A_R2      30842
    ## ATF3_R1         3848
    ## ATF3_R2        35018
    ## ATF3_R3        49209
    ## ATF3_R4        60518
    ## BHLHE40_R1     18183
    ## BHLHE40_R2      2618
    ## BHLHE40_R3     10119
    ## BHLHE40_R4      8074
    ## BRCA1_R1        2212
    ## BRCA1_R2        4073
    ## BRCA1_R3       44978
    ## BRCA1_R4       44173
    ## CEBPB_R1       47888
    ## CEBPB_R2       26625
    ## CEBPB_R3       53316
    ## CEBPB_R4      187109

# Now I am going to create consensus peaks for each protein

``` r
dbps <- unique(sapply(names(peak_list), function(x) {
   unlist(strsplit(x, "_"))[1]
}))
# Let's run it !
consensus_list <- lapply(dbps, consensus_from_reduced, peak_list)
names(consensus_list) <- dbps

# export consensus peaks to results folder
if (dir.exists("results") == FALSE) {
  dir.create("results") }

for(i in 1:length(consensus_list)) {
rtracklayer::export(consensus_list[[i]], paste0("results/", names(consensus_list)[i],"_consensus_peaks.bed")) }
```

# Now I am going to make my consensus peaks compatable with UCSC genome browser

``` r
# paths for UCSC output
consensus_peak_path <- "CLASS_2023/CLASSES/05_R_analyses/class_exercise/results"
consensusFilePath <- file.path(basepath, consensus_peak_path)
export_path <- file.path(consensusFilePath, "UCSC")

# print out consensus peak files in a results/UCSC directory
if (dir.exists("results/UCSC") == FALSE) {
  dir.create("results/UCSC") }

ucsc_formating(consensusFilePath = consensusFilePath, export_path = export_path)
```

    ## [1] 5

# I am curious if my proteins are transcription factors so I will use the annotations

# in a cell paper I found and see

``` r
annot_path <- file.path(basepath, "CLASS_2023/CLASSES/05_R_analyses/01_peak_features/results/peak_features.RData")
load(annot_path, verbose = T)
```

    ## Loading objects:
    ##   consensus_peaks
    ##   gencode_genes
    ##   lncrna_gene_ids
    ##   mrna_gene_ids
    ##   num_peaks_df
    ##   peak_occurence_df
    ##   lncrna_mrna_promoters
    ##   mrna_lncrna_genes

``` r
options(warn=-1)
human_tfs <- readxl::read_excel("results/TF_annotations.xlsx",
                                sheet = 2, skip = 1)
```

    ## New names:
    ## • `` -> `...4`

``` r
# let's rename the 4th column to indicate if it is a TF.
names(human_tfs)[4] <- "is_tf"

# now let's intersect gene names that are in our ChIP data and has TF identity.
length(which(tolower(num_peaks_df$dbp) %in% tolower(human_tfs$Name)))
```

    ## [1] 5

``` r
# first let's filter and grab the first 4 columns that match DBPs in num_peaks_df
human_tfs <- human_tfs[tolower(human_tfs$Name) %in% tolower(dbps), 1:4]
# if you leave the object name you just created in the environment
# it will print out in the knit. For example :
human_tfs
```

    ## # A tibble: 5 × 4
    ##   ID              Name    DBD         is_tf
    ##   <chr>           <chr>   <chr>       <chr>
    ## 1 ENSG00000116017 ARID3A  ARID/BRIGHT Yes  
    ## 2 ENSG00000134107 BHLHE40 bHLH        Yes  
    ## 3 ENSG00000162772 ATF3    bZIP        Yes  
    ## 4 ENSG00000172216 CEBPB   bZIP        Yes  
    ## 5 ENSG00000012048 BRCA1   Unknown     No

# Now I want to compare a protein with a previous analysis

``` r
# go to UCSC genome browser and load in a peak file for a given protein
# load in the data for the same protein from the previous analysis
# compare how your consensus peaks are similar or different to previous analyses

#Consensus peaks are NOT similar to previous analysis
knitr::include_graphics("results/UCSC/BRCA1_compare_2.png")
```

<img src="results/UCSC/BRCA1_compare_2.png" width="1713" />

# Now I am going to determine how my peaks for each protein overlap annotations of the genome

# First I will find the overlaps between my consensus peaks with promoters of lncRNA and mRNA promoters

``` r
cat("number of consensus peaks\n")
```

    ## number of consensus peaks

``` r
sapply(consensus_list, length)
```

    ##  ARID3A    ATF3 BHLHE40   BRCA1   CEBPB 
    ##   23669    2876    1234    1160   12519

``` r
# counting promoter overlaps for mrna and lncrna
promoter_peak_counts <- count_peaks_per_feature(lncrna_mrna_promoters, consensus_peaks, type = "counts")

# find overlaps of promoters for each protein
cat("\nlncrna and mRNA promoter overlaps\n")
```

    ## 
    ## lncrna and mRNA promoter overlaps

``` r
rowSums(promoter_peak_counts)
```

    ##  ARID3A    ATF3 BHLHE40   BRCA1   CEBPB 
    ##    4601    2659     628    1481    1742

## results:

\#1) What can you determine from these overlaps? From these overlaps, I
can determine that a large proportion of ATF3 and BRCA1 consensus peaks
overlap promoters. About half of BHLHE40 consensus peaks overlap
promoters and few ARID3A and CEBPB consensus peaks overlap promoters.

# Now I want to compare the overlaps with lncRNA and mRNA promoters seperately

``` r
# Now let's break these promoters into two groups "lncrna" and "mrna"
# lncrna promoter overlaps
cat("lncrna promoter overlaps\n")
```

    ## lncrna promoter overlaps

``` r
rowSums(promoter_peak_counts[,lncrna_gene_ids])
```

    ##  ARID3A    ATF3 BHLHE40   BRCA1   CEBPB 
    ##    1259     532     159     286     529

``` r
cat("\nmRNA promoter overlaps\n")
```

    ## 
    ## mRNA promoter overlaps

``` r
# mrna promoter overlaps
rowSums(promoter_peak_counts[,mrna_gene_ids])
```

    ##  ARID3A    ATF3 BHLHE40   BRCA1   CEBPB 
    ##    3342    2127     469    1195    1213

## results:

# 1) What is the difference in overlaps between mRNA and lncRNA promoters

For every protein, there are many more overlaps with mRNA compared to
lncRNA, around 3-4x more in mRNA than lncRNA

# Now I am going to test if there is more binding over gene bodies than promoters

# I will seperate lncRNA and mRNA gene bodies to find the overlaps

``` r
genebody_peak_counts <- count_peaks_per_feature(mrna_lncrna_genes, 
                                                consensus_peaks, 
                                                type = "counts")

# All gene bodies
cat("lncrna and mRNA genebody overlaps\n")
```

    ## lncrna and mRNA genebody overlaps

``` r
rowSums(genebody_peak_counts)
```

    ##  ARID3A    ATF3 BHLHE40   BRCA1   CEBPB 
    ##   19847    3719    1162    1731   10065

``` r
# lncRNA gene bodies 
cat("\nlncrna genebody overlaps\n")
```

    ## 
    ## lncrna genebody overlaps

``` r
rowSums(genebody_peak_counts[,lncrna_gene_ids])
```

    ##  ARID3A    ATF3 BHLHE40   BRCA1   CEBPB 
    ##    4527     722     233     306    2291

``` r
# mRNA gene bodies
cat("\nmRNA genebody overlaps\n")
```

    ## 
    ## mRNA genebody overlaps

``` r
rowSums(genebody_peak_counts[,mrna_gene_ids])
```

    ##  ARID3A    ATF3 BHLHE40   BRCA1   CEBPB 
    ##   15320    2997     929    1425    7774

## results:

# 1) Do my proteins have more overlaps with promoters or genebodies?

All of my proteins have more overlaps with genebodies, especially mRNA
genebodies

# It is nice and all to find overlaps, but I am interested in how many proteins

# bind a specific promoter. I will use my handy “occurence” parameter in

# " count peaks per feature"

``` r
promoter_peak_occurence <- count_peaks_per_feature(lncrna_mrna_promoters, consensus_peaks, 
                                               type = "occurrence")

max(colSums(promoter_peak_occurence))
```

    ## [1] 5

## results: I find the max number of proteins on a promoter to be X

The maximum number of protein on a promoter is 5, meaning that all my
proteins bind to a specific promoter

# Now I want to start plotting my results

# First I will see if there is a realtionship between peak number and total DNA covered

``` r
ggplot(num_peaks_df, aes(x = num_peaks, y = total_peak_length)) +
  geom_point() + 

  ylab("BP covered") +
  xlab("Number of peaks") +
  ggtitle("Peak count vs. total bases covered")
```

![](class_exercise_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

# Now I want to color my plot by wether the protein is a TF or not.

``` r
ggplot(num_peaks_df, aes(x = num_peaks, y = total_peak_length,
                         color = tf == "Yes")) +
  geom_point() + 

  ylab("BP covered") +
  xlab("Number of peaks") +
  ggtitle("Peak count vs. total bases covered")
```

![](class_exercise_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

# I want to make a histogram of the number of peaks for each of my proteins

``` r
ggplot(num_peaks_df, aes(x = num_peaks)) + 
  geom_histogram(bins = 20)
```

![](class_exercise_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

# Now I want to facet this by the type of DNA binding domain my protein has.

``` r
ggplot(num_peaks_df, aes(x = num_peaks)) + 
  facet_wrap(dbd ~ .) +
  geom_histogram(bins = 20)
```

![](class_exercise_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

# Cool now I am ready to send my result to my collaborator as a

# Knitted document
