# Load Gencode-v32: for genome features.
gencode_gr <- rtracklayer::import("/scratch/Shares/rinnclass/CLASS_2023/data/data/genomes/gencode.v32.annotation.gtf")

# let's add 1Kb upstream and downstream from the TSS to define "promoters"
gencode_promoters <- promoters(gencode_gr[gencode_gr$type == "gene"], 
                               upstream = 1e3, 
                               downstream = 1e3)

kurt_basepath <- "/scratch/Shares/rinnclass/CLASS_2023/kurt"
kurt_path <- "CLASS_2023/group_chip/kurt/results/bwa/mergedLibrary/macs/broadPeak"

# first we read the peak files in as gRanges object with rtracklayer function.
peaks1 <- rtracklayer::import(file.path(kurt_basepath,kurt_path, "BRCA1_R1_peaks.broadPeak"))
peaks2 <- rtracklayer::import(file.path(kurt_basepath,kurt_path, "BRCA1_R2_peaks.broadPeak"))
peaks3 <- rtracklayer::import(file.path(kurt_basepath,kurt_path, "BRCA1_R3_peaks.broadPeak"))
peaks4 <- rtracklayer::import(file.path(kurt_basepath,kurt_path, "BRCA1_R4_peaks.broadPeak"))

# get length of each peak file
lp1 <- length(peaks1)
lp2 <- length(peaks2)
lp3 <- length(peaks3)
lp4 <- length(peaks4)

# get overlaps for each combination of peak files
ovf12 <- findOverlaps(peaks1, peaks2)
ovf13 <- findOverlaps(peaks1, peaks3)
ovf14 <- findOverlaps(peaks1, peaks4)
ovf23 <- findOverlaps(peaks2, peaks3)
ovf24 <- findOverlaps(peaks2, peaks4)
ovf34 <- findOverlaps(peaks3, peaks4)

# get percent overlap for each peak file
p12_1 <- length(ovf12)/lp1
p12_2 <- length(ovf12)/lp2
p13_1 <- length(ovf13)/lp1
p13_3 <- length(ovf13)/lp3
p14_1 <- length(ovf14)/lp1
p14_4 <- length(ovf14)/lp4
p23_2 <- length(ovf23)/lp2
p23_3 <- length(ovf23)/lp3
p24_2 <- length(ovf24)/lp2
p24_4 <- length(ovf24)/lp4
p34_3 <- length(ovf34)/lp3
p34_4 <- length(ovf34)/lp4

pavg = mean(c(p12_1,p12_2,
              p13_1,p13_3,
              p14_1,p14_4,
              p23_2,p23_3,
              p24_2,p24_4,
              p34_3,p34_4))

min_p12 <- length(ovf12)/min(lp1,lp2)
min_p13 <- length(ovf13)/min(lp1,lp3)
min_p14 <- length(ovf14)/min(lp1,lp4)
min_p23 <- length(ovf23)/min(lp2,lp3)
min_p24 <- length(ovf24)/min(lp2,lp4)
min_p34 <- length(ovf34)/min(lp3,lp4)

min_pavg <- mean(c(min_p12,min_p13,min_p14,
                   min_p23,min_p24,min_p34))

# get promoter overlaps for each peak
po_peaks1 <- subsetByOverlaps(peaks1, gencode_promoters) #promoter overlaps
po_peaks2 <- subsetByOverlaps(peaks2, gencode_promoters)
po_peaks3 <- subsetByOverlaps(peaks3, gencode_promoters)
po_peaks4 <- subsetByOverlaps(peaks4, gencode_promoters)

# now overlapping the two overlaps by findOverlaps
spo_peaks12 <- findOverlaps(po_peaks1, po_peaks2) #subset promoter overlaps
spo_peaks34 <- findOverlaps(po_peaks3, po_peaks4)
min_p_spo_peaks12 <- length(spo_peaks12)/min(lp1,lp2)
min_p_spo_peaks34 <- length(spo_peaks34)/min(lp3,lp4)

# Consensus peaks of all
con_peaks_12 <- subsetByOverlaps(peaks1, peaks2)
con_peaks_34 <- subsetByOverlaps(peaks3, peaks4)
con_peaks <- subsetByOverlaps(con_peaks_12, con_peaks_34)
spo_con_peaks <- subsetByOverlaps(con_peaks, gencode_promoters)
# percent of all consensus peaks overlapping promoters
p_spo_con_peaks <- length(spo_con_peaks)/length(con_peaks) #0.85


#> length(subsetByOverlaps(con_peaks_12, gencode_promoters))/length(con_peaks_12)
#[1] 0.8342629
#> length(subsetByOverlaps(con_peaks_34, gencode_promoters))/length(con_peaks_34)
#[1] 0.4140567

### Insights
# BRCA1 binds promoters 85% of the time
# Replicates 1 and 2 have way fewer peaks
# Reps 1/2 show more promoter binding than reps 3/4 (0.83 vs 0.41)
# Reps 3/4 noisy?


# get percent promoter overlap of higher peak scores
quantile(peaks4$score, 0.99) # find number
peaks4_90th <- peaks4[peaks4$score > 89] # use number from above
po_peaks4_90th <- subsetByOverlaps(peaks4_90th, gencode_promoters)
po_peaks4_90th_percent = length(po_peaks4_90th)/length(peaks4_90th)
po_peaks4_90th_percent

