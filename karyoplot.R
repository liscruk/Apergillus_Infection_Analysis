library(tidyverse)
library(biomaRt)
library(GenomicRanges)
library(rtracklayer)
library(Rsamtools)
library(karyoploteR)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

setwd("/Volumes/Work_Hot2/Johannes_ChIP")
resout = "/Volumes/Work_Hot2/Johannes_ChIP/analysis/"

# Gene list
clist = c("IRF1", "GBP2", "GBP5", "RSAD2", 
          "OASL", "MX1", "MX2", "USP18",
          "CXCL10", "TNFSF10","IL1B","IL12B", 
          "IL18", "IL10", "IL4","IFI16",
          "IFIT3", "IFIT2", "IFI44L",
          "ISG15", "IDO1", "RIGI")

# Causes an error
# IRGB10 (mouse gene) substituted for IFI16
# DDX5 not correct gene symbol substituted for RIGI
# IL1b not correct gene symbol substituted for IL1B


# Load peaks
infected_peaks = import("data/IRF1_Chip_Infected/analysis/macs3/IRF1_Chip_Infected_peaks_mod_clean.narrowPeak")
strand(infected_peaks) = "*"
uninfected_peaks = import("data/IRF1_Chip_Uninfectd/analysis/macs3/IRF1_Chip_Uninfected_peaks_mod_clean.narrowPeak")
strand(uninfected_peaks) = "*"

# Load BAM files
bamInfected = "data/IRF1_Chip_Infected/data/bam/ChIP_06/ChIP_06.dedup.bam"
# bamUninfected = "data/IRF1_Chip_Uninfectd/data/bam/ChIP_03/ChIP_03.dedup.bam"
bamIgg = "data/IRF1_Chip_Infected/data/bam/ChIP_05/ChIP_05.dedup.bam"

# Map gene symbols to Entrez Gene IDs
entrez_ids = mapIds(org.Hs.eg.db, 
                     keys = clist,
                     column = "ENTREZID",
                     keytype = "SYMBOL",
                     multiVals = "first")

txdb = TxDb.Hsapiens.UCSC.hg38.knownGene

### Get genes of interest
for (i in seq_along(clist)) {

  gr = genes(TxDb.Hsapiens.UCSC.hg38.knownGene)[entrez_ids[[i]]]
  start(gr) = start(gr)-10000
  end(gr) = end(gr)+10000
  
  ### Infected
  ### Compute infected max to set scaling corretly
  x = plotKaryotype(zoom = gr, genome = "hg38")
  f = kpPlotBAMCoverage(x, data = bamInfected)
  s = kpPlotBAMCoverage(x, data = bamIgg)
  
  u = max(c(f$latest.plot$computed.values$max.coverage,
            s$latest.plot$computed.values$max.coverage))
  
  dev.off()
  
  name = paste0(clist[i],"_karyoplot.png")
  print(name)
  
  png(filename = paste0(resout,name),
      width = 1000,
      height = 600)
  
  ### Set up actual plot
  kp = plotKaryotype(zoom = gr, genome = "hg38", cex = 1, plot.type = 2)
  kpAddMainTitle(kp,clist[i])
  
  kp = kpPlotBAMCoverage(kp, data = bamInfected, ymax = u, col = "red", r0 = 0.6)
  kpAxis(kp, ymax = u, r0 = 0.6)
  kpAddLabels(kp, "Infected", r0 = 1.3, r1 = 0.5, label.margin = 0.05, srt = 90)
  
  kp = kpPlotBAMCoverage(kp, data = bamIgg, col = "blue", r1 = 0.4)
  kpAxis(kp, ymax = u, r1 = 0.4)
  
  kpAddLabels(kp, "IgG Control", r0 = 0.3, r1 = 0.5, label.margin = 0.05, srt = 90)
  
  genes.data = makeGenesDataFromTxDb(TxDb.Hsapiens.UCSC.hg38.knownGene, karyoplot = kp,
                                     plot.transcripts = TRUE, plot.transcripts.structure = TRUE)
  
  genes.data = mergeTranscripts(genes.data)
  
  kpPlotGenes(kp,add.gene.names = FALSE, data = genes.data, r0 = 0.1, r1 = 0.5, data.panel = 2,add.strand.marks = TRUE)
  
  kpPlotRegions(kp,infected_peaks,r0=0.5,r1=0.6, col = "grey",data.panel = 2,border = "black")
  kpAddLabels(kp,"Sig. Peaks \n Infected", r0=0.4,r1=0.7, data.panel = 2)
  kpAddBaseNumbers(kp, tick.dist = 15000, add.units = TRUE)
  
  dev.off()
  
  
}



