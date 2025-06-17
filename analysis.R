library(edgeR);library(tximport);library(readr);
library(GenomicFeatures);library(EnsDb.Hsapiens.v86);
library(org.Hs.eg.db);library(clusterProfiler);
library(EnhancedVolcano);library(gridExtra);
library(DOSE);library(RColorBrewer);library(msigdbr)
library(viridis);library(ggfortify);library(clusterProfiler);
library(RColorBrewer);library(ComplexHeatmap);library(circlize)
library(GeneOverlap);library(jsonlite)


### Latest Heiko list
clist = c("IRF1","GBP2","GBP5","RSAD2","OASL","MX1","MX2","USP18",
          "IFIT3","IFIT2","IFI44L","IRGB10","ISG15","DDX58","IDO1",
          "CXCL10","TNFSF10","IL12B","IL18","IL1b","IL10","IL4")

set.seed(1)

#### Some functions ####

coolBlueHotRed = function(n, alpha = 1) {
  rainbow(n, end=4/6, alpha=alpha)[n:1]
}
pheatmap_colors = function(n = 100) {
  colorRampPalette(c("navy", "white", "firebrick3"))(n)
}

symmetric_breaks = function(mat, n = 100) {
  max_val = max(abs(mat), na.rm = TRUE)
  seq(min(mat), max_val, length.out = n)
}


get_gene_lengths = function(gene_ids) {
  # Load biomaRt if not already loaded
  library(biomaRt)
  
  # Connect to Ensembl
  ensembl = useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  
  # Query gene coordinates
  gene_info = getBM(
    attributes = c("ensembl_gene_id", "start_position", "end_position"),
    filters = "ensembl_gene_id",
    values = gene_ids,
    mart = ensembl
  )
  
  # Calculate gene length
  gene_info$gene_length = gene_info$end_position - gene_info$start_position + 1
  
  return(gene_info)
}

### Load and Process Data ####
### Define ins and out
datadir = "~/../../Volumes/Work_Cold/Data/AG_Bruns_Novogen/Aspergillus/johannes_aspergillus/analysis/counts/"
resdir  = "Aspergilluis/analysis/"

### Read ssamples
samples = read.table("~/../../Volumes/Work_Cold/Data/AG_Bruns_Novogen/Aspergillus/Samples_RNAseq_JB_asp.csv", header = T, sep = ";")
files = file.path(datadir,samples$Sample,"quant.sf")
all(file.exists(files))
names(files) = samples$Sample

### read annotation and convert to txiMeta
txdb = EnsDb.Hsapiens.v86
k = keys(txdb, keytype = "TXNAME")
tx2gene = AnnotationDbi::select(txdb,k,"GENEID","TXNAME")

k = keys(txdb,keytype = "GENEID")

gene2sym = AnnotationDbi::select(txdb,k,keytype = "GENEID",columns = c("GENENAME","ENTREZID"))
txi = tximport(files, type = "salmon",tx2gene = tx2gene,ignoreTxVersion = T)
gene2sym = gene2sym[match(rownames(txi$counts),gene2sym$GENEID),]

samples$Treatment_Time = paste(samples$Treatment,samples$Time,sep="_")

### DGE
dge = DGEList(counts = txi$counts,
              group = samples$Treatment_time,
              genes = gene2sym)

Donor = factor(samples$Donor)
Time = factor(samples$Time)
Treatment = factor(samples$Treatment)

keep = filterByExpr(dge,group = samples$Treatment)
dge = dge[keep,]

### Model Setup  ####
design = model.matrix(~0+Treatment_Time,data = samples)
keep = filterByExpr(dge,group = samples$Treatment)
dge = dge[keep,]
dge = calcNormFactors(dge)
dge = normLibSizes(dge)
dge = estimateDisp(dge,design,robust = T)
fit = glmQLFit(dge,design)

dge$genes$geneid = rownames(dge$counts)
m = match(dge$genes$geneid,gene2sym$GENEID)
dge$genes$Symbol = list(gene2sym$GENENAME[m])

colnames(design)

### Contrasts ####
con1 = makeContrasts(Treatment_Timegermlings_6h-Treatment_Timeuntreated_6h,levels = colnames(design))
con2 = makeContrasts(Treatment_Timegermlings_3h-Treatment_Timeuntreated_3h,levels = colnames(design))
con3 = makeContrasts(
  ((Treatment_Timegermlings_3h+Treatment_Timegermlings_6h)/2)-
  ((Treatment_Timeuntreated_3h+Treatment_Timeuntreated_6h)/2),
levels = colnames(design))

x1 = glmQLFTest(fit,contrast = con1)
x2 = glmQLFTest(fit,contrast = con2)
x3 = glmQLFTest(fit,contrast = con3)
  
### Contrast 1 - 6h ####
x = x1
res = topTags(x,n = Inf)
res = as.data.frame(res)
slist = as.data.frame(topTags(x,n=10)[,2])

condition_name = "Treated_versus_Control_(6h)"

vp = EnhancedVolcano(res,x = "logFC",y = "FDR",
                     title = condition_name,
                     subtitle = "",
                     lab = res$GENENAME,
                     selectLab = clist,
                     pCutoff  = 0.05,
                     FCcutoff = 1.0,
                     boxedLabels = TRUE,
                     parseLabels = TRUE,
                     colAlpha = 4/5,
                     drawConnectors = TRUE,
                     legendPosition = 'bottom',
                     max.overlaps = 20)


png(filename = paste0(paste0("Aspergilluis/analysis//",condition_name),"_volcano_newList.png"),
    width = 600,height = 600,pointsize = 12)
plot(vp)

dev.off()

counts = cpm(x, prior.count = 1, log = TRUE)
profile = res[abs(res$logFC)>=0.5 & res$FDR<0.05,]
full_res = cbind(profile[,c("GENENAME","logFC","PValue","FDR")],counts[rownames(profile),])

write.table(full_res,paste0(paste0("Aspergilluis/analysis/",condition_name),".csv"),
            sep = ",",
            row.names = F,
            col.names = T)

gene_list = res[,"logFC"]
names(gene_list) = res$GENEID
gene_list = na.omit(gene_list)
gene_list = gene_list[unique(names(gene_list))]
gene_list = sort(gene_list,decreasing = T)

#### Heatmap - Split in two for abundances ####
tpm = as.data.frame(txi$abundance)
tpm = tpm[rownames(res),]
tpm$GENENAMES = res[rownames(tpm),"GENENAME"]
tpm = tpm[tpm$GENENAMES %in% clist,]
rownames(tpm) = tpm$GENENAMES
tpm = subset(tpm,select = -GENENAMES)

#### Setup annotation
ano = samples[samples$Treatment_Time == "germlings_6h" | samples$Treatment_Time == "untreated_6h",]
rownames(ano) = ano$Sample
ch = tpm[,rownames(ano)]
colnames(ch) = rownames(ano)
ch = log2(ch+1)

d = RColorBrewer::brewer.pal(5,"Set1")
names(d) = unique(ano$Donor)

col_ano = ComplexHeatmap::HeatmapAnnotation(group = ano$Treatment,donor = ano$Donor,
                                            col = list(group = c( "untreated" = "orange",
                                                                  "germlings" = "darkgreen"),
                                                       donor = d))

ano_row = full_res[full_res$GENENAME %in% clist,]
rownames(ano_row) = ano_row$GENENAME

col_fun = circlize::colorRamp2(c(-2, 3.5, 6), c("blue", "white", "red"))

row_ano = ComplexHeatmap::rowAnnotation(logFC = ano_row$logFC,
                                        col = list(logFC = col_fun))

col_fun3 = colorRamp2(
  c(min(ch),0,1.5,3),
  c("#1B6AB9", "#F7F7F7", "orange","red")
)

#### Setup Heatmap
png(filename = paste0(condition_name, "_heatmap_genes_of_interest.png"),
    width = 800,height = 600,pointsize = 12)

ComplexHeatmap::Heatmap(as.matrix(ch), 
                        name = "log2(TPM+1)",
                        border = T,
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          grid.rect(x = x, y = y, width = width, height = height,
                                    gp = gpar(col = "grey", fill = NA, lwd = 0.5))
                        },
                        column_title = paste0(condition_name," - Genes of Interest"),
                        top_annotation = col_ano,
                        # left_annotation = row_ano,
                        show_column_names = FALSE,
                        show_column_dend = TRUE)

dev.off()


#### SETUP NEW GSEA Analysis with IRF1 GeneSet
# https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/IRF1_01.html
gene_list = res
gene_list = gene_list[gene_list$FDR<0.05,]
gene_list = gene_list[,c("GENEID","logFC")]
n = gene_list$GENEID
gene_list = gene_list[,"logFC"]
names(gene_list) = n
gene_list = na.omit(gene_list)
gene_list = gene_list[unique(names(gene_list))]
gene_list = sort(gene_list,decreasing = T)

mt2irf = msigdbr(species = "Homo sapiens", 
                collection = "C3",
                subcollection = "TFT:TFT_LEGACY") %>% dplyr::select(gs_name,ensembl_gene)


h = clusterProfiler::GSEA(gene_list , 
                          TERM2GENE  = mt2irf, 
                          scoreType = "pos", 
                          pAdjustMethod = "fdr",
                          pvalueCutoff = 0.5)

h = h@result
h$Description = gsub("_.*$","",h$Description)

count_ensembl_ids = function(strings) {
  # strings: a character vector, each element is a string of IDs separated by "/"
  sapply(strings, function(x) {
    ids = unlist(strsplit(x, split = "/"))
    length(ids)
  })
}


f = count_ensembl_ids(h$core_enrichment)
h$count = f

png(filename = paste0(condition_name, "_NES_IRF1.png"),
    width = 800,height = 600,pointsize = 16)

ggplot(h[1:20,],aes(reorder(Description,NES),NES))+
  geom_point(aes(size = count ,fill = p.adjust),
             shape = 21,
             colour = "black",
             stroke = 0.6)+
  scale_fill_gradient(low = "#a2bffe", high="#FF746C")+
  coord_flip()+
  theme_bw()+
  theme(axis.text = element_text(color = "black", size = 16))
  
dev.off()


###
library(httr);library(jsonlite)
irf1_content = GET("https://maayanlab.cloud/Harmonizome/api/1.0/gene_set/IRF1/ChEA+Transcription+Factor+Targets")

### Get raw list of targets
### IRF1
irf1_content = content(irf1_content,as = "text", encoding = "UTF-8")
irf1_targets = fromJSON(irf1_content)
irf1_targets = irf1_targets$associations$gene[1]
irf1_targets = irf1_targets$symbol

### IRF7
### IRF7 targets
irf7_content = read_tsv("Aspergilluis/data//IRF7_01.v2025.1.Hs.tsv")
irf7_content = unlist (irf7_content[17,])
irf7_content = unlist(strsplit(irf7_content[2],","))

### SP1
### SP1 targets
sp1_content = read_tsv("Aspergilluis/data/SP1_01.v2025.1.Hs.tsv")
sp1_content = unlist(sp1_content[17,])
sp1_content = unlist(strsplit(sp1_content[2],","))

### TBP
### TBP Targets
tbp_content = GET("https://maayanlab.cloud/Harmonizome/api/1.0/gene_set/TBP/ENCODE+Transcription+Factor+Targets")
tbp_targets = content(tbp_content, as = "text", encoding = "UTF-8")
tbp_targets = fromJSON(tbp_targets)

### TBP has more targets then is expressed in data. Keep < Targets



tbp_targets = tbp_targets$associations$gene[1]
tbp_targets = tbp_targets$symbol

### IRF1 Enrichment
fcs = seq(1,7, by = 0.1)
odds = c()
ins = c()
pval = c()
lfc = c()
j = 1

for (i in fcs) {
  
  a = res[with(res,logFC>i & FDR<0.1),]
  b = irf1_targets
  a = a$GENENAME
  
  
  
  f = newGeneOverlap(a,b, genome.size = length(keep))
  f = testGeneOverlap(f)
  
  print(i)
  
  lfc[j] = i
  odds[j] = f@odds.ratio
  ins[j] = length(f@intersection)
  pval[j] = f@pval
  
  j = j+1
  
}

irf1 = data.frame(lfc,odds,ins,pval)
irf1$fdr = p.adjust(irf1$pval,method = "fdr")
irf1$sig = ifelse(irf1$pval<0.05,1,0)

ggplot(irf1,aes(-log10(pval),odds))+
  geom_point(aes(size = ins,color = lfc))+
  geom_vline(xintercept = -log10(0.05),linetype = "longdash", color = "grey")+
  geom_smooth()+
  scale_colour_gradient(low = "blue", high="red")+
  ggtitle("IRF1")+
  theme(axis.text = element_text(color = "black", size = 16))+
  theme_bw()


### IRF7 Enrichment
fcs = seq(1,7,by=0.1)
odds = c()
ins = c()
pval = c()
lfc = c()
j = 1

for (i in fcs) {
  
  a = res[with(res,logFC>i & FDR<0.1),]
  b = irf7_content
  a = a$GENENAME
  
  
  
  f = newGeneOverlap(a,b, genome.size = length(keep))
  f = testGeneOverlap(f)
  
  print(i)
  
  lfc[j] = i
  odds[j] = f@odds.ratio
  ins[j] = length(f@intersection)
  pval[j] = f@pval
  
  j = j+1
  
}

irf7 = data.frame(lfc,odds,ins,pval)
irf7$fdr = p.adjust(irf7$pval,method = "fdr")
irf7$sig = ifelse(irf7$pval<0.05,1,0)

ggplot(irf7,aes(-log10(pval),odds))+
  geom_point(aes(size = ins,color = lfc))+
  geom_vline(xintercept = -log10(0.05),linetype = "longdash", color = "grey")+
  geom_smooth()+
  scale_colour_gradient(low = "blue", high="red")+
  ggtitle("IRF7")+
  theme(axis.text = element_text(color = "black", size = 16))+
  theme_bw()


### SP1 Enrichment
fcs = seq(1, 7, by=0.1)
odds = c()
ins = c()
pval = c()
lfc = c()
j = 1

### SP1 Enrichment
for (i in fcs) {
  
  a = res[with(res,logFC>i & FDR<0.1),]
  b = sp1_content
  a = a$GENENAME
  
  
  
  f = newGeneOverlap(a,b, genome.size = length(keep))
  f = testGeneOverlap(f)
  
  print(i)
  
  lfc[j] = i
  odds[j] = f@odds.ratio
  ins[j] = length(f@intersection)
  pval[j] = f@pval
  
  j = j+1
  
}

sp1 = data.frame(lfc,odds,ins,pval)
sp1$fdr = p.adjust(sp1$pval,method = "fdr")
sp1$sig = ifelse(sp1$pval<0.05,1,0)

ggplot(sp1,aes(-log10(pval),odds))+
  geom_point(aes(size = ins,color = lfc))+
  geom_vline(xintercept = -log10(0.05),linetype = "longdash", color = "grey")+
  geom_smooth()+
  scale_colour_gradient(low = "blue", high="red")+
  ggtitle("SP1")+
  theme(axis.text = element_text(color = "black", size = 16))+
  theme_bw()



irf1$tf = "IRF1"
irf7$tf = "IRF7"
sp1$tf = "SP1"

tfres = rbind(irf1,irf7)
tfres = rbind(tfres,sp1)
tfres = tfres[tfres$ins != 0,]

### Smooth 
ggplot(tfres,aes(-log10(pval),odds))+
  geom_point(aes(size = ins, color = tf))+
  geom_vline(xintercept = -log10(0.05),linetype = "longdash")+
  geom_smooth(aes(group = tf, linetype = tf, colour = tf), se = FALSE)+
  ggtitle("All Targets")+
  theme(axis.text = element_text(color = "black", size = 16))+
  facet_wrap(~tf)+
  theme_bw()

### TBP Enrichment
fcs = seq(1,7, by = 0.1)
odds = c()
ins = c()
pval = c()
lfc = c()
j = 1

for (i in fcs) {
  
  a = res[with(res,logFC>i & FDR<0.1),]
  b = tbp_targets
  a = a$GENENAME
  
  
  
  
  f = newGeneOverlap(a,b, genome.size = length(keep))
  f = testGeneOverlap(f)
  
  print(i)
  
  lfc[j] = i
  odds[j] = f@odds.ratio
  ins[j] = length(f@intersection)
  pval[j] = f@pval
  
  j = j+1
  
}

irf1 = data.frame(lfc,odds,ins,pval)
irf1$fdr = p.adjust(irf1$pval,method = "fdr")
irf1$sig = ifelse(irf1$pval<0.05,1,0)

ggplot(irf1,aes(-log10(pval),odds))+
  geom_point(aes(size = ins,color = lfc))+
  geom_vline(xintercept = -log10(0.05),linetype = "longdash", color = "grey")+
  geom_smooth()+
  scale_colour_gradient(low = "blue", high="red")+
  ggtitle("IRF1")+
  theme(axis.text = element_text(color = "black", size = 16))+
  theme_bw()


 ### Density 
ggplot(tfres,aes(ins,odds, group = tf))+
  geom_point()+
  facet_wrap(~tf)+
  theme_bw()



### GSEA DO / KEGG First contrast ####
mt2gh = msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name,ensembl_gene)

# Only use highly changed stuff
gene_list = gene_list[gene_list>1]

h = enricher(names(gene_list), TERM2GENE  = mt2gh)

png(filename = paste0(paste0("analysis/",condition_name),"_hallmark_pathways_upregulated.png"),
    width = 600,height = 600,pointsize = 12)
dotplot(h, showCategory = 20, title = "Enriched Up-regulated Hallmark Pathways 6h")+scale_fill_viridis()
dev.off()

do_list = res[,"logFC"]
names(do_list) = res$ENTREZID

do_list = na.omit(do_list)
do_list = do_list[unique(names(do_list))]
do_list = sort(do_list,decreasing = T)

gse =  gseKEGG(do_list)

gse@result$Description
library(enrichplot)

# "Cytokine-cytokine receptor interaction"
gsea_res = "analysis/gesea/"

png(filename = paste0(paste0(gsea_res,condition_name),gsub(" ","_",gse@result$Description[1]),".png"),width = 400,height = 300,pointsize = 12)
gseaplot2(gse, geneSetID = 1, title = paste(gse@result$Description[1],"6h"))
dev.off()

# NF-kappa B signaling pathway 
png(filename = paste0(paste0(gsea_res,condition_name),gsub(" ","_",gse@result$Description[5]),".png"),width = 400,height = 300,pointsize = 12)
gseaplot2(gse, geneSetID = 5, title = paste(gse@result$Description[5],"6h"))
dev.off()

# Chemokine signaling pathway
png(filename = paste0(paste0(gsea_res,condition_name),gsub(" ","_",gse@result$Description[7]),".png"),width = 400,height = 300,pointsize = 12)
gseaplot2(gse, geneSetID = 7, title = paste(gse@result$Description[7],"6h"))
dev.off()

# "NOD-like receptor signaling pathway" 
png(filename = paste0(paste0(gsea_res,condition_name),gsub(" ","_",gse@result$Description[13]),".png"),width = 400,height = 300,pointsize = 12)
gseaplot2(gse, geneSetID = 13, title = paste(gse@result$Description[13],"6h"))
dev.off()

# Cytosolic DNA-sensing pathway
png(filename = paste0(paste0(gsea_res,condition_name),gsub(" ","_",gse@result$Description[18]),".png"),width = 400,height = 300,pointsize = 12)
gseaplot2(gse, geneSetID = 18, title = paste(gse@result$Description[18],"6h"))
dev.off()

# JAK-STAT signaling pathway
png(filename = paste0(paste0(gsea_res,condition_name),gsub(" ","_",gse@result$Description[26]),".png"),width = 400,height = 300,pointsize = 12)
gseaplot2(gse, geneSetID = 26, title = paste(gse@result$Description[26],"6h"))
dev.off()

# Proteasome
png(filename = paste0(paste0(gsea_res,condition_name),gsub(" ","_",gse@result$Description[27]),".png"),width = 400,height = 300,pointsize = 12)
gseaplot2(gse, geneSetID = 27, title = paste(gse@result$Description[27],"6h"))
dev.off()

# RIG-I-like
png(filename = paste0(paste0(gsea_res,condition_name),gsub(" ","_",gse@result$Description[24]),".png"),width = 400,height = 300,pointsize = 12)
gseaplot2(gse, geneSetID = 24, title = paste(gse@result$Description[24],"6h"))
dev.off()

# C-type
png(filename = paste0(paste0(gsea_res,condition_name),gsub(" ","_",gse@result$Description[40]),".png"),width = 400,height = 300,pointsize = 12)
gseaplot2(gse, geneSetID = 40, title = paste(gse@result$Description[40],"6h"))
dev.off()

###







#### FOR THE PAPER NOT RELEVANT ####

### Group 3h vs control ####
x = x2
counts = cpm(x,normalized.lib.sizes =T,log = T)
res = topTags(x,n = Inf)
res = as.data.frame(res)
slist = as.data.frame(topTags(x,n=10)[,2])

condition_name = "Treated_versus_Control_(3h)"

vp = EnhancedVolcano(res,x = "logFC",y = "FDR",
                     title = condition_name,
                     subtitle = "",
                     lab = res$GENENAME,
                     selectLab = dlist,
                     pCutoff  = 0.05,
                     FCcutoff = 1.0,
                     boxedLabels = TRUE,
                     parseLabels = TRUE,
                     colAlpha = 4/5,
                     drawConnectors = TRUE,
                     legendPosition = 'bottom',
                     max.overlaps = 20)

png(filename = paste0(paste0("analysis/",condition_name),"_volcano.png"),
    width = 600,height = 600,pointsize = 12)

plot(vp)
dev.off()

profile = res[abs(res$logFC)>=1 & res$FDR<0.05,]
full_res = cbind(profile[,c("GENENAME","logFC","PValue","FDR")],counts[rownames(profile),])

write.table(full_res,"analysis/Treated_versus_Control_(3h).csv",
            sep = ",",
            row.names = F,
            col.names = T)

gene_list = res[,"logFC"]
names(gene_list) = res$GENEID

gene_list = na.omit(gene_list)
gene_list = gene_list[unique(names(gene_list))]
gene_list = sort(gene_list,decreasing = T)

gseBP = gseGO(geneList = gene_list,
              keyType = "ENSEMBL",
              pvalueCutoff = 0.05,
              minGSSize = 50,
              maxGSSize = 500,
              seed = TRUE,
              verbose = T,
              eps = 0,
              ont          = "BP",
              OrgDb = org.Hs.eg.db)

gseMF = gseGO(geneList = gene_list,
              keyType = "ENSEMBL",
              pvalueCutoff = 0.05,
              verbose = T,
              eps = 0,
              ont          = "CC",
              OrgDb = org.Hs.eg.db)

p1 = dotplot(gseBP,showCategory = 20, 
             font.size = 10)+ggtitle("ONT:Biological Process")+scale_fill_viridis_c()

p2 = dotplot(gseMF,showCategory = 20, 
             font.size = 10)+ggtitle("ONT:Cellular Component")+scale_fill_viridis_c()

cw = cowplot::plot_grid(p1, p2)

png(filename = paste0(paste0("analysis/",condition_name),"_GO-Enrichment.png"), 
    width = 1000, height = 480,pointsize = 10)

plot(cw)
dev.off()

### Total treated versus untreated #####
condition_name = "Total_Treatment_against_Control"

x = x3
counts = cpm(x,normalized.lib.sizes =T,log = T)
res = topTags(x,n = Inf)
res = as.data.frame(res)
slist = as.data.frame(topTags(x,n=10)[,2])

vp = EnhancedVolcano(res,x = "logFC",y = "FDR",
                     title = condition_name,
                     subtitle = "",
                     lab = res$GENENAME,
                     selectLab = dlist,
                     pCutoff  = 0.05,
                     FCcutoff = 1.0,
                     boxedLabels = TRUE,
                     parseLabels = TRUE,
                     colAlpha = 4/5,
                     drawConnectors = TRUE,
                     legendPosition = 'bottom',
                     max.overlaps = 20)

png(filename = paste0(paste0("analysis/",condition_name),"_volcano.png"),
    width = 600,height = 600,pointsize = 12)

plot(vp)
dev.off()

profile = res[abs(res$logFC)>=1 & res$FDR<0.05,]
full_res = cbind(profile[,c("GENENAME","logFC","PValue","FDR")],counts[rownames(profile),])

write.table(full_res,"analysis/Total_Treated_versus_Control.csv",
            sep = ",",
            row.names = F,
            col.names = T)

gene_list = res[,"logFC"]
names(gene_list) = res$GENEID

gene_list = na.omit(gene_list)
gene_list = gene_list[unique(names(gene_list))]
gene_list = sort(gene_list,decreasing = T)

gseBP = gseGO(geneList = gene_list,
              keyType = "ENSEMBL",
              pvalueCutoff = 0.05,
              minGSSize = 50,
              maxGSSize = 500,
              seed = TRUE,
              verbose = T,
              eps = 0,
              ont          = "BP",
              OrgDb = org.Hs.eg.db)

gseMF = gseGO(geneList = gene_list,
              keyType = "ENSEMBL",
              pvalueCutoff = 0.05,
              verbose = T,
              eps = 0,
              ont          = "CC",
              OrgDb = org.Hs.eg.db)

p1 = dotplot(gseBP,showCategory = 20,
             font.size = 10)+ggtitle("ONT:Biological Process")+scale_fill_viridis_c()

p2 = dotplot(gseMF,showCategory = 20, 
             font.size = 10)+ggtitle("ONT:Cellular Component")+scale_fill_viridis_c()

cw = cowplot::plot_grid(p1, p2)

png(filename = paste0(paste0("analysis/",condition_name),"_GO-Enrichment.png"), 
    width = 1000, height = 480,pointsize = 10)

plot(cw)
dev.off()

r = gseBP@result["GO:0009615",]
q = gseBP@result
write.table(q,"GO_Genes.tsv",sep = "\t")

# Heatmap ####

f = full_res[,5:ncol(full_res)]

Time = RColorBrewer::brewer.pal(3,"Set1")
Treatment = RColorBrewer::brewer.pal(3,"Dark2")

ann_col = list(Treatment = c(untreated = Treatment[1], germlings_1to1 = Treatment[2]),
               Time = c("3h" = Time[1], "6h" = Time[2]))

png(filename = paste0("analysis/","Heatmap_log2_normCounts_allSigGenes.png"),width = 500, height = 480,pointsize = 10)

pheatmap::pheatmap(f,
                   annotation_col = samples[,c("Time","Treatment"),drop = F],
                   main = "log2Normalized Counts",
                   show_rownames = F,
                   show_colnames = F,
                   annotation_colors = ann_col,
                   gaps_col = 4)

dev.off()


cls = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)

cols = ComplexHeatmap::HeatmapAnnotation(Time = samples$Time,col = list(Time = c("3h"= "#800080","6h" = "#8bc34a")))

ComplexHeatmap::Heatmap(as.matrix(f), 
                        name = "log2Normalized Counts",
                        show_row_names = F,
                        top_annotation = cols,
                        show_column_names = F,
                        column_split = samples$Treatment,
                        column_gap = unit(10, "mm"),
                        col = cls)





png(filename = paste0("analysis/","Heatmap_log2_normCounts_GenesOfIntereset.png"),width = 600, height = 480,pointsize = 10)

rownames(full_res) = full_res$GENENAME

f = full_res[full_res$GENENAME %in% dlist,5:ncol(full_res)]

pheatmap::pheatmap(f,
                   annotation_col = samples[,c("Time","Treatment"),drop = F],
                   main = "log2Normalized Counts",
                   show_rownames = T,
                   show_colnames = F,
                   cluster_cols = T,
                   annotation_colors = ann_col)

dev.off()

cls = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)

cols = ComplexHeatmap::HeatmapAnnotation(Time = samples$Time,col = list(Time = c("3h"= "#800080","6h" = "#8bc34a")))

ComplexHeatmap::Heatmap(as.matrix(f), 
                        name = "log2Normalized Counts",
                        show_row_names = T,
                        top_annotation = cols,
                        show_column_names = F,
                        column_split = samples$Treatment,
                        column_gap = unit(10, "mm"),
                        col = cls)



# Second GOI Heatmap ####

profile = res
full_res = cbind(profile[,c("GENENAME","logFC","PValue","FDR")],counts[rownames(profile),])


full_res = full_res[!duplicated(full_res$GENENAME),]

rownames(full_res) = full_res$GENENAME

t = full_res[full_res$GENENAME %in% blist,5:ncol(full_res)]
colnames(t)

pheatmap::pheatmap(t,
                   annotation_col = samples[,c("Time","Treatment"),drop = F],
                   main = "log2Normalized Counts - Genes of Interest",
                   show_rownames = T,
                   show_colnames = F,
                   cluster_cols = T,
                   annotation_colors = ann_col,
                   border_color = NA)

cls = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)


samples$num = rownames(samples)

ordered_samples = samples %>% arrange(Treatment, Time)
ordered_samples$Treatment = factor(ordered_samples$Treatment, levels = c("untreated", "germlings_1to1"))
cols = ComplexHeatmap::HeatmapAnnotation(Time = ordered_samples$Time,col = list(Time = c("3h"= "#800080","6h" = "#8bc34a")))

colOrder = ordered_samples$num

t = as.matrix(t)
t = t[,colOrder]



ComplexHeatmap::Heatmap(t,
                        cluster_columns = F,
                        name = "log2Normalized Counts",
                        top_annotation = cols,
                        show_row_names = T,
                        show_column_names = F,
                        column_split = ordered_samples$Treatment,
                        column_gap = unit(10, "mm"),
                        col = cls)



# PCA ####
tpm = txi$abundance
tpm = tpm[apply(tpm,1,var)>0,]
tpm = t(log2(tpm+1))

pca =  prcomp(tpm, scale. = T, center = T)

trt = samples$Treatment
trt = ifelse(trt == "untreated",Treatment[1],Treatment[2])

names(trt) = rownames(tpm)
sp = samples$Sample

tpm = as.data.frame(tpm)

tpm[,"Group"] = samples$Treatment

p = autoplot(pca,data = tpm,colour = "Group",frame = TRUE,frame.type = 'norm',size = 3)+theme_light()+theme(legend.position = "right", axis.text= element_text(colour="black"),plot.title = element_text(hjust = 0.5))+ggtitle("PCA on all TPMs")

png(filename = paste0("analysis/","PCA_AllGenes_TPM.png"),width = 1000, height = 480,pointsize = 10)
plot(p)
dev.off()

# Meta-Transcriptomics ####
### Read ssamples

isolate_enst <- function(input_string) {
  # Split the string by the pipe (|) symbol
  split_string <- unlist(strsplit(input_string, "\\|"))
  
  # Extract the ENST number (first element) and remove version (e.g., .5)
  enst <- unlist(strsplit(split_string[1], "\\."))[1]
  
  return(enst)
}


dat = read.table("~/Projects/Collab_Heiko/Meta/quant.tsv", header = T)
rownames(dat) = dat$Name
hs = dat[,c(2,4,6)]
hs = hs[complete.cases(hs),]
hs = hs[rowSums(hs)>1,]

rownames(hs) = unlist(lapply(row.names(hs), isolate_enst))
hs$ENST = rownames(hs)

merged = merge(hs, tx2gene, by.x = "ENST", by.y = "TXNAME")

head(merged)

aggregated_expression_ensg <- merged %>%
  group_by(GENEID) %>%
  summarise(
    total_gfp_pos_hs = sum(gfp_pos_hs, na.rm = TRUE),
    total_K57_hs = sum(K57_hs, na.rm = TRUE),
    total_K41_hs = sum(K41_hs, na.rm = TRUE)
  )

aggregated_expression_ensg <- aggregated_expression_ensg %>%
  dplyr::left_join(gene2sym, by = "GENEID") %>%
  rename(sym = GENENAME)

aggregated_expression_ensg[,c("q75","q90")] = quantile(aggregated_expression_ensg[])

genes = c("IRF1","ISG15","MX2"," OASL","IL1B","CXCL10")
geneoi = as.data.frame(aggregated_expression_ensg[aggregated_expression_ensg$sym %in% genes,])

geneoi = geneoi[,c("sym","total_gfp_pos_hs","total_K57_hs","total_K41_hs")]
row.names(geneoi) = geneoi$sym
geneoi = geneoi[,-1]
colnames()

pheatmap::pheatmap(log2(geneoi),show_colnames = F)

af = dat[,c(3,5,7)]
af = af[complete.cases(af),]
af = af[rowSums(af)>1,]

quantile(af)
#af[,c("q90","q99")] = apply(af,1,quantile,c(0.9,0.99))
af$median = apply(af,1,median)

af = af[af$median>quantile(af$median,0.90),]

rownames(af) = gsub("\\-T","",rownames(af),)

rownames(af)

af_sm = c("ncRNA_gene-Afu4g02350","ncRNA_gene-Afu4g01850","hsp70-Afu1g07440",
          "Ubiquitin-Afu3g11260","Polyubiquitin-Afu4g10350","act1-Afu6g04740","ncRNA_gene-Afu4g02100",
          "Tubulin_alpha-1_subunit-Afu1g02550","Beta-tubulin-Afu1g10910")

rownames(af) = af_sm
pheatmap::pheatmap(log2(af[,-4]+1),
                   show_colnames = F)



