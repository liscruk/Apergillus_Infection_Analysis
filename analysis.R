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

library(httr);library(jsonlite)
### Get raw list of targets
### IRF1
irf1_content = GET("https://maayanlab.cloud/Harmonizome/api/1.0/gene_set/IRF1/ChEA+Transcription+Factor+Targets")

irf1_targets = fromJSON(
  content(
    irf1_content,
    as = "text",
    encoding = "UTF-8"))

irf1_targets = irf1_targets$associations$gene[1]$symbol

### IRF7
irf7_content = read_tsv("Aspergilluis/data//IRF7_01.v2025.1.Hs.tsv")
irf7_content = unlist (irf7_content[17,])
irf7_content = unlist(strsplit(irf7_content[2],","))

### SP1
sp1_content = read_tsv("Aspergilluis/data/SP1_01.v2025.1.Hs.tsv")
sp1_content = unlist(sp1_content[17,])
sp1_content = unlist(strsplit(sp1_content[2],","))

### PPARG
pparg_content = GET("https://maayanlab.cloud/Harmonizome/api/1.0/gene_set/PPARG/ChEA+Transcription+Factor+Targets")

pparg_targets = fromJSON(
  content(
    pparg_content, 
    as = "text", 
    encoding = "UTF-8"))

pparg_targets = pparg_targets$associations$gene[1]$symbol


### Overlapp tests
tfs = list(irf1_targets,irf7_content,sp1_content,pparg_targets)
names(tfs) = c("IRF1","IRF7","SP1","PPARG")

ovres = list()

for (z in 1:length(tfs)) {

  print(z)
  fcs = seq(1,max(res$logFC), by = 0.05)
  odds = c()
  ins = c()
  pval = c()
  lfc = c()
  j = 1
  
  for (i in fcs) {
    
    a = res[with(res,logFC>i & FDR<0.05),]
    b = tfs[[z]]
    a = a$GENENAME
    
    print(table(a %in% b))
    
    
    f = newGeneOverlap(a,b, genome.size = length(keep))
    f = testGeneOverlap(f)
    
    
    lfc[j] = i
    odds[j] = f@odds.ratio
    ins[j] = length(f@intersection)
    pval[j] = f@pval
    
    j = j+1
    
  }
  
  
  tfdf = data.frame(lfc,odds,ins,pval)
  tfdf$tf = names(tfs[z])
  tfdf$padjust = p.adjust(tfdf$pval,method = "fdr")
  tfdf$sig = ifelse(tfdf$padjust<0.05,1,0)
  tfdf = tfdf[tfdf$ins>0,]
  
  
  
  ovres[[z]] = tfdf
  
    
}

tfdf = do.call(rbind,ovres)


ggplot(tfdf,aes(-log10(padjust),odds))+
  geom_point(aes(size = ins,color = lfc))+
  geom_vline(xintercept = -log10(0.05),linetype = "longdash", color = "black")+
  geom_hline(yintercept = 1, linetype = "longdash", color = "black")+
  geom_smooth(se=FALSE)+
  scale_colour_gradient(low = "blue", high="red")+
  ggtitle("Overlapps")+
  theme(axis.text = element_text(color = "black", size = 16))+
  theme_bw()+
  facet_wrap(~tf)












