# Script to annotate structural variants based on the position of their breakpoints

library("optparse")
# get options

option_list = list(
  make_option(c("-r", "--reference"), type="character", default=".", help="Path to gtf file with gene annotations [default= %default]", metavar="character"),
  make_option(c("-i", "--input"), type="character", default=".", help="Folder with SV results [default= %default]", metavar="character"),
  make_option(c("-l", "--individual_label"), type="character", default="[-]*[B]*_T[U0-9]*$", help="Suffix indicating individual samples from a same group (e.g., regions of a same tumor, temporal samples) [default= %default]", metavar="character"),
  make_option(c("-m", "--multi"), type="logical", default=F, help="Trigger multi-sample mode where discarded low-confidence SVs are recovered if found as high-confidence in one region  [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# load libraries
library(data.table)
library(tidyverse)
library(GenomicAlignments)
library(readxl)

# load ref
ref = rtracklayer::import(opt$reference, format="gtf")
ref = as.data.frame(ref)

### Prepare ref exon and ref intron
## Exon
ref = as.data.table(ref)
ref.exon = ref[which(ref$type == "exon"),]

ref.exon = ref.exon[which(ref.exon$transcript_type %in% c("protein_coding","lncRNA")),]

ref.exon[, exonNum := ifelse(strand=="+", seq(.N), rev(seq(.N))), by=gene_id]
ref.exon[, exonCount := max(exonNum), by=gene_id]
ref.exon[, exonCountTrans := max(exon_number), by=list(gene_id, transcript_name)]

## Intron
ends = ref.exon[exonCount > 1 & exonCountTrans > 1, ifelse(strand=="+", start[seq(2, .N)], start[seq(1, .N-1)]) , by=list(gene_id, transcript_id)] 
starts = ref.exon[exonCount > 1 & exonCountTrans > 1, ifelse(strand=="+", end[seq(1, .N-1)], end[seq(2, .N)]) , by=list(gene_id, transcript_id)] 
stopifnot(nrow(ends) == nrow(starts))
stopifnot(all(ends$id == starts$id))
stopifnot(min(as.numeric(names(table(ref.exon[exonCount>1, exonCount])))) > 1)

ref.intron = data.table(id=ends$gene_id, start=starts$V1, end=ends$V1, transcript_id = ends$transcript_id)
ref.intron = unique(ref.intron[,list(id,transcript_id,start,end)])
short = ref.exon[!duplicated(gene_id), .(gene_id, gene_name, seqnames, strand, gene_type)]
setkey(ref.intron, id)
setkey(short, gene_id)
ref.intron = short[ref.intron]
ref.intron[, intronNum := ifelse(strand == "+", seq(.N), rev(seq(.N))), by=gene_id]
ref.intron[, intronCount := max(intronNum), by=gene_id]
ref.intron = ref.intron[which(ref.intron$gene_type %in% c("protein_coding","lncRNA")),]

## CDS
ref.cds = ref[which(ref$type == "CDS"),]
ref.cds = ref.cds[which(ref.cds$transcript_type %in% c("protein_coding","lncRNA")),]

ref.cds[, cdsNum := ifelse(strand=="+", seq(.N), rev(seq(.N))), by=gene_id]
ref.cds[, cdsCount := max(cdsNum), by=gene_id]

### Load raw data
vcf.SV.files = list.files(opt$input,pattern = "cns.integration.vcf",full.names = T)
vcf.SV.names = str_remove(list.files(opt$input,pattern = "cns.integration.vcf",full.names = F), ".cns.integration.vcf")
vcf.SVl      = lapply(vcf.SV.files, read_tsv,comment = "#", col_names = c("CHROM",	"POS",	"ID",	"REF",	"ALT",	"QUAL",	"FILTER",	"INFO",	"FORMAT","normal",	"tumor","file"),col_types="cicccicccccc" )

for(i in 1:length(vcf.SVl)){
  vcf.SVl[[i]]$Sample = vcf.SV.names[i]
  vcf.SVl[[i]]$Group  = str_remove(vcf.SV.names,opt$individual_label)[i]
}
dataset.sv = bind_rows(vcf.SVl)

## remove SV in chromosomes
dataset.sv = dataset.sv %>% filter(CHROM%in%paste0("chr",c(1:22,"X","Y") ) )

## Start with Breakpoint 1
dataset.sv1 = dataset.sv
colnames(dataset.sv1)[which(colnames(dataset.sv1) == "CHROM")] = "seqnames"
dataset.sv1$seqnames = gsub("^chr", "", dataset.sv1$seqnames)
all(str_remove(str_extract(dataset.sv1$INFO,"CIPOS=[0-9\\-]+"),"CIPOS=") <= 0)
dataset.sv1$start = as.numeric(dataset.sv1$POS) + as.numeric(str_remove(str_extract(dataset.sv1$INFO,"CIPOS=[0-9\\-]+"),"CIPOS="))
all(str_remove(str_extract(dataset.sv1$INFO,"CIPOS=[0-9\\-]+,[0-9\\-]+"),"CIPOS=[0-9\\-]+,") >= 0)
dataset.sv1$end = as.numeric(dataset.sv1$POS) + as.numeric(str_remove(str_extract(dataset.sv1$INFO,"CIPOS=[0-9\\-]+,[0-9\\-]+"),"CIPOS=[0-9\\-]+,"))

# Annotating overlapping exon
ref.exon$seqnames = gsub("^chr", "", ref.exon$seqnames)
ref.exon = as.data.table(ref.exon)
setnames(ref.exon,c("gene_name","gene_id"),c("geneSymbol","id"))

dataset.sv1 = as.data.table(dataset.sv1)
ref.exon = as.data.table(ref.exon)

rm(ref)
gc()
#


##
ro = findOverlaps(GRanges(dataset.sv1), GRanges(ref.exon)) 
dataset.sv1$exon_id     = dataset.sv1$exon_geneid = dataset.sv1$exon_type   = NA
dataset.sv1$exon_id[queryHits(ro)] = ref.exon$exon_id[subjectHits(ro)]
dataset.sv1$exon_geneid[queryHits(ro)] = ref.exon$id[subjectHits(ro)] 
dataset.sv1$exon_type[queryHits(ro)] = ref.exon$gene_type[subjectHits(ro)]

# Annotating overlapping coding exon
ref.cds$seqnames = gsub("^chr", "", ref.cds$seqnames)
ref.cds = as.data.table(ref.cds)
setnames(ref.cds,c("gene_name","gene_id"),c("geneSymbol","id"))

dataset.sv1 = as.data.table(dataset.sv1)

ro = findOverlaps(GRanges(dataset.sv1), GRanges(ref.cds)) 
dataset.sv1$cds_id     = dataset.sv1$cds_geneid = dataset.sv1$cds_type   = NA
dataset.sv1$cds_id[queryHits(ro)] = ref.cds$exon_id[subjectHits(ro)]
dataset.sv1$cds_geneid[queryHits(ro)] = ref.cds$id[subjectHits(ro)] 
dataset.sv1$cds_type[queryHits(ro)] = ref.cds$gene_type[subjectHits(ro)]

# Annotating overlapping intron
ref.intron$seqnames = gsub("^chr", "", ref.intron$seqnames)
ref.intron = as.data.table(ref.intron)
setnames(ref.intron, c("gene_name","gene_id"),c("geneSymbol","id"))

ref.intron = as.data.table(ref.intron)

ro = findOverlaps(GRanges(dataset.sv1), GRanges(ref.intron)) 
dataset.sv1$intron_id     = dataset.sv1$intron_geneid = dataset.sv1$intron_type   = NA
dataset.sv1$intron_id[queryHits(ro)] = ref.intron$intronNum[subjectHits(ro)]
dataset.sv1$intron_geneid[queryHits(ro)] = ref.intron$id[subjectHits(ro)] 
dataset.sv1$intron_type[queryHits(ro)] = ref.intron$gene_type[subjectHits(ro)]

## Breakpoint 2
dataset.sv2 = dataset.sv1
dataset.sv2$seqnames = str_remove(str_extract(dataset.sv2$INFO,"CHR2=chr[0-9MXY]+"),"CHR2=chr")

all(str_remove(str_extract(dataset.sv2$INFO,"CIEND=[0-9\\-]+"),"CIEND=") <= 0)
dataset.sv2$start = as.numeric(str_remove(str_extract(dataset.sv2$INFO,";END=[0-9\\-]+"),";END=")) + as.numeric(str_remove(str_extract(dataset.sv2$INFO,"CIEND=[0-9\\-]+"),"CIEND="))
all(str_remove(str_extract(dataset.sv2$INFO,"CIEND=[0-9\\-]+,[0-9\\-]+"),"CIEND=[0-9\\-]+,") >= 0)
dataset.sv2$end = as.numeric(str_remove(str_extract(dataset.sv2$INFO,";END=[0-9\\-]+"),";END=")) + as.numeric(str_remove(str_extract(dataset.sv2$INFO,"CIEND=[0-9\\-]+,[0-9\\-]+"),"CIEND=[0-9\\-]+,"))

# Annotating overlapping exon
ro = findOverlaps(GRanges(dataset.sv2), GRanges(ref.exon)) 
dataset.sv2$exon_id     = dataset.sv2$exon_geneid = dataset.sv2$exon_type   = NA
dataset.sv2$exon_id[queryHits(ro)] = ref.exon$exon_id[subjectHits(ro)]
dataset.sv2$exon_geneid[queryHits(ro)] = ref.exon$id[subjectHits(ro)] 
dataset.sv2$exon_type[queryHits(ro)] = ref.exon$gene_type[subjectHits(ro)]

# Annotating overlapping coding exon
ro = findOverlaps(GRanges(dataset.sv2), GRanges(ref.cds)) 
dataset.sv2$cds_id     = dataset.sv2$cds_geneid = dataset.sv2$cds_type   = NA
dataset.sv2$cds_id[queryHits(ro)] = ref.cds$exon_id[subjectHits(ro)]
dataset.sv2$cds_geneid[queryHits(ro)] = ref.cds$id[subjectHits(ro)] 
dataset.sv2$cds_type[queryHits(ro)] = ref.cds$gene_type[subjectHits(ro)]

# Annotating overlapping intron
ro = findOverlaps(GRanges(dataset.sv2), GRanges(ref.intron)) 
dataset.sv2$intron_id     = dataset.sv2$intron_geneid = dataset.sv2$intron_type   = NA
dataset.sv2$intron_id[queryHits(ro)] = ref.intron$intronNum[subjectHits(ro)]
dataset.sv2$intron_geneid[queryHits(ro)] = ref.intron$id[subjectHits(ro)] 
dataset.sv2$intron_type[queryHits(ro)] = ref.intron$gene_type[subjectHits(ro)]


### Merging B1 + B2
all(dim(dataset.sv1) == dim(dataset.sv2))
all(dataset.sv1$SV_ID == dataset.sv2$SV_ID)

colnames(dataset.sv1)[which(colnames(dataset.sv1) == "seqnames")] = "CHROM"
colnames(dataset.sv1)[which(colnames(dataset.sv1) == "start")] = "start.B1"
colnames(dataset.sv1)[which(colnames(dataset.sv1) == "end")] = "end.B1"
colnames(dataset.sv1)[which(colnames(dataset.sv1) == "exon_id")] = "exon_id.B1"
colnames(dataset.sv1)[which(colnames(dataset.sv1) == "intron_id")] = "intron_id.B1"
colnames(dataset.sv1)[which(colnames(dataset.sv1) == "cds_id")] = "cds_id.B1"
colnames(dataset.sv1)[which(colnames(dataset.sv1) == "exon_geneid")] = "exon_geneid.B1"
colnames(dataset.sv1)[which(colnames(dataset.sv1) == "intron_geneid")] = "intron_geneid.B1"
colnames(dataset.sv1)[which(colnames(dataset.sv1) == "cds_geneid")] = "cds_geneid.B1"
colnames(dataset.sv1)[which(colnames(dataset.sv1) == "exon_type")] = "exon_type.B1"
colnames(dataset.sv1)[which(colnames(dataset.sv1) == "cds_type")] = "cds_type.B1"
colnames(dataset.sv1)[which(colnames(dataset.sv1) == "intron_type")] = "intron_type.B1"

colnames(dataset.sv2)[which(colnames(dataset.sv2) == "seqnames")] = "CHROM2"
colnames(dataset.sv2)[which(colnames(dataset.sv2) == "start")] = "start.B2"
colnames(dataset.sv2)[which(colnames(dataset.sv2) == "end")] = "end.B2"
colnames(dataset.sv2)[which(colnames(dataset.sv2) == "exon_id")] = "exon_id.B2"
colnames(dataset.sv2)[which(colnames(dataset.sv2) == "intron_id")] = "intron_id.B2"
colnames(dataset.sv2)[which(colnames(dataset.sv2) == "cds_id")] = "cds_id.B2"
colnames(dataset.sv2)[which(colnames(dataset.sv2) == "exon_geneid")] = "exon_geneid.B2"
colnames(dataset.sv2)[which(colnames(dataset.sv2) == "intron_geneid")] = "intron_geneid.B2"
colnames(dataset.sv2)[which(colnames(dataset.sv2) == "cds_geneid")] = "cds_geneid.B2"
colnames(dataset.sv2)[which(colnames(dataset.sv2) == "exon_type")] = "exon_type.B2"
colnames(dataset.sv2)[which(colnames(dataset.sv2) == "cds_type")] = "cds_type.B2"
colnames(dataset.sv2)[which(colnames(dataset.sv2) == "intron_type")] = "intron_type.B2"

dataset.sv = cbind(dataset.sv1, dataset.sv2[,c("CHROM2","start.B2","end.B2","exon_id.B2","exon_geneid.B2","exon_type.B2","cds_id.B2","cds_geneid.B2","cds_type.B2","intron_id.B2","intron_geneid.B2","intron_type.B2")])
head(dataset.sv)
dataset.sv$CHROM = paste0("chr",dataset.sv$CHROM)
dataset.sv$CHROM2 = paste0("chr",dataset.sv$CHROM2)
save(dataset.sv, file = "dataset.sv.RData")

### find damaging alterations ###
## ID
sort(unique(dataset.sv$Sample))
dim(dataset.sv)
dataset.sv$STRANDS = str_remove(str_extract(dataset.sv$INFO,";STRANDS=[+-]+"),";STRANDS=")
dataset.sv$SVTYPE = str_remove(str_extract(dataset.sv$INFO,";SVTYPE=[A-Z]+"),";SVTYPE=")

dataset.sv$type = sapply(1:nrow(dataset.sv), function(i) ifelse(dataset.sv$CHROM[i] != dataset.sv$CHROM2[i],"inter","intra"))
all(dataset.sv$POS[which(dataset.sv$type == "intra")] <= dataset.sv$END[which(dataset.sv$type == "intra")])
table(dataset.sv$STRANDS,dataset.sv$SVTYPE,dataset.sv$type) # all inter should be TRA
any(dataset.sv$POS[which(dataset.sv$STRANDS == "++" & dataset.sv$SVTYPE == "INV")] == dataset.sv$END[which(dataset.sv$STRANDS == "++" & dataset.sv$SVTYPE == "INV")])

all(which(dataset.sv$cds_id.B1 != "") %in% which(dataset.sv$exon_id.B1 != "")) # matching genes in BP1
all(which(dataset.sv$cds_id.B2 != "") %in% which(dataset.sv$exon_id.B2 != "")) # matching genes in BP2

## Remove the non-coding SVs
length(which(is.na(dataset.sv$exon_geneid.B1) & is.na(dataset.sv$exon_geneid.B2)  & is.na(dataset.sv$intron_geneid.B1) & is.na(dataset.sv$intron_geneid.B2) ) )
dataset.sv.coding = dataset.sv[which(!is.na(dataset.sv$exon_geneid.B1) | !is.na(dataset.sv$exon_geneid.B2)  | !is.na(dataset.sv$intron_geneid.B1) | !is.na(dataset.sv$intron_geneid.B2) ),] #At least one breakpoint in a gene body

dataset.sv.coding$Gene.id1 = sapply(1:nrow(dataset.sv.coding), function(i) paste0(unique(unlist(strsplit(unlist(dataset.sv.coding[i,c("exon_geneid.B1","cds_geneid.B1","intron_geneid.B1")])[which(!is.na(unlist(dataset.sv.coding[i,c("exon_geneid.B1","cds_geneid.B1","intron_geneid.B1")])))],","))), collapse = "_"))
all(unlist(strsplit(dataset.sv.coding$Gene.id1[which(dataset.sv.coding$Gene.id1 != "")],"_")) %in% ref.exon$id)

dataset.sv.coding$Gene.id2 = sapply(1:nrow(dataset.sv.coding), function(i) paste0(unique(unlist(strsplit(unlist(dataset.sv.coding[i,c("exon_geneid.B2","cds_geneid.B2","intron_geneid.B2")])[which(!is.na(unlist(dataset.sv.coding[i,c("exon_geneid.B2","cds_geneid.B2","intron_geneid.B2")])))],","))), collapse = "_"))
all(unlist(strsplit(dataset.sv.coding$Gene.id2[which(dataset.sv.coding$Gene.id2 != "")],"_")) %in% ref.exon$id)
any(dataset.sv.coding$Gene.id1 == "" & dataset.sv.coding$Gene.id1 == dataset.sv.coding$Gene.id2) #non-coding SV removed; should be FALSE

all(grepl("lncRNA",dataset.sv.coding$exon_type.B1[which(grepl("protein_coding",dataset.sv.coding$exon_type.B1) & dataset.sv.coding$exon_type.B1 != "protein_coding")]))
all(grepl("lncRNA",dataset.sv.coding$intron_type.B1[which(grepl("protein_coding",dataset.sv.coding$intron_type.B1) & dataset.sv.coding$intron_type.B1 != "protein_coding")]))
all(grepl("lncRNA",dataset.sv.coding$exon_type.B2[which(grepl("protein_coding",dataset.sv.coding$exon_type.B2) & dataset.sv.coding$exon_type.B2 != "protein_coding")]))
all(grepl("lncRNA",dataset.sv.coding$intron_type.B2[which(grepl("protein_coding",dataset.sv.coding$intron_type.B2) & dataset.sv.coding$intron_type.B2 != "protein_coding")]))

### If (++;--) kept if it involves at least one lncRNA or protein coding gene
d1 = dataset.sv.coding[which(dataset.sv.coding$STRANDS %in% c("++","--") & (grepl("protein_coding|lncRNA", dataset.sv.coding$exon_type.B1) | grepl("protein_coding|lncRNA", dataset.sv.coding$intron_type.B1))),c("Gene.id1","Sample","type","SVTYPE","CHROM","start.B1","end.B1","CHROM2","start.B2","end.B2","tumor","STRANDS")]
d2 = dataset.sv.coding[which(dataset.sv.coding$STRANDS %in% c("++","--") & (grepl("protein_coding|lncRNA", dataset.sv.coding$exon_type.B2) | grepl("protein_coding|lncRNA", dataset.sv.coding$intron_type.B2))),c("Gene.id2","Sample","type","SVTYPE","CHROM","start.B1","end.B1","CHROM2","start.B2","end.B2","tumor","STRANDS")]
d1$gene_breakpoint = "g1"
d2$gene_breakpoint = "g2"
colnames(d1)[1:2] = c("ID","Sample_name")
colnames(d2)[1:2] = c("ID","Sample_name")

inter = intersect(d1$SV_ID, d2$SV_ID)
dbis = rbind(d1[which(d1$SV_ID %in% inter),], d2[which(d2$SV_ID %in% inter),])
dbis$gene_breakpoint = "g1&g2"
d1 = d1[which(!d1$SV_ID %in% inter),]
d2 = d2[which(!d2$SV_ID %in% inter),]

### If (+-;-+) need to damage coding part
## lncRNA at least one exon involved
# Single gene
d3 = dataset.sv.coding[which(dataset.sv.coding$STRANDS %in% c("+-","-+") & dataset.sv.coding$Gene.id1 == dataset.sv.coding$Gene.id2 & (grepl("lncRNA",dataset.sv.coding$exon_type.B1) | grepl("lncRNA",dataset.sv.coding$intron_type.B1)) & (dataset.sv.coding$exon_id.B1 != "" | dataset.sv.coding$exon_id.B2 != "")),c("Gene.id1","Sample","type","SVTYPE","CHROM","start.B1","end.B1","CHROM2","start.B2","end.B2","tumor","STRANDS")] #At least one exon is involved because cut
d4 = dataset.sv.coding[which(dataset.sv.coding$STRANDS %in% c("+-","-+") & dataset.sv.coding$Gene.id1 == dataset.sv.coding$Gene.id2 & 
                               (grepl("lncRNA",dataset.sv.coding$exon_type.B1) | grepl("lncRNA",dataset.sv.coding$intron_type.B1)) & 
                               is.na(dataset.sv.coding$exon_id.B1) & is.na(dataset.sv.coding$exon_id.B2) )[which(sapply(which(dataset.sv.coding$STRANDS %in% c("+-","-+") & dataset.sv.coding$Gene.id1 == dataset.sv.coding$Gene.id2 & 
                                                                                                                                (grepl("lncRNA",dataset.sv.coding$exon_type.B1) | grepl("lncRNA",dataset.sv.coding$intron_type.B1)) & 
                                                                                                                                is.na(dataset.sv.coding$exon_id.B1) & is.na(dataset.sv.coding$exon_id.B2) ), function(i) any(ref.exon$start[which(ref.exon$id %in% unlist(strsplit(dataset.sv.coding$Gene.id1[i],"_")))] > dataset.sv.coding$end.B1[i] & ref.exon$end[which(ref.exon$id %in% unlist(strsplit(dataset.sv.coding$Gene.id1[i],"_")))] < dataset.sv.coding$start.B2[i])))],c("Gene.id1","Sample","type","SVTYPE","CHROM","start.B1","end.B1","CHROM2","start.B2","end.B2","tumor","STRANDS")] #At least one exon is involved because overlapping
d3$gene_breakpoint = "g1&g2"
d4$gene_breakpoint = "g1&g2"

# Not same regions
d5 = dataset.sv.coding[which(dataset.sv.coding$STRANDS %in% c("+-","-+") & !is.na(dataset.sv.coding$Gene.id1) & dataset.sv.coding$Gene.id1 != dataset.sv.coding$Gene.id2 & !is.na(dataset.sv.coding$exon_id.B1) & grepl("lncRNA",dataset.sv.coding$exon_type.B1)),c("Gene.id1","Sample","type","SVTYPE","CHROM","start.B1","end.B1","CHROM2","start.B2","end.B2","tumor","STRANDS")] #At least one exon is involved because cut
d6 = dataset.sv.coding[which(dataset.sv.coding$STRANDS %in% c("+-","-+") & !is.na(dataset.sv.coding$Gene.id2) & dataset.sv.coding$Gene.id1 != dataset.sv.coding$Gene.id2 & !is.na(dataset.sv.coding$exon_id.B2) & grepl("lncRNA",dataset.sv.coding$exon_type.B2)),c("Gene.id2","Sample","type","SVTYPE","CHROM","start.B1","end.B1","CHROM2","start.B2","end.B2","tumor","STRANDS")] #At least one exon is involved because cut
d5$gene_breakpoint = "g1"
d6$gene_breakpoint = "g2"

d7 = dataset.sv.coding[which(dataset.sv.coding$STRANDS %in% c("+-","-+") & !is.na(dataset.sv.coding$Gene.id1) & dataset.sv.coding$Gene.id1 != dataset.sv.coding$Gene.id2 & (grepl("lncRNA",dataset.sv.coding$intron_type.B1) | grepl("lncRNA",dataset.sv.coding$exon_type.B1)))[which(sapply(which(dataset.sv.coding$STRANDS %in% c("+-","-+") & !is.na(dataset.sv.coding$Gene.id1) & dataset.sv.coding$Gene.id1 != dataset.sv.coding$Gene.id2 & (grepl("lncRNA",dataset.sv.coding$intron_type.B1) | grepl("lncRNA",dataset.sv.coding$exon_type.B1))), 
                                                                                                                                                                                                                                                                                             function(i) any(ref.exon$start[which(ref.exon$id %in% unlist(strsplit(dataset.sv.coding$Gene.id1[i],"_")))] > dataset.sv.coding$end.B1[i])))],c("Gene.id1","Sample","type","SVTYPE","CHROM","start.B1","end.B1","CHROM2","start.B2","end.B2","tumor","STRANDS")] #At least one exon is involved because overlapping
d8 = dataset.sv.coding[which(dataset.sv.coding$STRANDS %in% c("+-","-+") & !is.na(dataset.sv.coding$Gene.id2) & dataset.sv.coding$Gene.id1 != dataset.sv.coding$Gene.id2 & (grepl("lncRNA",dataset.sv.coding$intron_type.B2) | grepl("lncRNA",dataset.sv.coding$exon_type.B2)))[which(sapply(which(dataset.sv.coding$STRANDS %in% c("+-","-+") & !is.na(dataset.sv.coding$Gene.id2) & dataset.sv.coding$Gene.id1 != dataset.sv.coding$Gene.id2 & (grepl("lncRNA",dataset.sv.coding$intron_type.B2) | grepl("lncRNA",dataset.sv.coding$exon_type.B2))), 
                                                                                                                                                                                                                                                                                             function(i) any(ref.exon$start[which(ref.exon$id %in% unlist(strsplit(dataset.sv.coding$Gene.id2[i],"_")))] < dataset.sv.coding$start.B2[i])))],c("Gene.id2","Sample","type","SVTYPE","CHROM","start.B1","end.B1","CHROM2","start.B2","end.B2","tumor","STRANDS")] #At least one exon is involved because overlapping
d7$gene_breakpoint = "g1"
d8$gene_breakpoint = "g2"
## protein_coding at least one cds involved
# Single gene
d9  = dataset.sv.coding[which(dataset.sv.coding$STRANDS %in% c("+-","-+") & dataset.sv.coding$Gene.id1 == dataset.sv.coding$Gene.id2 & 
                                (dataset.sv.coding$exon_type.B1 == "protein_coding" | dataset.sv.coding$intron_type.B1 == "protein_coding") & 
                                (!is.na(dataset.sv.coding$cds_id.B1) | !is.na(dataset.sv.coding$cds_id.B2))),c("Gene.id1","Sample","type","SVTYPE","CHROM","start.B1","end.B1","CHROM2","start.B2","end.B2","tumor","STRANDS")] #At least one exon is involved because cut
d10 = dataset.sv.coding[which(dataset.sv.coding$STRANDS %in% c("+-","-+") & dataset.sv.coding$Gene.id1 == dataset.sv.coding$Gene.id2 & 
                                (dataset.sv.coding$exon_type.B1 == "protein_coding" | dataset.sv.coding$intron_type.B1 == "protein_coding") & 
                                is.na(dataset.sv.coding$cds_id.B1) & is.na(dataset.sv.coding$cds_id.B2) )[which(sapply(which(dataset.sv.coding$STRANDS %in% c("+-","-+") & dataset.sv.coding$Gene.id1 == dataset.sv.coding$Gene.id2 & (dataset.sv.coding$exon_type.B1 == "protein_coding" | dataset.sv.coding$intron_type.B1 == "protein_coding") & is.na(dataset.sv.coding$cds_id.B1) & is.na(dataset.sv.coding$cds_id.B2) ),
                                                                                                                       function(i) any(ref.cds$start[which(ref.cds$id %in% unlist(strsplit(dataset.sv.coding$Gene.id1[i],"_")))] > dataset.sv.coding$end.B1[i] & ref.cds$end[which(ref.cds$id %in% unlist(strsplit(dataset.sv.coding$Gene.id1[i],"_")))] < dataset.sv.coding$start.B2[i])))],c("Gene.id1","Sample","type","SVTYPE","CHROM","start.B1","end.B1","CHROM2","start.B2","end.B2","tumor","STRANDS")] #At least one exon is involved because overlapping
d9$gene_breakpoint = "g1&g2"
d10$gene_breakpoint = "g1&g2"
# Not same regions
d11 = dataset.sv.coding[which(dataset.sv.coding$STRANDS %in% c("+-","-+") & !is.na(dataset.sv.coding$Gene.id1) & dataset.sv.coding$Gene.id1 != dataset.sv.coding$Gene.id2 & !is.na(dataset.sv.coding$cds_id.B1) & dataset.sv.coding$cds_type.B1 == "protein_coding"),c("Gene.id1","Sample","type","SVTYPE","CHROM","start.B1","end.B1","CHROM2","start.B2","end.B2","tumor","STRANDS")] #At least one exon is involved because cut
d12 = dataset.sv.coding[which(dataset.sv.coding$STRANDS %in% c("+-","-+") & !is.na(dataset.sv.coding$Gene.id2) & dataset.sv.coding$Gene.id1 != dataset.sv.coding$Gene.id2 & !is.na(dataset.sv.coding$cds_id.B2) & dataset.sv.coding$cds_type.B2 == "protein_coding"),c("Gene.id2","Sample","type","SVTYPE","CHROM","start.B1","end.B1","CHROM2","start.B2","end.B2","tumor","STRANDS")] #At least one exon is involved because cut

d11$gene_breakpoint = "g1"
d12$gene_breakpoint = "g2"

d13 = dataset.sv.coding[which(dataset.sv.coding$STRANDS %in% c("+-","-+") & !is.na(dataset.sv.coding$Gene.id1) & dataset.sv.coding$Gene.id1 != dataset.sv.coding$Gene.id2 & 
                                is.na(dataset.sv.coding$cds_id.B1) & (dataset.sv.coding$exon_type.B1 == "protein_coding" | dataset.sv.coding$intron_type.B1 == "protein_coding"))[which(sapply(which(dataset.sv.coding$STRANDS %in% c("+-","-+") & !is.na(dataset.sv.coding$Gene.id1) & dataset.sv.coding$Gene.id1 != dataset.sv.coding$Gene.id2 & 
                                                                                                                                                                                                       is.na(dataset.sv.coding$cds_id.B1) & (dataset.sv.coding$exon_type.B1 == "protein_coding" | dataset.sv.coding$intron_type.B1 == "protein_coding")), 
                                                                                                                                                                                               function(i) any(ref.cds$start[which(ref.cds$id %in% unlist(strsplit(dataset.sv.coding$Gene.id1[i],"_")))] > dataset.sv.coding$end.B1[i])))],c("Gene.id1","Sample","type","SVTYPE","CHROM","start.B1","end.B1","CHROM2","start.B2","end.B2","tumor","STRANDS")] #At least one exon is involved because overlapping
d14 = dataset.sv.coding[which(dataset.sv.coding$STRANDS %in% c("+-","-+") & !is.na(dataset.sv.coding$Gene.id2) & dataset.sv.coding$Gene.id1 != dataset.sv.coding$Gene.id2 & 
                                is.na(dataset.sv.coding$cds_id.B2) & (dataset.sv.coding$exon_type.B2 == "protein_coding" | dataset.sv.coding$intron_type.B2 == "protein_coding"))[which(sapply(which(dataset.sv.coding$STRANDS %in% c("+-","-+") & !is.na(dataset.sv.coding$Gene.id2) & dataset.sv.coding$Gene.id1 != dataset.sv.coding$Gene.id2 & 
                                                                                                                                                                                                       is.na(dataset.sv.coding$cds_id.B2) & (dataset.sv.coding$exon_type.B2 == "protein_coding" | dataset.sv.coding$intron_type.B2 == "protein_coding")), 
                                                                                                                                                                                               function(i) any(ref.cds$start[which(ref.cds$id %in% unlist(strsplit(dataset.sv.coding$Gene.id2[i],"_")))] < dataset.sv.coding$start.B2[i])))],c("Gene.id2","Sample","type","SVTYPE","CHROM","start.B1","end.B1","CHROM2","start.B2","end.B2","tumor","STRANDS")] #At least one exon is involved because overlapping
d13$gene_breakpoint = "g1"
d14$gene_breakpoint = "g2"

### Merging
colnames(d3)[1:2] = c("ID","Sample_name")
colnames(d4)[1:2] = c("ID","Sample_name")
colnames(d5)[1:2] = c("ID","Sample_name")
colnames(d6)[1:2] = c("ID","Sample_name")
colnames(d7)[1:2] = c("ID","Sample_name")
colnames(d8)[1:2] = c("ID","Sample_name")
colnames(d9)[1:2] = c("ID","Sample_name")
colnames(d10)[1:2] = c("ID","Sample_name")
colnames(d11)[1:2] = c("ID","Sample_name")
colnames(d12)[1:2] = c("ID","Sample_name")
colnames(d13)[1:2] = c("ID","Sample_name")
colnames(d14)[1:2] = c("ID","Sample_name")

# To distinguish genes altered in B1 and B2 to get the list of damaged genes
#dbis$SV_ID = paste0(dbis$SV_ID,"-B1_B2")
#d1$SV_ID = paste0(d1$SV_ID,"-B1") 
#d2$SV_ID = paste0(d2$SV_ID,"-B2")
#d3$SV_ID = paste0(d3$SV_ID,"-B1_B2")
#d4$SV_ID = paste0(d4$SV_ID,"-B1_B2")
#d5$SV_ID = paste0(d5$SV_ID,"-B1") 
#d6$SV_ID = paste0(d6$SV_ID,"-B2")
#d7$SV_ID = paste0(d7$SV_ID,"-B1")
#d8$SV_ID = paste0(d8$SV_ID,"-B2")
#d9$SV_ID = paste0(d9$SV_ID,"-B1_B2")
#d10$SV_ID = paste0(d10$SV_ID,"-B1_B2")
#d11$SV_ID = paste0(d11$SV_ID,"-B1")
#d12$SV_ID = paste0(d12$SV_ID,"-B2")
#d13$SV_ID = paste0(d13$SV_ID,"-B1")
#d14$SV_ID = paste0(d14$SV_ID,"-B2")

# Type of SVs
unique(c(d1$type, d2$type))
d1$CN = ifelse(d1$type == "intra","TRANS_intra","TRANS_inter")
d2$CN = ifelse(d2$type == "intra","TRANS_intra","TRANS_inter")

d3$CN = "large_indel"
d4$CN = "large_indel"

d5$CN = ifelse(d5$type == "intra","TRANS_intra","TRANS_inter")
d6$CN = ifelse(d6$type == "intra","TRANS_intra","TRANS_inter")
d7$CN = ifelse(d7$type == "intra","TRANS_intra","TRANS_inter")
d8$CN = ifelse(d8$type == "intra","TRANS_intra","TRANS_inter")

d9$CN = "large_indel"
d10$CN = "large_indel"

d11$CN = ifelse(d11$type == "intra","TRANS_intra","TRANS_inter")
d12$CN = ifelse(d12$type == "intra","TRANS_intra","TRANS_inter")
d13$CN = ifelse(d13$type == "intra","TRANS_intra","TRANS_inter")
d14$CN = ifelse(d14$type == "intra","TRANS_intra","TRANS_inter")

d = rbind(d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13, d14)
all(unlist(strsplit(d$ID,"_")) %in% ref.exon$id)

#data = data.frame(ID = unlist(strsplit(d$ID,"_")),
#                  Sample_name = unlist(sapply(1:nrow(d), function(i) rep(d$Sample_name[i], length(strsplit(d$ID[i],"_")[[1]])))),
#                  #SV_ID = unlist(sapply(1:nrow(d), function(i) rep(d$SV_ID[i], length(strsplit(d$ID[i],"_")[[1]])))),
#                  SVTYPE = d$SVTYPE,
#                  CN = unlist(sapply(1:nrow(d), function(i) rep(d$CN[i], length(strsplit(d$ID[i],"_")[[1]])))),
                  #)
d = d[which(!duplicated(d)),]
dim(d)

d$Gene = sapply(d$ID, function(x) unique(ref.exon$geneSymbol[which(ref.exon$id == x)]))
d$Group = str_remove(d$Sample_name,opt$individual_label)
d$Gene = sapply(d$Gene,function(x){res=x;if(length(res)==0){res=NA};return(res)})

write_tsv(file = "SVs_annotated.tsv",d)


# Code for multi-region samples
if(opt$multi){
## recover alterations from filtered calls
### list SVs without match at 1kb
dataset.sv$time = "Parental"
dataset.sv$time[str_detect(dataset.sv$Sample,"p")] = "Organoid1"
#dataset.sv = dataset.sv[dataset.sv$Sample!="LNET2Np12",]
dataset.sv$time[dataset.sv$Sample == "LCNEC4Tp24"] = "Organoid2"
dataset.sv$time[dataset.sv$Sample == "PANEC1Tp14"] = "Organoid2"

d.uniq = c()
iexp = 1
for( exp in sort(unique(d.driver$Experiment))){#dataset.sv$Experiment)) ){
  d.tmp = d.driver %>% filter(Experiment==exp)  #dataset.sv %>% filter(Experiment==exp) %>% dplyr::select(-(INFO:file))
  colnames(d.tmp)[2] = "Sample"
  matchs.tmp = lapply(1:nrow(d.tmp) , function(i) which( (d.tmp$CHROM[i]==d.tmp$CHROM & d.tmp$CHROM2[i]==d.tmp$CHROM2 & 
                                                     abs(d.tmp$start.B1[i]-d.tmp$start.B1)<=1000 & abs(d.tmp$start.B2[i]-d.tmp$start.B2)<=1000 ) & d.tmp$Gene[i]==as.character(d.tmp$Gene) & d.tmp$SVTYPE[i] ==d.tmp$SVTYPE| 
                        ( d.tmp$CHROM[i]==d.tmp$CHROM2 & d.tmp$CHROM2[i]==d.tmp$CHROM & 
                            abs(d.tmp$start.B1[i]-d.tmp$start.B2)<=1000 & abs(d.tmp$start.B2[i]-d.tmp$start.B1)<=1000 & d.tmp$Gene[i]==as.character(d.tmp$Gene) & d.tmp$SVTYPE[i] ==d.tmp$SVTYPE )  )) # find matches at 1kb distance
  ## table of unique alterations
  d.tmp.firsts = sapply(matchs.tmp, function(x) x[1])
  d.tmp$Parental = c("NO","YES")[sapply(matchs.tmp , function(x) "Parental"%in% unique(d.tmp[x]$time) )+1]
  d.tmp$PDTO1    = c("NO","YES")[sapply(matchs.tmp , function(x) "Organoid1"%in% unique(d.tmp[x]$time) )+1]
  if(any(d.tmp$time=="Organoid2")){
    d.tmp$PDTO2    = c("NO","YES")[sapply(matchs.tmp , function(x) "Organoid2"%in% unique(d.tmp[x]$time) )+1]
  }else{
    d.tmp$PDTO2    = NA
  }
  d.uniq = rbind(d.uniq, d.tmp[unique(d.tmp.firsts),] %>% dplyr::select(-Sample,-time) )
 
  # increment
  iexp = iexp +1
}



## recover from low-quality calls
missing = which(!apply(d.uniq[,13:15]=="YES",1,all,na.rm=T ))

### load intermediate SVaba data
vcf.SV.intermSVaba.files.indels        = list.files("/data/lungNENomics/work/organoids/WGS/SV_calling/release2_sv-somatic-cns-nf_02112021/SVABA/",pattern = "svaba.unfiltered.somatic.indel.vcf$",full.names = T,recursive = T)
vcf.SV.intermSVaba.files.sv        = list.files("/data/lungNENomics/work/organoids/WGS/SV_calling/release2_sv-somatic-cns-nf_02112021/SVABA/",pattern = "svaba.unfiltered.somatic.sv.vcf$",full.names = T,recursive = T)
vcf.SV.intermSVaba.names        = str_remove(list.files("/data/lungNENomics/work/organoids/WGS/SV_calling/release2_sv-somatic-cns-nf_02112021/SVABA/",pattern = "svaba.unfiltered.somatic.sv.vcf$", full.names = F,recursive = F),
                                             ".svaba.unfiltered.somatic.sv.vcf$")
vcf.SV.intermSVaba.indell = lapply(vcf.SV.intermSVaba.files.indels, read_tsv,comment = "#", col_names = c("CHROM",	"POS",	"ID",	"REF",	"ALT",	"QUAL",	"FILTER",	"INFO",	"FORMAT","normal",	"tumor"))
vcf.SV.intermSVaba.svl    = lapply(vcf.SV.intermSVaba.files.sv, read_tsv,comment = "#", col_names = c("CHROM",	"POS",	"ID",	"REF",	"ALT",	"QUAL",	"FILTER",	"INFO",	"FORMAT","normal",	"tumor"))

### load intermediate Manta data
vcf.SV.intermManta.files        = list.files("/data/lungNENomics/work/organoids/WGS/SV_calling/release2_sv-somatic-cns-nf_02112021/MANTA/",pattern = "somaticSV.vcf.gz$",full.names = T,recursive = T)
vcf.SV.intermManta.names        = str_remove(list.dirs("/data/lungNENomics/work/organoids/WGS/SV_calling/release2_sv-somatic-cns-nf_02112021/MANTA/",full.names = F,recursive = F),
                                             "_results")
vcf.SV.intermMantal = lapply(vcf.SV.intermManta.files, read_tsv,comment = "#", col_names = c("CHROM",	"POS",	"ID",	"REF",	"ALT",	"QUAL",	"FILTER",	"INFO",	"FORMAT","normal",	"tumor"))

vcf.SV.intermDelly.files        = list.files("/data/lungNENomics/work/organoids/WGS/SV_calling/release2_sv-somatic-cns-nf_02112021/DELLY/",pattern = "somatic.vcf",full.names = T)
vcf.SV.intermDelly.names        = str_remove(list.files("/data/lungNENomics/work/organoids/WGS/SV_calling/release2_sv-somatic-cns-nf_02112021/DELLY/",pattern = ".delly_somatic.vcf",full.names = F),
                                             ".delly_somatic.vcf")
vcf.SV.intermDellyl = lapply(vcf.SV.intermDelly.files, read_tsv,comment = "#", col_names = c("CHROM",	"POS",	"ID",	"REF",	"ALT",	"QUAL",	"FILTER",	"INFO",	"FORMAT","tumor",	"normal"))

for(i in 1:length(vcf.SV.intermDellyl)){
  vcf.SV.intermSVaba.svl[[i]]$Sample            = vcf.SV.intermSVaba.names[i]
  vcf.SV.intermSVaba.svl[[i]]$Experiment        = str_remove(vcf.SV.intermSVaba.names,"Tp[0-9]*$|Mp[0-9]$|T$|M$")[i]
  vcf.SV.intermSVaba.svl[[i]]$QUAL              = as.numeric(vcf.SV.intermSVaba.svl[[i]]$QUAL)
  
  vcf.SV.intermSVaba.indell[[i]]$Sample            = vcf.SV.intermSVaba.names[i]
  vcf.SV.intermSVaba.indell[[i]]$Experiment        = str_remove(vcf.SV.intermSVaba.names,"Tp[0-9]*$|Mp[0-9]$|T$|M$")[i]
  vcf.SV.intermSVaba.indell[[i]]$QUAL              = as.numeric(vcf.SV.intermSVaba.indell[[i]]$QUAL)
  vcf.SV.intermSVaba.indell[[i]]$ID              = as.character(vcf.SV.intermSVaba.indell[[i]]$ID)
  
  vcf.SV.intermMantal[[i]]$Sample            = vcf.SV.intermManta.names[i]
  vcf.SV.intermMantal[[i]]$Experiment        = str_remove(vcf.SV.intermManta.names,"Tp[0-9]*$|Mp[0-9]$|T$|M$")[i]
  vcf.SV.intermMantal[[i]]$QUAL              = as.numeric(vcf.SV.intermMantal[[i]]$QUAL)
  
  vcf.SV.intermDellyl[[i]]$Sample            = vcf.SV.intermDelly.names[i]
  vcf.SV.intermDellyl[[i]]$Experiment        = str_remove(vcf.SV.intermDelly.names,"Tp[0-9]*$|Mp[0-9]$|T$|M$")[i]
  vcf.SV.intermDellyl[[i]]$QUAL              = as.numeric(vcf.SV.intermDellyl[[i]]$QUAL)
}
dataset.intermSVaba.sv = bind_rows(c(vcf.SV.intermSVaba.svl,vcf.SV.intermSVaba.indell ) )
dataset.intermManta.sv = bind_rows(vcf.SV.intermMantal)
dataset.intermDelly.sv = bind_rows(vcf.SV.intermDellyl)

rm(vcf.SV.intermMantal)
rm(vcf.SV.intermDellyl)
rm(vcf.SV.intermSVaba.svl)
rm(vcf.SV.intermSVaba.indell)
gc()

samplesTime = as.data.frame(matrix(c("LNET2T","LNET2Tp12",NA,
                   "LNET5T",	"LNET5Tp4",NA,
                   "LNET6T",	"LNET6Tp1",NA,
                   "LNET10T",	"LNET10Tp4",NA,
                   "LCNEC3T","LCNEC3Tp17",NA,
                   "LCNEC4T",	"LCNEC4Tp7", "LCNEC4Tp24",	
                   "PANEC1T","PANEC1Tp4","PANEC1Tp14",
                   "SINET7M",	"SINET7Mp2",NA,
                   "SINET8M",	"SINET8Mp2",NA,
                   "SINET9M",	"SINET9Mp1",NA) , ncol=3,byrow=T,dimnames = list(c("LNET2" , "LNET5", "LNET6", "LNET10","LCNEC3","LCNEC4","PANEC1","SINET7","SINET8","SINET9"),
                                                                               c("Parental","PDTO1","PDTO2"))) )

dataset.intermSVaba.sv = dataset.intermSVaba.sv %>% filter(!str_detect(ALT,"_decoy|_alt|chrUn|_random"),!str_detect(CHROM,"_decoy|_alt|chrUn|_random"))

for( i in missing ){
  d.uniq.tmp = d.uniq[missing[i],]
  missing.samp.tmp = samplesTime[d.uniq.tmp$Experiment,which(as.matrix(d.uniq[missing[i],15:17])=="NO")]#35:37])=="NO")]
  # check original direction
  tmpD = dataset.intermDelly.sv   %>% filter( Sample %in% missing.samp.tmp, CHROM==d.uniq.tmp$CHROM & abs(d.uniq.tmp$start.B1-POS)<=1000 &
                                                (d.uniq.tmp$CHROM==d.uniq.tmp$CHROM2 | str_detect(ALT,paste0(d.uniq.tmp$CHROM2,":") )) | 
                                                CHROM==d.uniq.tmp$CHROM2 & abs(d.uniq.tmp$start.B2-POS)<=1000 & # reverse
                                                (d.uniq.tmp$CHROM==d.uniq.tmp$CHROM2 | str_detect(ALT,paste0(d.uniq.tmp$CHROM,":") ))
  )
  tmpM = dataset.intermManta.sv   %>% filter( Sample %in% missing.samp.tmp, CHROM==d.uniq.tmp$CHROM & abs(d.uniq.tmp$start.B1-POS)<=1000 &
                                                (d.uniq.tmp$CHROM==d.uniq.tmp$CHROM2 | str_detect(ALT,paste0(d.uniq.tmp$CHROM2,":") )) | 
                                                CHROM==d.uniq.tmp$CHROM2 & abs(d.uniq.tmp$start.B2-POS)<=1000 & # reverse
                                                (d.uniq.tmp$CHROM==d.uniq.tmp$CHROM2 | str_detect(ALT,paste0(d.uniq.tmp$CHROM,":") ))
  )
  tmpS = dataset.intermSVaba.sv   %>% filter( Sample %in% missing.samp.tmp, CHROM==d.uniq.tmp$CHROM & abs(d.uniq.tmp$start.B1-POS)<=1000 &
                                                (d.uniq.tmp$CHROM==d.uniq.tmp$CHROM2 | str_detect(ALT,paste0(d.uniq.tmp$CHROM2,":") )) | 
                                                CHROM==d.uniq.tmp$CHROM2 & abs(d.uniq.tmp$start.B2-POS)<=1000 & # reverse
                                                (d.uniq.tmp$CHROM==d.uniq.tmp$CHROM2 | str_detect(ALT,paste0(d.uniq.tmp$CHROM,":") ))
  )
  tmp = bind_rows(tmpD[,-c(10:11)],tmpM[,-c(10:11)],tmpS[,-c(10:11)])
  tmp$CHROM2 = str_remove(str_extract(tmp$ALT,"[\\[,\\]]chr[0-9MXY]+"),"\\[|\\]")
  tmp[is.na(tmp$CHROM2) & !str_detect(tmp$ALT,"chr"),]$CHROM2 = tmp[is.na(tmp$CHROM2) & !str_detect(tmp$ALT,"chr"),]$CHROM #same chr
  tmp$POS2   = as.numeric(str_remove(str_extract(tmp$ALT,":[0-9]+"),":"))
  tmp[is.na(tmp$POS2) & !str_detect(tmp$ALT,"chr"),]$POS2 = tmp[is.na(tmp$POS2) & !str_detect(tmp$ALT,"chr"),]$POS #same chr
  tmp$start.B1 = as.numeric(tmp$POS) + as.numeric(str_remove(str_extract(tmp$INFO,"CIPOS=[0-9\\-]+"),"CIPOS="))
  tmp$end.B1 = as.numeric(tmp$POS) + as.numeric(str_remove(str_extract(tmp$INFO,"CIPOS=[0-9\\-]+,[0-9\\-]+"),"CIPOS=[0-9\\-]+,"))
  tmp$start.B2 = as.numeric(str_remove(str_extract(tmp$INFO,";END=[0-9\\-]+"),";END=")) + as.numeric(str_remove(str_extract(tmp$INFO,"CIEND=[0-9\\-]+"),"CIEND="))
  tmp$end.B2 = as.numeric(str_remove(str_extract(tmp$INFO,";END=[0-9\\-]+"),";END=")) + as.numeric(str_remove(str_extract(tmp$INFO,"CIEND=[0-9\\-]+,[0-9\\-]+"),"CIEND=[0-9\\-]+,"))
  # svaba SVs
  tmp[is.na(tmp$start.B1) & str_detect(tmp$ALT,"chr"),]$start.B1 =  tmp[is.na(tmp$start.B1) & str_detect(tmp$ALT,"chr"),]$POS
  tmp[is.na(tmp$start.B2) & str_detect(tmp$ALT,"chr"),]$start.B2 =  tmp[is.na(tmp$start.B2) & str_detect(tmp$ALT,"chr"),]$POS2
  # svaba indels
  tmp[is.na(tmp$start.B1) & !str_detect(tmp$ALT,"chr"),]$start.B1 = tmp[is.na(tmp$start.B1) & !str_detect(tmp$ALT,"chr"),]$POS
  tmp[is.na(tmp$start.B2) & !str_detect(tmp$ALT,"chr"),]$start.B2 = tmp[is.na(tmp$start.B2) & !str_detect(tmp$ALT,"chr"),]$POS + as.numeric(str_remove(str_extract(tmp[is.na(tmp$start.B2) & !str_detect(tmp$ALT,"chr"),]$INFO,";SPAN=[0-9\\-]+"),";SPAN="))
  
  tmp = tmp %>% filter( (d.uniq.tmp$CHROM==CHROM & abs(d.uniq.tmp$start.B1-start.B1)<=1000 & d.uniq.tmp$CHROM2==CHROM2  & abs(d.uniq.tmp$start.B2-start.B2)<=1000 ) | 
                          (d.uniq.tmp$CHROM==CHROM2 & abs(d.uniq.tmp$start.B1-start.B2)<=1000 & d.uniq.tmp$CHROM2==CHROM  & abs(d.uniq.tmp$start.B2-start.B1)<=1000)  )
  if(nrow(tmp)>0) { 
    colID = which(samplesTime[d.uniq.tmp$Experiment,] %in% unique(tmp$Sample))
    if(any(colID==1) ) d.uniq[missing[i],15] = list("Recovered") #"Recovered" Tumor
    if(any(colID==2) ) d.uniq[missing[i],16] = list("Recovered") #"Recovered" PDTO1
    if(any(colID==3) ) d.uniq[missing[i],17] = list("Recovered") #"Recovered" PDTO2
  }
}

#SV_uniq_recovered = d.uniq
driver_uniq_recovered = d.uniq
#save(SV_uniq_recovered , file = "/data/lungNENomics/work/organoids/WGS/SV_calling/release2_sv-somatic-cns-nf_02112021/SURVIVOR/SV_uniq_recovered.Rdata")
save(driver_uniq_recovered , file = "/data/lungNENomics/work/organoids/WGS/SV_calling/release2_sv-somatic-cns-nf_02112021/SURVIVOR/driver_uniq_recovered.Rdata")

write_tsv(SV_uniq_recovered, "/data/lungNENomics/work/organoids/figures/TableS4_SVs.tsv")


## add which break includes driver
load("/data/lungNENomics/work/organoids/WGS/SV_calling/release2_sv-somatic-cns-nf_02112021/SURVIVOR/SV_uniq_recovered.Rdata")
load("/data/lungNENomics/work/organoids/WGS/SV_calling/release2_sv-somatic-cns-nf_02112021/SURVIVOR/driver_uniq_recovered.Rdata")

driver_uniq_recovered$Gene = sapply( driver_uniq_recovered$Gene, function(x){res=x;if(length(x)==0){res=NA};return(res)})
SV_uniq_recovered$Gene = sapply(SV_uniq_recovered$Gene,function(x){res=x;if(length(x)==0){res=NA};return(res)})

# add correct drivers
SVs2 = full_join(SV_uniq_recovered,driver_uniq_recovered, by=c("type","SVTYPE", "CHROM", "start.B1", "end.B1", "CHROM2", "start.B2", "end.B2", 
                                             "Experiment","CN"))

# fill common columns
SVs2[which(SVs2$Gene.x!=SVs2$Gene.y),] %>% dplyr::select(ID.x:PDTO2.y,Gene.y)
SVs2$Gene = NA
SVs2$Gene[!is.na(SVs2$Gene.x)] = SVs2$Gene.x[!is.na(SVs2$Gene.x)]
SVs2$Gene[!is.na(SVs2$Gene.y)] = SVs2$Gene.y[!is.na(SVs2$Gene.y)]
SVs2$ID = NA
SVs2$ID[!is.na(SVs2$ID.x)] = SVs2$ID.x[!is.na(SVs2$ID.x)]
SVs2$ID[!is.na(SVs2$ID.y)] = SVs2$ID.y[!is.na(SVs2$ID.y)]

SVs2$Parental = NA
SVs2$Parental[!is.na(SVs2$Parental.x)] = SVs2$Parental.x[!is.na(SVs2$Parental.x)]
SVs2$Parental[!is.na(SVs2$Parental.y)] = SVs2$Parental.y[!is.na(SVs2$Parental.y)]
SVs2$PDTO1 = NA
SVs2$PDTO1[!is.na(SVs2$PDTO1.x)] = SVs2$PDTO1.x[!is.na(SVs2$PDTO1.x)]
SVs2$PDTO1[!is.na(SVs2$PDTO1.y)] = SVs2$PDTO1.y[!is.na(SVs2$PDTO1.y)]
SVs2$PDTO2 = NA
SVs2$PDTO2[!is.na(SVs2$PDTO2.x)] = SVs2$PDTO2.x[!is.na(SVs2$PDTO2.x)]
SVs2$PDTO2[!is.na(SVs2$PDTO2.y)] = SVs2$PDTO2.y[!is.na(SVs2$PDTO2.y)]


whichBl = sapply(1:nrow(SVs2) , function(i){
  geneinfo = ref.exon[ref.exon$geneSymbol==SVs2$Gene[i],]
  if(nrow(geneinfo)==0) return(NA)
  whichB = which.min( c( c(NA,1)[max(SVs2$CHROM[i]==paste0("chr",geneinfo$seqnames))+1]*min( c(abs(SVs2$start.B1[i]-geneinfo$start) , abs(SVs2$start.B1[i]-geneinfo$end) , abs(SVs2$end.B1[i]-geneinfo$start) , abs(SVs2$end.B1[i]-geneinfo$end) ) ),
                         c(NA,1)[max(SVs2$CHROM2[i]==paste0("chr",geneinfo$seqnames))+1]*min( c(abs(SVs2$start.B2[i]-geneinfo$start) , abs(SVs2$start.B2[i]-geneinfo$end) , abs(SVs2$end.B2[i]-geneinfo$start) , abs(SVs2$end.B2[i]-geneinfo$end) ) ) ) )
  return(whichB)} )
SVs2$Gene_breakpoint = whichBl

SVs2$Driver = !is.na(SVs2$Gene.y)

write_tsv(SVs2 %>% dplyr::select(-ID.x,-ID.y),"SVs_wbp_wdriv.tsv")
}