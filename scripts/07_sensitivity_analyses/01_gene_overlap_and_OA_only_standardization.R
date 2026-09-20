#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

project_dir <- normalizePath(getwd(), mustWork = TRUE)
out_dir <- file.path(project_dir, "results/frozen_summary_tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
s2c <- file.path(project_dir, "04_结果/Stage2C独立验证")

frozen <- read.csv(file.path(project_dir, "结果/Stage1B增量验证/frozen_program_manifest.csv"), check.names = FALSE)
program_names <- frozen$program
program_genes <- setNames(lapply(frozen$genes, function(x) unique(trimws(strsplit(x, ";", fixed = TRUE)[[1]]))), program_names)
ra_union <- unique(c(program_genes$frozen_ra_up, program_genes$frozen_ra_down))
comparators <- c("general_inflammation", "mhc_ii_apc_primary", "mhc_ii_apc_generic", "myeloid_context", "interferon_state", "t_cell_context", "b_cell_context", "fibroblast_ecm", "osteoclast_context")
primary5 <- c("general_inflammation", "mhc_ii_apc_primary", "mhc_ii_apc_generic", "myeloid_context", "interferon_state")

tpm <- read.delim(gzfile(file.path(project_dir, "04_结果/Stage2B基因矩阵/GSE283079_gene_tpm.tsv.gz")), check.names = FALSE, row.names = 1)
tpm <- as.matrix(tpm); storage.mode(tpm) <- "numeric"
sm <- read.csv(file.path(project_dir, "04_结果/Stage2B技术QC/GSE283079_analysis_sample_manifest_frozen.csv"), check.names = FALSE)
sm <- sm[match(colnames(tpm), sm$SRR), ]
sm$group <- ifelse(sm$group %in% c("non_OA", "non-OA"), "non-OA", "OA")
keep <- sm$final_quant_valid & sm$final_decision %in% c("KEEP", "QC-FLAG-BUT-KEEP")
oa <- sm$group == "OA" & keep
stopifnot(sum(oa) == 36L, sum(keep) == 41L)

tx <- read.delim(gzfile(file.path(project_dir, "04_结果/Stage2B基因矩阵/GENCODE_v47_tx2gene.tsv.gz")), check.names = FALSE)
tx$gene_id_clean <- sub("\\..*$", "", tx$gene_id); tx$gene_symbol <- trimws(tx$gene_symbol)
tx <- tx[!duplicated(tx$gene_id_clean), c("gene_id_clean", "gene_symbol")]
gid <- sub("\\..*$", "", rownames(tpm)); symbol <- tx$gene_symbol[match(gid, tx$gene_id_clean)]
symbol[is.na(symbol) | symbol == ""] <- gid[is.na(symbol) | symbol == ""]

make_z <- function(x, cols) { x <- x[, cols, drop=FALSE]; mu <- rowMeans(x, na.rm=TRUE); s <- apply(x,1,sd,na.rm=TRUE); s[!is.finite(s)|s==0] <- NA_real_; (x-mu)/s }
score <- function(z, sets) {
  rownames(z) <- symbol; out <- as.data.frame(matrix(NA_real_, nrow=ncol(z), ncol=length(sets)), check.names=FALSE); names(out) <- names(sets)
  for (p in names(sets)) { ii <- which(rownames(z) %in% sets[[p]]); if(length(ii)) out[[p]] <- colMeans(z[ii,,drop=FALSE], na.rm=TRUE) }
  out$ra_projection <- out$frozen_ra_up - out$frozen_ra_down; out
}
log_tpm <- log2(tpm+1)
primary_scores <- score(make_z(log_tpm, seq_len(ncol(tpm))), program_genes)
oa_scores <- primary_scores[oa,,drop=FALSE]
oa_only_scores <- score(make_z(log_tpm, which(oa)), program_genes)
safe_cor <- function(x,y) suppressWarnings(cor(x,y,method="spearman",use="complete.obs"))

# Analysis 1: remove RA_union from comparators only.
excluded_sets <- program_genes
overlap_info <- do.call(rbind, lapply(comparators, function(p) {
  ov <- intersect(program_genes[[p]], ra_union); excluded_sets[[p]] <<- setdiff(program_genes[[p]], ra_union)
  data.frame(comparator=p, original_gene_count=length(program_genes[[p]]), overlap_n=length(ov), overlap_percent=100*length(ov)/length(program_genes[[p]]), post_exclusion_gene_count=length(excluded_sets[[p]]), too_few_genes=length(excluded_sets[[p]])<3)
}))
excluded_scores <- score(make_z(log_tpm, seq_len(ncol(tpm))), excluded_sets)
old_gse <- read.csv(file.path(s2c,"GSE283079_program_correlations.csv"), check.names=FALSE)
old_gse <- old_gse[match(comparators, old_gse$program),]
new_gse <- vapply(comparators, function(p) safe_cor(oa_scores$ra_projection, excluded_scores[oa,p]), numeric(1))
overlap_info$original_GSE283079_rho <- old_gse$spearman_rho; overlap_info$overlap_excluded_GSE283079_rho <- new_gse

rels <- read.csv(file.path(s2c,"OA_crosscohort_program_relationships.csv"), check.names=FALSE)
old_meta <- read.csv(file.path(s2c,"OA_crosscohort_meta_analysis.csv"), check.names=FALSE)
rels_excl <- rels
for(p in primary5) rels_excl$spearman_rho[rels_excl$cohort=="GSE283079" & rels_excl$program==p] <- new_gse[match(p,comparators)]
dl <- function(d) { d<-d[is.finite(d$spearman_rho)&d$n>3,,drop=FALSE]; z<-atanh(pmax(pmin(d$spearman_rho,.999999),-.999999)); v<-1/(d$n-3); w<-1/v; m0<-sum(w*z)/sum(w); q<-sum(w*(z-m0)^2); c0<-sum(w)-sum(w^2)/sum(w); tau<-max(0,(q-(nrow(d)-1))/c0); wr<-1/(v+tau); m<-sum(wr*z)/sum(wr); se<-sqrt(1/sum(wr)); data.frame(k=nrow(d),pooled_rho=tanh(m),ci_low=tanh(m-1.96*se),ci_high=tanh(m+1.96*se),tau2=tau,i2=ifelse(q>0,max(0,(q-(nrow(d)-1))/q)*100,0)) }
reml_tau <- function(y,v) { nll<-function(tau) { w<-1/(v+tau); m<-sum(w*y)/sum(w); .5*(sum(log(v+tau))+log(sum(w))+sum(w*(y-m)^2)) }; max(0,optimize(nll,c(0,max(10,max(y^2)+1)),tol=1e-12)$minimum) }
hk <- function(d) { d<-d[is.finite(d$spearman_rho)&d$n>3,,drop=FALSE]; y<-atanh(pmax(pmin(d$spearman_rho,.999999),-.999999)); v<-1/(d$n-3); k<-length(y); tau<-reml_tau(y,v); w<-1/(v+tau); m<-sum(w*y)/sum(w); q<-sum(w*(y-m)^2); hv<-q/((k-1)*sum(w)); cr<-qt(.975,k-1); data.frame(REML_tau2=tau,REML_HK_pooled_rho=tanh(m),REML_HK_ci_low=tanh(m-cr*sqrt(hv)),REML_HK_ci_high=tanh(m+cr*sqrt(hv)),HK_df=k-1,HK_q=q) }
meta_excl <- do.call(rbind,lapply(primary5,function(p)cbind(comparator=p,dl(rels_excl[rels_excl$program==p,]))))
names(meta_excl)[names(meta_excl)=="pooled_rho"] <- "overlap_excluded_pooled_rho"; names(meta_excl)[names(meta_excl)=="ci_low"] <- "overlap_excluded_pooled_ci_low"; names(meta_excl)[names(meta_excl)=="ci_high"] <- "overlap_excluded_pooled_ci_high"; names(meta_excl)[names(meta_excl)=="tau2"] <- "overlap_excluded_tau2"; names(meta_excl)[names(meta_excl)=="i2"] <- "overlap_excluded_i2"
overlap_info$original_pooled_rho <- old_meta$pooled_rho[match(comparators,old_meta$program)]
overlap_info <- merge(overlap_info,meta_excl,by="comparator",all.x=TRUE,sort=FALSE); overlap_info <- overlap_info[match(comparators,overlap_info$comparator),]
overlap_info$direction_preserved <- sign(overlap_info$original_pooled_rho)==sign(overlap_info$overlap_excluded_pooled_rho)
overlap_info$interpretive_hierarchy_preserved <- overlap_info$direction_preserved & abs(overlap_info$overlap_excluded_pooled_rho-overlap_info$original_pooled_rho)<.15
write.table(overlap_info,file.path(out_dir,"Supplementary_Table_S10_Gene_Overlap_Sensitivity.tsv"),sep="\t",quote=FALSE,row.names=FALSE,na="")
hk_excl <- do.call(rbind,lapply(primary5,function(p)cbind(comparator=p,hk(rels_excl[rels_excl$program==p,]))))
write.table(hk_excl,file.path(out_dir,"Supplementary_Table_S10_primary5_REML_HK.tsv"),sep="\t",quote=FALSE,row.names=FALSE,na="")

# Analysis 2: OA-only gene-wise standardization.
oa_rel <- vapply(comparators,function(p)safe_cor(oa_only_scores$ra_projection[oa],oa_only_scores[oa,p]),numeric(1))
old9 <- old_gse
s11 <- data.frame(programme_relationship=c("RA projection concordance",comparators),primary_estimate=c(NA,old9$spearman_rho),oa_only_estimate=c(safe_cor(primary_scores$ra_projection[oa],oa_only_scores$ra_projection[oa]),oa_rel),stringsAsFactors=FALSE)
s11$absolute_change <- c(NA,abs(oa_rel-old9$spearman_rho)); s11$direction_preserved <- c(NA,sign(old9$spearman_rho)==sign(oa_rel)); near <- function(x) ifelse(abs(x)<.1,"near-zero",ifelse(x>0,"positive","negative")); s11$qualitative_interpretation_preserved <- c(TRUE,near(old9$spearman_rho)==near(oa_rel))
write.table(s11,file.path(out_dir,"Supplementary_Table_S11_OA_Only_Standardization.tsv"),sep="\t",quote=FALSE,row.names=FALSE,na="")
write.table(data.frame(primary_41_sample_ra_projection=primary_scores$ra_projection[oa],oa_only_ra_projection=oa_only_scores$ra_projection[oa]),file.path(out_dir,"Supplementary_Table_S11_RA_projection_pairs.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
summary <- data.frame(metric=c("ra_projection_41_vs_oa_only_spearman","primary5_direction_preserved","primary5_hierarchy_preserved"),value=c(safe_cor(primary_scores$ra_projection[oa],oa_only_scores$ra_projection[oa]),all(overlap_info$direction_preserved[overlap_info$comparator%in%primary5]),all(overlap_info$interpretive_hierarchy_preserved[overlap_info$comparator%in%primary5])))
write.table(summary,file.path(out_dir,"v36_sensitivity_summary.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
print(overlap_info); print(s11); print(summary)
