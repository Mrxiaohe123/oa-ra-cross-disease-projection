suppressPackageStartupMessages({library(readxl); library(data.table); library(ggplot2); library(patchwork); library(ragg)})
root <- normalizePath(getwd(), mustWork = TRUE)
pkg <- root
src <- file.path(root,'data','frozen_supplementary')
out <- file.path(root,'results','figure_exports')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
pal <- list(red='#C7473A',blue='#247BA0',amber='#C98314',purple='#6C54A3',teal='#238B8D',charcoal='#4B5563',grey='#B7C1CC')
theme_pub <- theme_minimal(base_family='Arial',base_size=9)+theme(panel.grid.minor=element_blank(),panel.grid.major=element_line(colour='#E5E7EB',linewidth=.25),strip.background=element_rect(fill='#F3F4F6',colour=NA),strip.text=element_text(face='bold',colour='#172033'),axis.title=element_text(face='bold',colour='#172033'),axis.text=element_text(colour='#374151'),plot.title=element_text(face='bold',colour='#172033'),legend.position='bottom')
save_fig <- function(p,stem,w=11,h=7){ggsave(file.path(out,paste0(stem,'.tif')),p,device=ragg::agg_tiff,width=w,height=h,units='in',dpi=600,compression='lzw',bg='white'); grDevices::pdf(file.path(out,paste0(stem,'.pdf')),width=w,height=h,onefile=TRUE); print(p); dev.off()}

# S1: expanded technical QC from frozen Table S2 only.
qc <- as.data.frame(read_excel(file.path(src,'Table_S2_GSE283079_technical_qc.xlsx'),sheet='table'))
qc$Run <- as.character(qc$Run); qc$Group <- as.character(qc$Group); qc$flagged <- qc$Run=='SRR31542944'
qc$`Mapping rate (%)` <- as.numeric(qc$`Mapping rate (%)`); qc$`Library size` <- as.numeric(qc$`Library size`); qc$`Detected genes` <- as.numeric(qc$`Detected genes`); qc$`Median sample correlation` <- as.numeric(qc$`Median sample correlation`)
p1 <- ggplot(qc,aes(x=reorder(Run,`Mapping rate (%)`),y=`Mapping rate (%)`,colour=Group))+geom_point(size=1.8)+geom_point(data=subset(qc,flagged),shape=21,fill='white',size=3,stroke=1.1,colour=pal$red)+coord_flip()+scale_colour_manual(values=c(OA=pal$blue,`non-OA`=pal$charcoal),drop=FALSE)+labs(x=NULL,y='Mapping rate (%)')+theme_pub+theme(axis.text.y=element_text(size=5.5),legend.position='none')
p2 <- ggplot(qc,aes(x=`Library size`,y=`Detected genes`,colour=Group))+geom_point(size=2)+geom_point(data=subset(qc,flagged),shape=21,fill='white',size=3,stroke=1.1,colour=pal$red)+scale_x_continuous(labels=scales::comma)+scale_colour_manual(values=c(OA=pal$blue,`non-OA`=pal$charcoal),drop=FALSE)+labs(x='Library size (counts)',y='Detected genes')+theme_pub
p3 <- ggplot(qc,aes(x=reorder(Run,`Median sample correlation`),y=`Median sample correlation`,colour=Group))+geom_point(size=1.8)+geom_point(data=subset(qc,flagged),shape=21,fill='white',size=3,stroke=1.1,colour=pal$red)+coord_flip()+scale_colour_manual(values=c(OA=pal$blue,`non-OA`=pal$charcoal),drop=FALSE)+labs(x=NULL,y='Median sample correlation')+theme_pub+theme(axis.text.y=element_text(size=5.5),legend.position='none')
sumdf <- data.frame(label=c('RNA-seq runs','OA runs','non-OA runs','QC-flag-but-keep'),value=c(nrow(qc),sum(qc$Group=='OA'),sum(qc$Group!='OA'),sum(qc$flagged)))
p4 <- ggplot(sumdf,aes(x=label,y=value,fill=label))+geom_col(width=.65,show.legend=FALSE)+geom_text(aes(label=value),vjust=-.3,fontface='bold')+scale_fill_manual(values=c(pal$blue,pal$blue,pal$charcoal,pal$red))+labs(x=NULL,y='Count')+theme_pub+theme(axis.text.x=element_text(angle=25,hjust=1))
s1 <- (p1|p2)/(p3|p4)+plot_annotation(title='Supplementary Figure S1. Expanded technical QC of the GSE283079 reconstruction',tag_levels='a')&theme(plot.title=element_text(face='bold',colour='#172033',size=12))
save_fig(s1,'Fig_S1',12,8)

# S2: frozen scoring sensitivity correlations; no underlying points are invented.
s2df <- data.frame(comparison=c('Rank-based score','Direction-control score','OA-only re-standardization'),rho=c(.944,.255,.9989704),label=c('rho = 0.944','rho = 0.255','rho = 0.9989704'))
p5 <- ggplot(s2df,aes(x=rho,y=reorder(comparison,rho)))+geom_vline(xintercept=0,colour='#9CA3AF',linewidth=.35)+geom_segment(aes(x=0,xend=rho,yend=reorder(comparison,rho)),colour=pal$blue,linewidth=2)+geom_point(size=3.5,colour=pal$blue)+geom_text(aes(label=label),hjust=-.15,size=3.2,colour='#172033')+coord_cartesian(xlim=c(-.05,1.08),clip='off')+labs(x='Spearman rho with primary RA projection',y=NULL)+theme_pub+theme(plot.margin=margin(5.5,35,5.5,5.5))
p6 <- ggplot(s2df,aes(x=seq_along(comparison),y=rho))+geom_hline(yintercept=0,colour='#9CA3AF',linewidth=.35)+geom_col(fill=pal$purple,width=.55)+geom_text(aes(label=label),vjust=-.35,size=3.2)+scale_x_continuous(breaks=1:3,labels=c('rank-based','direction-control','OA-only\nstandardization'))+coord_cartesian(ylim=c(0,1.12))+labs(x=NULL,y='Spearman rho')+theme_pub
p7 <- ggplot(data.frame(x=0:1,y=0:1),aes(x,y))+geom_blank()+annotate('text',x=.05,y=.78,hjust=0,label='Interpretation of the frozen comparisons',fontface='bold',size=4,colour='#172033')+annotate('text',x=.05,y=.58,hjust=0,label='The OA-only result is a scoring-standardization\nsensitivity result, not cross-platform reproducibility.',size=3.4,colour='#374151')+annotate('text',x=.05,y=.25,hjust=0,label='n = 36 OA samples in GSE283079\nNo new scores or correlations were generated.',size=3.4,colour='#374151')+theme_void()
s2 <- (p5|p6)/p7+plot_annotation(title='Supplementary Figure S2. Scoring robustness of the RA-derived projection',tag_levels='a')&theme(plot.title=element_text(face='bold',colour='#172033',size=12))
save_fig(s2,'Fig_S2',12,7)

# S3: secondary programme relationships from frozen Table S3.
s3 <- as.data.frame(read_excel(file.path(src,'Table_S3_programme_correlations.xlsx'),sheet='table'))
s3$`Spearman rho (rho)` <- as.numeric(s3$`Spearman rho (rho)`); s3$N <- as.numeric(s3$N)
s3$Cohort <- factor(as.character(s3$Cohort),levels=c('GSE89408','GSE55235','GSE55457','GSE55584','GSE82107','GSE206848','GSE283079'))
keep <- c('interferon_state','t_cell_context','b_cell_context','fibroblast_ecm','osteoclast_context')
s3m <- subset(s3,Programme %in% keep)
labs <- c(interferon_state='IFN',t_cell_context='T-cell',b_cell_context='B-cell',fibroblast_ecm='Fibroblast/ECM',osteoclast_context='Osteoclast')
s3m$Programme <- factor(s3m$Programme,levels=keep,labels=labs[keep])
p8 <- ggplot(s3m,aes(x=Cohort,y=Programme,fill=`Spearman rho (rho)`))+geom_tile(colour='white',linewidth=.35)+geom_text(aes(label=sprintf('%.3f',`Spearman rho (rho)`)),size=2.6)+scale_fill_gradient2(low='#3B6EA8',mid='white',high='#C45A4A',midpoint=0,limits=c(-1,1),name='rho')+labs(x=NULL,y=NULL)+theme_pub+theme(axis.text.x=element_text(angle=35,hjust=1),axis.text.y=element_text(face='bold'))
ns <- unique(s3m[,c('Cohort','N')]); ns$Cohort <- factor(ns$Cohort,levels=levels(s3m$Cohort))
p9 <- ggplot(ns,aes(x=Cohort,y=N))+geom_col(width=.65,fill=pal$teal)+geom_text(aes(label=N),vjust=-.25,fontface='bold')+labs(x=NULL,y='Samples per cohort')+theme_pub+theme(axis.text.x=element_text(angle=35,hjust=1))
s3fig <- p8/p9+plot_annotation(title='Supplementary Figure S3. Secondary programme relationships across OA cohorts',tag_levels='a')&theme(plot.title=element_text(face='bold',colour='#172033',size=12))
save_fig(s3fig,'Fig_S3',12,8)

# S4: frozen leave-one-cohort-out outputs for the five primary relationships.
s4 <- as.data.frame(read_excel(file.path(src,'Table_S5_leave_one_out.xlsx'),sheet='table'))
s4$`Pooled rho (rho)` <- as.numeric(s4$`Pooled rho (rho)`); s4$tau2 <- as.numeric(s4$tau2); s4$`I2 (%)` <- as.numeric(s4$`I2 (%)`)
s4$Programme <- as.character(s4$Programme); s4$`Omitted cohort` <- as.character(s4$`Omitted cohort`)
s4 <- subset(s4,Programme %in% c('general_inflammation','mhc_ii_apc_primary','mhc_ii_apc_generic','myeloid_context','interferon_state'))
s4 <- s4[!duplicated(s4[,c('Programme','Omitted cohort')]),]
s4$lo <- as.numeric(sub('–.*','',s4$`95% CI`)); s4$hi <- as.numeric(sub('.*–','',s4$`95% CI`)); s4$hi[is.na(s4$hi)] <- s4$`Pooled rho (rho)`[is.na(s4$hi)]
s4$rel <- factor(s4$Programme,levels=c('general_inflammation','mhc_ii_apc_primary','mhc_ii_apc_generic','myeloid_context','interferon_state'),labels=c('a  RA–general inflammation','b  RA–APC','c  RA–MHC-II','d  RA–myeloid','e  RA–IFN'))
s4$omit <- factor(s4$`Omitted cohort`,levels=c('none','GSE89408','GSE283079','GSE206848','GSE55235','GSE55457','GSE55584','GSE82107'),labels=c('All cohorts','GSE89408','GSE283079','GSE206848','GSE55235','GSE55457','GSE55584','GSE82107'))
s4$col <- ifelse(s4$omit=='GSE89408',pal$red,ifelse(s4$omit=='GSE283079',pal$blue,pal$charcoal))
s4fig <- ggplot(s4,aes(y=omit,x=`Pooled rho (rho)`))+geom_vline(xintercept=0,colour='#9CA3AF',linewidth=.35)+geom_errorbarh(aes(xmin=lo,xmax=hi),height=.18,colour='#9CA3AF',linewidth=.5)+geom_point(aes(colour=col),size=2.1)+facet_wrap(~rel,ncol=1,scales='free_y',strip.position='top')+scale_colour_identity()+labs(x='Pooled Spearman rho with 95% CI',y=NULL)+theme_pub+theme(strip.text=element_text(hjust=0),axis.text.y=element_text(size=7),panel.spacing=unit(.35,'lines'))+coord_cartesian(xlim=c(-.9,.9))
s4fig <- s4fig+plot_annotation(title='Supplementary Figure S4. Leave-one-cohort-out robustness of primary relationships')&theme(plot.title=element_text(face='bold',colour='#172033',size=12))
save_fig(s4fig,'Fig_S4',10,12)

fwrite(data.table(Figure=paste0('Fig_S',1:4),TIFF=file.path(out,paste0('Fig_S',1:4,'.tif')),VectorPDF=file.path(out,paste0('Fig_S',1:4,'.pdf'))),file.path(out,'supplementary_figure_manifest.tsv'),sep='\t')
