options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(e1071);library(GSVA);library(msigdbr);library(ggplot2);library(limma)})
root <- "C:/temp/ibi_tool"
input_dir <- file.path(root,"revision_supp_data")
output_dir <- file.path(root,"Result","revision_GSE51404_validation_pm1")
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)

read_series <- function(path){
 lines <- readLines(path,warn=FALSE)
 b <- grep("^!series_matrix_table_begin",lines); e <- grep("^!series_matrix_table_end",lines)
 stopifnot(length(b)==1L,length(e)==1L)
 x <- read.delim(path,skip=b,nrows=e-b-2L,row.names=1,check.names=FALSE,quote="\"",comment.char="")
 x <- as.matrix(x);storage.mode(x)<-"double"
 h <- lines[seq_len(min(100L,length(lines)))]
 field <- function(prefix){
  z <- strsplit(h[grep(prefix,h)[1L]],"\t",fixed=TRUE)[[1L]][-1L]
  sub('^"(.*)"$',"\\1",z)
 }
 ch <- h[grep("^!Sample_characteristics_ch1",h)]
 characteristic <- function(key){
  z <- ch[grep(paste0('^!Sample_characteristics_ch1\\t"',key,':'),ch)]
  z <- strsplit(z,"\t",fixed=TRUE)[[1L]][-1L]
  trimws(sub(paste0("^",key,":"),"",sub('^"(.*)"$',"\\1",z),ignore.case=TRUE))
 }
 meta <- data.frame(sample_id=field("^!Sample_geo_accession"),sample_title=field("^!Sample_title"),
  illness=characteristic("illness"),gender=characteristic("gender"),
  age=as.numeric(characteristic("age")))
 stopifnot(identical(colnames(x),meta$sample_id))
 list(expr=x,meta=meta)
}
map_genes <- function(x,a,probe_col){
 probe <- trimws(a[[probe_col]]);gene <- toupper(trimws(a[["Symbol"]]))
 gene <- gene[match(rownames(x),probe)];keep <- !is.na(gene)&gene!=""
 z <- avereps(x[keep,,drop=FALSE],ID=gene[keep])
 z[apply(z,1,function(v)all(is.finite(v))),,drop=FALSE]
}
make_pairs <- function(x,pairs){
 z <- matrix(-1L,ncol(x),length(pairs),dimnames=list(colnames(x),pairs))
 present <- logical(length(pairs))
 for(i in seq_along(pairs)){
  g <- strsplit(pairs[i],"_vs_",fixed=TRUE)[[1L]]
  present[i] <- all(g%in%rownames(x))
  if(present[i])z[,i] <- ifelse(x[g[1],]>x[g[2],],1L,-1L)
 }
 list(x=z,present=present)
}
kernel_name <- function(x)c("linear","polynomial","radial","sigmoid")[as.integer(x)+1L]

load(file.path(root,"best_svr0731.rda"));load(file.path(root,"X0723.rda"));load(file.path(root,"y0723.rda"))
old_model <- best_svr;features <- colnames(old_model$SV)
stopifnot(length(features)==390L,all(features%in%colnames(X)))
stopifnot(all(sort(unique(as.vector(as.matrix(X[,features]))))%in%c(-1,1)))
X_original_pm1 <- as.matrix(X[,features,drop=FALSE])
train <- read_series(file.path(input_dir,"GSE51405_series_matrix.txt"))
valid <- read_series(file.path(input_dir,"GSE51404_series_matrix.txt"))
a3 <- read.delim(file.path(input_dir,"GPL6947_HumanHT-12_V3_0_R1_11283641_A.bgx"),
                 skip=8L,quote="",check.names=FALSE,fill=TRUE)
a4 <- read.delim(file.path(input_dir,"GPL10558-50081.txt"),
                 comment.char="#",quote="",check.names=FALSE,fill=TRUE)
train_expr <- map_genes(train$expr,a3,"Probe_Id");valid_expr <- map_genes(valid$expr,a4,"ID")
train_pairs <- make_pairs(train_expr,features);valid_pairs <- make_pairs(valid_expr,features)
coverage <- data.frame(dataset=c("GSE51405 training set","GSE51404 validation set"),
 samples=c(ncol(train_expr),ncol(valid_expr)),mapped_genes=c(nrow(train_expr),nrow(valid_expr)),
 gene_pairs_present=c(sum(train_pairs$present),sum(valid_pairs$present)),
 gene_pairs_missing=c(sum(!train_pairs$present),sum(!valid_pairs$present)))
write.csv(coverage,file.path(output_dir,"platform_and_feature_coverage.csv"),row.names=FALSE)
if(any(!train_pairs$present)||any(!valid_pairs$present))warning("Missing model gene pairs; see coverage file")

train_label <- ifelse(tolower(train$meta$illness)=="healthy",0,-1)
x_aug <- rbind(X_original_pm1,train_pairs$x[,features,drop=FALSE])
y_aug <- c(as.numeric(y),train_label)
updated_svr <- svm(x=x_aug,y=y_aug,type="eps-regression",kernel=kernel_name(old_model$kernel),
 degree=old_model$degree,gamma=old_model$gamma,coef0=old_model$coef0,cost=old_model$cost,
 nu=old_model$nu,epsilon=old_model$epsilon,scale=TRUE)
save(updated_svr,features,file=file.path(output_dir,"updated_svr_with_GSE51405_pm1.rda"))

valid$meta$clinical_group <- ifelse(grepl("_XLA$",valid$meta$sample_title),"XLA",
 ifelse(valid$meta$illness=="CVID no complications","CVID without complications",
 ifelse(valid$meta$illness=="CVID with complications","CVID with complications","Healthy")))
valid$meta$clinical_group <- factor(valid$meta$clinical_group,
 levels=c("XLA","CVID without complications","CVID with complications","Healthy"))
if(anyNA(valid$meta$clinical_group))stop("Unrecognized illness label")
valid$meta$IBI <- as.numeric(predict(updated_svr,valid_pairs$x[,features,drop=FALSE]))
write.csv(valid$meta,file.path(output_dir,"GSE51404_sample_level_IBI.csv"),row.names=FALSE)
summary_table <- do.call(rbind,lapply(split(valid$meta,valid$meta$clinical_group),function(d)
 data.frame(clinical_group=as.character(d$clinical_group[1]),n=nrow(d),mean=mean(d$IBI),
 sd=sd(d$IBI),median=median(d$IBI),q1=unname(quantile(d$IBI,.25)),
 q3=unname(quantile(d$IBI,.75)),minimum=min(d$IBI),maximum=max(d$IBI))))
rownames(summary_table)<-NULL
write.csv(summary_table,file.path(output_dir,"GSE51404_IBI_group_summary.csv"),row.names=FALSE)
kw <- kruskal.test(IBI~clinical_group,data=valid$meta)
healthy <- valid$meta$IBI[valid$meta$clinical_group=="Healthy"]
tests <- do.call(rbind,lapply(levels(valid$meta$clinical_group)[1:3],function(g){
 disease <- valid$meta$IBI[valid$meta$clinical_group==g]
 w <- suppressWarnings(wilcox.test(disease,healthy,exact=FALSE))
 data.frame(comparison=paste(g,"vs Healthy"),n_disease=length(disease),n_healthy=length(healthy),
 median_disease=median(disease),median_healthy=median(healthy),
 median_difference=median(disease)-median(healthy),statistic=unname(w$statistic),p_value=w$p.value)
}))
tests <- rbind(data.frame(comparison="Overall Kruskal-Wallis",n_disease=NA,n_healthy=NA,
 median_disease=NA,median_healthy=NA,median_difference=NA,statistic=unname(kw$statistic),
 p_value=kw$p.value),tests)
write.csv(tests,file.path(output_dir,"GSE51404_IBI_group_tests.csv"),row.names=FALSE)

pal <- c("XLA"="#315A7D","CVID without complications"="#70A9A1",
         "CVID with complications"="#D79053","Healthy"="#777777")
p1 <- ggplot(valid$meta,aes(clinical_group,IBI,colour=clinical_group))+
 geom_hline(yintercept=0,colour="#B8B8B8",linewidth=.35,linetype="dashed")+
 geom_boxplot(width=.58,outlier.shape=NA,fill="white",colour="black",linewidth=.4)+
 geom_jitter(width=.12,size=1.5,alpha=.78)+scale_colour_manual(values=pal,guide="none")+
 labs(title="IBI distribution in the independent GSE51404 test set",
 subtitle=sprintf("Kruskal-Wallis P = %.3g; n = 53 independent participants",kw$p.value),
 x=NULL,y="Immune Balance Index (IBI)")+
 theme_classic(base_size=8,base_family="Arial")+
 theme(axis.line=element_line(linewidth=.35),axis.ticks=element_line(linewidth=.35),
 axis.text.x=element_text(angle=25,hjust=1),plot.title=element_text(face="bold",size=9))
cairo_pdf(file.path(output_dir,"GSE51404_IBI_distribution.pdf"),120/25.4,95/25.4,family="Arial")
print(p1);dev.off()

selected <- c("HALLMARK_ALLOGRAFT_REJECTION","HALLMARK_COMPLEMENT","HALLMARK_IL2_STAT5_SIGNALING",
 "HALLMARK_IL6_JAK_STAT3_SIGNALING","HALLMARK_INFLAMMATORY_RESPONSE",
 "HALLMARK_INTERFERON_ALPHA_RESPONSE","HALLMARK_INTERFERON_GAMMA_RESPONSE",
 "HALLMARK_TNFA_SIGNALING_VIA_NFKB")
pretty <- c(HALLMARK_ALLOGRAFT_REJECTION="Allograft rejection",HALLMARK_COMPLEMENT="Complement",
 HALLMARK_IL2_STAT5_SIGNALING="IL2-STAT5 signaling",
 HALLMARK_IL6_JAK_STAT3_SIGNALING="IL6-JAK-STAT3 signaling",
 HALLMARK_INFLAMMATORY_RESPONSE="Inflammatory response",
 HALLMARK_INTERFERON_ALPHA_RESPONSE="IFN-alpha response",
 HALLMARK_INTERFERON_GAMMA_RESPONSE="IFN-gamma response",
 HALLMARK_TNFA_SIGNALING_VIA_NFKB="TNF-alpha/NF-kB signaling")
model_genes <- unique(toupper(unlist(strsplit(features,"_vs_",fixed=TRUE))))
msig <- msigdbr(species="Homo sapiens",collection="H");msig <- msig[msig$gs_name%in%selected,]
sets <- lapply(selected,function(s)setdiff(unique(toupper(msig$gene_symbol[msig$gs_name==s])),model_genes))
names(sets)<-pretty[selected]
ss <- gsva(ssgseaParam(valid_expr,sets,minSize=10,maxSize=500,normalize=TRUE),verbose=FALSE,
 BPPARAM=BiocParallel::SerialParam(progressbar=FALSE))

partial_spearman <- function(x,y,covariates){
 ok <- complete.cases(x,y,covariates);x<-x[ok];y<-y[ok];covariates<-droplevels(covariates[ok,,drop=FALSE])
 design<-model.matrix(~.,covariates);rx<-residuals(lm.fit(design,rank(x)));ry<-residuals(lm.fit(design,rank(y)))
 rho<-cor(rx,ry);k<-qr(design)$rank-1L;df<-length(x)-k-2L
 p<-2*pt(-abs(rho*sqrt(df/pmax(1-rho^2,.Machine$double.eps))),df)
 se<-1/sqrt(length(x)-k-3L);z<-atanh(pmax(pmin(rho,.999999),-.999999))
 c(n=length(x),rho=rho,ci_low=tanh(z-qnorm(.975)*se),ci_high=tanh(z+qnorm(.975)*se),p_value=p)
}
rows<-list()
for(sig in rownames(ss)){
 score<-as.numeric(ss[sig,valid$meta$sample_id])
 ct<-suppressWarnings(cor.test(valid$meta$IBI,score,method="spearman",exact=FALSE))
 rho<-unname(ct$estimate);se<-1/sqrt(nrow(valid$meta)-3L);z<-atanh(pmax(pmin(rho,.999999),-.999999))
 rows[[length(rows)+1L]]<-data.frame(signature=sig,analysis="Overall Spearman",subgroup="All",
 n=nrow(valid$meta),rho=rho,ci_low=tanh(z-qnorm(.975)*se),ci_high=tanh(z+qnorm(.975)*se),p_value=ct$p.value)
 a<-partial_spearman(valid$meta$IBI,score,data.frame(clinical_group=valid$meta$clinical_group))
 rows[[length(rows)+1L]]<-data.frame(signature=sig,analysis="Clinical-group adjusted",subgroup="All",
 n=a["n"],rho=a["rho"],ci_low=a["ci_low"],ci_high=a["ci_high"],p_value=a["p_value"])
 b<-partial_spearman(valid$meta$IBI,score,data.frame(clinical_group=valid$meta$clinical_group,
 age=valid$meta$age,gender=factor(valid$meta$gender)))
 rows[[length(rows)+1L]]<-data.frame(signature=sig,analysis="Group, age and sex adjusted",subgroup="All",
 n=b["n"],rho=b["rho"],ci_low=b["ci_low"],ci_high=b["ci_high"],p_value=b["p_value"])
 for(g in levels(valid$meta$clinical_group)){
  idx<-valid$meta$clinical_group==g
  ctg<-suppressWarnings(cor.test(valid$meta$IBI[idx],score[idx],method="spearman",exact=FALSE))
  r<-unname(ctg$estimate);se2<-1/sqrt(sum(idx)-3L);z2<-atanh(pmax(pmin(r,.999999),-.999999))
  rows[[length(rows)+1L]]<-data.frame(signature=sig,analysis="Within clinical group",subgroup=g,
  n=sum(idx),rho=r,ci_low=tanh(z2-qnorm(.975)*se2),ci_high=tanh(z2+qnorm(.975)*se2),p_value=ctg$p.value)
 }
}
cors<-do.call(rbind,rows)
cors$FDR<-ave(cors$p_value,interaction(cors$analysis,cors$subgroup,drop=TRUE),
 FUN=function(p)p.adjust(p,method="BH"))
write.csv(cors,file.path(output_dir,"GSE51404_IBI_signature_correlations.csv"),row.names=FALSE)
write.csv(data.frame(signature=names(sets),genes_after_model_gene_exclusion=lengths(sets),
 genes_present_in_GSE51404=vapply(sets,function(g)sum(g%in%rownames(valid_expr)),integer(1L))),
 file.path(output_dir,"GSE51404_signature_gene_coverage.csv"),row.names=FALSE)
pc<-cors[cors$analysis%in%c("Overall Spearman","Clinical-group adjusted","Group, age and sex adjusted"),]
pc$analysis<-factor(pc$analysis,levels=c("Overall Spearman","Clinical-group adjusted","Group, age and sex adjusted"))
pc$signature<-factor(pc$signature,levels=rev(unname(pretty[selected])))
pc$label<-sprintf("%.2f%s",pc$rho,ifelse(pc$p_value<.05,"*",""))
p2<-ggplot(pc,aes(analysis,signature,fill=rho))+geom_tile(colour="white",linewidth=.55)+
 geom_text(aes(label=label),size=2.5)+scale_fill_gradient2(low="#4C78A8",mid="white",
 high="#D65F5F",midpoint=0,limits=c(-1,1),name="Spearman rho")+
 labs(title="Association of IBI with immune signatures in GSE51404",
 subtitle="Signature scores exclude all model genes; *P < 0.05",x=NULL,y=NULL)+
 theme_classic(base_size=7.5,base_family="Arial")+
 theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_text(angle=25,hjust=1),
 plot.title=element_text(face="bold",size=9),legend.position="right")
cairo_pdf(file.path(output_dir,"GSE51404_IBI_immune_signature_correlations.pdf"),
 145/25.4,100/25.4,family="Arial");print(p2);dev.off()

counts<-as.data.frame(table(source=c(rep("Original development data",length(y)),rep("GSE51405",length(train_label))),
 label=c(as.character(y),as.character(train_label))))
write.csv(counts[counts$Freq>0,],file.path(output_dir,"updated_model_training_counts.csv"),row.names=FALSE)
write.csv(data.frame(item=c("Feature selection","Gene-pair encoding","Updated model training",
 "Hyperparameter policy","Independent validation","Primary correlation analysis",
 "Sensitivity correlation analysis","Statistical threshold"),
 value=c("Existing 390 gene pairs retained without re-screening",
 "1 when gene A expression > gene B expression; otherwise -1",
 "Original model-development matrix retained as -1/1, plus all GSE51405 samples encoded as -1/1",
 "Original SVR hyperparameters retained; GSE51404 not used for tuning","All 53 GSE51404 samples",
 "Partial Spearman correlation controlling clinical group",
 "Partial Spearman correlation controlling clinical group, age and sex","Two-sided P < 0.05")),
 file.path(output_dir,"analysis_provenance.csv"),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(output_dir,"R_sessionInfo.txt"))
cat("Completed:",output_dir,"\n")
