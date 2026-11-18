options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(e1071); library(glmnet); library(ranger); library(xgboost)
  library(pROC); library(ggplot2); library(cowplot); library(svglite)
})

root <- Sys.getenv("IBI_TOOL_ROOT", unset = "C:/temp/ibi_tool")
input_file <- file.path(root, "Result", "revision_781_fixed_top390_final", "fixed_top390_all_samples_781_development.rda")
output_dir <- file.path(root, "Result", "revision_781_six_model_nested_comparison_20261105")
if (dir.exists(output_dir)) stop("Output directory already exists: ", output_dir)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20261105)

cohort_order <- c("GTEx", "GSE232381", "GSE297198", "GSE264125", "GSE222408", "GSE175839", "GSE181228", "GSE199911", "GSE140713", "GSE59216_batch1", "GSE51405")
model_order <- c("SVR", "Random Forest", "Lasso", "Elastic Net", "XGBoost", "Ridge")
model_cols <- c("SVR"="#D55E00", "Random Forest"="#4C78A8", "Lasso"="#59A14F", "Elastic Net"="#B279A2", "XGBoost"="#E15759", "Ridge"="#76B7B2")
lambda_grid <- 10^seq(-4, 3, length.out = 25)
svr_grid <- expand.grid(cost=c(0.1,1,10), gamma=c(0.0001,0.001,0.01), epsilon=c(0.1,0.2,0.3))
rf_grid <- expand.grid(mtry_fraction=c(0.1,0.25,0.5), min_node_size=c(5,15))
xgb_grid <- expand.grid(max_depth=c(2,4), eta=c(0.03,0.1), min_child_weight=c(1,5))

e <- new.env(parent=emptyenv()); load(input_file, envir=e)
keep <- e$sample_manifest$cohort %in% cohort_order
manifest <- e$sample_manifest[keep,,drop=FALSE]
manifest <- manifest[order(match(manifest$cohort, cohort_order)),,drop=FALSE]
x <- as.matrix(e$top390_feature_matrix[manifest$sample_key,,drop=FALSE]); y <- as.numeric(manifest$label)
if (all(unique(as.numeric(x)) %in% c(-1,1))) x <- (x+1)/2
stopifnot(nrow(x)==781, ncol(x)==390, all(y %in% c(-1,0,1)))

rmse <- function(a,b) sqrt(mean((a-b)^2)); mae <- function(a,b) mean(abs(a-b))
rho <- function(a,b) if(length(unique(b))<2) NA_real_ else suppressWarnings(cor(a,b,method="spearman"))
prepare <- function(train,test) {
  modes <- apply(train,2,function(v) as.numeric(mean(v[is.finite(v)]==1)>=0.5))
  imp <- function(z){for(j in seq_len(ncol(z)))z[!is.finite(z[,j]),j]<-modes[j];z}
  train<-imp(train);test<-imp(test); s<-apply(train,2,sd); k<-is.finite(s)&s>0
  train<-train[,k,drop=FALSE];test<-test[,k,drop=FALSE];mu<-colMeans(train);ss<-apply(train,2,sd)
  list(train=sweep(sweep(train,2,mu,"-"),2,ss,"/"),test=sweep(sweep(test,2,mu,"-"),2,ss,"/"))
}
group_folds <- function(indices, cohorts, k=5) {
  g<-unique(cohorts[indices]);sz<-sapply(g,function(z)sum(cohorts[indices]==z));g<-g[order(-sz,g)]
  total<-rep(0,min(k,length(g)));a<-setNames(integer(length(g)),g)
  for(z in g){j<-which.min(total);a[z]<-j;total[j]<-total[j]+sum(cohorts[indices]==z)}
  lapply(seq_along(total),function(j)indices[cohorts[indices]%in%names(a)[a==j]])
}

tune_fold <- function(indices) {
  folds<-group_folds(indices,manifest$cohort,5)
  loss_glm<-array(NA,c(length(folds),3,length(lambda_grid)),dimnames=list(NULL,c("Ridge","Elastic Net","Lasso"),NULL))
  loss_svr<-matrix(NA,length(folds),nrow(svr_grid));loss_rf<-matrix(NA,length(folds),nrow(rf_grid));loss_xgb<-matrix(NA,length(folds),nrow(xgb_grid))
  for(f in seq_along(folds)){
    va<-folds[[f]];tr<-setdiff(indices,va);p<-prepare(x[tr,,drop=FALSE],x[va,,drop=FALSE])
    for(ai in seq_along(c(0,0.5,1))){fit<-glmnet(p$train,y[tr],alpha=c(0,0.5,1)[ai],lambda=lambda_grid,standardize=FALSE);pr<-predict(fit,p$test,s=lambda_grid);loss_glm[f,ai,]<-apply(pr,2,function(z)rmse(y[va],z))}
    for(j in seq_len(nrow(svr_grid))){q<-svr_grid[j,];fit<-svm(p$train,y[tr],type="eps-regression",kernel="radial",cost=q$cost,gamma=q$gamma,epsilon=q$epsilon,scale=FALSE);loss_svr[f,j]<-rmse(y[va],predict(fit,p$test))}
    for(j in seq_len(nrow(rf_grid))){q<-rf_grid[j,];fit<-ranger(x=p$train,y=y[tr],num.trees=300,mtry=max(1,round(ncol(p$train)*q$mtry_fraction)),min.node.size=q$min_node_size,seed=20261105+f+j,num.threads=1);loss_rf[f,j]<-rmse(y[va],predict(fit,p$test)$predictions)}
    dtr<-xgb.DMatrix(p$train,label=y[tr]);dva<-xgb.DMatrix(p$test)
    for(j in seq_len(nrow(xgb_grid))){q<-xgb_grid[j,];fit<-xgb.train(params=list(objective="reg:squarederror",max_depth=q$max_depth,eta=q$eta,min_child_weight=q$min_child_weight,subsample=0.8,colsample_bytree=0.8,nthread=1),data=dtr,nrounds=200,verbose=0);loss_xgb[f,j]<-rmse(y[va],predict(fit,dva))}
  }
  glm_mean<-apply(loss_glm,c(2,3),mean); gi<-apply(glm_mean,1,which.min)
  list(lambda=setNames(lambda_grid[gi],rownames(glm_mean)),svr=svr_grid[which.min(colMeans(loss_svr)),,drop=FALSE],rf=rf_grid[which.min(colMeans(loss_rf)),,drop=FALSE],xgb=xgb_grid[which.min(colMeans(loss_xgb)),,drop=FALSE])
}

fit_outer <- function(cohort) {
  te<-which(manifest$cohort==cohort);tr<-which(manifest$cohort!=cohort);tune<-tune_fold(tr);p<-prepare(x[tr,,drop=FALSE],x[te,,drop=FALSE]);out<-data.frame(sample_key=manifest$sample_key[te],cohort=cohort,observed=y[te])
  for(nm in c("Ridge","Elastic Net","Lasso")){alpha<-c(Ridge=0,`Elastic Net`=0.5,Lasso=1)[nm];fit<-glmnet(p$train,y[tr],alpha=alpha,lambda=tune$lambda[nm],standardize=FALSE);out[[nm]]<-as.numeric(predict(fit,p$test,s=tune$lambda[nm]))}
  q<-tune$svr;fit<-svm(p$train,y[tr],type="eps-regression",kernel="radial",cost=q$cost,gamma=q$gamma,epsilon=q$epsilon,scale=FALSE);out$SVR<-as.numeric(predict(fit,p$test))
  q<-tune$rf;fit<-ranger(x=p$train,y=y[tr],num.trees=500,mtry=max(1,round(ncol(p$train)*q$mtry_fraction)),min.node.size=q$min_node_size,seed=20261105,num.threads=1);out$`Random Forest`<-predict(fit,p$test)$predictions
  q<-tune$xgb;fit<-xgb.train(params=list(objective="reg:squarederror",max_depth=q$max_depth,eta=q$eta,min_child_weight=q$min_child_weight,subsample=0.8,colsample_bytree=0.8,nthread=1),data=xgb.DMatrix(p$train,label=y[tr]),nrounds=200,verbose=0);out$XGBoost<-predict(fit,xgb.DMatrix(p$test))
  list(pred=out,tuning=tune)
}

results<-list();tunings<-list()
for(i in seq_along(cohort_order)){cat(sprintf("Outer fold %d/%d: %s\n",i,length(cohort_order),cohort_order[i]));z<-fit_outer(cohort_order[i]);results[[i]]<-z$pred;tunings[[cohort_order[i]]]<-z$tuning;saveRDS(z,file.path(output_dir,paste0("checkpoint_",cohort_order[i],".rds")))}
pred<-do.call(rbind,results)

auc_one <- function(obs,score){labs<-sort(unique(obs));dis<-setdiff(labs,0);if(!(0%in%labs)||length(dis)!=1)return(NA_real_);s<-if(dis==-1)-score else score;as.numeric(auc(roc(as.integer(obs==dis),s,levels=c(0,1),direction="<",quiet=TRUE)))}
cohort_metrics<-do.call(rbind,lapply(cohort_order,function(cc){d<-pred[pred$cohort==cc,];do.call(rbind,lapply(model_order,function(m)data.frame(cohort=cc,model=m,n=nrow(d),RMSE=rmse(d$observed,d[[m]]),MAE=mae(d$observed,d[[m]]),Spearman_rho=rho(d$observed,d[[m]]),AUC=auc_one(d$observed,d[[m]]))))}))

set.seed(20261105);B<-2000;overall<-list()
for(m in model_order){point<-c(RMSE=rmse(pred$observed,pred[[m]]),MAE=mae(pred$observed,pred[[m]]),Spearman_rho=rho(pred$observed,pred[[m]]),Macro_AUC=mean(cohort_metrics$AUC[cohort_metrics$model==m],na.rm=TRUE));boot<-matrix(NA,B,4)
  for(b in 1:B){cs<-sample(cohort_order,length(cohort_order),replace=TRUE);d<-do.call(rbind,lapply(cs,function(cc)pred[pred$cohort==cc,]));av<-sapply(cs,function(cc)cohort_metrics$AUC[cohort_metrics$model==m&cohort_metrics$cohort==cc][1]);boot[b,]<-c(rmse(d$observed,d[[m]]),mae(d$observed,d[[m]]),rho(d$observed,d[[m]]),mean(av,na.rm=TRUE))}
  for(j in seq_along(point)){ci<-quantile(boot[,j],c(.025,.975),na.rm=TRUE);overall[[length(overall)+1]]<-data.frame(model=m,metric=names(point)[j],estimate=point[j],CI_low=ci[1],CI_high=ci[2])}}
overall<-do.call(rbind,overall)

# Main figure: raw values in cells; fill is direction-aware within-metric scaling only.
main<-overall;main$performance<-ave(main$estimate,main$metric,FUN=function(v)(v-min(v))/(max(v)-min(v)));main$performance[main$metric%in%c("RMSE","MAE")]<-1-main$performance[main$metric%in%c("RMSE","MAE")]
main$model<-factor(main$model,levels=rev(model_order));main$metric<-factor(main$metric,levels=c("RMSE","MAE","Spearman_rho","Macro_AUC"),labels=c("RMSE\n(lower better)","MAE\n(lower better)","Spearman rho\n(higher better)","Macro AUC\n(higher better)"))
p_main<-ggplot(main,aes(metric,model,fill=performance))+geom_tile(colour="white",linewidth=1)+geom_text(aes(label=sprintf("%.3f",estimate)),family="Arial",size=3)+scale_fill_gradientn(colours=c("#D9E2E8","#F4F4F4","#D7A84B"),limits=c(0,1),name="Relative performance\nwithin each metric")+labs(x=NULL,y=NULL,title="Nested cohort-grouped comparison of six IBI models",subtitle="Raw out-of-fold metrics; colour is scaled within metrics without aggregation")+theme_minimal(base_size=8,base_family="Arial")+theme(panel.grid=element_blank(),axis.text=element_text(colour="black"),axis.text.x=element_text(face="bold"),plot.title=element_text(face="bold",size=10),legend.title=element_text(size=7),legend.text=element_text(size=6.5))

sup<-overall;sup$model<-factor(sup$model,levels=rev(model_order));sup$metric<-factor(sup$metric,levels=c("RMSE","MAE","Spearman_rho","Macro_AUC"),labels=c("RMSE (lower is better)","MAE (lower is better)","Spearman rho (higher is better)","Macro AUC (higher is better)"))
p_sup<-ggplot(sup,aes(estimate,model,xmin=CI_low,xmax=CI_high,colour=model))+geom_errorbar(width=.18,orientation="y")+geom_point(size=2)+facet_wrap(~metric,scales="free_x",nrow=1)+scale_colour_manual(values=model_cols,guide="none")+labs(x="Estimate with 95% cohort-bootstrap CI",y=NULL,title="Separate performance metrics from nested cohort-grouped cross-validation")+theme_classic(base_size=7,base_family="Arial")+theme(strip.background=element_blank(),strip.text=element_text(face="bold",size=6),plot.title=element_text(face="bold"))

save_plot<-function(p,base,w=183,h=95){
  wi<-w/25.4; he<-h/25.4
  grDevices::cairo_pdf(paste0(base,".pdf"),width=wi,height=he,family="Arial");print(p);dev.off()
  svglite::svglite(paste0(base,".svg"),width=wi,height=he);print(p);dev.off()
  png(paste0(base,".png"),width=wi,height=he,units="in",res=300,type="cairo",family="Arial",bg="white");print(p);dev.off()
  tiff(paste0(base,".tiff"),width=wi,height=he,units="in",res=600,compression="lzw",type="cairo",family="Arial",bg="white");print(p);dev.off()
}
save_plot(p_main,file.path(output_dir,"main_six_model_all_metrics_heatmap"),150,100);save_plot(p_sup,file.path(output_dir,"supplement_six_model_metrics_with_ci"),183,85)
write.csv(pred,file.path(output_dir,"nested_oof_predictions.csv"),row.names=FALSE);write.csv(cohort_metrics,file.path(output_dir,"cohort_performance_metrics.csv"),row.names=FALSE);write.csv(overall,file.path(output_dir,"overall_performance_with_ci.csv"),row.names=FALSE);saveRDS(tunings,file.path(output_dir,"outer_fold_tuning.rds"));capture.output(sessionInfo(),file=file.path(output_dir,"R_sessionInfo.txt"))
cat("Complete:",output_dir,"\n")
