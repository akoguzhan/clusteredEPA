#' Clustered EPA Test Procedure for Balanced Panels
#'
#' @description Implements the clustered EPA test on a balanced panel dataset. Supports known and unknown clusters, and selective inference using homogeneity-adjusted p-values.
#'
#' @param df A data frame containing panel data.
#' @param DL A string; the name of the loss differential variable.
#' @param H A character vector; names of conditioning variables.
#' @param id A string; name of the panel unit identifier.
#' @param time A string; name of the time identifier.
#' @param gamma A string; name of the column indicating predetermined cluster identifiers.
#' @param G An integer; number of clusters for k-means.
#' @param Gmax An integer; maximum number of clusters for BIC selection.
#' @param maxlag An integer; maximum lag length for Newey-West estimation.
#' @param autocorr Character string; either "no" or "yes", determines long-run variance adjustment.
#' @param prop A scalar between 0 and 1; proportion of data used for R sample in sample-splitting.
#' @param r A numeric scalar; parameter used in p-value combination (e.g., r=20).
#' @param Ninit An integer; number of random initializations for k-means.
#' @param iter.max An integer; maximum iterations allowed for Lloyd's algorithm.
#'
#' @return A list with components:
#' \describe{
#'   \item{epa_test_stats}{Named list of Wald or PC statistics.}
#'   \item{epa_test_pvals}{Named list of corresponding p-values.}
#'   \item{comp_test_stats}{Overall EPA and homogeneity statistics.}
#'   \item{comp_test_pvals}{Overall EPA and homogeneity p-values.}
#'   \item{clustering_stats}{Additional information about estimated clusters and BIC-selected groups.}
#' }
#' @export
clusteredEPA_balanced_panel = function(df,DL=NULL,H=NULL,
                                       id=NULL,time=NULL,gamma=NULL,
                                       G=NULL,Gmax=NULL,
                                       maxlag=NULL,autocorr="no",
                                       prop=0.5,r=20,Ninit=1000,iter.max=100) {
  require(plm)
  require(dplyr)
  
  if (!is.data.frame(df)) {
    stop("df has to be a data frame")
  }
  if (is.null(DL)) {
    stop("Loss differential should be specified")
  }
  if (is.null(H)) {
    # warning("H is not specified, calculating unconditional test statistics")
    Z_names = DL
    K = 1
  } else {
    Z_names = c(DL,H)
    K = length(H)+1
  }
  if (is.null(id)|is.null(time)) {
    stop("Unit and time identifiers should be specified")
  }
  if (!(id %in% colnames(df))|!(time %in% colnames(df))) {
    stop("Unit and time identifiers are not found in the data set")
  }
  if (!is.null(gamma)) {
    if (!(gamma %in% colnames(df))) {
      stop("Predetermined cluster identifier is not found in the data set")
    }
  }
  if (is.null(G)) {
    stop("Number of estimated clusters should be specified")
  }
  
  data = pdata.frame(df[,c(Z_names,id,time,gamma)],index=c(id,time))
  if (!is.pbalanced(data)) {
    stop("Panel has to balanced")
  }
  if (!is.null(H)) {
    for (k in 1:length(H)) {
      data[,H[k]] = data[,DL]*data[,H[k]]
    }
  }

  colnames(data)[which(names(data) == id)] = "id"
  colnames(data)[which(names(data) == time)] = "time"
  Tobs = length(unique(data[,"time"]))
  
  if (is.null(maxlag)) {
    maxlag = 1.3*Tobs^(1/2)
  }

  ################################################################################
  # Full sample
  ################################################################################
  # Taking cross-sectional averages
  Zbar_Nt_df = data.frame(aggregate(data[,Z_names],list(data$time),FUN=mean)[,-1])
  colnames(Zbar_Nt_df) = c(Z_names)
  Zbar_Nt = as.matrix(Zbar_Nt_df)

  # k-means
  km_Z_T = panel_kmeans_estimation(data,id="id",time="time",Z_names=Z_names,G=G,Ninit=Ninit,iter.max=iter.max)
  gid.df = data.frame(cbind(names(km_Z_T$final_cluster),km_Z_T$final_cluster))
  colnames(gid.df) = c("id","gamma_hat_kmeans_T")
  data = merge(data,gid.df,by="id",sort=F)

  # k-means - IC
  km_IC_Z_T = panel_kmeans_estimation_bic(data,id="id",time="time",Z_names=Z_names,Gmax=Gmax,Ninit=Ninit,iter.max=iter.max)
  gid.df = data.frame(cbind(names(km_IC_Z_T$final_cluster),km_IC_Z_T$final_cluster))
  colnames(gid.df) = c("id","gamma_hat_kmeans_IC_T")
  data = merge(data,gid.df,by="id",sort=F)
  Ghat_T = length(unique(km_IC_Z_T$final_cluster))

  # Taking cross-sectional averages by cluster
  if (!is.null(gamma)) {
    Zbar_ngt_gamma_df = data.frame(aggregate(data[,Z_names],list(data$time,data[,gamma]),FUN=mean))
    colnames(Zbar_ngt_gamma_df) = c("time",gamma,Z_names)
    Zbar_ngt_gamma = as.matrix(reshape(Zbar_ngt_gamma_df,idvar="time",timevar=gamma,direction="wide")[,-1])
  }
  
  Zbar_ngt_gamma_kmeans_df = data.frame(aggregate(data[,Z_names],list(data$time,data$gamma_hat_kmeans_T),FUN=mean))
  colnames(Zbar_ngt_gamma_kmeans_df) = c("time","gamma_hat_kmeans_T",Z_names)
  Zbar_ngt_gamma_kmeans = as.matrix(reshape(Zbar_ngt_gamma_kmeans_df,idvar="time",timevar="gamma_hat_kmeans_T",direction="wide")[,-1])

  Zbar_ngt_gamma_kmeans_IC_df = data.frame(aggregate(data[,Z_names],list(data$time,data$gamma_hat_kmeans_IC_T),FUN=mean))
  colnames(Zbar_ngt_gamma_kmeans_IC_df) = c("time","gamma_hat_kmeans_IC_T",Z_names)
  Zbar_ngt_gamma_kmeans_IC = as.matrix(reshape(Zbar_ngt_gamma_kmeans_IC_df,idvar="time",timevar="gamma_hat_kmeans_IC_T",direction="wide")[,-1])

  ################################################################################
  # Split sample
  ################################################################################
  data_R = data %>% group_by(id) %>% slice_head(prop=prop)
  data_P = anti_join(data,data_R,by=c("id","time"))
  data_R = data_R %>% group_by(id) %>% slice(1:(n()-maxlag))
  P = length(unique(data_P$time))
  
  # Taking cross-sectional averages
  Zbar_Nt_inR_df = data.frame(aggregate(data_R[,Z_names],list(data_R$time),FUN=mean)[,-1])
  colnames(Zbar_Nt_inR_df) = c(Z_names)
  Zbar_Nt_inR = as.matrix(Zbar_Nt_inR_df)

  # k-means
  km_Z_R = panel_kmeans_estimation(as.data.frame(data_R),id="id",time="time",Z_names=Z_names,G=G,Ninit=Ninit,iter.max=iter.max)
  gid.df = data.frame(cbind(names(km_Z_R$final_cluster),km_Z_R$final_cluster))
  colnames(gid.df) = c("id","gamma_hat_kmeans_R")
  data_R = merge(data_R,gid.df,by="id",sort=F)
  data_P = merge(data_P,gid.df,by="id",sort=F)
  
  # k-means - IC
  km_IC_Z_R = panel_kmeans_estimation_bic(as.data.frame(data_R),id="id",time="time",Z_names=Z_names,Gmax=Gmax,Ninit=Ninit,iter.max=iter.max)
  gid.df = data.frame(cbind(names(km_IC_Z_R$final_cluster),km_IC_Z_R$final_cluster))
  colnames(gid.df) = c("id","gamma_hat_kmeans_IC_R")
  data_R = merge(data_R,gid.df,by="id",sort=F)
  data_P = merge(data_P,gid.df,by="id",sort=F)
  Ghat_R = length(unique(km_IC_Z_R$final_cluster))

  # Taking cross-sectional averages by cluster
  Zbar_ngt_inP_gamma_kmeans_df = data.frame(aggregate(data_P[,Z_names],list(data_P$time,data_P$gamma_hat_kmeans_R),FUN=mean))
  colnames(Zbar_ngt_inP_gamma_kmeans_df)[1:2] = c("time","gamma_hat_kmeans_R")
  Zbar_ngt_inP_gamma_kmeans = as.matrix(reshape(Zbar_ngt_inP_gamma_kmeans_df,idvar="time",timevar="gamma_hat_kmeans_R",direction="wide")[,-1])

  Zbar_ngt_inP_gamma_kmeans_IC_df = data.frame(aggregate(data_P[,Z_names],list(data_P$time,data_P$gamma_hat_kmeans_IC_R),FUN=mean))
  colnames(Zbar_ngt_inP_gamma_kmeans_IC_df)[1:2] = c("time","gamma_hat_kmeans_IC_R")
  Zbar_ngt_inP_gamma_kmeans_IC = as.matrix(reshape(Zbar_ngt_inP_gamma_kmeans_IC_df,idvar="time",timevar="gamma_hat_kmeans_IC_R",direction="wide")[,-1])

  ################################################################################
  # Test statistics
  ################################################################################
  if (autocorr=="no") {
    B.overall = Tobs
    B = Tobs
    B.split = Tobs*prop
    B_IC = Tobs*prop
    B.split_IC = Tobs*prop
  } else {
    B.overall = 2*max(floor((K+4)/2),1)
    B = 2*max(floor((G*K+4)/2),1)
    B.split = 2*max(floor((G*K+4)/2),1)
    B_IC = 2*max(floor((Ghat_T*K+4)/2),1)
    B.split_IC = 2*max(floor((Ghat_R*K+4)/2),1)
  }
  
  adj.term_overall = ((B.overall-K+1)/B.overall)/(K)
  adj.term = ((B-G*K+1)/B)/(G*K)
  adj.term_split = ((B.split-G*K+1)/B.split)/(G*K)
  adj.term_IC = ((B_IC-Ghat_T*K+1)/B_IC)/(Ghat_T*K)
  adj.term_split_IC = ((B.split_IC-Ghat_R*K+1)/B.split_IC)/(Ghat_R*K)
  
  # Overall EPA Test
  # mu.hat.oepa = colMeans(Zbar_Nt)
  # S.oepa = NeweyWest(Zbar_Nt,maxlag)
  # W.oepa = Tobs*t(mu.hat.oepa)%*%solve(S.oepa)%*%(mu.hat.oepa)
  # p.oepa = pchisq(q=W.oepa,df=K,lower.tail=FALSE)
  
  mu.hat.oepa = colMeans(Zbar_Nt)
  S.oepa = EWC(Zbar_Nt,B.overall)
  W.oepa = adj.term_overall*Tobs*t(mu.hat.oepa)%*%solve(S.oepa)%*%(mu.hat.oepa)
  p.oepa = pf(W.oepa,K,B.overall-K+1,lower.tail=FALSE)

  # Clusters known
  if (!is.null(gamma)) {
    # mu.hat = colMeans(Zbar_ngt_gamma)
    # S = NeweyWest(Zbar_ngt_gamma,maxlag)
    # W = Tobs*t(mu.hat)%*%solve(S)%*%(mu.hat)
    # p = pchisq(q=W,df=G*K,lower.tail=FALSE)
    
    mu.hat = colMeans(Zbar_ngt_gamma)
    S = EWC(Zbar_ngt_gamma,B)
    W = adj.term*Tobs*t(mu.hat)%*%solve(S)%*%(mu.hat)
    p = pf(W,G*K,B-G*K+1,lower.tail=FALSE)
  } else {
    W = NULL
    p = NULL
  }

  # Naive test statistics - clusters unknown - kmeans
  # mu.hat.naive.kmeans = colMeans(Zbar_ngt_gamma_kmeans)
  # S.naive.kmeans = NeweyWest(Zbar_ngt_gamma_kmeans,maxlag)
  # W.naive.kmeans = Tobs*t(mu.hat.naive.kmeans)%*%solve(S.naive.kmeans)%*%(mu.hat.naive.kmeans)
  # p.naive.kmeans = pchisq(q=W.naive.kmeans,df=G*K,lower.tail=FALSE)
  
  mu.hat.naive.kmeans = colMeans(Zbar_ngt_gamma_kmeans)
  S.naive.kmeans = EWC(Zbar_ngt_gamma_kmeans,B)
  W.naive.kmeans = adj.term*Tobs*t(mu.hat.naive.kmeans)%*%solve(S.naive.kmeans)%*%(mu.hat.naive.kmeans)
  p.naive.kmeans = pf(W.naive.kmeans,G*K,B-G*K+1,lower.tail=FALSE)

  # Split sample test statistics - clusters unknown - kmeans
  # mu.hat.split.kmeans = colMeans(Zbar_ngt_inP_gamma_kmeans)
  # S.split.kmeans = NeweyWest(Zbar_ngt_inP_gamma_kmeans,maxlag)
  # W.split.kmeans = P*t(mu.hat.split.kmeans)%*%solve(S.split.kmeans)%*%(mu.hat.split.kmeans)
  # p.split.kmeans = pchisq(q=W.split.kmeans,df=G*K,lower.tail=FALSE)
  
  mu.hat.split.kmeans = colMeans(Zbar_ngt_inP_gamma_kmeans)
  S.split.kmeans = EWC(Zbar_ngt_inP_gamma_kmeans,B.split)
  W.split.kmeans = adj.term_split*P*t(mu.hat.split.kmeans)%*%solve(S.split.kmeans)%*%(mu.hat.split.kmeans)
  p.split.kmeans = pf(W.split.kmeans,G*K,B.split-G*K+1,lower.tail=FALSE)
  
  # mu.hat.split.kmeans_IC = colMeans(Zbar_ngt_inP_gamma_kmeans_IC)
  # S.split.kmeans_IC = NeweyWest(Zbar_ngt_inP_gamma_kmeans_IC,maxlag)
  # W.split.kmeans_IC = P*t(mu.hat.split.kmeans_IC)%*%solve(S.split.kmeans_IC)%*%(mu.hat.split.kmeans_IC)
  # p.split.kmeans_IC = pchisq(q=W.split.kmeans_IC,df=Ghat_R*K,lower.tail=FALSE)
  
  mu.hat.split.kmeans_IC = colMeans(Zbar_ngt_inP_gamma_kmeans_IC)
  S.split.kmeans_IC = EWC(Zbar_ngt_inP_gamma_kmeans_IC,B.split_IC)
  W.split.kmeans_IC = adj.term_split_IC*P*t(mu.hat.split.kmeans_IC)%*%solve(S.split.kmeans_IC)%*%(mu.hat.split.kmeans_IC)
  p.split.kmeans_IC = pf(W.split.kmeans_IC,Ghat_R*K,B.split_IC-Ghat_R*K+1,lower.tail=FALSE)

  # Selective p-value - kmeans
  # S.naive.kmeans_selective = NeweyWest(Zbar_ngt_gamma_kmeans,maxlag)
  S.naive.kmeans_selective = S.naive.kmeans
  pairwise_pvalues_kmeans = panel_homogeneity_test(data,id="id",time="time",
                                                   Z_names=Z_names,
                                                   Omega=S.naive.kmeans_selective,
                                                   estimated_k_means=km_Z_T)
  pvalues.selective = c(p.oepa,pairwise_pvalues_kmeans$pairwise_pvalues)
  PC.kmeans = (1/length(pvalues.selective))*sum(pvalues.selective^(-r))^(1/r)
  p.selective.kmeans = min((r/(r-1))*(1/PC.kmeans),1)
  
  # S.naive.kmeans_selective_IC = NeweyWest(Zbar_ngt_gamma_kmeans_IC,maxlag)
  S.naive.kmeans_selective_IC = EWC(Zbar_ngt_gamma_kmeans_IC,B_IC)
  pairwise_pvalues_kmeans_IC = panel_homogeneity_test(data,id="id",time="time",
                                                      Z_names=Z_names,
                                                      Omega=S.naive.kmeans_selective_IC,
                                                      estimated_k_means=km_IC_Z_T)
  pvalues.selective_IC = c(p.oepa,pairwise_pvalues_kmeans_IC$pairwise_pvalues)
  PC.kmeans_IC = (1/length(pvalues.selective_IC))*sum(pvalues.selective_IC^(-r))^(1/r)
  p.selective.kmeans_IC = min((r/(r-1))*(1/PC.kmeans_IC),1)

  ################################################################################
  # Results
  ################################################################################
  # Test statistics
  r1 = list("W.known" = W,
            "W.naive.kmeans" = W.naive.kmeans,
            "W.split.kmeans" = W.split.kmeans,
            "PC.kmeans" = PC.kmeans,
            "W.split.kmeans_IC" = W.split.kmeans_IC,
            "PC.kmeans_IC" = PC.kmeans_IC)
  
  # p-values
  r2 = list("p.known" = p,
            "p.naive.kmeans" = p.naive.kmeans,
            "p.split.kmeans" = p.split.kmeans,
            "p.selective.kmeans" = p.selective.kmeans,
            "p.split.kmeans_IC" = p.split.kmeans_IC,
            "p.selective.kmeans_IC" = p.selective.kmeans_IC)
  
  # Overall EPA and heterogeneity test statistics
  r3 = list("W.oepa" = W.oepa,
            "PC.homogeneity.kmeans" = pairwise_pvalues_kmeans$PC.kmeans,
            "PC.homogeneity.kmeans_IC" = pairwise_pvalues_kmeans_IC$PC.kmeans)
  
  # Overall EPA and heterogeneity p-values
  r4 = list("p.oepa" = p.oepa,
            "p.homogeneity.kmeans" = pairwise_pvalues_kmeans$pvalue_combination,
            "p.homogeneity.kmeans_IC" = pairwise_pvalues_kmeans_IC$pvalue_combination)
  
  # Additional statistics
  r5 = list("Est_Cl_FS" = km_Z_T$final_cluster,
            "Est_Cl_SS" = km_Z_R$final_cluster,
            "Ghat_FS" = Ghat_T,
            "Ghat_SS" = Ghat_R,
            "Est_Cl_FS_Ghat" = km_IC_Z_T$final_cluster,
            "Est_Cl_SS_Ghat" = km_IC_Z_R$final_cluster)
  
  RESULTS = list("epa_test_stats" = r1,
                 "epa_test_pvals" = r2,
                 "comp_test_stats" = r3,
                 "comp_test_pvals" = r4,
                 "clustering_stats" = r5)
  return(RESULTS)
}