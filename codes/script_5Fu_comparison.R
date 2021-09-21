library("survminer");
library(survival)

library(ggplot2)
require(gdata)


# --------------------------------------------------------------------------
#   "Development and Validation of a Prognostic and Predictive 
# 32-Gene Signature for Gastric Cancer"
# 
#  
# This script is written to reproduce the result shown in Figure 4.
# 
# Please contanct dubuck@gmail.com if you have any concerns or comments about the implementation. 
# --------------------------------------------------------------------------
  
  
#--- data preprocessing 
# Please set the root directory of the scripts to the working directory
setwd('/home/cs-com/Downloads/Yonsei_gastric_cancer_32genes-main')

load(file.path(getwd(), "data/data_Yonsei_5Fu_comparison.R"))

#--------------------------------------------------------------
mstr_save_path = file.path(getwd(), "Yonsei_5_Fu_comparison_");

mc_colors = c('black', 'green', 'blue', 'red');
mstr_Chemo = c('None', '5-FU alone', '5-FU + platinum', '5-FU + others');

for (mn_seldat in 1:length(mc_colors)){
  # selection
  mv_ColorIDX = data_f$group_inof == mc_colors[mn_seldat];  
  data_f_sub = data_f[mv_ColorIDX,];
  
  # stage 
  xst_old = data_f[mv_ColorIDX,]$Stage;
  
  mstr_stage_cnt = sprintf('%d|| stage 1: %d, 2: %d, 3: %d, 4: %d', mn_seldat,
                           sum(xst_old==1), sum(xst_old==2), sum(xst_old==3), sum(xst_old==4))
  
  if (mc_colors[mn_seldat] == 'blue'){
    mstr_stageBlue = c('1', '2', '3', '4');
  } else {
    mstr_stageBlue = c('2', '1', '3', '4');
  } 
  
  data_f_sub$Stage = factor(xst_old, levels = mstr_stageBlue);
  
  # regimen code
  xctx_old = data_f[mv_ColorIDX,]$CTx;
  xrc_old = data_f[mv_ColorIDX,]$Regimen.Code;
  
  xtmp = rep('NA', rep=nrow(data_f_sub));
  
  chemo_idx_none = which(xctx_old == 0)
  xtmp[chemo_idx_none] = 'None'
  
  chemo_idx_5FuA = which((xrc_old == 1) | (xrc_old == 3) | (xrc_old == 5)) 
  xtmp[chemo_idx_5FuA] = '5-FU alone';
  
  chemo_idx_5FuP = which((xrc_old == 2) | (xrc_old == 6))
  xtmp[chemo_idx_5FuP] = '5-FU + platinum';
  
  chemo_idx_5FuO = which((xrc_old == 7) | (xrc_old == 9))
  xtmp[chemo_idx_5FuO] = '5-FU + others';
  
  data_f_sub$Regimen.Code = factor(xtmp, levels = mstr_Chemo);
  
  mstr_stage_rc = sprintf('%d|| None: %d, 5-Fu alone: %d, 5-Fu + Platimum: %d, 5-Fu  + others: %d', 
                          mn_seldat,
                          length(chemo_idx_none), length(chemo_idx_5FuA), 
                          length(chemo_idx_5FuP), length(chemo_idx_5FuO))
  print(mstr_stage_cnt)
  print(mstr_stage_rc) 
  
  # Cox 
  fit <- coxph(formula = Surv(surv_time_m, vital_sign_1_d) 
               ~ Age + Stage + LaurenType + Perineural_Invasion + Regimen.Code, data = data_f_sub) 
  
  # Prepare the columns
  beta <- exp(coef(fit))
  CI   <- round(exp(confint(fit)), 5)
  
  coeffs <- coef(summary(fit))
  p    <- as.matrix(coeffs[,5])
  
  # Bind columns together, and select desired rows
  res <- cbind(beta, CI, p)
  
  mn_i <- 10; # Regimen.Code5-FU + platinum vs No chemo
  mstr_line = sprintf("5FU + platinum>> HR: %.2f (%.2f, %.2f) p=%.2f",
                      beta[mn_i], CI[mn_i,1], CI[mn_i,2], p[mn_i])
  print(mstr_line) 
  
  
  #-----
  # https://github.com/kassambara/survminer/issues/67
  library(data.table)
  
  sel_idx = c(which(data_f_sub$Regimen.Code == '5-FU + platinum'), which(data_f_sub$Regimen.Code == 'None'))
  new_df <- copy(data_f_sub[sel_idx, ]) %>% setDT() 
  # new_df2 <- new_df[complete.cases(new_df)] #
  
  pred <- survfit(fit, new_df) 
  
  timepoints <- c(0, as.numeric(pred$time))
  timenames <- paste0("time", timepoints)
  
  adj_surv <- t(pred$surv) %>% as.data.table() # taking the predicted survival curves for each person in the study based on the cox model
  adj_surv <- cbind(time0 = 1, adj_surv) # adding survival probability = 1 at time = 0
  setnames(adj_surv, timenames) 
  
  new_df2 <- cbind(new_df, adj_surv) 
  new_df3 <- new_df2[, lapply(.SD, mean), .SDcols = timenames, by = Regimen.Code] # key step -- calculate mean survival at each time point, by group of interest (e.g., sex)
  
  graph <- cbind(timepoints, t(new_df3[, -1])) %>% as.data.table() # make data for graph
  
  strata_ = data.matrix(c(rep("5-FU + platinum", times=dim(graph[,1])[1]), 
                          rep("None", times=dim(graph[,2])[1])))
  surv_df = cbind(rbind(data.matrix(graph[,1]), data.matrix(graph[,1])), 
                  rbind(data.matrix(graph[,2]), data.matrix(graph[,3])), strata_)
  
  surv_df = cbind(surv_df, c(rep("", time=dim(surv_df)[1])))
  
  colnames(surv_df) <- c("time", "surv", "strata", "label") # variable names "5-FU + platinum", "None"
  surv_df = data.frame(surv_df)
  surv_df$time <- as.numeric(as.character(surv_df$time))
  surv_df$surv <- as.numeric(as.character(surv_df$surv))
  
  surv_df$label <- as.character(surv_df$label)
  
  #
  sel_idx1 = which((data_f_sub$Regimen.Code == '5-FU + platinum') & (data_f_sub$vital_sign_1_d == 0))
  timepoint_1 = data_f_sub$surv_time_m[sel_idx1]
  for (time_p in timepoint_1)
  {
    diff_ = abs(graph$timepoints - time_p)
    idx = which(diff_ == min(diff_))
    
    surv_df$label[idx] = '|'
  }
  
  sel_idx2 = which((data_f_sub$Regimen.Code == 'None') & (data_f_sub$vital_sign_1_d == 0))
  
  n_offset = length(graph$timepoints)
  timepoint_2 = data_f_sub$surv_time_m[sel_idx2]
  for (time_p in timepoint_2)
  {
    diff_ = abs(graph$timepoints - time_p)
    idx = which(diff_ == min(diff_))
    
    surv_df$label[n_offset + idx] = '|'  
  }
  # 
  
  # draw a survival curve
  ggsurv <- ggplot(data=surv_df, aes(x=time, y=surv, group=strata)) + 
    geom_line(aes(color=strata),size=1) + 
    scale_color_manual(values=c("black","red")) + 
    scale_x_continuous(breaks = seq(0, max(timepoints), by = 10)) +
    # scale_y_continuous(breaks = seq(-0.5, max(surv_df$surv), by = 0.1)) +
    ylim(0.0, 1.0) + 
    ylab("Survival probability") +
    xlab("Time in months") + 
    theme(axis.line = element_line(colour = "black"),
          legend.position="top",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    geom_text(aes(x=time, y=surv, label=label, colour=strata), vjust=0.20, size=3, fontface=3)
  
  mstr_filename = paste(mstr_save_path, '_', mc_colors[mn_seldat], "_KMplot.pdf", sep="", collapse="")  
  
  pdf(mstr_filename)
  print(ggsurv, newpage = FALSE)
  dev.off()
  #-----
}






