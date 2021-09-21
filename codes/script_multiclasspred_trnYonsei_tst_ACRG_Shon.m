clc; clear; close all
%--------------------------------------------------------------------------
% "Development and Validation of a Prognostic and Predictive 
%       32-Gene Signature for Gastric Cancer"
%
%
% This script is written to reproduce the result shown in Supplementary Figure 7.
%
% A multiclass classifier (implemented based on multiple SVMs in all pairwise
% comparisons) is trained using the Yonsei cohort (the cluster labels, from 1
% to 4, are used as multi-class labels). 

% Please contanct dubuck@gmail.com if you have any concerns or comments about the implementation. 
%--------------------------------------------------------------------------

%---
mc_datasets = {'ACRG', 'Shon'};
mn_datanum = 2;

mstr_tst_data = mc_datasets{mn_datanum};

disp(['Data set: Yonsei + ', mstr_tst_data]);

%- read the combat file
load(['../data/Combat_data_Yonsei_', mstr_tst_data,'.mat'])

mm_Xtrn = mc_data.Xtrn;
% gene names: mc_data.Xfea
mv_Ytrn= mc_data.Ytrn;

mm_Xtst = mc_data.Xtst;

mm_Surv = mc_data.Surv;

%- Multiclass SVM
K = 4;

%- please add the path where your LIBSVM files are located in
% addpath('../utils') % for Windows
addpath('/home/cs-com/Downloads/libsvm-3.25/matlab');

mstr_code = 'APs';
mm_Codemat = myfunc_ConstCode((1:K)', mstr_code);
mn_M = size(mm_Codemat, 2);

mstr_svmsetting = '-s 0 -t 0 -b 1 -q -c 1';

mm_ProbMat_trn = zeros(size(mm_Xtrn,1), mn_M);
mm_ProbMat_tst = zeros(size(mm_Xtst,1), mn_M);
for mn_sub = 1:mn_M
    mv_Targets = mm_Codemat(mv_Ytrn, mn_sub);
    mv_Targets(mv_Targets==0) = -1;
    
    %- training
    mv_trnIdx = ~isnan(mv_Targets);
    model_ = svmtrain...
        (mv_Targets(mv_trnIdx), sparse(mm_Xtrn(mv_trnIdx, :)), mstr_svmsetting);
    
    %- evaluate    
    [~, ~, m_mTrnOuts] = svmpredict(ones(size(mm_Xtrn,1),1), sparse(mm_Xtrn), model_, '-b 1 -q');
    [~, ~, m_mTstOuts] = svmpredict(ones(size(mm_Xtst,1), 1), sparse(mm_Xtst), model_, '-b 1 -q');
        
    mv_idx = model_.Label == 1;    
    mv_trnVals = m_mTrnOuts(:, mv_idx);
    mv_tstVals = m_mTstOuts(:, mv_idx);
    
    mm_ProbMat_trn(:, mn_sub) = mv_trnVals;
    mm_ProbMat_tst(:, mn_sub) = mv_tstVals;
end
        
%- Ours: Multiclass SVM
m_vWeights = ones(mn_M, 1)/mn_M;

% [mv_Pred_trn, ~] = myfunc_Prediction_W(mm_ProbMat_trn, m_vWeights, mm_Codemat);
% fprintf('(Acc = %.5f) \n', sum(mv_Ytrn==mv_Pred_trn)/length(mv_Ytrn));

% test
[mv_Pred_tst, m_mOuts_tst] = myfunc_Prediction_W(mm_ProbMat_tst, m_vWeights, mm_Codemat);

% draw a KM plot 
addpath('../utils')

mm_line_colors = [0 0 0;  % black
                  0 1 0;  % green
                  0 0 1;  % blue
                  1 0 0]; % red
              
[p, fh] = MatSurv(mm_Surv(:,1), mm_Surv(:,2), cellstr(string(mv_Pred_tst)),...
                  'LineColor',mm_line_colors,'Print',false);




