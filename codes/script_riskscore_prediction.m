
clc; clear; close all;
%--------------------------------------------------------------------------
% "Development and Validation of a Prognostic and Predictive 
%       32-Gene Signature for Gastric Cancer"
%
% 
%
% This script is written to reproduce the risk prediction shown in Figure 3-B 

% The risk score predictor (an SVM with a linear kernel, implemented by libSVM) 
% is trained using the best outcome group (black, set to the negative class) and 
% the worst outcome group (red, set to the positive class). 
% To reproduce the exact same result in the manucript, we provide the
% trained SVM model (User also can train his/her own SVM model by setting mflag_model_reuse to false.

% Please contanct dubuck@gmail.com if you have any concerns or comments about the implementation. 
%--------------------------------------------------------------------------

%- loading the data
load('../data/Combat_data_All.mat');   

%- please add the path where your LIBSVM files are located in
% addpath('../utils') % for Windows
addpath('/home/cs-com/Downloads/libsvm-3.25/matlab');

% model reuse
mflag_model_reuse = true;

if mflag_model_reuse
    tmp = load('../data/RiskPrediction/RiskscoreSVMPred_trn_Yonsei_tst_ALL_thres_25_75__predSVMres.mat');
    model_SVM = tmp.mc_SVMreslt.model;
else
    m_strSVM = '-s 0 -t 0 -c 1';
    
    % mc_data.Ytrn>> #(Y=-1): 114 and #(Y=1): 162 
    model_SVM = svmtrain(mc_data.Ytrn, mc_data.Xtrn, m_strSVM);
end

mn_Ntst = size(mc_data.Xtst, 1);

[mv_label_new, b_new, m_vScores_new] = ...
                svmpredict(ones(mn_Ntst,1), mc_data.Xtst, model_SVM, '-q');

if model_SVM.Label(1) == -1
    disp('converting the sign of the scores')
    m_vXvals = -m_vScores_new;
else
    m_vXvals = m_vScores_new;
end

y_min = quantile(m_vXvals,0.25);
y_max = quantile(m_vXvals,0.75);

m_vgroupIDX = zeros(length(m_vXvals), 1);
m_vgroupIDX(m_vXvals <= y_min) = 1;
m_vgroupIDX((m_vXvals > y_min) & (m_vXvals <= y_max)) = 2;
m_vgroupIDX(m_vXvals > y_max) = 3;

mm_line_colors = [0 0 0;  % black
                  0 1 0;  % green
                  1 0 0]; % red

% draw a KM plot              
addpath('../utils')
[p, fh] = MatSurv(mc_data.Surv_month, mc_data.Surv_status, ...
                  cellstr(string(m_vgroupIDX)),...
                  'LineColor',mm_line_colors,'Print',false);
              
