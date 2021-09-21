clc; clear; close all;
%--------------------------------------------------------------------------
% "Development and Validation of a Prognostic and Predictive 
%       32-Gene Signature for Gastric Cancer"
%
% 
%
% This script is written to reproduce the consensus clustering resutls shown in
% Figure 2 B-C-D. 
%
% Input data (X): samples X features (=32 genes)
% Output variable (Y): samples X 1
%
% We run a consensus clustering method [Stefano Monti et. al, Machine Learning 
% (2003)] with the fixed number of clusters K=4. We have used a version of NMF [Jingu Kim   
% and Haesun Park, ICDM (2008)] as a base clustering method. 
% To produce the exact clustering results in the main paper, we have provided 
% the ramdom indices of the samples for the bootstrapping and all the inital matrices 
% of the NMF. User also can run the consensus clustering method with random
% initialization by setting mflag_random to True;

% Please contanct dubuck@gmail.com if you have any concerns or comments about the implementation. 
%--------------------------------------------------------------------------

%--- Loading the input data 
mc_datasets = {'Yonsei', 'ACRG', 'Shon_etal'}; % for figure 1 B, C, and D respectively 
mn_datanum = 1;

mstr_dataset = mc_datasets{mn_datanum};

load(['../data/data_', mstr_dataset, '.mat']);

disp(['Data set: ', mstr_dataset])
mstr_data_org = mc_data.discription;

OX = mc_data.X;

mc_genes = mc_data.gene;
mc_sampleids = mc_data.pID;
mm_Surv = mc_data.surv;

mc_init_W = mc_data.init_W;
mc_init_H = mc_data.init_H;
mc_init_idx = mc_data.init_idx;

mflag_random = false;
%--- 


%--- Consensus clustering [Stefano Monti et. al, Machine Learning (2003)] 
mn_MaxIters = 1000; % the number of repetition 
mr_ResRate = 0.8; % the resampling rate

mn_Patients = length(mc_sampleids);
mn_SubSamples = round(mn_Patients*mr_ResRate);

%- Add the path of the NMF [Jingu Kim and Haesun Park, ICDM (2008)] to MATLAB
% We have used the nmf package available at https://www.cc.gatech.edu/~hpark/nmfbpas.html
addpath('../utils/nmf_bpas/')

mstr_savePath = '../figures/';
if ~exist(mstr_savePath, 'dir')
    mkdir(mstr_savePath)
end
    
%- input 
if (sum(sum(OX<0)) ==  0)
    X = OX;        
else
    X = [max(OX,0) max(-OX,0)];
    
    disp('Tranformation: Duplication');
end

% consens clustering
K = 4;
    
fprintf('K=%d case is being processed \n', K);

%- Run the NMF 1000 (m_nMaxIters) times with resampling
mm_Group = inf*ones(mn_Patients, mn_MaxIters);
for mn_iter = 1:mn_MaxIters
    %- resampling: randomly selected 80% of the samples
    if mflag_random
        mv_RandIDX = randperm(mn_Patients);
        
        W_init = rand(mn_SubSamples, K);
        H_init = rand(K, length(mc_genes));
    else
        mv_RandIDX = mc_init_idx{mn_iter};
        
        W_init = mc_init_W{mn_iter, 1};
        H_init = mc_init_H{mn_iter, 1};
    end
    
    mv_DataIDX = mv_RandIDX(1:mn_SubSamples);    
    mm_X = X(mv_DataIDX,:);    
    
    %- run the base clustering method
    [W, H] = nmf(mm_X, K, 'W_INIT',W_init, 'H_INIT',H_init, ...
        'type','plain', 'tol',1e-4, 'NNLS_SOLVER','bp');
    
    %- clustering results from each NMF run
    [val, mv_idx] = max(W, [], 2);
    
    mm_Group(mv_DataIDX, mn_iter) = mv_idx;
end

%- Construct a consensus matrix (samples X samples): each element of the
% matrix is interpreted as a probability that 
% the corresponding two samples are in the same cluster   
mm_Sim = zeros(mn_Patients, mn_Patients);
for mn_i = 1:mn_Patients
    mv_CurIDs = mm_Group(mn_i,:);
    
    mm_Chk = (mm_Group)./repmat(mv_CurIDs,[mn_Patients,1]);
    mm_Chk = mm_Chk == 1.0;
    
    mm_Inc = (mm_Group<Inf) & repmat(mv_CurIDs<Inf,[mn_Patients,1]);
    
    mm_Sim(mn_i,:) = sum(mm_Chk,2)./sum(mm_Inc,2);
end

mm_Sim(1:mn_Patients+1:mn_Patients^2) = 1;

%- Obtain the final clusters from the consensus matrix by the hierarchical clustering
mv_pdist = squareform(1-mm_Sim);
Z = linkage(mv_pdist, 'average');
mv_groupIDX = cluster(Z, 'maxclust', K); % cluster assignment

%
[m_vdummy, mv_ReorderIDX] = sort(mv_groupIDX);

fig_open = figure();
colormap('jet');   % set colormap
h = imagesc((mm_Sim(mv_ReorderIDX, mv_ReorderIDX)));        % draw image and scale colormap to values range
colorbar;
    

%--- Final outcome 
Yest = mv_groupIDX;

% when mflag_random is false
T_ref = readtable(['../data/ConsensusClustering/', mstr_data_org,...
                   '_4_groupInfo.txt'], 'Delimiter', '\t');

Yref = cell(length(mc_sampleids), 1);
[mc_dummy, mv_Aidx, mv_Bidx] = intersect(mc_sampleids, T_ref.P_ID);
Yref(mv_Aidx) = T_ref.GroupColor(mv_Bidx);

mc_groups = {'black', 'green', 'blue', 'red'};
YrefNum = zeros(length(Yref), 1);
for mn_i = 1:length(mc_groups)
    YrefNum(strcmpi(Yref, mc_groups{mn_i})) = mn_i;
end

% The following code (bestMap) matches the cluster labels obtained by the above method with the ones we already reported in the paper.  
addpath('../utils/')
Yest_adj = bestMap(YrefNum, Yest);

disp(confusionmat(YrefNum, Yest_adj))

mm_line_colors = [0 0 0;  % black
                  0 1 0;  % green
                  0 0 1;  % blue
                  1 0 0]; % red
              
[p, fh] = MatSurv(mm_Surv(:,1), mm_Surv(:,2), cellstr(string(Yest_adj)),...
                  'LineColor',mm_line_colors,'Print',false);
