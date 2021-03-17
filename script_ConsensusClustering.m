clc; clear; close all;
%--------------------------------------------------------------------------
% "Development and Validation of a Prognostic and Predictive 
%       32-Gene Signature for Gastric Cancer"
%
% Submitted to Nature Medicine
% 
%
% This script is written to run consensus clustering on Yonsei 32 gene data
% Input data (X): 567 (# samples) X 32 (# genes)
% Output variable (Y): 567 X 1 (cluster indices)
%
% We run a consensus clustering method [Stefano Monti et. al, Machine Learning (2003)] 
% with increasing the number of clusters K=[2,...,8] and select the best K based on 
% several metrics, including the proportion change in the area under the consensus CDF
% Cophenetic correlation coefficient, and visual inspection.
%
% For each K, we run a NMF [Jingu Kim and Haesun Park, ICDM (2008)] multiple 
% times (1000) with the resampling on the input data (for each run, 
% we randomly select 80% of the samples without replacement and 
% obtain clustering results by the NMF with the given K on the selected samples). 
% To aggregate all the clustering results, we construct a consensus matrix 
% (samples X samples), where each element represents the probability of 
% the corresponding pair of two samples being in the same cluster.     
% The consensus matrix can be used as new input data to any clustering methods 
% (e.g., the hierarchical clustering) in order to produce more stable clustering results.
%
% Please contanct parks@ccf.org if you have any concerns or comments about the implementation. 
%--------------------------------------------------------------------------

%--- Loading the input data   
fid = fopen('./data/Yonsei_gastric_cancer_32genes_data_norm.txt', 'r');

tline = fgets(fid);
m_cheaders_en = strsplit(tline,'\t');
m_nLen = length(m_cheaders_en);

m_strpattern = ['%s', repmat('%f',[1,m_nLen]), '%[^\n\r]'];
m_mData = textscan(fid,m_strpattern,'Delimiter','\t');

fclose(fid);

% X (samples X genes)
mc_genes = m_cheaders_en(2:end)';
mc_sampleids = m_mData{1, 1};

OX = zeros(length(mc_sampleids), length(mc_genes));
for mn_sub = 1:length(mc_genes)
    OX(:,mn_sub) = m_mData{1, mn_sub+1};
end
%--- 


%--- Consensus clustering [Stefano Monti et. al, Machine Learning (2003)] 
m_nMaxIters = 1000; % the number of repetition 
m_rResRate = 0.8; % the resampling rate

m_nPatients = length(mc_sampleids);

%- Add the path of the NMF [Jingu Kim and Haesun Park, ICDM (2008)] to the system
% We use a nmf package available at https://www.cc.gatech.edu/~hpark/nmfbpas.html
addpath('E:\Research_Mar092017\Code\tools\nmf_bpas');

m_str_savePath = './figures/';
if ~exist(m_str_savePath, 'dir')
    mkdir(m_str_savePath)
end
    
%- input 
if (sum(sum(OX<0)) ==  0)
    X = OX;        
else
    X = [max(OX,0) max(-OX,0)];
    
    disp('Tranformation: Duplication');
end

% resample
m_nSubSamples = round(m_nPatients*m_rResRate);

m_vKs = [2, 3, 4, 5, 6, 7, 8];

m_vsihouette = zeros(length(m_vKs), 1);
m_vAk = zeros(length(m_vKs), 1);
m_vCCr = zeros(length(m_vKs), 1);

% random number controll
rng shuffle;

for m_nselK = 1:length(m_vKs)
    K = m_vKs(m_nselK);
    
    fprintf('K=%d case is being processed \n', K);
    
    %- Run the NMF 1000 (m_nMaxIters) times with resampling
    m_mGroup = inf*ones(m_nPatients, m_nMaxIters);
    for m_niter = 1:m_nMaxIters
        %- resampling: randomly 80% of the samples
        m_vRandIDX = randperm(m_nPatients);    
        m_vDataIDX = m_vRandIDX(1:m_nSubSamples);

        m_mX = X(m_vDataIDX,:);
        [n, m1] = size(m_mX);
        
        %- initialization
        W_init = rand(n, K);
        H_init = rand(K, m1);
        
        [W, H] = nmf(m_mX, K, 'W_INIT',W_init, 'H_INIT',H_init, ...
            'type','plain', 'tol',1e-4, 'NNLS_SOLVER','bp');
        
        %- clustering results obtained from the NMF
        [val, m_vidx] = max(W, [], 2);

        m_mGroup(m_vDataIDX, m_niter) = m_vidx;
    end
    
    %- Construct a consensus matrix (samples X samples)
    m_mSim = zeros(m_nPatients, m_nPatients);
    for m_ni = 1:m_nPatients
        m_vCurIDs = m_mGroup(m_ni,:);
        
        m_mChk = (m_mGroup)./repmat(m_vCurIDs,[m_nPatients,1]);
        m_mChk = m_mChk == 1.0;
        
        m_mInc = (m_mGroup<Inf) & repmat(m_vCurIDs<Inf,[m_nPatients,1]);
        
        m_mSim(m_ni,:) = sum(m_mChk,2)./sum(m_mInc,2);
    end
    
    m_mSim(1:m_nPatients+1:m_nPatients^2) = 1;
        
    m_cResults.Group{m_nselK,1} = m_mGroup;
    m_cResults.SimMat{m_nselK,1} = m_mSim;
        
    %- Obtain the final clusters from the consensus matrix by the hierarchical clustering
    mv_pdist = squareform(1-m_mSim);
    Z = linkage(mv_pdist, 'average');
    m_vgroupIDX = cluster(Z, 'maxclust', K); % cluster assignment
    
    m_cResults.Clusters{m_nselK,1} = m_vgroupIDX;
    
    %- Draw the consensus matrix
    [m_vdummy, m_vReorderIDX] = sort(m_vgroupIDX);
        
    fig_open = figure('visible','off');
    colormap('jet');   % set colormap
    h = imagesc((m_mSim(m_vReorderIDX, m_vReorderIDX)));        % draw image and scale colormap to values range
    colorbar; 

    m_strname = [m_str_savePath, 'Cluster_K_', num2str(K), '.pdf'];
    saveas(fig_open, m_strname)
    close all hidden
    
    %- Calculate the coefficients
    m_vupperidx = find(~tril(ones(size(m_mSim))));
    m_vUpperVals = sort(m_mSim(m_vupperidx), 'ascend');
    
    % CDF
    m_vuniquevals = unique(m_vUpperVals);
    
    N = hist(m_vUpperVals, m_vuniquevals);
    m_vdf = N'/length(m_vUpperVals);
    m_vCDF = cumsum(m_vdf);
    
    m_vdiff = m_vuniquevals(2:end) - m_vuniquevals(1:end-1);
    
    % A(k)
    m_vAk(m_nselK) = sum(m_vdiff.*m_vCDF(2:end));
    
    %- Cophenetic correlation coefficient
    Y = 1 - m_mSim;
    Y = Y(~triu(ones(size(m_mSim))))';
    
    [c, D] = cophenet(Z,Y);
    r = corr(Y',D','type','Pearson');
    
    m_vCCr(m_nselK) = r;
end
    
%- Calcluate the proportion change in the area under the consensus CDF
m_vDeltaK_num1 = zeros(length(m_vKs),1);
m_vDeltaK_num1(1) = m_vAk(1);

m_vDeltaK_num1(2:end-1) = ( m_vAk(3:end)-m_vAk(2:end-1) )./m_vAk(3:end);
m_vDeltaK_num1(end) = NaN;

m_vDeltaK_num2 = zeros(length(m_vKs),1);

m_vAhatk = zeros(length(m_vKs),1);
for m_ni = 1:length(m_vKs)
    m_vAhatk(m_ni) = max(m_vAk(1:m_ni));
end

m_vDeltaK_num2(1) = m_vAhatk(1);
m_vDeltaK_num2(2:end-1) = ( m_vAhatk(3:end)-m_vAhatk(2:end-1) )./m_vAhatk(3:end);
m_vDeltaK_num2(end) = NaN;

disp('>>>>> K ----- CC ----- Ak ----- Delta')
disp([m_vKs', m_vCCr, m_vAk, m_vDeltaK_num2])

fig_open = figure('visible','off');
plot(m_vKs(1:end-1), m_vDeltaK_num2(1:end-1), '-bs',...
    'LineWidth',2,'MarkerEdgeColor','k',...
    'MarkerFaceColor','g','MarkerSize',10)
title('Relative Change in Area (\Delta) under the CDF plot');
ylabel('\Delta(K)');
xlabel('Number of subtypes(K)');

m_strname = [m_str_savePath, 'Delta_Ks.pdf'];
saveas(fig_open, m_strname)
close all hidden

%--- Final outcome (K is set to 4)
K = 4;
m_nselK = K - 1;

Y = m_cResults.Clusters{m_nselK,1}; 
