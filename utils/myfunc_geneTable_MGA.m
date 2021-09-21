function T = myfunc_geneTable_MGA...
                        (mv_Pred_tst, mc_Categories, mc_UniqLabels, mc_groups)
    
    M = length(unique(mv_Pred_tst));
    K = length(mc_UniqLabels);
    
    mm_ResMat = zeros(K, length(mc_UniqLabels));
       
    for mn_i = 1:M
        mv_chki = mv_Pred_tst == mn_i;

        for mn_sub = 1:K
            mv_chkj = strcmpi(mc_Categories, mc_UniqLabels{mn_sub});
            mm_ResMat(mn_i, mn_sub) = sum(mv_chki & mv_chkj);        
        end
    end

    mstr_cmd1 = 'T = table(';
    mstr_cmd2 = [];
    for mn_i = 1:K
        mstr_cmd2 = [mstr_cmd2, 'mm_ResMat(:,', num2str(mn_i), '), '];
    end
    mstr_cmd2(end-1:end) = [];
    
    eval([mstr_cmd1,mstr_cmd2, ');']);
    T.Properties.RowNames = mc_groups;
    
    for mn_i = 1:K
        mstr_pos = sprintf('Var%d',mn_i);
        mstr_cmd = ['T.Properties.VariableNames{', char(39), mstr_pos,  char(39), '} = mc_UniqLabels{mn_i};'];
        eval(mstr_cmd);
    end    
end