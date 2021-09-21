function [m_vPreds, m_mOuts] = myfunc_Prediction_W(m_mTstVals, m_vWeights, m_mCodeMat)
    %--------------------------------------------------------------------------
    %- Precition with Aggregtaion Weights
    %--------------------------------------------------------------------------
    m_nN = size(m_mTstVals, 1);  
    [m_nK, m_nM] = size(m_mCodeMat);
    
    %--- Proc: calculation of difference() 
    %- Outputs    
    m_mCurInX = zeros(m_nN, m_nM);    
    m_mDist = zeros(m_nN, m_nK);
    
    for m_ni = 1:m_nK
        m_vidxzero = m_mCodeMat(m_ni, :)==0;
        m_vidxnans = isnan(m_mCodeMat(m_ni, :));

        m_mCurInX(:,:) = m_mTstVals;
        m_mCurInX(:, m_vidxzero) = 1 - m_mTstVals(:,m_vidxzero);
        
        m_mCurInX(:,~m_vidxnans) = -log(m_mCurInX(:,~m_vidxnans) + eps);
        m_mCurInX(:,m_vidxnans) = 0;
        
        m_mDist(:, m_ni) = m_mCurInX*m_vWeights;
    end
    %--- end of Proc() 
    
    
    %- converting distance to prob.:
    m_mOuts = exp(-m_mDist);
    m_mOuts = m_mOuts./repmat(sum(m_mOuts, 2),[1, m_nK]);
           
    [~, m_vPreds] = max(m_mOuts,[],2);
end  %--- end of function


















