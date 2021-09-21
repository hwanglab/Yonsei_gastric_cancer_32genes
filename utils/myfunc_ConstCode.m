function m_mCodmat = myfunc_ConstCode(m_vTrnY, m_strCode)   
    m_nC = length(unique(m_vTrnY));
                
    switch upper(m_strCode)
        % all pairwise comparison
        case 'APS'
            m_mind = nchoosek(1:m_nC, 2);
            
            m_nN = size(m_mind, 1);
            m_mCodmat = NaN*ones(m_nC,m_nN);
            
            for m_ni = 1:m_nN
                m_mCodmat(m_mind(m_ni,1),m_ni) = 1;
                m_mCodmat(m_mind(m_ni,2),m_ni) = 0;
            end
            
        % one versus rest
        case 'OVA'
            m_mCodmat = eye(m_nC, m_nC);            
        otherwise
    end

end