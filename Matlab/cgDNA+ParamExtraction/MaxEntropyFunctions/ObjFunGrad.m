function GradVec = ObjFunGrad(reData,mdData,ParamSetInfo,G)
% MAX ENTROPY

M = length(mdData) ;

GradMat = InitParamStructure(ParamSetInfo.NA);
kT  = ParamSetInfo.kT;
Wts = ParamSetInfo.Wts;

for mu = 1:M
    reW = reData(mu).shape ;
    reK = reData(mu).stiff ;
    
    mdW = mdData(mu).shape ;
    mdK = mdData(mu).stiff ;
    
    [col,row] = size(reK);
    I = eye(col,row);
    reKinv = reK\I ;
    dFdreK =  0.5*(reKinv - reK\mdK/reK) ;
    dFdreW =  mdK*(reW-mdW)/kT ;
    
    GradSigma = reK\dFdreW ;
    GradStiff = dFdreK - 0.5*(GradSigma*reW' + reW*GradSigma') ;
    
    seqInfo = mdData(mu).seqInfo;
    for d = 1:size(seqInfo,1)
        dimerInfo = seqInfo{d,:};
        dimer = dimerInfo{1};
        type  = dimerInfo{2};
        id    = dimerInfo{3};
        sym   = dimerInfo{4};
        
        tmpGradStiff = Wts(mu)*sym*GradStiff(id,id)*sym';
        tmpGradSigma = Wts(mu)*sym*GradSigma(id);
        
        stiff_type = [ 'stiff_' type ];
        sigma_type = [ 'sigma_' type ];
        
        GradMat.(stiff_type).(dimer) = GradMat.(stiff_type).(dimer) + tmpGradStiff ;
        GradMat.(sigma_type).(dimer) = GradMat.(sigma_type).(dimer) + tmpGradSigma ;
        
    end
    
end

GradVec = ModelMatToVec(GradMat,ParamSetInfo) ;
GradVec = G.*GradVec ;

end
