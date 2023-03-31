function GradVec = ObjFunGrad(reData,mdData,ParamSetInfo,G)
% MAX LIKE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = length(mdData) ;

GradMat = InitParamStructure(ParamSetInfo.NA);

kT  = ParamSetInfo.kT;
Wts = ParamSetInfo.Wts;

for mu = 1:M
    
    reW = reData(mu).shape ;
    reK = reData(mu).stiff ;
    
    mdW = mdData(mu).shape ;
    mdKinv = mdData(mu).stiff_inv ;
    
    [col,row] = size(reK);
    I = eye(col,row);
    reKinv = reK\I ;
    
    GradSigma = (reW-mdW)/kT ;
    GradStiff = 0.5*(mdKinv - reKinv - (reW*reW')/kT + (mdW*mdW')/kT) ;
    
    seqInfo = mdData(mu).seqInfo;
    for d = 1:size(seqInfo,1)
        dimerInfo = seqInfo{d,:};
        dimer = dimerInfo{1};
        type  = dimerInfo{2};
        id    = dimerInfo{3};
        sym   = dimerInfo{4};
        
%        tmpGradStiff = Wts(mu)*sym'*GradStiff(id,id)*sym;
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

