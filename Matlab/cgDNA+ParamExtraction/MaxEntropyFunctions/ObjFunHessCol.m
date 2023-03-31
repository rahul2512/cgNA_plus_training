function ColVec = ObjFunHessCol(reData,mdData,ArrayDir,ParamSetInfo,Gmat)
% MAX ENTROPY
M = length(mdData) ;
                       
kT  = ParamSetInfo.kT;
Wts = ParamSetInfo.Wts;

ColMat = InitParamStructure();

for mu = 1:M

  reW = reData(mu).shape ;
  reK = reData(mu).stiff ;
  
  mdW = mdData(mu).shape ;
  mdK = mdData(mu).stiff ;
  
  B = mdK ;
  G = (reK\B)/reK ;
  qvec = (1/kT)*reK\B*(reW-mdW) ;

  Kprime = ArrayDir(mu).stiff ;
  Sprime = ArrayDir(mu).sigma ;

  HessStiff = ...
         - 0.5*(reK\Kprime)/reK ...
         + 0.5*(reK\Kprime)*G ...
         + 0.5*(G*Kprime)/reK ...
         + (reK\Kprime*qvec)*reW' ...
         + (1/kT)*(G*Kprime*reW)*reW' ...
         + qvec*(reK\Kprime*reW)' ...
         - (1/kT)*(G*Sprime)*reW' ...
         - qvec*(reK\Sprime)' ;
         
  HessStiff = 0.5*HessStiff + 0.5*(HessStiff)' ;

  HessSigma = - reK\Kprime*qvec - (1/kT)*G*Kprime*reW + (1/kT)*G*Sprime ;
  
  seqInfo = mdData(mu).seqInfo; 
  for d = 1:size(seqInfo,1)
      dimer = seqInfo{d,1};
      type  = seqInfo{d,2};
      id    = seqInfo{d,3};
      sym   = seqInfo{d,4};
      
      tmpHessStiff = Wts(mu)*sym'*HessStiff(id,id)*sym;
      tmpHessSigma = Wts(mu)*sym'*HessSigma(id)*sym;
      
      stiff_type = [ 'stiff_' type ];
      sigma_type = [ 'sigma_' type ];
      
      ColMat.(stiff_type).(dimer) = ColMat.(stiff_type).(dimer) + tmpHessStiff ;
      ColMat.(sigma_type).(dimer) = ColMat.(sigma_type).(dimer) + tmpHessSigma ;
  
  end
  
  
end

ColVec = ModelMatToVec(ColMat,ParamSetInfo) ;

ColVec = Gmat*ColVec ;


end
