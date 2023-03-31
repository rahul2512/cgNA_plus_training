function [H, Rhs] = ObjFunHess_Fisher(mdData,ParamSetInfo,G)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dimVec = sum(ParamSetInfo.dimVec) ;

fprintf('Computing hessian ... ') ;

myCluster = parcluster('local') ;
myCluster.NumWorkers = 20;
parpool(20) ;

H = zeros(dimVec,dimVec);
%for index=1:dimVec
parfor index=1:dimVec
index
  UnitParamVec = zeros(dimVec,1) ;
  UnitParamVec(index) = 1 ;

  UnitParamSet = ModelVecToMat(UnitParamVec,ParamSetInfo,'unit');
  UnitData = ModelReconstr(mdData,UnitParamSet,1);

  [ColVec, RhsVal] = ObjFunHessCol_Fisher(mdData,UnitData,ParamSetInfo,G) ;

  H(:,index) = ColVec ;

  Rhs(index,1) = RhsVal ;
end

fprintf('done. \n') ;

delete(gcp)

end

