function H = ObjFunHess(reData,mdData,ParamSetInfo,Gmat)

dimVec = ParamSetInfo.dimVec ;

fprintf('Computing hessian ... ') ;

myCluster = parcluster('local') ;
myCluster.NumWorkers = 24;

parpool(24) ;

H = zeros(dimVec,dimVec);
a = 'uninowebdkjwhdqj.h'
for index=1:dimVec
  tic
  UnitParamVec = zeros(dimVec,1) ;
  UnitParamVec(index) = 1 ;

  UnitParamSet = ModelVecToMat(UnitParamVec,ParamSetInfo,'unit');
  
  UnitData = ModelReconstr(mdData,UnitParamSet,1);
  
  ColVec = ObjFunHessCol(reData,mdData,UnitData,ParamSetInfo,Gmat) ;

  H(:,index) = ColVec ;
  toc
end

fprintf('done. \n') ;

delete(gcp)

end



