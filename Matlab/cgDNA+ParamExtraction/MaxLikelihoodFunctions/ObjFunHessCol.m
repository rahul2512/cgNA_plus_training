function ColVec = ObjFunHessCol(reStiff,reShape,mdStiff,mdShape,...
                  StiffDir,SigmaDir,kT,Wts,AssemMap,...
                  ParamSetInfo,ParamVec,Gmat)
% MAX LIKE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ParamSet = ModelVecToMat(ParamVec,ParamSetInfo) ;

[n N] = size(ParamSet) ; % N = number of indep param arrays
                         % n = 1 (unused dimension)

[n M] = size(reStiff) ; % M = number of sequences
                        % n = 1 (unused dimension)

for mu = 1:M
  HQ(mu).array = zeros(size(reStiff(mu).array));
  Hq(mu).array = zeros(size(reShape(mu).array));

  reW = reShape(mu).array ;
  mdW = mdShape(mu).array ;
  reK = reStiff(mu).array ;
  mdK = mdStiff(mu).array ;

  Kprime = StiffDir(mu).array ;
  Sprime = SigmaDir(mu).array ;

  X = reK\eye(558) ;
  B = mdK ;

  g = X*Sprime - X*Kprime*reW ;

  HQ(mu).array = 0.5*X*Kprime*X ...
                 - 0.5*(g*reW')/kT - 0.5*(reW*g')/kT ;

  Hq(mu).array = g/kT ;
end

for pnum = 1:N
  Column(pnum).stiff = zeros(size(ParamSet(pnum).stiff));
  Column(pnum).sigma = zeros(size(ParamSet(pnum).sigma));
  seqlist = AssemMap(pnum).seqlist ;
  rangelist = AssemMap(pnum).rangelist ;
  cwmatlist = AssemMap(pnum).cwmatlist ;
  for k=1:length(seqlist)
    mu = seqlist(k,1) ;
    range = rangelist(k,:) ;
    [sx sy sz] = size(cwmatlist(k,:,:)) ;
    CWmat = reshape(cwmatlist(k,:,:),[sy sz]) ;
    weight = Wts(mu) ;
    Column(pnum).stiff = ...
        Column(pnum).stiff + weight*CWmat*(HQ(mu).array(range,range))*CWmat ;
    Column(pnum).sigma = ...
        Column(pnum).sigma + weight*CWmat*(Hq(mu).array(range)) ;
  end
end

ColVec = ModelMatToVec(Column,ParamSetInfo) ;

ColVec = Gmat*ColVec ;

