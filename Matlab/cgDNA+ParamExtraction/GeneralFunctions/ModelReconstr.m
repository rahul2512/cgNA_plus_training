function reData = ModelReconstr(mdData,ParamSet,flagSigma)

if nargin < 3 
    flagSigma = 0;
end

M = length(mdData) ; % M = number of sequences

if flagSigma == 1
    type = 'sigma';
else
    type = 'shape';
end

for m = 1:M
    [tmpVec,tmpStiff] = constructSeqParms(mdData(m).seq, ParamSet,flagSigma);
    reData(m).(type) = tmpVec; 
    reData(m).stiff = tmpStiff;
end


end
