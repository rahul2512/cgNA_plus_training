function G = ModelGmat(ParamSetInfo)

TwoDim = ParamSetInfo.twomerdim ;

dimVec = ParamSetInfo.dimVec;
G = ones(dimVec(1),1); 
tmpG = [];
for i = 1:length(TwoDim);
   
    A = ones(TwoDim(i));
    A = triu(2*A - diag(diag(A))) ;
    A = A(:);
    A = A(A~=0);
    
    tmpG = [tmpG ; A];
    
end

G = [G ; tmpG];

end
