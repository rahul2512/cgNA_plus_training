function Pmat = ModelPmat(ParamSetInfo)

activelist = ParamSetInfo.activelist_int;
if  ParamSetInfo.NA =='DNA'
SymDimer = {'AT','GC','TA','CG','MN','NM' };
end
if  ParamSetInfo.NA =='MDNA'
SymDimer = {'AT','GC','TA','CG','MN','NM' };
end
if ParamSetInfo.NA =='DRNA'
SymDimer = { };
end
% Good thing to do but not working on lcvmmsrv7
% This is just for testing dioporco
s = contains(SymDimer,activelist);
SymDimer = SymDimer(s);
[ E_42 , ~ ] = Etransform();

SymParamSetInfo = InitSymParamSetInfo(SymDimer);

dimVec = ParamSetInfo.dimVec;
dimVecSym = SymParamSetInfo.dimVec;

diff = dimVec-dimVecSym;

tmpPsym = zeros(sum(dimVecSym));
ZeroVec = zeros(sum(dimVecSym),1);

% Sigma Part
for j = 1:length(SymDimer)
    dimer = SymDimer{j};
    
    first = 1 + 42*(j-1);
    last  = 42 + 42*(j-1);
    for index = first:last
        UnitParamVec = ZeroVec ;
        UnitParamVec(index) = 1 ; 
        UnitParamSet = ModelVecToMat(UnitParamVec,SymParamSetInfo) ;
        UnitParamSet
        TempVec = UnitParamSet.sigma_int.(dimer) ;
        UnitParamSet.sigma_int.(dimer) = 0.5*(TempVec + E_42*TempVec) ;
        Column = ModelMatToVec(UnitParamSet,SymParamSetInfo) ;
        tmpPsym(:,index) = Column ;
        
    end
     
end

% Stiff Part
for j = 1:length(SymDimer)
    dimer = SymDimer{j};
    
    first = 42*length(SymDimer) + 1 + 903*(j-1);
    last  = 42*length(SymDimer) + 903 + 903*(j-1);
    for index = first:last
	numel(ZeroVec)
        UnitParamVec = ZeroVec ;
        UnitParamVec(index) = 1 ; 
        UnitParamSet = ModelVecToMat(UnitParamVec,SymParamSetInfo) ;
                      
        TempMat = UnitParamSet.stiff_int.(dimer) ;
        UnitParamSet.stiff_int.(dimer) = 0.5*(TempMat + E_42*TempMat*E_42) ;
        
        Column = ModelMatToVec(UnitParamSet,SymParamSetInfo) ;
        
        tmpPsym(:,index) = Column ;
    end
     
end

dSig = dimVecSym(1);
SigmaSym = tmpPsym(1:dSig,1:dSig);
SigmaPart = blkdiag(SigmaSym,eye(diff(1)));
StiffSym = tmpPsym(42*length(SymDimer)+1:end,42*length(SymDimer)+1:end);
StiffPart = blkdiag(StiffSym,eye(diff(2)));

Pmat = sparse(blkdiag(SigmaPart,StiffPart));

end

function ParamSetInfo = InitSymParamSetInfo(SymDimer)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize lists of active dimers and their complementaries 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ParamSetInfo.activelist_int = SymDimer;
ParamSetInfo.activelist_complementary_int = SymDimer;

ParamSetInfo.activelist_end = [];
ParamSetInfo.activelist_complementary_end = [];

Nint = length(ParamSetInfo.activelist_int);
Nend = length(ParamSetInfo.activelist_end);

ParamSetInfo.twomerdim = [42*ones(1,Nint) , 36*ones(1,Nend) ] ;

Dint = 42;
Dend = 36;

dimSigma =  Nint*Dint+Nend*Dend ;
dimStiff = (Nint*Dint*(Dint+1))/2 + (Nend*Dend*(Dend+1))/2;
ParamSetInfo.dimVec = [dimSigma dimStiff]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load and store the initial parameter set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ParamSetInfo.initvalues = InitParamStructure() ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize weights and kT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ParamSetInfo.Wts = 0;
ParamSetInfo.kT  = 0;

end

