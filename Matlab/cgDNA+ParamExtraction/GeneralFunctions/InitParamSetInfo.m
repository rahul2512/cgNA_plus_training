function ParamSetInfo = InitParamSetInfo(Opt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialise the type of nucleic acid 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ParamSetInfo.NA = Opt.NA ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Provide palindromic dimers for given type of nucleic acid 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(Opt.NA,'DNA')
	ParamSetInfo.SymDimer = {'AT','GC','TA','CG'};

elseif strcmp(Opt.NA,'RNA')
	ParamSetInfo.SymDimer = {'AU','GC','UA','CG'};

elseif strcmp(Opt.NA,'MDNA')
	ParamSetInfo.SymDimer = {'AT','GC','TA','CG','MN','NM'};

elseif strcmp(Opt.NA,'HDNA')
	ParamSetInfo.SymDimer = {'AT','GC','TA','CG','HK','KH'};

elseif strcmp(Opt.NA,'DRNA')
	ParamSetInfo.SymDimer = {};

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize lists of dimers 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[d_list,dc_list,ind,tot] = DimerLists(Opt.NA) ;

[int,ends]=TestInputDimerInd(ind,tot,Opt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize lists of active dimers and their complementaries 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ParamSetInfo.activelist_int = d_list(int);
ParamSetInfo.activelist_complementary_int = dc_list(int);
ParamSetInfo.activelist_end = d_list(ends);
ParamSetInfo.activelist_complementary_end = dc_list(ends);

Nint = length(ParamSetInfo.activelist_int);
Nend = length(ParamSetInfo.activelist_end);

ParamSetInfo.twomerdim = [42*ones(1,Nint) , 36*ones(1,Nend) ] ;

Dint = 42;
Dend = 36;

dimSigma =  Nint*Dint+Nend*Dend ; 
dimStiff = (Nint*Dint*(Dint+1))/2 + (Nend*Dend*(Dend+1))/2 ;
ParamSetInfo.dimVec = [dimSigma dimStiff] ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load and store the initial parameter set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ParamSet = load(Opt.InitGuess);
ParamSetInfo.initvalues = ParamSet ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize weights and kT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ParamSetInfo.Wts = Opt.Wts;
ParamSetInfo.kT  = Opt.kT;

end

function [int,ends]=TestInputDimerInd(ind,tot,Opt)

int = Opt.activelist_int_id ;
ends = Opt.activelist_end_id ;

if max(int) > ind
    error('Invalid number in activelist_int_id')
elseif max(ends) > tot
    error('Invalid number in activelist_end_id')
end

end
