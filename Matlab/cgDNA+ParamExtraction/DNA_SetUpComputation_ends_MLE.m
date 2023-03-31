%--------------------------------------------------------------------------
% cgDNA+ParamExtr function = SetUpComputation()
%--------------------------------------------------------------------------
% This function is used to initialise the computation of a cgDNA+ parameter
% set 
%--------------------------------------------------------------------------

%% Enter a name of the project and the folder where the results will be
% saved.

Opt.ProjectName = 'DNA_CGF_MLE';
Opt.RunName = '12mus_ends_CGF_BSC1';
%% Enter the name of the Run you want to use. The value will be used to name
% the resulting parameter set. (Give the number in a string format.
Opt.PSNbr = '1.2';

%% Give the type of nucleic acid data you are using 
% 'DNA', 'RNA', 'DRNA' -- for hybrid DNA-RNA  
% 'MDNA' -- for methylyated DNA and 'HDNA' -- for Hydroxymethylyated DNA
Opt.NA = 'DNA'; 

%% Give the path to the training data set you want use to train model.
Opt.PathData = './Data/palin.bscl.tip3p.jc.ends.comb.stats.3mus.cgF.mat';

%% Give a description of the MD training data used to trained the model. 
Opt.Description = 'Which MD protocol?';

%% Give the path to the initial parameter set guess. If empty the initial
% guess will be computed using the Fisher system. 
Opt.InitGuess = 'DNA_CGF_MLE/cgDNA+_MLE_ends12mus_ends_CGF_BSC1.mat';

%% See README for details about the numbering of the dimer steps.

% Give the number corresponding to the independent dimer parameters that
% should be computed. Read the README file to find out index
% 1:10 for interior blocks
Opt.activelist_int_id = []; %1:10 %[]

% Give the number corresponding to the 5' end dimer parameters that should
% be computed.
End_Map = [9 1 7 11 3 6 14 15 10 8 5 4 16 13 12] ;
Opt.activelist_end_id = End_Map(j); %[2]
% Specify which sequences to use in the data set. If empty all the
% sequences in the data set will be used. 
Opt.whichseqs = 4*(j-1)+1:4*j; % [];  %49:52 ;

%% Provide weights for the sum of KLd. Default values Wts = 1. 
% !! The number of weights should correspond to the number of sequences in
% the training data set.  %16,1
Opt.Wts = ones(4,1);

%% Set kT value (=1 in non-dimensional case).
Opt.kT = 1 ;

%% Two possible choiches for the ordering in the KLd
% Position 1) Classic cgDNA. Flag_Ordering = 1
% Position 2) Maximum Likelihood Flag_Ordering = 2
Opt.Flag_Ordering = 2

%% Which stiffness matrix and shape vector to use in the oligomer-based statistic:
% Default:
%stiff_name = 's1b' (s1b_sym if palindromic);
%shape_stiff = 'shape' (shape_sym if palindromic);
Opt.stiff_name = 's1b';
Opt.shape_name = 'shape';

%==========================================================================
%                     DO NOT MODIFY AFTER THESE LINES
%==========================================================================

%% Prepare folder for results
if ~exist(Opt.ProjectName, 'dir')
    mkdir(Opt.ProjectName)
end

%%%% Save the Opt to keep a track of what was done
filename = [ './' Opt.ProjectName '/' Opt.RunName '_Initial_opt.mat' ];
save(filename,'Opt') ;

%% Add to path the folder for the gradient function corresponding to the KLd 
% odrering.
if Opt.Flag_Ordering == 1
    addpath(genpath('./MaxEntropyFunctions'),genpath('./GeneralFunctions'),genpath('./FisherFunctions'))
elseif Opt.Flag_Ordering == 2
    addpath(genpath('./MaxLikelihoodFunctions'),genpath('./GeneralFunctions'),genpath('./FisherFunctions'))
end

%% Initialise all the usefull information of the parameter set using the 
% given options.
ParamSetInfo = InitParamSetInfo(Opt);
%% Initialise the training data set.
mdData = InitmdData(Opt,ParamSetInfo);

%% Initialise the model related matrices.
% The G matrix compensate the change between Frobenius inner product for
% matrices and Euclidean inner product for vectors. See ModelGMat.mat function 
% for details.
G = ModelGmat(ParamSetInfo) ;
% The P matrix will ensure that the palindromic symmetries are respected
% during the computation. See ModelPMat.mat function for details.
P = ModelPmat(ParamSetInfo) ;
