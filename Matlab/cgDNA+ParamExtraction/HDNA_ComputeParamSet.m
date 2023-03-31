%--------------------------------------------------------------------------
% cgDNA+_ParamExtr function = ComputeParamSet()
%--------------------------------------------------------------------------
% Computes a cgDNA+ parameter set using the data set by SetUpCompInterior()
% using a preconditioned gradient descent which preconditioner is the
% Fisher information matrix computed as the second derivative of the
% sum of the Kullback-Liebler divergences between oligomer-based MD
% statistics and cgDNA+ reconstructions.
%--------------------------------------------------------------------------


%% Initialize variables
HDNA_SetUpComputation();

%% Assemble the Fisher matrix and the Rhs of the Fisher System
[H, Rhs] = ObjFunHess_Fisher(mdData,ParamSetInfo,G);

% Compute pseudo inverse of the Fisher matrix which will be used as
% preconditioner of the gradient flow (very slow in high dimension)
fprintf('computing Hinv') ; 
Hinv = P*pinv(P*H*P)*P ;
fprintf('Done computing Hinv') ;
%load ./DNA_BSTJ_cgF/4mus_endsFisherSystem.mat

% Set initial guess for the gradient flow.
% If no initial guess is given the solution of the Fisher system is used

if isempty(ParamSetInfo.initvalues)
  
    ParamVec = Hinv*Rhs;
    ParamSetInfo.initvalues = ModelVecToMat(ParamVec,ParamSetInfo);
    fprintf('Using Fisher solution');

else   
    ParamVec = ModelMatToVec(ParamSetInfo.initvalues,ParamSetInfo);
    fprintf('Using Initial solution'); 
end


% Save Fisher matrix, Fisher matrix pseudoinverse and Rhs
filename = [ './' Opt.ProjectName '/' Opt.RunName 'FisherSystem.mat' ];
save(filename,'H','Hinv','Rhs') ;

% Clear useless variable
%clear H Ths

%% Fisher-preconditioned gradient flow

% Set interation counter and max number of iterations
count = 1 ;
countmax = 10000000;

% Set the setpsize to use and the tolerance to reach with the given
% stepsize
stepsize =  [5e-6, 1e-5, 5e-5, 1e-4, 5e-3, 5e-2,     1e-1];
stepsize = 10* [5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 5e-3,     1e-2];
tolerance = [300 ,  100,   50,   10,    1, 1e-6, 1e-12];

% Compute the parameter set
Descent_Fisher();

%% Save final results
ParamSet = ModelVecToMat(ParamVec,ParamSetInfo) ;

filename = [ './' Opt.ProjectName  '/test.mat' ];

save(filename,'-struct','ParamSet') ;
