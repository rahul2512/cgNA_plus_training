%--------------------------------------------------------------------------
% cgDNA+_ParamExtr function = Descent_Fisher()
%--------------------------------------------------------------------------
% This function implements the preconditioned gradient flow algorithm.
%--------------------------------------------------------------------------

%% Compute value of gradient of the initial guess
% Convert the initial guess in vectorial for to initial guess in parameter
% set form 
ParamSet = ModelVecToMat(ParamVec,ParamSetInfo) ;

% Reconstruct the model prodecitions for initial guess
reData = ModelReconstr(mdData,ParamSet);

% Compute gradient 
GradVec = ObjFunGrad(reData,mdData,ParamSetInfo,G);

% The P matrix ensure that the palindromic symmetries are respected
GradVec = P*GradVec ;
fprintf('Step number : %g,  ', count) ;
fprintf('GradF-norm : %g  \n', norm(GradVec,inf)) ;

nstep = length(stepsize);
% For loop over the number of computation steps
for s = 1:nstep
    
    alpha = stepsize(s);
    tol = tolerance(s);
    
    % While loop over the value of GradVec norm
    while(norm(GradVec,inf)>tol && count<=countmax)
        
        % Preconditioned Gradien flow step
        ParamVec = ParamVec - alpha*Hinv*GradVec ;
        
        % Convert current vectorial solution to parameter set type
        ParamSet = ModelVecToMat(ParamVec,ParamSetInfo) ;
        
        % Reconstruct the model prodecitions for current solution
        reData = ModelReconstr(mdData,ParamSet);
        
        % Compute the gradient 
        GradVec = ObjFunGrad(reData,mdData,ParamSetInfo,G);
        
        % Ensure the palindromic symmetries
        GradVec = P*GradVec ;
        
        fprintf('Step number : %g,  ', count) ;
        fprintf('GradF-norm : %g  \n', norm(GradVec,inf)) ;
        
        count = count + 1 ;
    end
    
end
