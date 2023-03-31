function ParamSet = ModelVecToMat(Vec,ParamSetInfo,varargin)

[E_42,E_36] = Etransform();

list_int = ParamSetInfo.activelist_int ;
list_int_c = ParamSetInfo.activelist_complementary_int ;
nbr_int = length(list_int);

list_end = ParamSetInfo.activelist_end ;
list_end_c = ParamSetInfo.activelist_complementary_end ;
nbr_end = length(list_end);

%%%%%%%%%%%%%%% this block is added by rahul
if nargin < 3
ParamSet = ParamSetInfo.initvalues ;
%fprintf('Using provided  initial prmset');
end
if nargin==3
ParamSet = InitParamStructure(ParamSetInfo.NA) ;
%fprintf('computing hessian and initiating zero prmset');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%ParamSet = ParamSetInfo.initvalues ;
[ Vec_sigma_int, Vec_sigma_end, Vec_stiff_int , Vec_stiff_end ] = DecomposeVec(Vec,nbr_int,nbr_end);


for i = 1:nbr_int
    
    ParamSet.sigma_int.(list_int{i}) = Vec_sigma_int(:,i);
    
    Tmp_Stiff = reshape(Vec_stiff_int(:,i),[42,42]);
    Tmp_Stiff = Tmp_Stiff + Tmp_Stiff' - diag(diag(Tmp_Stiff));
    ParamSet.stiff_int.(list_int{i}) = Tmp_Stiff;
   
    switch list_int{i}
        
        case {'AT','GC','TA','CG','MN','NM'}
            continue 
        otherwise
            
            ParamSet.sigma_int.(list_int_c{i}) = E_42*Vec_sigma_int(:,i);
            ParamSet.stiff_int.(list_int_c{i}) = E_42*Tmp_Stiff*E_42; 
            
    end
    
end

for i = 1:nbr_end
    
    ParamSet.sigma_end5.(list_end{i}) = Vec_sigma_end(:,i);
    ParamSet.sigma_end3.(list_end_c{i}) = E_36'*Vec_sigma_end(:,i);
    
    Tmp_Stiff = reshape(Vec_stiff_end(:,i),[36,36]);
    Tmp_Stiff = Tmp_Stiff + Tmp_Stiff' - diag(diag(Tmp_Stiff));
    
    ParamSet.stiff_end5.(list_end{i}) = Tmp_Stiff;
    ParamSet.stiff_end3.(list_end_c{i}) = E_36'*Tmp_Stiff*E_36;
    
end


end

function [ Vec_sigma_int, Vec_sigma_end, Vec_stiff_int , Vec_stiff_end ] = DecomposeVec(Vec,nbr_int,nbr_end)
Vec_sigma_int = reshape(Vec(1:nbr_int*42), [42 nbr_int]);
Vec(1:nbr_int*42) = [];

Vec_sigma_end = reshape(Vec(1:nbr_end*36),[36 nbr_end]); 
Vec(1:nbr_end*36) = [];

[ id_int, id_end ] = getIndexSymMatrix();

Vec_stiff_int = zeros(42*42,nbr_int);
Vec_stiff_end = zeros(36*36,nbr_end);

Vec_stiff_int(id_int,:) = reshape(Vec(1:nbr_int*(42*43)/2), [ (42*43)/2 nbr_int ] )  ;
Vec(1:nbr_int*(42*43)/2) = [];
Vec_stiff_end(id_end,:) = reshape(Vec, [ (36*37)/2 nbr_end ] ) ;
end

function [ id_int, id_end ] = getIndexSymMatrix()

A = triu(ones(42,42));
B = triu(ones(36,36));

id_int = find(A);
id_end = find(B);


end
