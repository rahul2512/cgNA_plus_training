function Vec = ModelMatToVec(ParamSet,ParamSetInfo)

[ id_int, id_end ] = getIndexSymMatrix();

list_int = ParamSetInfo.activelist_int ;
nbr_int = length(list_int);

list_end = ParamSetInfo.activelist_end ;
nbr_end = length(list_end);

Vec_sigma = nan(42,nbr_int + nbr_end); 
Vec_stiff = nan(903,nbr_int + nbr_end);

for i = 1:nbr_int
   
    Vec_sigma(:,i) = ParamSet.sigma_int.(list_int{i});
    Tmp_Vec_Stiff = triu(ParamSet.stiff_int.(list_int{i}));
    Tmp_Vec_Stiff = Tmp_Vec_Stiff(:);
    
    Vec_stiff(:,i) = Tmp_Vec_Stiff(id_int);
        
end

for i = 1:nbr_end
   
    Vec_sigma(1:36,i+nbr_int) = ParamSet.sigma_end5.(list_end{i});
    Tmp_Vec_Stiff = triu(ParamSet.stiff_end5.(list_end{i}));
    Tmp_Vec_Stiff = Tmp_Vec_Stiff(:);
    
    Vec_stiff(1:666,i+nbr_int) = Tmp_Vec_Stiff(id_end);
       
end

Vec = [Vec_sigma(:) ; Vec_stiff(:)];
Vec = Vec(isnan(Vec)~=1);
end

function [ id_int, id_end ] = getIndexSymMatrix()

A = triu(ones(42,42)); 
B = triu(ones(36,36)); 

id_int = find(A);
id_end = find(B);


end