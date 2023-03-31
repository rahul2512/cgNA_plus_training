function [ColVec, RhsVal] = ObjFunHessCol_Fisher(mdData,UnitData,ParamSetInfo,G)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = length(mdData) ;

ColMat = InitParamStructure(ParamSetInfo.NA);

kT  = ParamSetInfo.kT;
Wts = ParamSetInfo.Wts;

RhsVal = 0 ;

for mu = 1:M
    
    mdW = mdData(mu).shape ;
    mdK = mdData(mu).stiff ;    
    mdSigma = mdData(mu).sigma ;
    mdKinv = mdData(mu).stiff_inv ;
    
    Kprime = UnitData(mu).stiff ;
    Sprime = UnitData(mu).sigma ;
    
    % mdK\Kprime
    n = length(Kprime);
    [row,col,~] = find(Kprime);
    row_S = find(Sprime);
   
    % Contribution of Sigma prime
    if isempty(row) && isempty(row_S)~=1
        % mdK\Sprime
        mdSprime = zeros(n,1);
        for i = 1:length(row_S)
            mdSprime = mdSprime + mdKinv(:,row_S(i))*Sprime(row_S(i));
        end
        mdSprime = (1/kT)*mdSprime;        
        HessStiff = -mdSprime*mdW' ;

        HessSigma = mdSprime ;
        HessStiff = 0.5*(HessStiff + HessStiff') ;
        RhsVal = RhsVal + sum(sum(HessStiff.*mdK)) + (HessSigma)'*mdSigma ;
        
    % Contribution of K prime
    elseif isempty(row)~=1
        row_id = repmat(1:n,[ length(row) 1 ])';
        col_id = repmat(col',[ n 1 ]);
        i = row_id(:);
        j = col_id(:);
        s = mdKinv(:,row);
        for id = 1:length(row)
            if Kprime(row(id),col(id))< 0 
                s(:,id) = -s(:,id);
            end    
        end
               
        s = s(:);
        % mdKprime = mdK\Kprime
        mdKprime = sparse(i,j,s,n,n);
        
        % mdKprime/mdK
        col = unique(col);
        A = mdKprime(:,col)*mdKinv(col,:);
        
        % mdKprime*mdW
        B = mdKprime(:,col)*mdW(col);
        HessStiff = 0.5*A + (1/kT)*B*mdW';        
        HessSigma = -(1/kT)*B ;
        HessStiff = 0.5*(HessStiff + HessStiff') ;        
        RhsVal = RhsVal + sum(sum(HessStiff.*mdK)) + (HessSigma)'*mdSigma ;
                
    else
        continue
    end
    seqInfo = mdData(mu).seqInfo;
    for d = 1:size(seqInfo,1)
        dimerInfo = seqInfo{d,:};
        dimer = dimerInfo{1};
        type  = dimerInfo{2};
        id    = dimerInfo{3};
        sym   = dimerInfo{4};
        tmpHessStiff = Wts(mu)*sym*HessStiff(id,id)*sym';
        tmpHessSigma = Wts(mu)*sym*HessSigma(id);
        stiff_type = [ 'stiff_' type ] ;
        sigma_type = [ 'sigma_' type ] ;
        ColMat.(stiff_type).(dimer) = ColMat.(stiff_type).(dimer) + tmpHessStiff ;
        ColMat.(sigma_type).(dimer) = ColMat.(sigma_type).(dimer) + tmpHessSigma ;
    end
end
ColVec = ModelMatToVec(ColMat,ParamSetInfo) ;
ColVec = G.*ColVec ;
end
