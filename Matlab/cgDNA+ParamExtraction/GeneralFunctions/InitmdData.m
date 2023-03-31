function mdData = InitmdData(Opt,ParamSetInfo)

load(Opt.PathData)

if isempty(Opt.whichseqs)
    mdData = olig;
else
    mdData = olig(Opt.whichseqs);
end

clear olig

a_list = ParamSetInfo.activelist_int;
c_list = ParamSetInfo.activelist_complementary_int;

a_e_list = ParamSetInfo.activelist_end;
c_e_list = ParamSetInfo.activelist_complementary_end;

stiff_name = Opt.stiff_name;
shape_name = Opt.shape_name;

[ E_42 , E_36 ] = Etransform();

M = length(mdData);

for mu = 1:M
    seq = mdData(mu).seq;
    nbp = mdData(mu).nbp;
    count = 1;
    for i = 1:length(seq)-1
        
        d = seq(i:i+1);
        
        if(i==1)
            
            res = contains(a_e_list, d);
            if sum(res) > 0
                type = 'end5';
                range = 1:36 ;
                sym = speye(36);
                
            else
                continue
            end
            
        elseif(i==(nbp-1))
            
            res = contains(c_e_list, d);
            if sum(res) > 0
                type = 'end5';
                range = 24*nbp-18-35:24*nbp-18 ;
                sym = E_36;
                d = a_e_list{res};
                
            else
                continue
            end
            
        else
            res = contains(c_list, d);
            if sum(contains(a_list, d))>0
                
                type = 'int';
                range = 24*(i-2) + 19: 24*(i-1) + 36 ;
                sym = speye(42);
                
            elseif sum(res)>0
                
                type = 'int';
                range = 24*(i-2) + 19: 24*(i-1) + 36 ;
                sym = E_42;
                d = a_list{res};
                
            else
                continue
            end
            
        end
        
        tmp{count,:} = {d, type, range, sym};
        count = count + 1;
    end
    
    mdData(mu).seqInfo= tmp;
    
    [m,n] = size(mdData(mu).(shape_name));
    
    if m<n
        mdData(mu).shape = mdData(mu).(shape_name);
    end
    
    mdData(mu).stiff = sparse(mdData(mu).(stiff_name));
    
    mdData(mu).stiff_inv = mdData(mu).stiff\eye(m+n-1);
    mdData(mu).sigma = mdData(mu).stiff*mdData(mu).shape;
end

mdData = rmfield(mdData,{'s1b','stiff_me','nsnap'});

end

