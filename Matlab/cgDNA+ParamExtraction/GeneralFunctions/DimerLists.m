function [d_list,dc_list,ind,tot] = DimerLists(NA_type)

rna_list  = { ...
            % 4 palindromic dimers 
            'AU','GC','UA','CG', ...
        ... % 6 independent dimers    
            'GU','UG','AG','GA','AA','GG' , ...
        ... % 6 dependent dimer
            'AC','CA','CU','UC','UU','CC' , ...
          } ;

dna_list  = { ...
            % 4 palindromic dimers 
            'AT','GC','TA','CG', ...
        ... % 6 independent dimers    
            'GT','TG','AG','GA','AA','GG' , ...
        ... % 6 dependent dimer
            'AC','CA','CT','TC','TT','CC' , ...
          } ;



mdna_list  = { ...
           % 2 methylated palindromic dimer
           'AT','GC','TA','CG','MN','NM', ...
         ... % 5 independent met dimers    
           'GT','TG','AG','GA','AA','GG','AM','TM','GM','CM','MA','MT','MG','MC', ...
         ... % 5 dependent met dimers    
           'AC','CA','CT','TC','TT','CC','NT','NA','NC','NG','TN','AN','CN','GN', ...
          } ;

hdna_list  = { ...
           % 2 hydroxymethylated palindromic dimer
           'AT','GC','TA','CG','HK','KH', ...
         ... % 5 independent met dimers    
           'GT','TG','AG','GA','AA','GG','HG','AH','TH','GH','CH', ...
         ... % 5 dependent met dimers    
           'AC','CA','CT','TC','TT','CC','CK','KT','KA','KC','KG', ...
          } ;
    if ismember(NA_type, {'DNA','RNA','DRNA'})
            %% In case of new dimer type modify the following dimer list
            % List of possible dimers:
            if NA_type == "DNA"
                d_list = dna_list ;
	        pal = 4;
	        % dc_list is list complementary dimers of d_list
        	[d_list,dc_list,ind,tot] = sets(d_list,pal) ;
    
    	    elseif NA_type == "RNA"
                d_list = rna_list ;
        	pal = 4;
	        % dc_list is list complementary dimers of d_list
        	[d_list,dc_list,ind,tot] = sets(d_list,pal) ;

            elseif NA_type == "DRNA"
                d_list = dna_list ;
	        pal = 4;
	        [d_list,dc_list,ind,tot] = sets(d_list,pal) ;
            end
    end

    if ismember(NA_type, {'MDNA','HDNA'})

       if NA_type == 'MDNA'
               d_list  = mdna_list ; 
       elseif NA_type == 'HDNA'
               d_list  = hdna_list ;
       end

        pal = 6;
        [d_list,dc_list,ind,tot] = sets(d_list,pal) ;

    end

end

function [d_list,dc_list,ind,tot] = sets(d_list,pal)
        tot = length(d_list);
        ind = pal + (tot-pal)/2;
        % List of complementary dimers of d_list
        dc_list = d_list([1:pal, ind+1:tot , pal+1:ind ]) ;
%	if NA_type == "DRNA"
%	dc_list ={ ...
            % 4 palindromic dimers 
%            'AT','GC','TA','CG', ...
        ... % 6 dependent dimer
%           'AC','CA','CT','TC','TT','CC' , ...
        ... % 6 independent dimers    
%            'GT','TG','AG','GA','AA','GG' , ...
%         } ;
%	ind=16
%	end
end
