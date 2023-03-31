function ParamSet = InitParamStructure(NA)

list = PrepList(zeros(42,42),NA);
ParamSet.stiff_int = struct(list{:});
list = PrepList(zeros(42,1),NA);
ParamSet.sigma_int = struct(list{:});

if ismember(NA, {'DNA','RNA','DRNA'})
        list1 = PrepList(zeros(36,36),NA);
        list2 = PrepList(zeros(36,1),NA);

elseif ismember(NA, {'MDNA','HDNA'})
        %% as there are no end-blocks for modefied basepair step
        list1 = PrepList(zeros(36,36),'DNA');
        list2 = PrepList(zeros(36,1),'DNA');
end

ParamSet.stiff_end5 = struct(list1{:});
ParamSet.stiff_end3 = struct(list1{:});
ParamSet.sigma_end5 = struct(list2{:});
ParamSet.sigma_end3 = struct(list2{:});

end

function list = PrepList(value,NA)

if strcmp(NA,'DNA')
        dimerlist = @(value){'AT',value,'GC',value,'TA',value,'CG',value,'GT',value,'TG',value,'AG',value,'GA',value,'AA',value,'GG',value,'AC',value,'CA',value,'CT',value,'TC',value,'TT',value,'CC',value};

elseif strcmp(NA,'RNA')
        dimerlist = @(value){'AU',value,'GC',value,'UA',value,'CG',value,'GU',value,'UG',value,'AG',value,'GA',value,'AA',value,'GG',value,'AC',value,'CA',value,'CU',value,'UC',value,'UU',value,'CC',value};

elseif strcmp(NA,'DRNA')
        dimerlist = @(value){'AT',value,'GC',value,'TA',value,'CG',value,'GT',value,'TG',value,'AG',value,'GA',value,'AA',value,'GG',value,'AC',value,'CA',value,'CT',value,'TC',value,'TT',value,'CC',value};

elseif strcmp(NA,'MDNA')
        dimerlist = @(value){'AT',value,'GC',value,'TA',value,'CG',value,'MN',value,'NM',value,'GT',value,'TG',value,'AG',value,'GA',value,'AA',value,'GG',value,'AM',value,'TM',value,'GM',value,'CM',value,'MA',value,'MT',value,'MG',value,'MC',value,'AC',value,'CA',value,'CT',value,'TC',value,'TT',value,'CC',value,'NT',value,'NA',value,'NC',value,'NG',value,'TN',value,'AN',value,'CN',value,'GN',value};

elseif strcmp(NA,'HDNA')
        dimerlist = @(value){'AT',value,'GC',value,'TA',value,'CG',value,'HK',value,'KH',value,'GT',value,'TG',value,'AG',value,'GA',value,'AA',value,'GG',value,'HG',value,'AH',value,'TH',value,'GH',value,'CH',value,'AC',value,'CA',value,'CT',value,'TC',value,'TT',value,'CC',value,'CK',value,'KT',value,'KA',value,'KC',value,'KG',value};

end

list = dimerlist(value);

end
