function [outName,animal,date] = fFormatBehFilename(rawName,format)
%FFORMATBEHFILENAME format .beh file name, since different session may have
%different filename format, e.g. for FP yyyy_mm_dd_name; for 2P,
%name_yyyymmdd. 
%   rawName- string, the original file name, usually contain both animal
%   and date information
%   format- char,e.g.  yyyy_mm_dd_name;, name_yyyymmdd. 
%Output
%   outName- output name
%   animal- name of animal
%   date- date of the session

info=split(rawName,'_');
ind_animal=contains(info,'pyx');%currently only support my animal
ind_animal=find(ind_animal);
if (ind_animal==1)% (length(info)==2)%animal_yymmdd/animal_yymmdd_2P/animal_yymmdd_FP
    date=info{2};
    animal=info{1};
elseif ind_animal==4
    date=strcat(info{1:3});
    animal=info{4};
end

formatchar={'animal','yyyy','mm','dd'};
datachar=cell(1,4);
datachar{1}=animal;
datachar{2}=date(1:4);
datachar{3}=date(5:6);
datachar{4}=date(7:8);
outName=format;
for i=1:length(formatchar)
    outName = strrep(outName,formatchar{i},datachar{i});
end

end

