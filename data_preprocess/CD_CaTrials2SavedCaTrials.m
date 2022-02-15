close all;
[num,txt,raw] =xlsread('C:\Users\PYX\Documents\DataSummary\imaging_data_summary.xlsx');%criteria to choose sessions come from this file
T=cell2table(raw(2:end,1:15));
T.Properties.VariableNames=strrep(raw(1,1:15),' ','_');%table variable name can't have ' ',so replace them
Tchoose=T([139:142],:);
patharray=Tchoose.file_path;
cellfun(@fCD_CaTrials2SavedCaTrials,patharray);