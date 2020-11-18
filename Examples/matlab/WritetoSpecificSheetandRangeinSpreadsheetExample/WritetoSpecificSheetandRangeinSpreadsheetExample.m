%% Write to Specific Sheet and Range in Spreadsheet  
% Write mixed text and numeric data to an Excel(R) file starting at cell
% |E1| of |Sheet2|.

%%
filename = 'testdata.xlsx';
A = {'Time','Temperature'; 12,98; 13,99; 14,97};
sheet = 2;
xlRange = 'E1';
xlswrite(filename,A,sheet,xlRange)   

