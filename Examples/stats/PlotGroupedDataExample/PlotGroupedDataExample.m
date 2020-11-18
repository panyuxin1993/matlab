%% Plot Grouped Data  

%% 
% Load the sample data. Create data vector |x| from the first column of
% the data matrix, which contains sepal length measurements from three species
% of iris flowers. Create data vector |y| from the second column of the
% data matrix, which contains sepal width measurements from the same flowers. 
load fisheriris.mat;
x = meas(:,1);
y = meas(:,2);  

%% 
% Create a scatter plot and six kernel density plots to visualize the relationship
% between sepal length and sepal width, grouped by species. 
scatterhist(x,y,'Group',species,'Kernel','on')    

%%
% The plot shows that the relationship between sepal length and width varies
% depending on the flower species.   
