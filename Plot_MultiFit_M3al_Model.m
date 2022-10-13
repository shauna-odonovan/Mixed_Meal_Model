function out=Plot_MultiFit_M3al_Model(fitting,input_data,time)
%Plot multiple fittings of the Mixed Meal Model to the measured meal challenge
%test data generated using Fit_M3al_Model_LatinHyperCube function. 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% fitting    - output of Fit_M3al_Model_MultipleInitial, a struct
%              containing .p_opt which is an array of parameter set estimates
% input_data - struct of measured meal challenge test data used for fitting
%              should contain time series of glucose, insulin,
%              triglyceride, NEFA. Should also specify body weight and meal
%              composition. 
% time       - specify simulation time for plotting. 
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for further information contact Shauna O'Donovan at
% s.d.odonovan@tue.nl
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p_opt_array = fitting.p_opt;

color_plot = colormap(parula(numel(p_opt_array(:,1))));

for j=1:numel(p_opt_array(:,1))
    Plot_M3al_Model(input_data,p_opt_array(j,:),time,color_plot(j,:));
end