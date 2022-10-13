% Script to evaluate parameter identifiability using profile likelihood
% analysis in accordance with the method proposed by Raue et al.
% Bioinformatics (2009). 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial_parameters - vector of opitmised parameter values
% input_data - structure of measured meal challenge test data to which the
%              Mixed Meal Model will be fit - should contain time series of
%              glucose, insulin, triglyceride, and NEFA as well as the meal
%              composition and body weight. 
% time       - vector specifying time span to be simulated for the error
%              calculation - should corrispond to value used in
%              optimisation. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% each parameter will be varied in turn, the remaining parameters will be
% re-estimated. Consquenly this code can run for a long time. 
% A parameter is deemed identifiable if the value of the
% cost function increases as you move away from the optimal value. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for further information please contact Shauna O'Donovan at
% s.d.odonovan.tue.nl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%specify inital parameter guess for lsqnonlin algorithm (only parameters
%being estimated)
initial_parameters = [0.0105,0.0424,2.2975,0.00045,0.0385,0.0713,208.88,0.0163,0.0119];
%specify upper and lower bound for each parameter
lower_bounds = [0.005,0,0,0,0,0,60,0.005,0];
upper_bounds = [0.1,1,5,1,10,1,600,0.1,1];

%Specify simulation time for error function
time = 0:1:720;

%define measure meal challenge test being fitted
input_data = sample_data;

% re-optimise parameter fitting to ensure a local minima has been found
% (this step can be removed)
fitting = Fit_M3al_Model(initial_parameters,sample_data,time,lower_bounds,upper_bounds);
p_opt=fitting.p_opt;

% specify parameter names for plotting purposes
par_name={'k_{1}','k_{5}','k_{6}','k_{11}','K_{ATL}','k_{12}','tau_{LPL}','k_{14}','k_{16}'}; 


h=figure();
for i= 1:9
    %define options for PLA algorithm
    name=par_name{i};
    out.name=zeros(numel(p_opt),2);
    threshold = chi2inv(0.05,8);
    lsq_options=optimset('Algorithm','trust-region-reflective','MaxFunEvals',500,'TolX',1e-8,'Display','iter','UseParallel',1);
    %run PLA for decreasing parameter values
    minStep=0.001;
    maxStep=0.01;
    minChange=0.001;
    maxChange=0.01;
    nr=25;
    [plPar_decreasing,plRes_decreasing]=PLA(@(p_opt)M3al_Model_ErrorFunc(p_opt,input_data,time),p_opt,i,lower_bounds(i),threshold,lower_bounds,upper_bounds,lsq_options,minStep,maxStep,minChange,maxChange,nr);
    %run PLA for increasing parameter values
    minStep=0.001;
    maxStep=0.01;
    minChange=0.001;
    maxChange=0.01;
    nr=25;
    [plPar_increasing,plRes_increasing]=PLA(@(p_opt)M3al_Model_ErrorFunc(p_opt,input_data,time),p_opt,i,upper_bounds(i),threshold,lower_bounds,upper_bounds,lsq_options,minStep,maxStep,minChange,maxChange,nr);
    
    p_opt_res=sum(M3al_Model_ErrorFunc(p_opt,input_data,time).^2);
    subplot(3,3,i)
    plot(plPar_decreasing,plRes_decreasing,'-','Color',[0 0.4470 0.7410],'LineWidth',1.5)
    hold on
    plot(plPar_increasing,plRes_increasing,'-','Color',[0 0.4470 0.7410],'LineWidth',1.5)
    plot(p_opt(i),p_opt_res,'x','Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
    title(par_name(i));
end
saveas(h,'PLA.fig');