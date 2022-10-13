function out = M3al_Model_ErrorFunc(p_opt,input_data,time)
%Cost function to fit Mixed Meal Model to measured data. The error between 
%the Mixed Meal Model (M3al Model) simulation for a given parameter set and 
%supplied measured meal challenge test data is calculated. Additional
%constraints are applied to ensure certain physioloigical constraints are
%met; 
%i)   All glucose contained in the meal appears in 4 hours following the meal. 
%ii)  All triglyceride contained in the meal appears in 12 hours following the meal.
%iii) At fasting (t=0 mins)NEFA uptake into tissues is in steady state (equal to measured value)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p_opt        - vector of parameter values for which error is being
%                calculated.
% input_data   - struct of measured challenge test data for caculating error
%              - mean and standard deviation values are required for
%                glucose and insulin (need to extend to TG and NEFA)
%              - vector of sampling time points are also required.
%              - any additional variables required for model simulation.
% time         - time span for model simulation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for further information contact Shauna O'Donovan at
% s.d.odonovan@tue.nl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% simulate model for given parameter set

%form full parameter vector for simulation
parameters = M3al_Model_Parameters(p_opt,input_data);
%define intial values and model constants needed for simulation of eDES model
[initial_values,constants]=M3al_Model_Initial(input_data,parameters);

%define global parameters for simulation
global t_saved G_PL_saved;
%initialise gloabl parameters
t_saved = 0;
G_PL_saved = input_data.glucose(1);

%specify options for ODE solver (Integrator function)
ODE_options = odeset('RelTol',1e-5,'OutputFcn',@integratorfunG);

%simulate model
[T,X] = ode45(@M3al_Model_ODE,time,initial_values,ODE_options,parameters,constants,input_data);

%% Calculate error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model fit error - data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%glucose error
measured_time_G=ismember(T,input_data.time_G);
G_err = (X(measured_time_G,2)' - input_data.glucose)./max(input_data.glucose);

%%insulin error
measured_time_I=ismember(T,input_data.time_I);
I_err = (X(measured_time_I,4)' - input_data.insulin)./max(input_data.insulin);

%NEFA error
measured_time_NEFA=ismember(T,input_data.time_NEFA);
NEFA_err = (X(measured_time_NEFA,9)' - input_data.NEFA)./max(input_data.NEFA);

%TG error
measured_time_TG=ismember(T,input_data.time_TG);
TG_err = (X(measured_time_TG,13)' - input_data.TG)./max(input_data.TG);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% regularisation error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%error if AUC of rate of appearence of glucose from gut < meal content. 
k2  = parameters(2);
BW  = input_data.BW;
f_G = constants.f_G;
V_G = constants.V_G;  

G_gut = k2.*(f_G/(V_G*BW)).*X(1:240,1);

AUC_G=trapz(G_gut);
AUC_G_norm = ((V_G*BW)/f_G)*AUC_G;

err_auc_G = abs((AUC_G_norm - input_data.meal.G)./10000);

%error if AUC of rate of appearnce of triglyceride in plasma < meal content.
k14  = parameters(23);
BW   = input_data.BW;
f_TG = constants.f_TG;
V_TG = constants.V_TG; 
 
TG_gut = k14.*(f_TG/(V_TG*BW)).*X(1:480,12);
 
AUC_TG=trapz(TG_gut);
AUC_TG = ((V_TG*BW)/f_TG)*AUC_TG;

err_auc_TG = abs((AUC_TG - input_data.meal.TG)./10000);

%constrain steady state TG to measured fasting value.
TG_12hours = input_data.TG(1)-X(720,13);

%constrain steady state NEFA to measured fasting value. 
k_12    = parameters(20);
spill   = parameters(16);
k_11    = parameters(17);
ATL_max = parameters(18);
K_ATL   = parameters(19);
I_b     = parameters(14);
NEFA_b  = input_data.NEFA(1);
TG_b    = input_data.TG(1);

model_fasting_NEFA = (3.*(spill/100).*k_11.*TG_b.*I_b + (ATL_max./(1+K_ATL.*(I_b^2))))./k_12;

NEFA_difference = NEFA_b - model_fasting_NEFA;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% total error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=[max(input_data.glucose).*G_err,max(input_data.glucose).*I_err,max(input_data.glucose).*NEFA_err,max(input_data.glucose).*TG_err,err_auc_G,8.*NEFA_difference,err_auc_TG,TG_12hours];

