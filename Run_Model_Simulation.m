%To simulate the Mixed Meal Model for specific parameters and other
%phenotypic traits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Specify phenotypic traits necessary for model  simulation
sample_person.glucose = 5;    %fasting glucose (mmol/l)
sample_person.insulin = 18;   %fasting insulin (uIU/ml)
sample_person.TG      = 1.3;  %fasting plasma triglyceride (mmol/l)
sample_person.NEFA    = 0.33; %fasting plasma NEFA(mmol/l)
sample_person.BW      = 84.2; %body weight (kg)

%% specify composition of meal being simulated
sample_person.meal.G  = 75000; %mass of glucose in meal (mg)
sample_person.meal.TG = 60000; %mass of triglyceride in meal (mg)

%% specify value for each model parameter
%glucose + insulin parameters (Rozendaal et al. (2018))
parameters(1) = 0.0105;  %k1 rate constant for glucose stomach emptying (fast)[1/min]
parameters(2) = 0.28;      %k2 rate constant for glucose appearence from gut [1/min]
parameters(3) = 6.07e-3;   %k3 rate constant for suppresstion of hepatic glucose release by change of plasma glucose
parameters(4) = 2.35e-4;   %k4 rate constant for suppression of hepatic glucose release by delayed (remote) insulin
parameters(5) = 0.0424;  %k5 rate constant for delayed insulin depedent uptake of glucose
parameters(6) = 2.2975;  %k6 rate constant for stimulation of insulin production by the change of plasma glucose concentration (beta cell funtion)
parameters(7) = 1.15;      %k7 rate constant for integral of glucose on insulin production (beta cell function)
parameters(8) = 7.27;      %k8 rate constant for the simulation of insulin production by the rate of change in plasma glucose concentration (beta cell function)
parameters(9) = 3.83e-2;   %k9 rate constant for outflow of insulin from plasma to interstitial space
parameters(10) = 2.84e-1;  %k10 rate constant for degredation of insulin in remote compartment
parameters(11) = 1.4;      %sigma shape factor (appearance of meal)
parameters(12) = 13.2;     %Km michaelis-menten coefficient for glucose uptake
parameters(13) = sample_person.glucose(1);%G_b basal plasma glucose [mmol/l]
parameters(14) = sample_person.insulin(1); %I_PL/_b basal plasma glucose [microU/ml]
parameters(15) = 0.043;    %EGP_bbasal hepatic glucose release

%triglyceride + NEFA parameters (new)
parameters(16) = 30;       %f_spill - fractional spillover of LPL derived NEFA
parameters(17) = 0.00045; %k11 - rate coeficient LPL lipolysis sips/jelic models
parameters(18) = 0.215;    %ATL_max maximum rate of ATL lipolysis in adipose tissues
parameters(19) = 0.0385; %K_ATL michealis menten coeficient for ATL lipolysis of store TG in adipose tissue
parameters(20) = 0.0713; %k12 rate constanst for uptake of NEFA into tissues (currently insulin indenpendent)
parameters(21) = 208.88; %tau_LPL time delay for insulin stimulation of LPL lipolysis
parameters(22) = 0.0088;   %k13 - rate constant for stomach emptying TG(very slow)
parameters(23) = 0.0163; %k14 - rate constant for rate of TG appearance from gut
parameters(24) = 1e-5;     %k15 coefficient for inhibition of TG secretion from liver by insulin
parameters(25) = 0.0119; %k16 basal secretion of TG from liver

%% visualise model simulation for the specifies parameter values

%specify plotting colour
plot_colour = [0, 0.4470, 0.7410]; 
%specify simulation time
time = 0:1:500; 
%generate plot of model simualtion. 
figure()
Simulate_M3al_Model(sample_person,parameters,time,plot_colour);