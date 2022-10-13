%Scrpt to perform local parameter sensitivity analysis for the Mixed Meal Model
%
% values for sepcific phenotypic traits and model parameter values can be
% specified below. 
%
% each model parameter will be varied through a range of values proportional to the
% specifed/optimal value for each parameter (provided as a percentage i.e.
% for 50% range = 0.5).
% The number of parameter values in the specified range that will be
% simulated is specified by the variable n. 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for further information contact Shauna O'Donovan at s.d.odonovan@tue.nl
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

%% a vector of opmitimsed parameters can be specified by uncommenting the below line of code. 
%parameters = M3al_Model_Parameters(p_opt,phenotypic_data);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define parameters for the sensitivity analysis

%specify range through which each parameter will be varied
range = 0.5; %50% 
%specify number of parameter values in range to be simulated. 
n=50; 
%specify the time span for the model simulation
t_span = 0:1:500;
%define parameter names for plotting and saving
par_name = {'k_1','k_2','k_3','k_4','k_5','k_6','k_7','k_8','k_9','k_{10}','sigma','K_M','G_b','I_b','EGP_b','f_{spill}','k_{11}','ATL_{max}','K_{ATL}','k_{12}','tau_{LPL}','k_{13}','k_{14}','k_{15}','k_{16}'};


for p=1:25
    %skip parameters with 0 values (cannot calculate range)
    if parameters(p)==0
        continue
    end
    %define range and step size for parameter p
    p_lb=parameters(p)-range*parameters(p);
    p_ub=parameters(p)+range*parameters(p);
    step=(p_ub-p_lb)/(n-1);
    p_values=p_lb:step:p_ub;
    %generate parameter set matrix
    parameter_values=repmat(parameters,n,1);
    parameter_values(:,p)=p_values';
    h=figure();
    %define colour map for plotting
    color_plot=colormap(parula(n));
    j=1;
    while j < n+1    
        try Simulate_M3al_Model(sample_person,parameter_values(j,:),t_span,color_plot(j,:));
            j=j+1;
        catch e
            j=j+1;
            break
        end
        
    end
    Simulate_M3al_Model(sample_person,parameters,t_span,[0,0,0]);
    subplot(2,5,6)
    colorbar('Ticks',[0,0.5,1],'TickLabels',[p_lb,parameters(p),p_ub]);
    title(par_name(p));
    %save plot
    name=string(par_name(p));
    saveas(h,sprintf('sensitivty_50_%s.fig',name));
end