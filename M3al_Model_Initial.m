function [x0,c] = M3al_Model_Initial(input_data,parameters)
%Specify intitial values and model constants for Mixed Meal Model (M3al
%Model)from input data and parameters
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%input_data - structure of measured meall challenge test data and
%             individidual information (BW, meal composition, ect).
%parameters - model parameters.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for further information contact Shauna O'Donovan at
% s.d.odonova@tue.nl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% define initial values for state variables
M_G_0     = 0; %intial mass of glucose in digestive tract (assume 0 as fasting)
G_PL_0    = input_data.glucose(1); %fasting glucose concentration
G_int_0   = 0; %integrated plasma glucose (assume 0 as fasting)
I_PL_0    = input_data.insulin(1);%fasting insulin concentration
I_d1_0    = 0; %insulin concentrtion in remote compartment
I_d2_0    = I_PL_0; %insulin delay 3 (LPL compartment 1);
I_d3_0    = I_PL_0; %insulin delay 3 (LPL compartment 2);
I_d4_0    = I_PL_0; %insulin delay 4 (LPL compartment 3);
NEFA_PL_0 = input_data.NEFA(1); %fasting plasma NEFA concentration 
M_TG_g1_0 = 0; %initial maas of triglyceride in the digestive tract (assume 0 as fasting);
M_TG_g2_0 = 0; %initial maas of triglyceride in the digestive tract/lymphatic system (assume 0 as fasting);
M_TG_g3_0 = 0; %initial maas of triglyceride in the lymphatic system (assume 0 as fasting);
TG_PL = input_data.TG(1); %intial concentration of plasma triglyceride

x0=[M_G_0,G_PL_0,G_int_0,I_PL_0,I_d1_0,I_d2_0,I_d3_0,I_d4_0,NEFA_PL_0,M_TG_g1_0,M_TG_g2_0,M_TG_g3_0,TG_PL];

%% define model consntants

% conversion factor glucose - convert glucose from mg/l to mmol/l
c.f_G = 0.005551; %if mg/l

%conversion factor triglyceride - convert from mg/l to mmol/l
c.f_TG = 0.00113; 

%convert insulin from uIU/ml to mmol/l
c.f_I = 1; 

c.V_G     = 17/70; %volume of distribution for glucose 
c.V_TG    = 0.06; %volume of distribution of triglycerides (volume of blood)
c.G_liv_b = parameters(15);%basal hepatic glucose concentration 
c.tau_i   = 31; %(min)
c.tau_d   = 3; %(min)
c.G_th_PL = 9; %threshold for renal extraction
c.t_integralwindow = 30; %30
c.c1      = 0.1;
%fixed to values specified in parameter vector (expected fasting values)
c.c2     = c.G_liv_b.*(parameters(12) + parameters(13))./parameters(13) - parameters(5).*c.f_I.*parameters(15);
c.c3     = parameters(7).*parameters(13)./(c.f_I*c.tau_i.*parameters(14)).*c.t_integralwindow;

