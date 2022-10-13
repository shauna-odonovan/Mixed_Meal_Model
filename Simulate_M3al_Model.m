function out = Simulate_M3al_Model(phenotypic_data,parameters,time,plot_col)
% Simulate M3al Model for a given parameter set
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% phenotypic_data - struct of information specifying fasting glucose,
%                   insulin, TG, and NEFA. Must also sepcify body weight (BW)
%                   and meal composition. 
% parameters      - parameter values to be simulated.
% time            - timespan for ode simulation
% plot_col        - Colour to plot line 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for further information contact Shauna O'Donovan at
%s.d.odonovan@tue.nl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% simulate model for given parameter set

%define intial values and model constants needed for simulation of M3al Model model
[initial_values,constants] = M3al_Model_Initial(phenotypic_data,parameters);

%define global parameters for simulation
global t_saved G_PL_saved;
%initialise gloabl parameters
t_saved = 0;
G_PL_saved = phenotypic_data.glucose(1);

%specify options for ODE solver (Integrator function)
ODE_options = odeset('RelTol',1e-5,'OutputFcn',@integratorfunG);

%simulate model
[T,X] = ode45(@M3al_Model_ODE,time,initial_values,ODE_options,parameters,constants,phenotypic_data);

%% Generate figure

 subplot(2,5,1)
 %glucose from gut
 k2 = parameters(2);
 BW = phenotypic_data.BW;
 f_G = constants.f_G;
 V_G = constants.V_G; 
 
 G_gut = k2.*(f_G/(V_G*BW)).*X(:,1);
 
 plot(T,G_gut,'Color',plot_col,'LineWidth',1.5);
 hold on
 xlabel('time (mins)')
 ylabel('glucose (mmol/l)');
 title('glucose from gut')


 subplot(2,5,2)
 %net hepatic glucose flux
 k3 = parameters(3);
 k4 = parameters(4);
 f_I = 1;
 G_liv_b = parameters(15);
 G_b = parameters(13);
 
 G_liv = G_liv_b - k4.*f_I.*X(:,5) - k3.*(X(:,2)-G_b);
 
 plot(T,G_liv,'Color',plot_col,'LineWidth',1.5);
 hold on
 xlabel('time (mins)')
 ylabel('glucose (mmol/l)');
 title('net hepatic glucose flux')
 
 
 subplot(2,5,3)
 %glucose uptake into tissues
 KM = parameters(12);
 k5=parameters(5);

 %insulin independent glucose utilisation (brain, erythrocyte,adipose,muscle)
 U_ii = G_liv_b*((KM + G_b)./G_b).*(X(:,2)./(KM + X(:,2)));
 %insulin dependent glucose utilisation (muscle,adipose)
 U_id = k5.*f_I.*X(:,5).*(X(:,2)./(KM + X(:,2)));
 %total glucose uptake
 G_U = U_ii + U_id;
 
 plot(T,G_U,'Color',plot_col,'LineWidth',1.5);
 hold on
 xlabel('time (mins)')
 ylabel('glucose (mmol/l)');
 title('glucose uptake')
 
 
 subplot(2,5,4)
 %plot plasma glucose

 plot(T,X(:,2),'Color',plot_col,'LineWidth',1.5);
 hold on
 ylabel('glucose (mmol/l)');
 title('plasma Glucose')
 xlabel('time (mins)')
 
 
 subplot(2,5,5)
 %plasma insulin
 plot(T,X(:,4),'Color',plot_col,'LineWidth',1.5);
 hold on
 ylabel('insulin (uIU/ml)');
 title('plasma insulin')
 
 subplot(2,5,7)
 % triglyceride from the gut
 k14  = parameters(23);
 f_TG = constants.f_TG;
 V_TG = constants.V_TG; 
 
 TG_gut = k14.*(f_TG/(V_TG*BW)).*X(:,12);
 
 plot(T,TG_gut,'Color',plot_col,'LineWidth',1.5);
 hold on
 xlabel('time (mins)')
 ylabel('TG (mmol/l)');
 title('TG from gut')
 
 subplot(2,5,8)
 %TG secretion from liver
 I_PL_b = parameters(14);
 k15    = parameters(24);
 k16    = parameters(25);

 VLDL = k16 - k15.*(X(:,8)-I_PL_b);
 
 plot(T,VLDL,'Color',plot_col,'LineWidth',1.5);
 hold on
 xlabel('time (mins)')
 ylabel('TG (mmol/l)');
 title('TG secretion liver (VLDL)') 
 
 subplot(2,5,9)
 %plasma triglyceride
 plot(T,X(:,13),'Color',plot_col,'LineWidth',1.5);
 hold on
 ylabel('TG (mmol/l)');
 title('plasma TG')
 xlabel('time (mins)')
 
 subplot(2,5,10)
 %plasma NEFA
 plot(T,X(:,9),'Color',plot_col,'LineWidth',1.5)
 hold on
 ylabel('NEFA (mmol/l)');
 title('plasma NEFA')
 xlabel('time (mins)')
 
 out=1;