function out = Plot_M3al_Model(input_data,p_opt,time,plot_col)
% Plot measured data and Mixed Meal Model simulation for a given parameter set
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input_data - struct of measured challenge test data must contain
%              vectors of time series of mean values and standard deviations
%              for glucose, insulin with that nomenclature.
%              Must also contain vector of time points corresponding to
%              sampling time points of measured data. 
% p_opt      - parameter values to be simulated.
% time       - timespan for ode simulation
% plot_col   - Colour to plot line 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for further information contact Shauna O'Donovan at
% s.d.odonovan@tue.nl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% simulate model for given parameter set
%generate complete set of model parameters from p_opt vector
parameters = M3al_Model_Parameters(p_opt,input_data);

%define intial values and model constants needed for simulation of eDES model
[initial_values,constants] = M3al_Model_Initial(input_data,parameters);

%define global parameters for simulation
global t_saved G_PL_saved;
%initialise gloabl parameters
t_saved = 0;
G_PL_saved = input_data.glucose(1);

%specify options for ODE solver (Integrator function)
ODE_options = odeset('RelTol',1e-5,'OutputFcn',@integratorfunG);

%simulate model
[T,X] = ode45(@M3al_Model_ODE,time,initial_values,ODE_options,parameters,constants,input_data);

%% Generate figure

 subplot(2,5,1)
 %glucose from gut
 k2 = parameters(2);
 BW = input_data.BW;
 f_G = constants.f_G;
 V_G = constants.V_G; 
 
 G_gut = k2.*(f_G/(V_G*BW)).*X(:,1);
 
 AUC_G=trapz(G_gut);
 AUC_G = ((V_G*BW)/f_G)*AUC_G*0.001;
 
 plot(T,G_gut,'Color',plot_col,'LineWidth',1.5);
 hold on
 message=['AUC = ',num2str(AUC_G),'g'];
 text(240,max(G_gut),message,'HorizontalAlignment','right','Color',plot_col)
 xlabel('time (mins)')
 ylabel('glucose (mmol/l)');
 title('glucose from gut')
 xticks(input_data.time_G)
 

 subplot(2,5,2)
 %net hepatic glucose flux
 k3 = parameters(3);
 k4 = parameters(4);
 f_I = 1;
 G_liv_b = parameters(15);
 G_b = parameters(13);%input_data.glucose(individual,1);
 
 G_liv = G_liv_b - k4.*f_I.*X(:,5) - k3.*(X(:,2)-G_b);
 
 plot(T,G_liv,'Color',plot_col,'LineWidth',1.5);
 hold on
 xlabel('time (mins)')
 ylabel('glucose (mmol/l)');
 title('net hepatic glucose flux')
 xticks(input_data.time_G)
 
 
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
 xticks(input_data.time_G)
 
 
 subplot(2,5,4)
 %plot plasma glucose

 plot(T,X(:,2),'Color',plot_col,'LineWidth',1.5);
 hold on
 plot(input_data.time_G,input_data.glucose,'kx','MarkerSize',10,'LineWidth',1.5)
 ylabel('glucose (mmol/l)');
 title('plasma Glucose')
 xlabel('time (mins)')
 xticks(input_data.time_G)
 
 
 subplot(2,5,5)
 %plasma insulin
 plot(T,X(:,4),'Color',plot_col,'LineWidth',1.5);
 hold on
 plot(input_data.time_I,input_data.insulin,'kx','MarkerSize',10,'LineWidth',1.5)
 %plot(T,zeros(size(T)),'k--');
 ylabel('insulin (uIU/ml)');
 title('plasma insulin')
 xticks(input_data.time_I)
 xlabel('time (mins)')
 
 subplot(2,5,7)
 % triglyceride from the gut
 k14 = parameters(23);
 BW = input_data.BW;
 f_TG = constants.f_TG;
 V_TG = constants.V_TG; 
 
 TG_gut = k14.*(f_TG/(V_TG*BW)).*X(:,12);
 
 AUC_TG=trapz(TG_gut);
 AUC_TG = ((V_TG*BW)/f_TG)*AUC_TG*0.001;
 
 plot(T,TG_gut,'Color',plot_col,'LineWidth',1.5);
 hold on
 message=['AUC = ',num2str(AUC_TG),'g'];
 text(240,max(TG_gut),message,'HorizontalAlignment','right','Color',plot_col)
 xlabel('time (mins)')
 ylabel('TG (mmol/l)');
 title('TG from gut')
 xticks(input_data.time_TG)
 
 subplot(2,5,8)
 %TG secretion from liver
 I_PL_b = parameters(14);
 k15 = parameters(24);
 k16 = parameters(25);

 VLDL = k16 - k15.*(X(:,8)-I_PL_b);
 
 plot(T,VLDL,'Color',plot_col,'LineWidth',1.5);
 hold on
 xlabel('time (mins)')
 ylabel('TG (mmol/l)');
 title('TG secretion liver (VLDL)')
 xticks(input_data.time_G)
 
 
 subplot(2,5,9)
 %plasma triglyceride
 plot(T,X(:,13),'Color',plot_col,'LineWidth',1.5);
 hold on
 plot(input_data.time_TG,input_data.TG,'kx','MarkerSize',10,'LineWidth',1.5)
 ylabel('TG (mmol/l)');
 title('plasma TG')
 xlabel('time (mins)')
 xticks(input_data.time_TG)
 
 subplot(2,5,10)
 %plasma NEFA
 plot(T,X(:,9),'Color',plot_col,'LineWidth',1.5)
 hold on
 plot(input_data.time_NEFA,input_data.NEFA,'kx','MarkerSize',10,'LineWidth',1.5)
 ylabel('NEFA (mmol/l)');
 title('plasma NEFA')
 xlabel('time (mins)')
 xticks(input_data.time_NEFA)
 
 out=1;