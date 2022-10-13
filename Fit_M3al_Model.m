function out = Fit_M3al_Model(initial_par,input_data,time,lb,ub)
%Fit Mixed Meal glucose-insulin-NEFA-TG model to  measured meal 
%challenge data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initial_par - initial guess for parameters to be optimised
%input_data  - structure of measured postprandial data
%              - must contain postprandial measurement of glucose, insulin, TG, NEFA
%              - must also contain body weight and meal composition
%time        - vector of time for simulation
%lb          - lower bounds for parameter fitting (same size as initial_par)
%ub          - upper bounds for parameter fitting (same size as initial_par)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for further information contact Shauna O'Donovan at
% s.d.odonovan@tue.nl
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define options for lsqnonlin algorithm
lsq_options=optimset('Algorithm','trust-region-reflective','MaxFunEvals',1000,'TolX',1e-8,'Display','iter','UseParallel',1);
%Fit model to data
[p_opt,resnorm,residual,exitflag,output,lambda,jacobian]=lsqnonlin(@M3al_Model_ErrorFunc,initial_par,lb,ub,lsq_options,input_data,time);
 
out.p_opt    = p_opt;
out.resnorm  = resnorm;
out.residual = residual;
out.exitflag = exitflag;
out.output   = output;
out.lambda   = lambda;
out.jacobian = jacobian;