function status = integratorfunG(t,xmodelG,flag,~,~,~,~,~,~)  
%intergraterfunG is used from Rozendaal et al.(2018)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%integratorfunG - function used to implement the integrator part of the
%insulin secretion model. 
%Current plasma glucose concentation is stored in global vector Gpl_saved. 
%
%   status = integratorfunG(t,xmodelG,flag)
%
%inputs:
% t:       current simulation time
% xmodelG: vector of current model state variables 
%

% all parameters that are passed to the ode solver need to be listed here
% using ~ if they are unused. First the outputs are listed (here t and
% xmodelG), then the corresponding flag, and then all inputs (in this case
% p,p0,c,input,constant_parameters,variable_parameters. Since these are all
% unused, we replace them with ~.

%History
% 5-Mar-2016 Help info added (Natal van Riel, TU/e)

global t_saved GPL_saved

if nargin < 3 || isempty(flag)
    t_saved = [t_saved t];
    GPL_saved = [GPL_saved xmodelG(2,:)];
end

status = 0;