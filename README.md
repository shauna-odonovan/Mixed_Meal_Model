# Mixed_Meal_Model
Code accompanying "The Mixed Meal Model: quantifying the contribution of 
triglyceride to metabolic resilience."  
The Mixed Meal Model MATLAB implementation including scripts to run local 
parameter sensitivity analysis and parameter identifiability analysis. 
Matlab version: R2019b
_________________________________________________________________________________________________
%% CONTENTS:
_________________________________________________________________________________________________
Run_Parameter_Estimation.m       - Script fitting the Mixed Meal Model to 
                                   measured challenge test data. 
Run_Model_Simulation.m           - Script simulating the Mixed Meal Model for user 
                                   specified parameters and phenotypic traits.
Run_Sensitivity_Analysis.m       - Script to run the local parameter sensitivity
                                   presented in supplementary file S2. 
Run_Identifiability_Analysis.m   - Script to perform parameter identifiability
                                   analysis using the profile liklihood approach. 
                                   This script will reproduce supplementary figure S3. 

M3al_Model_ODE.m                - Mixed Meal Model equations.  
M3al_Model_Initial.m            - Specify initial values and model constants 
                                  necessary for the Mixed Meal Model to simulated. 
M3al_Model_Parameters.m         - Defines parameter values for Mixed Meal Model. 
integratorfunG.m                - Calculate the integral of plasma glucose for 
                                  PID controller. 
M3al_Model_ErrorFunc.m          - Calculate the error used in the model 
                                  fitting/optimisation procedure. 
Fit_M3al_Model.m                - Function fitting the Mixed Meal model to 
                                  supplied challenge test data.(single fitting
                                  using lsqnonlin)
Fit_M3al_Model_LatinHyperCube.m - Function fitting the Mixed Meal model to 
                                  supplied challenge test data using multiple 
                                  initial values defined using a latin-hypercube 
                                  design. 
Plot_M3al_Model.m               - Generate figure showing Mixed Meal Model 
                                  simulation for a specified parameter set. 
Plot_MultiFit_M3al_Model.m      - Generate figure showing multiple Mixed Meal Model 
                                  simulations for each fitting generated  using 
                                  multiple initial values from the 
                                  Fit_M3al_Model_LatinHyperCube.m function.
Simulate_M3al_Model.m           - Generate figure showing Mixed Meal Model
                                  simulation without measured data.  
PLA.m                           - Function to perform profile likelihood analysis
                                  as outlined by Raue et al. (2009) Bioinformatics.

sample_data.mat                 - sample structure of measured challenge test data 
__________________________________________________________________________________________________________
%% USE:
__________________________________________________________________________________________________________
To fit the Mixed Model to the provided sample challenge test data run the Run_Parameter_Estimate.m script.
This script will load the provided sample data and estimate values for 9 parameters from the Mixed Meal 
Model using 5 initial values defined using latin-hypercube sampling. The resulting model fits will be 
visualised. This script was used to generate the results presented in Figures 2,3,4,5, and Tables 1 and 2. 

To simulate the Mixed Meal Model for specific parameters and phenotypic traits open the Run_Model_Simulation.m 
script and adjust the parameters/traits as desired. Save and run the script. The resulting model simulation will
be visualised. 

To perform local parameter sensitivity analysis for the Mixed Meal model open the Run_Sensitivty_Analysis.m script
and specify the range over which the parameters will be vairied, save the changes, and run. An individual plot will
generated and saved for each parameters. 

To perform profile likelihood analysis to determine parameter identifiability open the Run_Indetifiability_Analysis.m 
script, sepcify the initial/optimial parameter values which are being evaluated, save and run the script. Note as this 
script runs the parameter estimation proceedure many times it can take several minutes to run. 

%% For further information or questions please contact Shauna O'Donovan at s.d.odonovan@tue.nl. 
