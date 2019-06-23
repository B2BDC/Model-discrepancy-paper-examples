% This script provides the workflow to reproduce B2BDC calculation on the
% hydrogen example

%  Wenyu Li, 06/21/2019

%% load data
clear;clc
load('H2_data.mat');

%% generate dataset with specified epsilon and model discrepancy
eps = 0.01;    % 0.01 or 0.005
n_poly = 1;  % discrepancy polynomial order (-1 for no model discrepancy case)
ds = supplementary_function.generateDataset(H2_data,eps,n_poly);

%% set B2BDC option
opt = generateOpt('Display',false,'AddFitError',false,'Prediction','inner');

%% dataset consistency
opt.OptimOption.RandomStart = 5;
ds.isConsistent(opt);

%% if dataset is feasible, distance of the nominal value from feasible set
d_lambda = calculate_function.Distance_from_feasible_set(ds,opt);

%% calculate prediction
t_pred = calculate_function.B2Bprediction(H2_data,ds,opt);

%% calculate posterior interval of parameters
posterior_interval = calculate_function.B2Bpred_posterior_interval_of_parameters(ds,opt);

%% posterior uncertainty of model discrepancy at specified scenario slice
s_fix.scenario = 'T';  % T, P or Phi
s_fix.value = 0;    % slice is defined by fix T at s1=0
ngrid = 51;
posterior_delta = calculate_function.B2Bpred_posterior_delta_in_2D_scenario_slice(ds,opt,s_fix,ngrid);
