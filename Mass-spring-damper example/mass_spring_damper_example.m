% This script provides the workflow to generate all the table contents
% and figures in the mass-spring-damper example

%  Wenyu Li, 06/21/2019

%% load data and set parameters
clear;clc
flag_usedata = true;  % whether to use the existing result or calculate everything from scratch
supplementary_function.generateData;

%% behavior of true and inadequate model, and its difference (figure 2)
plot_function.plot_True_and_inadequate_models(envir_data);

%% prediction with true model (figure 3)
% calculate result
z_pred_trueModel = calculate_function.B2Bpredict_with_true_model(flag_usedata,spring_data,envir_data,t_data,z_measure);
% plot result
plot_function.plot_B2Bprediction_with_true_model(z_pred_trueModel,envir_data.t_pred,envir_data.z_pred_true);

%% consistency and k* \in F with inadequate model (table 1)
CMresult = calculate_function.Dataset_consistency_and_kstart_in_feasible_set(flag_usedata,spring_data,envir_data,t_data,z_measure);

%% prediction with inadequate model (figure 4)
% calculate result
z_pred_inadequateModel = calculate_function.B2Bpredict_with_inadequate_model(flag_usedata,spring_data,envir_data,t_data,z_measure,CMresult);
% plot result
plot_function.plot_B2Bprediction_with_inadequate_model(z_pred_inadequateModel,envir_data.t_pred,envir_data.z_pred_true);

%% posterior uncertainy interval of k using the true model
parameter_posterior_interval1 = calculate_function.B2Bposterior_interval_with_true_model(flag_usedata,spring_data,envir_data,t_data,z_measure);

%% posterior uncertainy interval of k and c's using the inadequate model (table 2)
% calculate result
parameter_posterior_interval2 = calculate_function.B2Bposterior_interval_with_inadequate_model(flag_usedata,spring_data,envir_data,t_data,z_measure,CMresult);

%% volume ratio
nsample = 1e6;
volume_ratio = calculate_function.Volume_ratio_computaion(flag_usedata,spring_data,envir_data,t_data,z_measure,nsample,parameter_posterior_interval2,CMresult);

%% model behavior over stiffness prior uncertainty (figure 5)
plot_function.plot_model_behavior_over_prior(envir_data,t_data,z_measure);

%% confounding (figure 6)
% calculate result
n_grid = 51;
confounding_mesh = calculate_function.B2Bconfounding(flag_usedata,spring_data,envir_data,t_data,z_measure,n_grid);
% plot result
plot_function.plot_Confounding(confounding_mesh);

%% posterior uncertainty of model discrepancy function (figure 7)
% calculate posterior delta
n_grid = 1000;
tRange = [0 4];
delta_posterior = calculate_function.B2Bposterior_delta(flag_usedata,spring_data,envir_data,t_data,z_measure,n_grid,CMresult,parameter_posterior_interval2,tRange);
% plot result
plot_function.plot_B2Bposterior_delta(delta_posterior,envir_data);

%% prior constraint on model discrepancy
% calculate result
ntest = 10;
pred_with_prior_constraint = calculate_function.B2Bpred_with_prior_constraint(flag_usedata,spring_data,envir_data,t_data,z_measure,parameter_posterior_interval2{1,5},ntest);
% plot result
plot_function.plot_prior_constrain_on_delta(pred_with_prior_constraint,z_pred_inadequateModel{5},envir_data,ntest);