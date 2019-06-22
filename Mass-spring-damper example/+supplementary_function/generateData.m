% This script sets up all relevant parameter values used in the example
envir_data.k = 0.25;
envir_data.b = 0.05;
envir_data.z0 = -1.5;
envir_data.v0 = 1;
envir_data.lambda = 0.5*sqrt(4*envir_data.k-envir_data.b^2);
envir_data.c1 = (envir_data.v0+0.5*envir_data.b*envir_data.z0)/envir_data.lambda;
envir_data.c2 = envir_data.z0;
z_true = @(t) exp(-0.5*envir_data.b*t).*(envir_data.c1*sin(envir_data.lambda*t)+envir_data.c2*cos(envir_data.lambda*t));
envir_data.H_t = [0 3];
envir_data.H_k = [0.2 0.3];
envir_data.ndata = 20;
envir_data.t_pred = [1.5; 3.2; 4];
envir_data.z_pred_true = z_true(envir_data.t_pred);
envir_data.n_pred = length(envir_data.t_pred);
envir_data.eps = [0.1 0.05];
envir_data.n_eps = length(envir_data.eps);
envir_data.n_start = 10;
envir_data.n_polytest = 4;
envir_data.seedID = 20;
if flag_usedata
   load('spring_data.mat');
   t_data = spring_data.t_data;
   z_data = spring_data.z_data;
   z_measure = spring_data.z_measure;
else
   rng(envir_data.seedID);
   t_data = 3*rand(ndata,1);
   t_data = sort(t_data);
   z_data = z_true(t_data);
   z_measure = cell(2,1);
   for i = 1:n_eps
      z_measure{i} = z_data + 2*eps(i)*(rand(ndata,1)-0.5);
   end
end
