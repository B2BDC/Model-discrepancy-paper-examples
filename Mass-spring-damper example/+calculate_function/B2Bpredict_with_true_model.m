function  y = B2Bpredict_with_true_model(flag_usedata,spring_data,envir_data,t_data,z_measure)
% Calculate B2BDC predicted interval using the true model

% Output is a cell array of length n+1, corresponding to cases with
% different delta functions. The 1st element corresponds to the case where
% no delta function is used. Each cell element is a n_eps-by-2-by-n_pred 
% matrix with its 1st, and 3rd dimensions
% correspond to the tested epsilon and prediction time scenarios,
% respectively. The 2nd dimension corresponds to the predicted interval.

if flag_usedata
   y = spring_data.z_pred_trueModel;
else
   eps = envir_data.eps;
   n_eps = envir_data.n_eps;
   n = envir_data.n_start;
   t_pred = envir_data.t_pred;
   n_pred = envir_data.n_pred;
   y = zeros(n_eps,2,n_pred);
   b = envir_data.b;
   c2 = envir_data.c2;
   v0 = envir_data.v0;
   for i = 1:n_eps
      for j = 1:n_pred
         y(i,:,j) = supplementary_function.makePrediction_true(t_data,z_measure{i},t_pred(j),eps(i),b,c2,v0,n);
      end
   end
end