function  z_pred_trueModel = B2Bpredict_with_true_model(flag_usedata,spring_data,t_data,z_measure,envir_data)

if flag_usedata
   z_pred_trueModel = spring_data.z_pred_trueModel;
else
   eps = envir_data.eps;
   n_eps = envir_data.n_eps;
   n = envir_data.n_start;
   t_pred = envir_data.t_pred;
   n_pred = envir_data.n_pred;
   z_pred_trueModel = zeros(n_eps,2,n_pred);
   b = envir_data.b;
   c2 = envir_data.c2;
   v0 = envir_data.v0;
   for i = 1:n_eps
      for j = 1:n_pred
         z_pred_trueModel(i,:,j) = subfun.makePrediction_true(t_data,z_measure{i},t_pred(j),eps(i),b,c2,v0,n);
      end
   end
end