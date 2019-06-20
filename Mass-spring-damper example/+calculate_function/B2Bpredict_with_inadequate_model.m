function z_pred_inadequateModel = B2Bpredict_with_inadequate_model(flag_usedata,spring_data,t_data,z_measure,envir_data)

if flag_usedata
   z_pred_inadequateModel = spring_data.z_pred_inadequateModel;
else
   eps = envir_data.eps;
   n_eps = envir_data.n_eps;
   n = envir_data.n_start;
   t_pred = envir_data.t_pred;
   n_pred = envir_data.n_pred;
   c2 = envir_data.c2;
   v0 = envir_data.v0;
   n_poly = envir_data.n_polytest;
   z_pred_inadequateModel = cell(n_poly+1,1);
   for i = 0:n_poly
      Hc = repmat([-inf inf],i,1);
      z_pred_inadequateModel{i+1} = zeros(n_eps,2,n_pred);
      for j = 1:n_eps
         for k = 1:n_pred
            z_pred_inadequateModel(j,:,k) = supplementary_function.makePrediction(t_data,z_measure{j},t_pred(k),eps(j),Hc,n,c2,v0);
         end
      end
   end
end