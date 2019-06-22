function y = Volume_ratio_computaion(flag_usedata,spring_data,envir_data,t_data,z_measure,nsample,Posterior_bound,CMresult)
% Calculate the volume ratio of the feasible set to the multidimensional
% posterior box.

% Output is a n_eps-by-(n+1) matrix. The 1st and 2nd dimensions correspond
% to the tested epsilon and used polynomial discrepancy functions,
% respectively

if flag_usedata
   y = spring_data.volume_ratio;
else
   eps = envir_data.eps;
   n_eps = length(eps);
   n_poly = envir_data.n_polytest;
   y = zeros(n_eps,n_poly+1);
   v0 = envir_data.v0;
   nD = length(t_data);
   c2 = envir_data.c2;
   for i = 1:2
      ydata = z_measure{i};
      for j = 0:4
         if strcmp(CMresult{i,j+1},'Dataset is inconsistent')
            continue
         else
            boxH = -ones(j+1,2);
            for k = 1:j+1
               boxH(k,:) = Posterior_bound{i,j+1}(k).Posterior_interval;
            end
            boxH = boxH';
            dH = diff(boxH);
            xCand = boxH(1,:)+rand(nsample,j+1).*dH;
            if j > 0
               Tfit = supplementary_function.makePower(t_data,j);
            else
               Tfit = 0;
            end
            xx = sqrt(xCand(:,1)');
            C1 = v0./xx;
            ySample = repmat(C1,nD,1).*sin(t_data*xx)+c2*cos(t_data*xx)+Tfit*xCand(:,2:end)'-...
               repmat(ydata,1,nsample);
            flags = all(abs(ySample)<=eps(i));
            y(i,j+1) = sum(flags)/nsample;
         end
      end
   end
end