function y = B2Bposterior_interval_with_inadequate_model(flag_usedata,spring_data,envir_data,tdata,z_measure,CMresult)
% Calculate the posterior interval of stiffness k. The B2BDC computation is
% conducted with the inadequate model.

% Output is a n_eps-by-(n+1) cell array, corresponding to the tested
% epsilon and used discrepancy function, respectively. The first column
% corresponds to the case where no delta is used. Each cell element is a 
% structure array, with its fields representing the name and posterior 
% interval of calculated parameters.

if flag_usedata
   y = spring_data.Posterior_with_inadequate_model;
else
   n = envir_data.n_start;
   c2 = envir_data.c2;
   v0 = envir_data.v0;
   eps = envir_data.eps;
   n_eps = length(eps);
   Hk = envir_data.H_k;
   n_poly = envir_data.n_polytest;
   y = cell(2,n_poly+1);
   nD = length(z_measure{1});
   opt = optimoptions('fmincon','Display','off','MaxFunctionEvaluations',1e4,...
      'ConstraintTolerance',1e-6,'OptimalityTolerance',1e-14,'MaxIterations',1e3,...
      'GradConstr','on','GradObj','off','DerivativeCheck','on','FiniteDifferenceType','central');
   for i = 1:n_eps
      ydata = z_measure{i};
      for j = 0:n_poly
         Hc = repmat([-inf inf],j,1);
         H = [Hk;Hc];
         aa = repmat(struct('Parameter_name','','Posterior_interval',[]),1,j+1);
         Tfit = supplementary_function.makePower(tdata,j);
         if strcmp(CMresult{i,j+1},'Dataset is inconsistent')
            for k = 1:j+1
               if k == 1
                  aa(k).Parameter_name = 'k';
               else
                  aa(k).Parameter_name = ['c' num2str(k-2)];
               end
               aa(k).Posterior_interval = [inf -inf];
            end
         else
            for k = 1:j+1
               x_start = [diff(Hk)*rand(2*n,1)+Hk(1) randn(2*n,j)];
               y_min = zeros(n,1);
               y_max = zeros(n,1);
               exitflag_min = zeros(n,1);
               exitflag_max = zeros(n,1);
               for kk = 1:n
                  [~,y_min(kk),exitflag_min(kk)] = fmincon(@fun_min,x_start(2*kk-1,:)',[],[],[],[],...
                     H(:,1),H(:,2),@neq,opt);
                  [~,y_max(kk),exitflag_max(kk)] = fmincon(@fun_max,x_start(2*kk,:)',[],[],[],[],...
                     H(:,1),H(:,2),@neq,opt);
               end
               y_min(exitflag_min<=0) = [];
               y_max(exitflag_max<=0) = [];
               if k == 1
                  aa(k).Parameter_name = 'k';
               else
                  aa(k).Parameter_name = ['c' num2str(k-2)];
               end
               if isempty(y_min)
                  aa(k).Posterior_interval(1) = inf;
               else
                  aa(k).Posterior_interval(1) = min(y_min);
               end
               if isempty(y_max)
                  aa(k).Posterior_interval(2) = -inf;
               else
                  aa(k).Posterior_interval(2) = -min(y_max);
               end
            end
         end
         y{i,j+1} = aa;
      end
   end 
end

  

   function [y,gy] = fun_min(x)
      y = x(k);
      gy = zeros(j+1,1);
      gy(k) = 1;
   end

   function [y,gy] = fun_max(x)
      y = -x(k);
      gy = zeros(j+1,1);
      gy(k) = -1;
   end

   function [c,ceq,g,geq] = neq(x)
      rr = sqrt(x(1));
      drr = 0.5/rr;
      c1 = v0/rr;
      dc1 = -c1*drr/rr;
      c = zeros(2*nD,1);
      g = zeros(j+1,2*nD);
      ceq = [];
      geq = [];
      if j == 0
         c1*sin(rr*tdata)+c2*cos(rr*tdata);
      else
         yy = c1*sin(rr*tdata)+c2*cos(rr*tdata)+Tfit*x(2:end);
      end
      c(1:nD) = yy-ydata-eps(i);
      c(nD+1:end) = ydata-eps(i)-yy;
      g(1,1:nD) = drr*tdata.*(c1*cos(rr*tdata)-c2*sin(rr*tdata))+dc1*sin(rr*tdata);
      g(2:end,1:nD) = Tfit';
      g(:,nD+1:end) = -g(:,1:nD);
   end
end
