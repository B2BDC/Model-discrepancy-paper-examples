function y = B2Bposterior_interval_with_true_model(flag_usedata,spring_data,envir_data,tdata,z_measure)
% Calculate the posterior interval of stiffness k. The B2BDC computation is
% conducted with the true model.

% Output is a cell array with length n_eps, corresponds to the tested
% epsilon values. Each cell element is a structure array, with its fields
% representing the name and posterior interval of calculated parameters.

if flag_usedata
   y = spring_data.Posterior_with_true_model;
else
   n = envir_data.n_start;
   c2 = envir_data.c2;
   v0 = envir_data.v0;
   b = envir_data.b;
   eps = envir_data.eps;
   n_eps = length(eps);
   H = envir_data.H_k;
   y = cell(2,1);
   nD = length(tdata);
   opt = optimoptions('fmincon','Display','off','MaxFunctionEvaluations',1e4,...
      'ConstraintTolerance',1e-6,'OptimalityTolerance',1e-14,'MaxIterations',1e3,...
      'GradConstr','on','GradObj','on','DerivativeCheck','off','FiniteDifferenceType','central');
   for i = 1:n_eps 
      ydata = z_measure{i};
      x_start = 0.1*rand(2*n,1)+0.2;
      y_min = zeros(n,1);
      y_max = zeros(n,1);
      exitflag_min = zeros(n,1);
      exitflag_max = zeros(n,1);
      aa = struct('Parameter_name','','Posterior_interval',[]);
      for j = 1:n
         [~,y_min(j),exitflag_min(j)] = fmincon(@fun_min,x_start(2*j-1,:)',[],[],[],[],...
            H(:,1),H(:,2),@neq,opt);
         [~,y_max(j),exitflag_max(j)] = fmincon(@fun_max,x_start(2*j,:)',[],[],[],[],...
            H(:,1),H(:,2),@neq,opt);
      end
      y_min(exitflag_min<=0) = [];
      y_max(exitflag_max<=0) = [];
      if isempty(y_min)
         kRange(1) = inf;
      else
         kRange(1) = min(y_min);
      end
      if isempty(y_max)
         kRange(2) = -inf;
      else
         kRange(2) = -min(y_max);
      end
      aa.Posterior_interval = kRange;
      y{i} = aa;
   end
end


   function [y,gy] = fun_min(x)
      y = x;
      gy = 1;
   end

   function [y,gy] = fun_max(x)
      y = -x;
      gy = -1;
   end

   function [c,ceq,g,geq] = neq(x)
      rr = sqrt(x);
      drr = 0.5/rr;
      c1 = (v0+0.5*b*c2)/rr;
      dc1 = -c1*drr/rr;
      c = zeros(2*nD,1);
      g = zeros(1,2*nD);
      ceq = [];
      geq = [];
      yy = exp(-0.5*b*tdata).*(c1*sin(rr*tdata)+c2*cos(rr*tdata));
      c(1:nD) = yy-ydata-eps(i);
      c(nD+1:end) = ydata-eps(i)-yy;
      g(1,1:nD) = exp(-0.5*b*tdata).*(drr*tdata.*(c1*cos(rr*tdata)-c2*sin(rr*tdata))+dc1*sin(rr*tdata));
      g(:,nD+1:end) = -g(:,1:nD);
   end

end