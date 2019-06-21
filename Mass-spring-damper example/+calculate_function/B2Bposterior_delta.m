function y = B2Bposterior_delta(flag_usedata,spring_data,envir_data,t_data,z_measure,ndata,CMresult,c_posterior,tRange)
% Calculate the posterior uncertainty interval of model discrepancy
% function, for the cases of a quadratic and cubic delta.

% Output is a structure variable with calculated minimum and maximum of the
% discrepancy function over the specified scenario domain

if flag_usedata
   y = spring_data.posterior_delta;
else
   tCal = linspace(tRange(1),tRange(2),ndata);
   v0 = envir_data.v0;
   c2 = envir_data.c2;
   nD = length(t_data);
   eps = envir_data.eps;
   n_eps = length(eps);
   n_poly = envir_data.n_polytest;
   n = envir_data.n_start;
   Hk = envir_data.H_k;
   opt = optimoptions('fmincon','Display','off','MaxFunctionEvaluations',1e4,...
      'ConstraintTolerance',1e-6,'OptimalityTolerance',1e-14,'MaxIterations',1e3,...
      'GradConstr','on','GradObj','on','DerivativeCheck','off','FiniteDifferenceType','central');
   y.min = cell(n_eps,n_poly+1);
   y.max = cell(n_eps,n_poly+1);
   y.tRange = tRange;
   for i = 1:n_eps
      ydata = z_measure{i};
      for j = 1:n_poly
         Hc = repmat([-inf inf],j,1);
         H = [Hk,Hc];
         pmin = zeros(1,ndata);
         pmax = zeros(1,ndata);
         if strcmp(CMresult{i,j},'Dataset is inconsistent')
            continue
         else
            c_pos = zeros(j,2);
            for k = 1:j
               c_pos(k,:) = c_posterior{i,j}(k+1).Posterior_interval;
            end
            for k = 1:ndata
               ymin = zeros(n,1);
               ymax = zeros(n,1);
               flag_min = zeros(n,1);
               flag_max = zeros(n,1);
               xStart = [Hk(1)+rand(2*n,1)*diff(Hk) mean(c_pos,2)'+0.1*diff(c_pos').*randn(2*n,j)];
               for kk = 1:n
                  [~,ymin(kk),flag_min(kk)] = fmincon(@fun_min,xStart(2*k-1,:)',[],[],[],[],...
                     H(:,1),H(:,2),@neq,opt);
                  [~,ymax(kk),flag_max(kk)] = fmincon(@fun_max,xStart(2*k,:)',[],[],[],[],...
                     H(:,1),H(:,2),@neq,opt);
               end
               ymin(flag_min<=0) = [];
               ymax(flag_max<=0) = [];
               if isempty(ymin)
                  pmin(k) = inf;
               else
                  pmin(k) = min(ymin);
               end
               if isempty(ymax)
                  pmax(k) = -inf;
               else
                  pmax(k) = -min(ymax);
               end
            end
         end
         y.min{i,j} = pmin;
         y.max{i,j} = pmax;
      end
   end
end

   function [y,gy] = fun_min(x)
      Tfit = supplementary_function.makePower(tCal(k),j);
      y = Tfit*x(2:end);
      gy = [0; Tfit'];
   end

   function [y,gy] = fun_max(x)
      Tfit = supplementary_function.makePower(tCal(k),j);
      y = -Tfit*x(2:end);
      gy = [0; -Tfit'];
   end

   function [c,ceq,g,geq] = neq(x)
      rr = sqrt(x(1));
      drr = 0.5/rr;
      c1 = v0/rr;
      dc1 = -c1*drr/rr;
      c = zeros(2*nD,1);
      g = zeros(nV,2*nD);
      ceq = [];
      geq = [];
      Tfit = supplementary_function.makePower(t_data,j);
      yy = c1*sin(rr*t_data)+c2*cos(rr*t_data)+Tfit*x(2:end);
      c(1:nD) = yy-ydata-eps(i);
      c(nD+1:end) = ydata-eps(i)-yy;
      g(1,1:nD) = drr*t_data.*(c1*cos(rr*t_data)-c2*sin(rr*t_data))+dc1*sin(rr*t_data);
      g(2:end,1:nD) = Tfit';
      g(:,nD+1:end) = -g(:,1:nD);
   end

end