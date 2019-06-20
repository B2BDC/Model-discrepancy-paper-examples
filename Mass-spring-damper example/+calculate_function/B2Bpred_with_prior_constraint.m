function y = B2Bpred_with_prior_constraint(flag_usedata,spring_data,envir_data,t_data,z_measure,post_interval,ntest)

if flag_usedata
   y = spring_data.prediction_with_prior_on_delta;
else
   c2 = envir_data.c2;
   v0 = envir_data.v0;
   ydata = z_measure{1};
   nD = length(t_data);
   eps = envir_data.eps(1);
   tpred = envir_data.t_pred;
   nPred = length(tpred);
   Hk = envir_data.H_k;
   alpha_opt = inf(nPred,1);
   Hc = zeros(4,1);
   for i = 1:4
      Hc(i,:) = post_interval(i+1).Posterior_interval;
   end
   H = [Hk;Hc];
   dH = diff(H,[],2);
   nV = size(H,1);
   n_poly = nV-1;
   xOpt = cell(nPred,1);
   opt = optimoptions('fmincon','Display','off','MaxFunctionEvaluations',1e4,...
      'ConstraintTolerance',1e-6,'OptimalityTolerance',1e-14,'MaxIterations',1e3,...
      'GradConstr','on','GradObj','on','DerivativeCheck','off','FiniteDifferenceType','central');
   for i = 1:nPred
      % minimum |B/M|
      tcond = [t_data;tpred(i)];
      xStart = [Hk(1)+diff(Hk)*rand(2*n,1) mean(Hc')+diff(Hc')*randn(2*n,n_poly)];
      xopt = cell(n,1);
      yopt = zeros(n,1);
      flag = zeros(n,1);
      for j = 1:n
         [xopt{j},yopt(j),flag(j)] = fmincon(@minB,xStart(j,:)',[],[],[],[],...
            H(:,1),H(:,2),@neq,opt);
      end
      xopt(flag<=0) = [];
      yopt(flag<=0) = [];
      if ~isempty(yopt)
         [alpha_opt(i),id] = min(yopt);
         xOpt{i} = xopt{id};
      else
         alpha_opt(i) = inf;
      end
   end
   if any(isinf(yopt))
      return
   else
      alpha_large = 10*alpha_opt;
   end
   yPred = zeros(ntest,2,nPred);
   for i = 1:nPred
      tmpy = zeros(ntest,2);
      tcond = [t_data;tpred(i)];
      cc = linspace(alpha_opt(i),alpha_large(i),ntest+1);
      cc = cc(2:end);
      for j = 1:ntest
         xStart = repmat(xOpt{i}',2*n,1)+2e-3*(rand(2*n,nV)-0.5).*repmat(dH',2*n,1);
         ymin = zeros(n,1);
         flagmin = zeros(n,1);
         ymax = zeros(n,1);
         flagmax = zeros(n,1);
         for k = 1:n
            [~,ymin(j),flagmin(j)] = fmincon(@minfun,xStart(2*k-1,:)',[],[],[],[],...
               H(:,1),H(:,2),@neq2,opt);
            [~,ymax(j),flagmax(j)] = fmincon(@maxfun,xStart(2*k,:)',[],[],[],[],...
               H(:,1),H(:,2),@neq2,opt);
         end
         ymin(flagmin<=0) = inf;
         ymax(flagmax<=0) = inf;
         tmpy(j,:) = [min(ymin) -min(ymax)];
      end
      yPred(:,:,i) = tmpy;
   end
   y.prediction = yPred;
   y.minimal_alpha = alpha_opt;
   y.maximal_alpha = alpha_large;
   
end
   
   function [y,gy] = minB(x)
      gy = zeros(length(x),1);
      Tfit = supplementary_function.makePower(tcond,n_poly);
      delta = Tfit*x(2:end);
      ss = sign(delta);
      y = mean(abs(delta));
      gy(2:end) = mean(Tfit'.*repmat(ss',n_poly,1),2);
   end
   
   function [y,gy] = minfun(x)
      rr = sqrt(x(1));
      drr = 0.5/rr;
      c1 = v0/rr;
      dc1 = -c1*drr/rr;
      Tfit = supplementary_function.makePower(tpred(i),n_poly);
      gy = zeros(nV,1);
      y = c1*sin(rr*tpred(i))+c2*cos(rr*tpred(i))+Tfit*x(2:end);
      gy(1) = drr*tpred(i)*(c1*cos(rr*tpred(i))-c2*sin(rr*tpred(i)))+dc1*sin(rr*tpred(i));
      gy(2:end) = Tfit;
   end
   
   function [y,gy] = maxfun(x)
      rr = sqrt(x(1));
      drr = 0.5/rr;
      c1 = v0/rr;
      dc1 = -c1*drr/rr;
      Tfit = supplementary_function.makePower(tpred(i),n_poly);
      gy = zeros(nV,1);
      y = -c1*sin(rr*tpred(i))-c2*cos(rr*tpred(i))-Tfit*x(2:end);
      gy(1) = -drr*tpred(i)*(c1*cos(rr*tpred(i))-c2*sin(rr*tpred(i)))-dc1*sin(rr*tpred(i));
      gy(2:end) = -Tfit;
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
      Tfit = supplementary_function.makePower(t_data,n_poly);
      yy = c1*sin(rr*t_data)+c2*cos(rr*t_data)+Tfit*x(2:end);
      c(1:nD) = yy-ydata-eps;
      c(nD+1:end) = ydata-eps-yy;
      g(1,1:nD) = drr*t_data.*(c1*cos(rr*t_data)-c2*sin(rr*t_data))+dc1*sin(rr*t_data);
      g(2:end,1:nD) = Tfit';
      g(:,nD+1:end) = -g(:,1:nD);
   end

   function [c,ceq,g,geq] = neq2(x)
      rr = sqrt(x(1));
      drr = 0.5/rr;
      c1 = v0/rr;
      dc1 = -c1*drr/rr;
      c = zeros(2*nD+1,1);
      g = zeros(nV,2*nD+1);
      ceq = [];
      geq = [];
      Tfit = supplementary_function.makePower(t_data,n_poly);
      yy = c1*sin(rr*t_data)+c2*cos(rr*t_data)+Tfit*x(2:end);
      c(1:nD) = yy-ydata-eps;
      c(nD+1:end-1) = ydata-eps-yy;
      g(1,1:nD) = drr*t_data.*(c1*cos(rr*t_data)-c2*sin(rr*t_data))+dc1*sin(rr*t_data);
      g(2:end,1:nD) = Tfit';
      g(:,nD+1:end-1) = -g(:,1:nD); 
      Tfit = subfun.makePower(tcond,n_poly);
      delta = Tfit*x(2:end);
      ss = sign(delta);
      c(end) = mean(abs(delta))-cc(j);
      g(2:end,end) = mean(Tfit'.*repmat(ss',n_poly,1),2);  
   end
end