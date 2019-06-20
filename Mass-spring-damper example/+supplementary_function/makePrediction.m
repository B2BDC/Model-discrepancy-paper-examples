function y = makePrediction(tdata,ydata,tpred,eps,Hc,n,c2,v0)

nD = length(ydata);
opt = optimoptions('fmincon','Display','off','MaxFunctionEvaluations',1e4,...
  'ConstraintTolerance',1e-6,'OptimalityTolerance',1e-14,'MaxIterations',1e3,...
  'GradConstr','on','GradObj','on','DerivativeCheck','off','FiniteDifferenceType','central');
Hx = [0.2 0.3];
nMI = size(Hc,1);
nV = nMI+1;
H = [Hx;Hc];
x_start = zeros(2*n,nV);
x_start(:,1) = 0.1*rand(2*n,1)+0.2;
x_start(:,2:end) = randn(2*n,nMI);
y_min = zeros(n,1);
y_max = zeros(n,1);
exitflag_min = zeros(n,1);
exitflag_max = zeros(n,1);
lambda_min = cell(n,1);
lambda_max = cell(n,1);
if nMI > 0
   for i = 1:n
      [~,y_min(i),exitflag_min(i),~,lambda_min{i}] = fmincon(@fun_min_MI,x_start(2*i-1,:)',[],[],[],[],...
         H(:,1),H(:,2),@neq_MI,opt);
      [~,y_max(i),exitflag_max(i),~,lambda_max{i}] = fmincon(@fun_max_MI,x_start(2*i,:)',[],[],[],[],...
         H(:,1),H(:,2),@neq_MI,opt);
   end
else
   for i = 1:n
      [~,y_min(i),exitflag_min(i),~,lambda_min{i}] = fmincon(@fun_min,x_start(2*i-1,:)',[],[],[],[],...
         H(:,1),H(:,2),@neq,opt);
      [~,y_max(i),exitflag_max(i),~,lambda_max{i}] = fmincon(@fun_max,x_start(2*i,:)',[],[],[],[],...
         H(:,1),H(:,2),@neq,opt);
   end
end
y_min(exitflag_min<=0) = [];
y_max(exitflag_max<=0) = [];

if isempty(y_min)
   y(1) = inf;
else
   y(1) = min(y_min);
end
if isempty(y_max)
   y(2) = -inf;
else
   y(2) = -min(y_max);
end


   function [y,gy] = fun_min(x)
      rr = sqrt(x);
      drr = 0.5/rr;
      c1 = v0/rr;
      dc1 = -c1*drr/rr;
      y = c1*sin(rr*tpred)+c2*cos(rr*tpred);
      gy = drr*tpred*(c1*cos(rr*tpred)-c2*sin(rr*tpred))+...
         dc1*sin(rr*tpred);
   end

   function [y,gy] = fun_max(x)
      rr = sqrt(x);
      drr = 0.5/rr;
      c1 = v0/rr;
      dc1 = -c1*drr/rr;
      y = -c1*sin(rr*tpred)-c2*cos(rr*tpred);
      gy = -drr*tpred*(c1*cos(rr*tpred)-c2*sin(rr*tpred))-...
         dc1*sin(rr*tpred);
   end

   function [y,gy] = fun_min_MI(x)
      rr = sqrt(x(1));
      drr = 0.5/rr;
      c1 = v0/rr;
      dc1 = -c1*drr/rr;
      Tfit = supplementary_function.makePower(tpred,nMI);
      gy = zeros(nV,1);
      y = c1*sin(rr*tpred)+c2*cos(rr*tpred)+Tfit*x(2:end);
      gy(1) = drr*tpred*(c1*cos(rr*tpred)-c2*sin(rr*tpred))+dc1*sin(rr*tpred);
      gy(2:end) = Tfit;
   end

   function [y,gy] = fun_max_MI(x)
      rr = sqrt(x(1));
      drr = 0.5/rr;
      c1 = v0/rr;
      dc1 = -c1*drr/rr;
      Tfit = supplementary_function.makePower(tpred,nMI);
      gy = zeros(nV,1);
      y = -c1*sin(rr*tpred)-c2*cos(rr*tpred)-Tfit*x(2:end);
      gy(1) = -drr*tpred*(c1*cos(rr*tpred)-c2*sin(rr*tpred))-dc1*sin(rr*tpred);
      gy(2:end) = -Tfit;
   end

   function [c,ceq,g,geq] = neq(x)
      rr = sqrt(x);
      drr = 0.5/rr;
      c1 = v0/rr;
      dc1 = -c1*drr/rr;
      c = zeros(2*nD,1);
      g = zeros(1,2*nD);
      ceq = [];
      geq = [];
      yy = c1*sin(rr*tdata)+c2*cos(rr*tdata);
      c(1:nD) = yy-ydata-eps;
      c(nD+1:end) = ydata-eps-yy;
      g(1:nD) = drr*tdata.*(c1*cos(rr*tdata)-c2*sin(rr*tdata))+dc1*sin(rr*tdata);
      g(nD+1:end) = -g(:,1:nD);
   end

   function [c,ceq,g,geq] = neq_MI(x)
      rr = sqrt(x(1));
      drr = 0.5/rr;
      c1 = v0/rr;
      dc1 = -c1*drr/rr;
      c = zeros(2*nD,1);
      g = zeros(nV,2*nD);
      ceq = [];
      geq = [];
      Tfit = supplementary_function.makePower(tdata,nMI);
      yy = c1*sin(rr*tdata)+c2*cos(rr*tdata)+Tfit*x(2:end);
      c(1:nD) = yy-ydata-eps;
      c(nD+1:end) = ydata-eps-yy;
      g(1,1:nD) = drr*tdata.*(c1*cos(rr*tdata)-c2*sin(rr*tdata))+dc1*sin(rr*tdata);
      g(2:end,1:nD) = Tfit';
      g(:,nD+1:end) = -g(:,1:nD);
   end

end