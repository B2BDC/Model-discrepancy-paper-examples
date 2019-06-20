function y = makePrediction_true(tdata,ydata,tpred,eps,b,c2,v0,n)


nD = size(ydata,1);
opt = optimoptions('fmincon','Display','off','MaxFunctionEvaluations',1e4,...
  'ConstraintTolerance',1e-6,'OptimalityTolerance',1e-14,'MaxIterations',1e3,...
  'GradConstr','on','GradObj','on','DerivativeCheck','off','FiniteDifferenceType','central');
Hx = [0.2 0.3];

% Displacement
x_start = 0.1*rand(2*n,1)+0.2;
y_min = zeros(n,1);
y_max = zeros(n,1);
exitflag_min = zeros(n,1);
exitflag_max = zeros(n,1);
for i = 1:n
   [~,y_min(i),exitflag_min(i)] = fmincon(@fun_min,x_start(2*i-1,:)',[],[],[],[],...
      Hx(:,1),Hx(:,2),@neq,opt);
   [~,y_max(i),exitflag_max(i)] = fmincon(@fun_max,x_start(2*i,:)',[],[],[],[],...
      Hx(:,1),Hx(:,2),@neq,opt);
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
      rr = 0.5*sqrt(4*x-b^2);
      drr = (4*x-b^2)^(-0.5);
      c1 = v0/rr;
      dc1 = -c1*drr/rr;
      y = exp(-0.5*b*tpred)*(c1*sin(rr*tpred)+c2*cos(rr*tpred));
      gy = tpred*drr*exp(-0.5*b*tpred)*(c1*cos(rr*tpred)-c2*sin(rr*tpred))+...
         exp(-0.5*b*tpred)*sin(rr*tpred)*dc1;
   end

   function [y,gy] = fun_max(x)
      rr = 0.5*sqrt(4*x-b^2);
      drr = (4*x-b^2)^(-0.5);
      c1 = v0/rr;
      dc1 = -c1*drr/rr;
      y = -exp(-0.5*b*tpred)*(c1*sin(rr*tpred)+c2*cos(rr*tpred));
      gy = -tpred*drr*exp(-0.5*b*tpred)*(c1*cos(rr*tpred)-c2*sin(rr*tpred))-...
         exp(-0.5*b*tpred)*sin(rr*tpred)*dc1;
   end

   function [c,ceq,g,geq] = neq(x)
      rr = 0.5*sqrt(4*x-b^2);
      drr = (4*x-b^2)^(-0.5);
      c1 = v0/rr;
      dc1 = -c1*drr/rr;
      c = zeros(2*nD,1);
      g = zeros(1,2*nD);
      ceq = [];
      geq = [];
      yy = exp(-0.5*b*tdata).*(c1*sin(rr*tdata)+c2*cos(rr*tdata));
      c(1:nD) = yy-ydata-eps;
      c(nD+1:2*nD) = ydata-eps-yy;
      g(1:nD) = drr*tdata.*exp(-0.5*b*tdata).*(c1*cos(rr*tdata)-c2*sin(rr*tdata))+...
         dc1*exp(-0.5*b*tdata).*sin(rr*tdata);
      g(nD+1:2*nD) = -g(1:nD);
   end
end