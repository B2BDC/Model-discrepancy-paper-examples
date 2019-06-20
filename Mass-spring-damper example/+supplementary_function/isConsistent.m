function y = isConsistent(tdata,zdata,eps,Hc,c2,v0,n)

nD = length(zdata);
opt = optimoptions('fmincon','Display','off','MaxFunctionEvaluations',1e4,...
  'ConstraintTolerance',1e-6,'OptimalityTolerance',1e-14,'MaxIterations',1e3,...
  'GradConstr','on','GradObj','on','DerivativeCheck','off','FiniteDifferenceType','central');
Hx = [0.2 0.3];
nMI = size(Hc,1);
nV = nMI+1;
H = [Hx;Hc];
x_start = zeros(n,nV);
x_start(:,1) = 0.1*rand(n,1)+0.2;
x_start(:,2:end) = randn(n,nMI);
y = false;
if nMI > 0
   for i = 1:n
      [~,~,exitflag] = fmincon(@fun_min,x_start(i,:)',[],[],[],[],...
         H(:,1),H(:,2),@neq_MI,opt);
      if exitflag > 0
         y = true;
         break
      end
   end
else
   for i = 1:n
      [~,~,exitflag] = fmincon(@fun_min,x_start(i,:)',[],[],[],[],...
         H(:,1),H(:,2),@neq,opt);
      if exitflag > 0
         y = true;
         break
      end
   end
end


   function [y,gy] = fun_min(x)
      y = 0;
      gy = zeros(nV,1);
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
      c(1:nD) = yy-zdata-eps;
      c(nD+1:end) = zdata-eps-yy;
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
      c(1:nD) = yy-zdata-eps;
      c(nD+1:end) = zdata-eps-yy;
      g(1,1:nD) = drr*tdata.*(c1*cos(rr*tdata)-c2*sin(rr*tdata))+dc1*sin(rr*tdata);
      g(2:end,1:nD) = Tfit';
      g(:,nD+1:end) = -g(:,1:nD);
   end

end