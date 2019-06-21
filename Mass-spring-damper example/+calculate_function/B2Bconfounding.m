function y = B2Bconfounding(flag_usedata,spring_data,envir_data,t_data,z_measure,ngrid)
% Calculate the 3D feasible set and its 2D projections. The linear
% delta and eps=0.01 case is used.

% Output is a structure variable

if flag_usedata
   y = spring_data.confounding_mesh;
else
   n = envir_data.n_start;
   c2 = envir_data.c2;
   v0 = envir_data.v0;
   ydata = z_measure{1};
   eps = envir_data.eps(1);
   Hk = envir_data.H_k;
   nD = length(t_data);
   opt = optimoptions('fmincon','Display','off','MaxFunctionEvaluations',1e4,...
      'ConstraintTolerance',1e-6,'OptimalityTolerance',1e-14,'MaxIterations',1e3,...
      'GradConstr','on','GradObj','on','DerivativeCheck','off','FiniteDifferenceType','central');
   n_poly = 2;
   Tfit = supplementary_function.makePower(tdata,n_poly);
   H = repmat([-inf,inf],n_poly,1);
   kdata = linspace(Hk(1),Hk(2),ngrid);
   kmesh = repmat(kdata,ngrid,1);
   c0_range = zeros(ngrid,2);
   c1_range1 = zeros(ngrid,2);
   c1_range2 = zeros(ngrid,2);
   c0mesh = zeros(ngrid);
   c1mesh = zeros(ngrid,ngrid,2);
   for i = 1:ngrid
      x_start = 0.01*randn(2*n,n_poly);
      y_min = zeros(n,1);
      y_max = zeros(n,1);
      x_min = cell(n,1);
      x_max = cell(n,1);
      exitflag_min = zeros(n,1);
      exitflag_max = zeros(n,1);
      for j = 1:n
         [x_min{j},y_min(j),exitflag_min(j)] = fmincon(@fun_min,x_start(2*j-1,:)',[],[],[],[],...
            H(:,1),H(:,2),@neq_MI,opt);
         [x_max{j},y_max(j),exitflag_max(j)] = fmincon(@fun_max,x_start(2*j,:)',[],[],[],[],...
            H(:,1),H(:,2),@neq_MI,opt);
      end
      y_min(exitflag_min<=0) = [];
      y_max(exitflag_max<=0) = [];
      x_min(exitflag_min<=0) = [];
      x_max(exitflag_max<=0) = [];
      if isempty(y_min)
         y_min = inf;
      end
      if isempty(y_max)
         y_max = inf;
      end
      [c0_range(i,1),id] = min(y_min);
      x_min_start = x_min{id}(2);
      [c0_range(i,2),id] = min(y_max);
      x_max_start = x_max{id}(2);
      c0_range(i,2) = -c0_range(i,2);
      c0data = linspace(c0_range(i,1),c0_range(i,2),ngrid);
      c0mesh(:,i) = c0data;
      for j = 1:ngrid
         y_min = zeros(n,1);
         y_max = zeros(n,1);
         x_min = zeros(n,1);
         x_max = zeros(n,1);
         exitflag_min = zeros(n,1);
         exitflag_max = zeros(n,1);
         for k = 1:n
            [x_min(k),y_min(k),exitflag_min(k)] = fmincon(@fun_min2,x_min_start+0.01*randn,[],[],[],[],...
               H(2:end,1),H(2:end,2),@neq_MI2,opt);
            [x_max(k),y_max(k),exitflag_max(k)] = fmincon(@fun_max2,x_max_start+0.01*randn,[],[],[],[],...
               H(2:end,1),H(2:end,2),@neq_MI2,opt);
         end
         y_min(exitflag_min<=0) = [];
         y_max(exitflag_max<=0) = [];
         x_min(exitflag_min<=0) = [];
         x_max(exitflag_max<=0) = [];
         if isempty(y_min)
            y_min = inf;
         end
         if isempty(y_max)
            y_max = inf;
         end
         [c1mesh(j,i,1),id] = min(y_min);
         x_min_start = x_min(id);
         [c1mesh(j,i,2),id] = min(y_max);
         x_max_start = x_max(id);
      end
      c1mesh(:,i,2) = -c1mesh(:,i,2);
   end
   for i = 1:ngrid
      x_start = 0.01*randn(2*n,n_poly);
      y_min = zeros(n,1);
      y_max = zeros(n,1);
      exitflag_min = zeros(n,1);
      exitflag_max = zeros(n,1);
      for j = 1:n
         [~,y_min(j),exitflag_min(j)] = fmincon(@fun_min3,x_start(2*j-1,:)',[],[],[],[],...
            H(:,1),H(:,2),@neq_MI,opt);
         [~,y_max(j),exitflag_max(j)] = fmincon(@fun_max3,x_start(2*j,:)',[],[],[],[],...
            H(:,1),H(:,2),@neq_MI,opt);
      end
      y_min(exitflag_min<=0) = [];
      y_max(exitflag_max<=0) = [];
      if isempty(y_min)
         y_min = inf;
      end
      if isempty(y_max)
         y_max = inf;
      end
      c1_range1(i,:) = [min(y_min) -min(y_max)];
   end
   c0data = linspace(min(c0mesh(:)),max(c0mesh(:)),ngrid);
   for i = 1:ngrid
      x_start = [0.2+0.1*rand(2*n,1) 0.01*randn(2*n,1)];
      y_min = zeros(n,1);
      y_max = zeros(n,1);
      exitflag_min = zeros(n,1);
      exitflag_max = zeros(n,1);
      H = [0.2 0.3;-mag mag];
      for j = 1:n
         [~,y_min(j),exitflag_min(j)] = fmincon(@fun_min4,x_start(2*j-1,:)',[],[],[],[],...
            H(:,1),H(:,2),@neq_MI3,opt);
         [~,y_max(j),exitflag_max(j)] = fmincon(@fun_max4,x_start(2*j,:)',[],[],[],[],...
            H(:,1),H(:,2),@neq_MI3,opt);
      end
      y_min(exitflag_min<=0) = [];
      y_max(exitflag_max<=0) = [];
      if isempty(y_min)
         y_min = inf;
      end
      if isempty(y_max)
         y_max = inf;
      end
      c1_range2(i,:) = [min(y_min) -min(y_max)];
   end
   y.xmesh = kmesh;
   y.ymesh = c0mesh;
   y.zmesh = c1mesh;
   y.xdata = kdata;
   y.yVSx = c0_range;
   y.zVSx = c1_range1;
   y.ydata = c0data;
   y.zVSy = c1_range2;
end


   function [y,gy] = fun_min(x)
      y = x(1);
      gy = zeros(n_poly,1);
      gy(1) = 1;
   end

   function [y,gy] = fun_max(x)
      y = -x(1);
      gy = zeros(n_poly,1);
      gy(1) = -1;
   end

   function [y,gy] = fun_min2(x)
      y = x;
      gy = 1;
   end

   function [y,gy] = fun_max2(x)
      y = -x;
      gy = -1;
   end

   function [y,gy] = fun_min3(x)
      y = x(2);
      gy = zeros(n_poly,1);
      gy(2) = 1;
   end

   function [y,gy] = fun_max3(x)
      y = -x(2);
      gy = zeros(n_poly,1);
      gy(2) = -1;
   end

   function [y,gy] = fun_min4(x)
      y = x(2);
      gy = zeros(n_poly,1);
      gy(2) = 1;
   end

   function [y,gy] = fun_max4(x)
      y = -x(2);
      gy = zeros(n_poly,1);
      gy(2) = -1;
   end

   function [c,ceq,g,geq] = neq_MI(x)
      rr = sqrt(kdata(i));
      c1 = v0/rr;
      c = zeros(2*nD,1);
      g = zeros(n_poly,2*nD);
      ceq = [];
      geq = [];
      yy = c1*sin(rr*tdata)+c2*cos(rr*tdata)+Tfit*x;
      c(1:nD) = yy-ydata-eps;
      c(nD+1:end) = ydata-eps-yy;
      g(:,1:nD) = Tfit';
      g(:,nD+1:end) = -g(:,1:nD);
   end

   function [c,ceq,g,geq] = neq_MI2(x)
      rr = sqrt(kdata(i));
      c1 = v0/rr;
      c = zeros(2*nD,1);
      g = zeros(n_poly-1,2*nD);
      ceq = [];
      geq = [];
      yy = c1*sin(rr*tdata)+c2*cos(rr*tdata)+Tfit*[c0data(j);x];
      c(1:nD) = yy-ydata-eps;
      c(nD+1:end) = ydata-eps-yy;
      g(:,1:nD) = Tfit(:,2:end)';
      g(:,nD+1:end) = -g(:,1:nD);
   end

   function [c,ceq,g,geq] = neq_MI3(x)
      rr = sqrt(x(1));
      drr = 0.5/rr;
      c1 = v0/rr;
      dc1 = -c1*drr/rr;
      c = zeros(2*nD,1);
      g = zeros(n_poly,2*nD);
      ceq = [];
      geq = [];
      yy = c1*sin(rr*tdata)+c2*cos(rr*tdata)+Tfit*[c0data(i);x(2)];
      c(1:nD) = yy-ydata-eps;
      c(nD+1:end) = ydata-eps-yy;
      g(1,1:nD) = drr*tdata.*(c1*cos(rr*tdata)-c2*sin(rr*tdata))+dc1*sin(rr*tdata);
      g(2,1:nD) = Tfit(:,2:end)';
      g(:,nD+1:end) = -g(:,1:nD);
   end

end