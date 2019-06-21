function y = Dataset_consistency_and_kstart_in_feasible_set(flag_usedata,spring_data,envir_data,tdata,z_measure)
% Calculate dataset consistency and whether k* is in the feasible set with
% the inadequate model

% Output is a n_eps-by-(n+1) cell array. The 1st dimension corresponds to
% the tested epsilon. The 2nd dimension corresponds to used model
% discrepancy function delta, with the first column specifying the case
% where no delta function is used.

if flag_usedata
   y = spring_data.CMresult;
else
   k_star = envir_data.k;
   eps = envir_data.eps;
   n_eps = envir_data.n_eps;
   n_poly = envir_data.n_polytest;
   n = envir_data.n_start;
   c2 = envir_data.c2;
   v0 = envir_data.v0;
   rr = sqrt(k_star);
   c1 = v0/rr;
   nD = size(z_measure{1},1);
   opt = optimoptions('fmincon','Display','off','MaxFunctionEvaluations',1e4,...
      'ConstraintTolerance',1e-6,'OptimalityTolerance',1e-14,'MaxIterations',1e3,...
      'GradConstr','on','GradObj','on','DerivativeCheck','off','FiniteDifferenceType','central');
   y = cell(n_eps,n_poly+1);
   for i = 1:n_eps
      zdata = z_measure{i};
      for j = 0:n_poly
         Hc = repmat([-inf inf],j,1);
         consisFlag = supplementary_function.isConsistent(tdata,zdata,eps(i),Hc,c2,v0,n);
         if ~consisFlag
            y{i,j+1} = 'Dataset is inconsistent';
            continue
         else
            y{i,j+1} = 'k* is not in F';
         end
         if j > 0
            x_start = 0.05*randn(n,j);
            yfun = zeros(n,1);
            for k = 1:n
               [~,yfun(k),exitflag] = fmincon(@fun,x_start(k,:)',[],[],[],[],...
                  Hc(:,1),Hc(:,2),@neq_MI,opt);
               if exitflag > 0
                  y{i,j+1} = 'k* is in F';
                  break
               end
            end
         else
            ypred = c1*sin(rr*tdata)+c2*cos(rr*tdata);
            dy = abs(ypred-zdata);
            if all(dy<=eps(i))
               y = true;
            else
               y = false;
            end
         end
      end
   end
end

   function [y,gy] = fun(x)
      y = 0;
      gy = zeros(j,1);
   end

   function [c,ceq,g,geq] = neq_MI(x)
      c = zeros(2*nD,1);
      g = zeros(j,2*nD);
      ceq = [];
      geq = [];
      Tfit = supplementary_function.makePower(tdata,j);      
      yy = c1*sin(rr*tdata)+c2*cos(rr*tdata)+Tfit*x;
      c(1:nD) = yy-zdata-eps;
      c(nD+1:end) = zdata-eps-yy;
      g(:,1:nD) = Tfit';
      g(:,nD+1:end) = -g(:,1:nD);
   end

end