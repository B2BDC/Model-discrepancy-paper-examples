function y = B2Bpred_posterior_delta_in_2D_scenario_slice(ds,opt,s_fix,ngrid)
% Calculate the posterior interval of model discrepancy function at a 2D
% scenario space slice, specififed by s_fix

% Output is a ngrid-by-ngrid-by-2 matrix. The 1st and 2nd dimension
% corresponds to the grid point in normalized scenarios. The 3rd dimension
% corresponds to the predicted lower and upper bounds of discrepancy
% function.

nVar = ds.Variables.Length;
if nVar == 5
   error('No model discrepancy function is used')
elseif nVar == 6
   n_poly = 0;
elseif nVar == 9
   n_poly = 1;
else
   n_poly = 2;
end
svalue = s_fix.value;
if strcmpi(s_fix.scenario,'T')
   fB=@(x,n_poly)supplementary_function.calBasis([svalue x],n_poly);
elseif strcmpi(s_fix.scenario,'P')
   fB=@(x,n_poly)supplementary_function.calBasis([x(1) svalue x(2)],n_poly);
elseif strcmpi(s_fix.scenario,'Phi')
   fB=@(x,n_poly)supplementary_function.calBasis([x svalue],n_poly);
else
   error('Invalid scenario parameter name');
end

if ds.isConsistent(opt)
   design_grid = linspace(-1,1,ngrid);
   y = zeros(ngrid,ngrid,2);
   vv = zeros(1,6);
   for i = 1:ngrid
      for j = 1:ngrid
         vecCoef = [vv fB(design_grid([i j]),n_poly)]';
         predQ = generateModel(vecCoef,ds.Variables);
         pp = ds.predictQOI(predQ,opt);
         y(i,j,:) = [max(pp.min) min(pp.max)];
      end
   end
else
   error('Dataset is inconsistent');
end
