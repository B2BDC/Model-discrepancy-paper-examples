function y = B2Bpred_posterior_interval_of_parameters(ds,opt)
% Calculate B2BDC predicted interval for model and discrepancy-function
% coefficients

% Output is a structure array, with fields Parameter_name and
% Posterior_interval specifying the name and prediction results for each
% parameter

if ds.isConsistent(opt)
   y = struct('Parameter_name','','Posterior_interval',[]);
   nVar = ds.Variables.Length;
   vv = [zeros(1,nVar); eye(nVar)];
   for i = 1:5
      tmpModel = generateModel(vv(:,i),ds.Variables);
      pp = ds.predictQOI(tmpModel,opt);
      y(i).Parameter_name = ['lambda ' num2str(i)];
      y(i).Posterior_interval = [max(pp.min) min(pp.max)];
   end
   vNames = ds.VarNames;
   if nVar > 5
      for i = 6:nVar
         tmpModel = generateModel(vv(:,i),ds.Variables);
         pp = ds.predictQOI(tmpModel,opt);
         y(i).Parameter_name = vNames{i};
         y(i).Posterior_interval = [max(pp.min) min(pp.max)];
      end
   end
else
   error('Dataset is inconsistent');
end
   