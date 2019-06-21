function y = Distance_from_feasible_set(ds,opt)
% Calculate the distance of nominal model parameter value from the feasible
% set

% Output is the squred distance, i.e. d_lambda

if ds.isConsistent(opt)
   nVar = ds.Variables.Length;
   tmpCoef = zeros(nVar+1);
   tmpCoef(2:6,2:6) = eye(5);
   tmpQ = generateModel(tmpCoef,ds.Variables);
   pp = ds.predictQOI(tmpQ,opt);
   y = sqrt(max(pp.min));
else
   error('Dataset is inconsistent')
end
