function y = B2Bprediction(H2_data,ds,opt)
% Calculate B2BDC predicted interval for the prediction QOIs

% Output is a npred-by-2 matrix, whose 1st dimension corresponds to
% preiction cases.

if ds.isConsistent(opt)
   Spred = H2_data.Surrogate_for_prediction;
   npred = length(Spred);
   y = zeros(npred,2);
   design = H2_data.scenario_prediction;
   prior_c = 10;
   nVar = ds.Variables.Length;
   if nVar == 5
      n_poly = -1;
   elseif nVar == 6
      n_poly = 0;
   elseif nVar == 9
      n_poly = 1;
   else
      n_poly = 2;
   end
   if n_poly < 0
      for i = 1:npred
         pp = ds.predictQOI(Spred{i},opt);
         y(i,:) = [max(pp.min) min(pp.max)];
      end
   else
      cNames = cell(10,1);
      for i = 0:3
         cNames{i+1} = ['c' num2str(i)];
      end
      cc = 5;
      for i = 1:3
         for j = i:3
            cNames{cc} = ['c' num2str(i) num2str(j)];
            cc = cc+1;
         end
      end
      ncoeff=nchoosek(n_poly+3,n_poly);
      fB=@(x)supplementary_function.calBasis(x,n_poly);
      Hc = repmat(prior_c*[-1,1],ncoeff,1);
      vCoeff = generateVar(cNames(1:ncoeff),Hc);
      for i = 1:npred
         tmpCoef = zeros(6+ncoeff);
         tmpCoef(1:6,1:6) = Spred{i}.CoefMatrix;
         T=design(i,1); P=design(i,2); Phi=design(i,3);
         basis = [T P Phi];
         tmpCoef(1,7:end) = 0.5*fB(basis);
         tmpCoef(7:end,1) = 0.5*fB(basis);
         vAll = Spred{i}.Variables.addList(vCoeff);
         tmpQ = generateModel(tmpCoef,vAll);
         tmpQ.ErrorStats = Spred{i}.ErrorStats;
         pp = ds.predictQOI(tmpQ,opt);
         y(i,:) = [max(pp.min) min(pp.max)];
      end
   end
else
   error('Dataset is inconsistent')
end
   