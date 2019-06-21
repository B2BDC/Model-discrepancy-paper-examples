function y = generateDataset(H2_data,eps,n_poly)
% Construct the dataset with specified epsilon and degree of model
% discrepancy polynomial

% Output is a B2BDC dataset object

if eps == 0.01
   t_measure = H2_data.t_measure{1};
elseif eps == 0.005
   t_measure = H2_data.t_measure{2};
else
   error('Invalid input epsilon value');
end

if n_poly > 2
   error('Invalid input polynomial degree for the discrepancy function')
end

prior_c = 10;

center_T = H2_data.scenario_center_parameter.T(1); dT = H2_data.scenario_center_parameter.T(2);
center_P = H2_data.scenario_center_parameter.P(1); dP = H2_data.scenario_center_parameter.P(2);
center_Phi = H2_data.scenario_center_parameter.Phi(1); dPhi = H2_data.scenario_center_parameter.Phi(2);
Sdata = H2_data.Surrogate_for_design;
ndata = length(Sdata);
design = H2_data.scenario_design;

if n_poly < 0 
   y = generateDataset('H2 combustion with reduced mechanism');
   for i = 1:ndata
      T=design(i,1); P=design(i,2); Phi=design(i,3);
      T=T*dT+center_T; P=P*dP+center_P; Phi=Phi*dPhi+center_Phi;
      T=1000/T; P=exp(P);
      UnitName = ['Experiment at T=' num2str(round(T)) 'K, P=' num2str(P,2) 'atm, \Phi=' num2str(Phi,3)];
      tmpQ = Sdata{i};
      tUnit = generateDSunit(UnitName,tmpQ,t_measure(i)+[log(1-eps) log(1+eps)],t_measure(i));
      y.addDSunit(tUnit);
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
   Sdata2 = cell(ndata,1);
   for i = 1:ndata
      tmpCoef = zeros(6+ncoeff);
      tmpCoef(1:6,1:6) = Sdata{i}.CoefMatrix;
      T=design(i,1); P=design(i,2); Phi=design(i,3);
      basis = [T P Phi];
      tmpCoef(1,7:end) = 0.5*fB(basis);
      tmpCoef(7:end,1) = 0.5*fB(basis);
      vAll = Sdata{i}.Variables.addList(vCoeff);
      tmpQ = generateModel(tmpCoef,vAll);
      tmpQ.ErrorStats = Sdata{i}.ErrorStats;
      Sdata2{i} = tmpQ;
   end
   y = generateDataset(['H2 combustion with reduced mechanism and discrepancy polynomial ' num2str(n_poly)]);
   for i = 1:ndata
      T=design(i,1); P=design(i,2); Phi=design(i,3);
      T=T*dT+center_T; P=P*dP+center_P; Phi=Phi*dPhi+center_Phi;
      T=1000/T; P=exp(P);
      UnitName = ['Experiment at T=' num2str(round(T)) 'K, P=' num2str(P,2) 'atm, \Phi=' num2str(Phi,3)];
      tmpQ = Sdata2{i};
      tUnit = generateDSunit(UnitName,tmpQ,t_measure(i)+[log(1-eps) log(1+eps)],t_measure(i));
      y.addDSunit(tUnit);
   end
end