function Tfit = makePower(T,n)
if size(T,1) == 1
   T = T';
end
nT = length(T);
T = T;
Tfit = repmat(T,1,n).^repmat(0:n-1,nT,1);
end