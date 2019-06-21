function y = calBasis(x,np)

nx = size(x,1);
n_coeff = nchoosek(3+np,np);
y = zeros(nx,n_coeff);
for ii = 1:nx
   yy = zeros(1,n_coeff);
   yy(1) = 1;
   cc = 2;
   cp = 1;
   while cp<=np
      nc = (cp+1)^3;
      bb = zeros(nc,3);
      b3 = cp:-1:0;
      b2 = zeros(1,(cp+1)^2);
      for i = 1:cp+1
         b2((1:cp+1)+(i-1)*(cp+1)) = cp+1-i;
      end
      bb(:,3) = repmat(b3,1,(cp+1)^2);
      bb(:,2) = repmat(b2,1,cp+1);
      for i = 0:cp
         bb((cp+1)^2*i+(1:(cp+1)^2),1) = cp-i;
      end
      flag = sum(bb,2)~=cp;
      bb(flag,:) = [];
      n1 = size(bb,1);
      yy(cc:cc+n1-1) = prod(repmat(x(ii,:),n1,1).^bb,2);
      cc = cc+n1;
      cp = cp+1;
   end
   y(ii,:) = yy;
end