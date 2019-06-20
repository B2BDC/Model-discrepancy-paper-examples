function [x,y,z] = createPatch3D(X,Y,Z)

nd = size(X,1)-1;
nlu = nd^2;
nlr = nd;
x = zeros(4,nlu+nlr);
y = zeros(4,nlu+nlr);
z = zeros(4,nlu+nlr);
% lower and upper surface
id = 1:nlu;
% v1
pX = X(1:end-1,1:end-1)'; 
pY = Y(1:end-1,1:end-1)';
pZ1 = Z(1:end-1,1:end-1,1)';
pZ2 = Z(1:end-1,1:end-1,2)';
x(1,id) = pX(:);
x(1,id+nlu) = pX(:);
y(1,id) = pY(:);
y(1,id+nlu) = pY(:);
z(1,id) = pZ1(:);
z(1,id+nlu) = pZ2(:);
%v2
pX = X(1:end-1,2:end)'; 
pY = Y(1:end-1,2:end)';
pZ1 = Z(1:end-1,2:end,1)';
pZ2 = Z(1:end-1,2:end,2)';
x(2,id) = pX(:);
x(2,id+nlu) = pX(:);
y(2,id) = pY(:);
y(2,id+nlu) = pY(:);
z(2,id) = pZ1(:);
z(2,id+nlu) = pZ2(:);
%v3
pX = X(2:end,2:end)';
pY = Y(2:end,2:end)';
pZ1 = Z(2:end,2:end,1)';
pZ2 = Z(2:end,2:end,2)';
x(3,id) = pX(:);
x(3,id+nlu) = pX(:);
y(3,id) = pY(:);
y(3,id+nlu) = pY(:);
z(3,id) = pZ1(:);
z(3,id+nlu) = pZ2(:);
%v4
pX = X(1:end-1,2:end)';
pY = Y(1:end-1,2:end)';
pZ1 = Z(1:end-1,2:end,1)';
pZ2 = Z(1:end-1,2:end,2)';
x(4,id) = pX(:);
x(4,id+nlu) = pX(:);
y(4,id) = pY(:);
y(4,id+nlu) = pY(:);
z(4,id) = pZ1(:);
z(4,id+nlu) = pZ2(:);

% left and right surface
id = 2*nlu+1:2*nlu+nlr;
x(:,id) = X(1,1);
x(:,id+nlr) = X(1,end);
% v1
y(1,id) = Y(1:end-1,1);
z(1,id) = Z(1:end-1,1,1);
y(1,id+nlr) = Y(1:end-1,end);
z(1,id+nlr) = Z(1:end-1,end,1);
% v2
y(2,id) = Y(2:end,1);
z(2,id) = Z(2:end,1,1);
y(2,id+nlr) = Y(2:end,end);
z(2,id+nlr) = Z(2:end,end,1);
% v3
y(3,id) = Y(2:end,1);
z(3,id) = Z(2:end,1,2);
y(3,id+nlr) = Y(2:end,end);
z(3,id+nlr) = Z(2:end,end,2);
% v1
y(4,id) = Y(1:end-1,1);
z(4,id) = Z(1:end-1,1,2);
y(4,id+nlr) = Y(1:end-1,end);
z(4,id+nlr) = Z(1:end-1,end,2);