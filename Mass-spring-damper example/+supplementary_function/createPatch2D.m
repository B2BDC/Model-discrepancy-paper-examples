function [x,y] = createPatch2D(X,Y)

nd = length(X)-1;
x = zeros(4,nd);
y = zeros(4,nd);
x(1,:) = X(1:end-1);
x(2,:) = X(2:end);
x(3,:) = X(2:end);
x(4,:) = X(1:end-1);
y(1,:) = Y(1:end-1,1);
y(2,:) = Y(2:end,1);
y(3,:) = Y(2:end,2);
y(4,:) = Y(1:end-1,2);