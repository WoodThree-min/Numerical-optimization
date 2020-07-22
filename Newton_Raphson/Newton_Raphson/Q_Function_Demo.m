
%-----------------------------------------------
%   Q_Function_Demo.m
%
%   Demonstrating the effects of Q(x,par)   
% as a general 2-dimensional function.
%
%----------------------------------------------

m = 90; 
n = 80;

Y   = zeros(m,n);
par = 0;

x  =  -8 + 16 * (0:(n-1))/(n-1);    %4 * pi * (0:(n-1))/(n-1);
y  =  -8 + 16 * (0:(m-1))/(m-1) ;        %3 * pi * (0:(m-1))/(m-1);
 
    for i = 1:m
          yy = y(i);
        for j= 1:n
          xx = x(j);
            X = [xx,yy]';
           Y(i,j) = Q(X,par);
        end
    end

figure (1); surfc(x,y,Y); shading flat; alpha(0.9); 
     colormap  hot;
  % axis off;
  