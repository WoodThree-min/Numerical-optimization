%=========================================================================
%|                                                                       |
%|           xs = NR2D(x,par).m                                          |
%|                                                                       |
%|  This program performs the Newton-Raphson optimization method upon    |
%| the 2-dimensional target function Q in a neighbourhood of [x].        |
%|                                                                       |
%|    x        :   point in Euclidean R^2;                               |
%|    par      :   parameters to be transmitted;                         |
%|    Q(x,par) :   target function.                                      |
%|                                                                       |
%=========================================================================

function xs = NR2D(x,par)

    N       =  1000;
   xs       =  zeros(2,N); 
  xs(:,1)   = x;              % records of minimization trajectory;
 
 
  dg    = 1E-2;                        % steepest gradient increment; 
  
  dx    =  1E-6;              % local mesh increments;
  dy    =  1E-6;
  
  d1    = [-dx,0,dx];
  d2    = [dy,dy,dy];
  ds    = [d1,d1,d1; d2, 0*d2, -d2];   % 9 ordered increment vectors

  
    xn    = x;                        %  Current Point;
  for i = 2:N
      X     =   xn * ones(1,9) + ds;          % local mesh (ordered);
      q     =   zeros(9,1);
    for k = 1:9
     q(k) = Q(X(:,k),par);     % target function values over local mesh; 
    end
      
 
      G = [(q(6)-q(4))/(2*dx); (q(2)-q(8))/(2*dy)];  % gradent vector;
  
      H11 = (q(6) -2 * q(5) + q(4)) / (dx * dx);
      H22 = (q(2) -2 * q(5) + q(8)) / (dy * dy);
      H12 = (q(3) + q(7) -q(1) -q(9))/( 4 *  dx * dy );
  
      H   =[H11, H12;  H12,  H22];                  % Hessian matrix;
  
      D1 = sqrt(G' * G);
      D2 = abs(det(H));
  
       if D1 < 1E-14
          xs(:,i) = xn; break;
       elseif D2 < 1E-10
          xs(:,i) = xn - dg * G;       % Steepest descent; 
       else
          xs(:,i) = xn - H \ G;          % Newon-Raphson; 
       end
   
       d3 = (xs(:,i)-xn);
       D3 = sqrt( d3' * d3);
       if  D3 < 1E-13 
             break; 
       end
       
      xn = xs(:,i);               % Updating current point;
  end
  
  xs = xs(:,1:i);                 % Output trajectory;
   
  return
  