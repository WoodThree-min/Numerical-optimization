%=========================================================================
%|                                                                       |
%|           [ZB, TB, BET, Yh, R, xs] = BTP_Search(ZB0,TB0,ys)           |
%|                                                                       |
%|  This program performs the Newton-Raphson minimization method upon    |
%| the 2-dimensional target function Q in a neighbourhood of [x].        |
%|                                                                       |
%|    x        :   point in Euclidean R^2;                               |
%|    ys       :   temperature records at boxes: [21]---[28];            |
%|    Q(x,par) :   SSE of the candidate [ZB, TB].                        |
%|                                                                       |
%|  The Outputs are                                                      |
%|     [ZB, TB] :  Optimized virtual BTP [ZB], and the temperature at it.|                                                                     |
%|        BET   :  [a1,b1,a2,b2]' LSE estimates of the parameters.       |
%|        Yh    :       Fitting values of [21]-[28] boxes' temperature;  |
%|        R     :       Sum of squares of the residuals.                 |
%|        xs    :   trajectory of minimization.
%=========================================================================

function [ZB, TB, BET, Yh, R, xs] = BTP_Search(ZB0,TB0,ys)

%% Parameter setting-------------------------------------------
    N       =  500;
    x       =  [ZB0, TB0]';
   xs       =  zeros(2,N); 
  xs(:,1)   = x;              % records of minimization trajectory;
 
  te    = 1E-5;               % tolerance error during evaluation;
  te2   = 1E-2;               % tolerance error of outputs;
  
  dg    = 5.0E-1;                        % steepest gradient increment; 
  
  dx    =  1E-5;              % local mesh increments;
  dy    =  1E-5;
  
  d1    = [-dx,0,dx];
  d2    = [dy,dy,dy];
  ds    = [d1,d1,d1; d2, 0*d2, -d2];   % 9 ordered increment vectors
   

  %% Minimization -----------------------------------------------
    xn    = x;                            %  Current Point;
 
  for i = 2:N,  
      
      X     =   xn * ones(1,9) + ds;          % local mesh (ordered);
      q     =   zeros(9,1); 

      for k = 1:9
           [BETk, Yhk, Rk] = Q(X(1,k), X(2,k), ys);
           if k == 5, 
               BET = BETk; 
               Yh = Yhk; 
               R  = Rk; 
           end
        q(k) = Rk;     % target function values over local mesh; 
    end 
 
      G   = [(q(6)-q(4))/(2*dx); (q(2)-q(8))/(2*dy)];  % gradent vector;
      D1  = sqrt(G' * G);  
        
              H11 = (q(6) -2 * q(5) + q(4)) / (dx * dx);
              H22 = (q(2) -2 * q(5) + q(8)) / (dy * dy);
              H12 = (q(3) + q(7) -q(1) -q(9))/( 4 *  dx * dy ); 
      H   =  [H11, H12;  H12,  H22];                  % Hessian matrix;  
      E  =   eig(H); 
      
     
     if D1 < te  && E(1) >= 0 && E(2) >=0
          xs(:,i) = xn; break;                     % Local MIN. reached. 
     elseif D1 < te  && (E(1) < 0 || E(2) <0 )    % Saddle point or Max;
          vi   =   - dg * randn(2,1);
         %  D3   =  sqrt(vi' * vi);
         % xs(:,i) = xn + vi;                            % Jump out;
     elseif D1 >= te && (E(1) > 0 && E(2) >0)      % Downward Convex 
          vi   =  - H\G;                           % Candidat NR increment;
           D3   =  sqrt(vi' * vi);
          zp   =  xn(1) + vi(1);                 % New candidate ZB   
        if D3 >= 15 ||  zp <  61 || (zp > 76)...
                          || abs(vi(1))> 0.8 || abs(vi(2)) > 10                    
               vi  =  - (dg / D1 ) * G;     %===== Steepest descent =======  
        end
         % xs(:,i) = xn + vi;                 % --- Potential Newton Raphson    
     elseif D1 >= te && (E(1) < 0 || E(2) <0 )  % Not downward convex;
           vi  = - ( D1^2 / (G' * H * G) ) * G;      % Modified NR;
           
           D3   =  sqrt(vi' * vi);
           zp   =  xn(1) + vi(1);    
        if D3 >= 15 ||  zp < 61 || (zp > 76)...
                    || abs(vi(1))> 0.8 || abs(vi(2)) > 10                     
               vi  =  - (dg / D1 ) * G;     %===== Steepest descent =======  
         end
         % xs(:,i) = xn + vi;      % --- Potential Modified Newton Raphson      
     end  
     
      
      vi = 2 * vi; 
      Rt = R + 1;
      while Rt > R
               vi   =      0.5 * vi;
               xt   =      xn  + vi; 
         [~, ~, Rt] =   Q(xt(1), xt(2), ys);
               D3   =      sqrt(vi' * vi);
             %  disp([' i = ', num2str(i),'  D3 = ', num2str(D3)]);
         if D3 < te2
            break; 
         end
      end 
 
         xs(:,i) =  xt;
         xn      =  xt;               % Updating current point;
  end 
  
  xs = xs(:,1:i);                 % Output trajectory;

  ZB = xs(1,i);
  TB = xs(2,i); 
  
  return
  