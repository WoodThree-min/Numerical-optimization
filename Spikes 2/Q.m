%=========================================================================
%|                                                                       |
%|         [BET, Yh, R] = Q(ZB, TB, ys)  [Q_form_lambda.m]               |
%|                                                                       |
%|    Target function to bo minimized by the Newton-Raphson optimization |
%|   method.                                                             |
%|                                                                       |
%|                                                                       |
%=========================================================================


function [BET, Yh, R] = Q(ZB, TB, ys)   

alfs   =     [55,58,61,64,67,70,72,74]';
bets   =        [58,61,64,67,70,72,74,76]';
delts  =     0.5 * ( bets - alfs );
z      =     alfs + delts; 

Y      =    ys(21:28)';

 
 
ind     =  (1:8)';
i_max   =  max(ind(alfs <= ZB) );     % alf of the box [ ZB ] falls in;
 
       alf    =  alfs(i_max);
       bet    =  bets(i_max);
       dd     =  0.25;
       d1     =  ZB - alf;
       d2     =  bet - ZB; 
       
        if d1 < dd && i_max >= 3                    % Polishing needed;
           NL     = 2;   
           lambda = 1 - (dd - d1)/( 2 * dd);
           ZB1    = alf + dd;
           ZB2    = alf - dd;
           ZBS    = [ZB1,ZB2]; 
       elseif d2 < dd  && i_max <= 7
           NL     = 2;
           lambda = 1 - (dd - d2)/( 2 * dd);
           ZB1    = bet - dd;
           ZB2    = bet + dd;
           ZBS    = [ZB1,ZB2];
       else
           NL     = 1;                    % No polishing;
           ZBS    = ZB;
       end
        
      ZB_bak  =  ZB; 
      RS      =  ones(NL,1);
       
for j = 1:NL 
    
    ZB  = ZBS(j);
    
  I_start = 0;                   %  Backup ;
  BETb    = ones(4,1);   
  Yhb     = ys;
  Rb      = 0;              
%   ib      = 0;
  
 for i = 1:5,                           % Choose the best of Q1---Q5;
             par   =   [i, ys];           
  [I_flag, BET, Yh, R] = QS(ZB, TB, par);  % evaluation over ordered layers; 
    
      if i == 1 && I_flag == 0
          BETb = BET;
          Yhb  = Yh;
          Rb   = R; 
          break; 
      elseif  I_flag == 0 && I_start == 0
          I_start = 1;
          BETb = BET;
          Yhb  = Yh;
          Rb   = R; 
      elseif I_flag == 0 && R < Rb
          BETb = BET;
          Yhb  = Yh;
          Rb   = R; 
      end
 end
      
    if j  ==  1
        BETj  =  BETb;               % Update with the minimized SSE.
        Yhj   =  Yhb;
    end
         RS(j)   =  Rb; 
end
     BET = BETj;
     Yh  = Yhj;
     
     if NL == 1
        R  = RS;
     else
        R  = lambda * RS(1) + (1 - lambda) * RS(2); 
     end
  
  return