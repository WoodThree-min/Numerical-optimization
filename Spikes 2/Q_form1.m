

%=========================================================================
%|                                                                       |
%|         [BET, Yh, R] = Q(ZB, TB, ys)  [Q_form1.m]                     |
%|                                                                       |
%|    Target function to bo minimized by the Newton-Raphson optimization |
%|   method.                                                             |
%|                                                                       |
%|                                                                       |
%=========================================================================


function [BET, Yh, R] = Q(ZB, TB, ys)     
 
  I_start = 0;                   %  Backup ;
  BETb    = ones(4,1);   
  Yhb     = ys;
  Rb      = 0;              
  ib      = 0;
  
for i = 1:5,                           % Choose the best of Q1---Q5;
             par   =   [i, ys];           
  [I_flag, BET, Yh, R] = QS(ZB, TB, par);  % evaluation over ordered layers; 
    
      if i == 1 && I_flag == 0
          BETb = BET;
          Yhb  = Yh;
          Rb   = R;
          ib   = 1;         break; 
      elseif  I_flag == 0 && I_start == 0
          I_start = 1;
          BETb = BET;
          Yhb  = Yh;
          Rb   = R;
          ib   = i;
      elseif I_flag == 0 && R < Rb
          BETb = BET;
          Yhb  = Yh;
          Rb   = R;
          ib   = i;
     end
end

        BET  =  BETb;               % Update with the minimized SSE.
        Yh   =  Yhb;
         R   =  Rb; 
  
  return