


%=========================================================================
%|                                                                       |
%|          [I_flag, BET, Yh, R] =  QS(ZB,TB,par).m                      |
%|                                                                       |
%|    Target function of the Features-Fitting.                           |
%|                                                                       |
%|    ZB   :   Candidate BTP;                                            |
%|    TB   :   Candidate background maximum temperature;                 |
%|    par  :   parameters to be transmitted.                             |
%|      par(1) : Optimization level number {I, II, III, IV, V};          |
%|      par(2:29) : ys(1:28), Flue gas temperatures from box 1 to 28;    |
%|                                                                       |
%|   The outputs are :                                                   |
%|      I_flag  :       [0] Constrains satisfied; [1] otherwise;         |
%|        BET   :  [a1,b1,a2,b2]' LSE estimates of the parameters.       |
%|      Yh      :       Fitting values of [21]-[28] boxes' temperature;  |
%|      R       :       Sum of squares of the residuals.                 |
%|                                                                       |
%=========================================================================

 function [I_flag, BET, Yh, R] = QS(ZB, TB, par)

     I_flag   = 1; 
     I_level  = par(1);
     ys       = par(2:29);
     
   if     I_level == 1
       [I_flag, BET, Yh, R] =  Q1(ZB,TB,ys);
   elseif I_level == 2 
       [I_flag, BET, Yh, R] =  Q2(ZB,TB,ys);
   elseif I_level == 3
       [I_flag, BET, Yh, R] =  Q3(ZB,TB,ys);
   elseif I_level == 4
       [I_flag, BET, Yh, R] =  Q4(ZB,TB,ys);     
   elseif I_level == 5
       [I_flag, BET, Yh, R] =  Q5(ZB,TB,ys);         
   end
   
  return