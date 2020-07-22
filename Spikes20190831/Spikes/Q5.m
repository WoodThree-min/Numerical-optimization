%=========================================================================
%|         Q5.m                                                          |
%|                                                                       |
%|           [I_flag, BET, Yh, R] =  Q5(ZB,TB,ys)                        |
%|                                                                       |
%|   This function returns the LSE estimates of the level [ V ]          |
%|optimization of the stable iron ore sintering process.                 |
%|                                                                       |
%|   The input consists of :                                             |
%|        ys:   Flue gas temperatre records from Box 1 to Box 28;        |
%|        TB:  assumed virtual background maximum flue gas temperatre.   |
%|   The output consists of:                                             |
%|        I_flag:  [0] if all the constrains are satisfied; [1] otherwise|
%|        BET   :  [a1,b1,a2,b2]' LSE estimates of the parameters.       |
%|                                                                       |
%|                       School of Math. Shandong University             |
%|                                2019.08.18.                            |
%|                                                                       |
%=========================================================================


function   [I_flag, BET, Yh, R] =  Q5(ZB,TB,ys)
 
alfs   =     [55,58,61,64,67,70,72,74]';
bets   =        [58,61,64,67,70,72,74,76]';
delts  =     0.5 * ( bets - alfs );
z      =     alfs + delts; 

Y      =    ys(21:28)';
 
ind     =  (1:8)';
i_max   =  max(ind(alfs <= ZB) );


n1     =   i_max - 1;
n2     =   8 - i_max; 

z1     =   z( 1 : n1 );             %  Seperating interpreating variables;
z2     =   z( (i_max+1) : 8);

alf    = alfs(i_max);            %  Interval of max flue gas temperature;
bet    = bets(i_max);
delt   = bet - alf;



%-----------------------   Design Matrix ------------------------------

  X1  = [ z1 - ZB * ones(n1,1)];
  X2  = [ z2 - ZB * ones(n2,1)];
       
  tmp  = (0.5 / delt);
   
  x2   = tmp * (  (ZB^2 - alf^2)/2 - ZB   * (ZB - alf )); 
  x4   = tmp * (  (bet^2 - ZB^2)/2 - ZB   * (bet - ZB ));
  
  if n2 > 0 
        X    = [X1,           zeros(n1,1);...
                x2,           x4;...
                zeros(n2,1),  X2            ];          % Design Matrix
  else
        X    = [X1,           zeros(n1,1);...
                x2,           x4            ];      % Design Matrix      
  end
     
      
%--------------------------LSE -------------------------------------

       Ym = Y - TB;                          % Modified responses;
      BET =  (X' * X) \ (X' * Ym);           % Parameter estimates;
      
    Ymh  =  X * BET;
      R  =  Ym - Ymh;
      R  =  R' * R;   
      Yh =  Ymh + TB; 
        
       I_flag   =   0;
       BET      =   [0,BET(1),0,BET(2)]';
       
       return
       