%=========================================================================
%|         Q1.m                                                          |
%|                                                                       |
%|           [I_flag, BET, Yh, R] =  Q1(ZB,TB,ys)                        |
%|                                                                       |
%|   This function returns the LSE estimates of the level [ I ]          |
%|optimization of the stable iron ore sintering process.                 |
%|                                                                       |
%|   The input consists of :                                             |
%|        ys:   Flue gas temperatre records from Box 1 to Box 28;        |
%|        ZB:  assumed virtual BTP;                                      |
%|        TB:  assumed virtual background maximum flue gas temperatre.   |
%|   The output consists of:                                             |
%|        I_flag:  [0] if all the constrains are satisfied; [1] otherwise|
%|        BET   :  [a1,b1,a2,b2]' LSE estimates of the parameters.       |
%|                                                                       |
%|                       School of Math. Shandong University             |
%|                                2019.08.18.                            |
%|                                                                       |
%=========================================================================


function   [I_flag, BET, Yh, R] =  Q1(ZB,TB,ys)
 
alfs   =     [55,58,61,64,67,70,72,74]';
bets   =        [58,61,64,67,70,72,74,76]';
delts  =     0.5 * ( bets - alfs );
z      =     alfs + delts; 

Y      =    ys(21:28)';
 
ind     =  (1:8)';
i_max   =  max(ind(alfs <= ZB) );


n1     =   i_max - 1;
n2     =   8 - i_max; 

    if  n2 == 0,                 % ZB falls in the last box;
       I_flag = 1;
       BET2   = ones(4,1);
       return;
    end

z1     =   z( 1 : n1 );             %  Seperating interpreating variables;
z2     =   z( (i_max+1) : 8);

alf    = alfs(i_max);            %  Interval of max flue gas temperature;
bet    = bets(i_max);
delt   = bet - alf;



%-----------------------   Design Matrix ------------------------------

  X1  = [z1.^2 + ( (1/3) * delts(1:n1).^2 - ZB^2 ) .* ones(n1,1),...
           z1 - ZB * ones(n1,1)];
  X2  = [z2.^2 + ( (1/3) * delts((i_max+1):8).^2 -ZB^2 ) .* ones(n2,1),...
           z2 - ZB * ones(n2,1)];
       
  tmp  = (0.5 / delt);
  x1   = tmp * (  (ZB^3 - alf^3)/3 - ZB^2 * (ZB - alf ));
  x2   = tmp * (  (ZB^2 - alf^2)/2 - ZB   * (ZB - alf ));
  x3   = tmp * (  (bet^3 - ZB^3)/3 - ZB^2 * (bet - ZB ));
  x4   = tmp * (  (bet^2 - ZB^2)/2 - ZB   * (bet - ZB ));
  
 X    = [X1,           zeros(n1,2);...
         x1,  x2,      x3,  x4;...
         zeros(n2,2),  X2            ];          % Design Matrix
      
%------------------------- LSE -------------------------------------

       Ym = Y - TB;                          % Modified responses;
      BET =  (X' * X) \ (X' * Ym);           % Parameter estimates;
      
%------------------ Constrains Test     ----------------------------
      
      a1  =  BET(1);
      b1  =  BET(2);
      a2  =  BET(3);
      b2  =  BET(4);
     Ymh  =  X * BET;
      R  =  Ym - Ymh;
      R  = R' * R;   
      Yh = Ymh + TB;
      
       par      = [a1, b1, a2, b2, ZB]';         
       I_flag   =   Constrains_Test(par);
       
       
       
       return
       