%=========================================================================
%|         Q3.m                                                          |
%|                                                                       |
%|           [I_flag, BET, Yh, R] =  Q3(ZB,TB,ys)                        |
%|                                                                       |
%|   This function returns the LSE estimates of the level [ III ]        |
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


function   [I_flag, BET, Yh, R] =  Q3(ZB,TB,ys)
 
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
 
%====================  Case (7)  ========================================
  X1  = [z1 - ZB * ones(n1,1)];
  X2  = [z2.^2 + ( (1/3) * delts((i_max+1):8).^2 -ZB^2 ) .* ones(n2,1),...
           z2 - ZB * ones(n2,1)];
       
  tmp  = (0.5 / delt);
  x2   = tmp * (  (ZB^2 - alf^2)/2 - ZB   * (ZB - alf ));
  x3   = tmp * (  (bet^3 - ZB^3)/3 - ZB^2 * (bet - ZB ));
  x4   = tmp * (  (bet^2 - ZB^2)/2 - ZB   * (bet - ZB ));
  
 X    = [X1,           zeros(n1,2);...
              x2,      x3,  x4;...
         zeros(n2,1),  X2            ];          % Design Matrix
      
%------------------------- LSE -------------------------------------

       Ym = Y - TB;                          % Modified responses;
      BET =  (X' * X) \ (X' * Ym);           % Parameter estimates;
  
    Ymh1  =  X * BET;
      R1  =  Ym - Ymh1;
      R1  = R1' * R1;    
      
      a1  =  0;
      b1  =  BET(1);
      a2  =  BET(2);
      b2  =  BET(3);
   
      BET1 = [a1, b1, a2, b2]';
      
       par1      = [a1, b1, a2, b2, ZB]';         
       I_flag1   =   Constrains_Test(par1);

%====================  Case (5)  ========================================
  X1  = [z1.^2 + ( (1/3) * delts(1:n1).^2 +ZB^2 ) .* ones(n1,1)...
         - 2 * ZB * z1];
  X2  = [z2.^2 + ( (1/3) * delts((i_max+1):8).^2 +ZB^2 ) .* ones(n2,1)...
         - 2 * ZB * z2];
       
  tmp  = (0.5 / delt);
  x1   = tmp * (  (ZB^3 - alf^3)/3 - ZB^2 * (ZB - alf )+...
          (-2 * ZB) * (  (ZB^2 - alf^2)/2 - ZB   * (ZB - alf ))  );  
   x3   = tmp * (  (bet^3 - ZB^3)/3 - ZB^2 * (bet - ZB )+...
          (-2 * ZB) *  (  (bet^2 - ZB^2)/2 - ZB   * (bet - ZB ))  );
  
  
  
 X    = [X1,           zeros(n1,1);...
         x1,           x3;...
         zeros(n2,1),  X2            ];          % Design Matrix
      
%------------------------- LSE -------------------------------------

       Ym = Y - TB;                          % Modified responses;
      BET =  (X' * X) \ (X' * Ym);           % Parameter estimates;
 
     Ymh2  =  X * BET;
      R2  =  Ym - Ymh2;
      R2  = R2' * R2;      
      
      a1  =  BET(1);
      b1  =  -2 * ZB * a1;
      a2  =  BET(2);
      b2  =  -2 * ZB * a2;

      BET2 = [a1, b1, a2, b2]';

       par2      = [a1, b1, a2, b2, ZB]';         
       I_flag2   =   Constrains_Test(par2);

%====================  Case (3)  ========================================
  X1  = [z1.^2 + ( (1/3) * delts(1:n1).^2 -ZB^2 ) .* ones(n1,1),...
           z1 - ZB * ones(n1,1)];
  X2  = [z2 - ZB * ones(n2,1)];
       
  tmp  = (0.5 / delt);
  x1   = tmp * (  (ZB^3 - alf^3)/3 - ZB^2 * (ZB - alf ));
  x2   = tmp * (  (ZB^2 - alf^2)/2 - ZB   * (ZB - alf )); 
  x4   = tmp * (  (bet^2 - ZB^2)/2 - ZB   * (bet - ZB ));
  
 X    = [X1,           zeros(n1,1);...
         x1,  x2,      x4;...
         zeros(n2,2),  X2            ];          % Design Matrix
      
%------------------------- LSE -------------------------------------

       Ym = Y - TB;                          % Modified responses;
      BET =  (X' * X) \ (X' * Ym);           % Parameter estimates;
    
    Ymh3  =  X * BET;
      R3  =  Ym - Ymh3;
      R3  = R3' * R3; 
            
      a1  =  BET(1);
      b1  =  BET(2);
      a2  =  0;
      b2  =  BET(3);

       BET3 = [a1, b1, a2, b2]';
            
       par3      = [a1, b1, a2, b2, ZB]';         
       I_flag3   =   Constrains_Test(par3);
       
       if I_flag1 == 1 && I_flag2 == 1 && I_flag3 == 1
           I_flag = 1;
           BET = BET1;
           Yh = Ymh1 + TB;  R = R1;
           return;
       elseif    I_flag1 == 0 && I_flag2 == 1 && I_flag3 == 1
           I_flag = 0;
           BET = BET1;
               Yh = Ymh1+TB;  
               R = R1;
          return;
       elseif    I_flag1 == 1 && I_flag2 == 0  && I_flag3 == 1
          I_flag = 0;
           BET = BET2;
              Yh = Ymh2 + TB;  
               R = R2;
          return;   
       elseif    I_flag1 == 1 && I_flag2 == 1  && I_flag3 == 0
          I_flag = 0;
           BET = BET3;
              Yh = Ymh3 + TB;  
               R = R3;
          return;    
       elseif   I_flag1 == 0 && I_flag2 == 0  && I_flag3 == 1
           I_flag = 0;
            if R1 <= R2
              BET = BET1;
              Yh = Ymh1 + TB;  
               R = R1;
            else
              BET = BET2;
              Yh = Ymh2 + TB;  
               R = R2;
            end
            return;
       elseif   I_flag1 == 0 && I_flag2 == 1  && I_flag3 == 0
           I_flag = 0;
            if R1 <= R3
              BET = BET1;
               Yh = Ymh1 + TB;  
               R = R1;
            else
              BET = BET3;
              Yh = Ymh3 + TB;  
               R = R3;
            end
            return;
       elseif   I_flag1 == 1 && I_flag2 == 0  && I_flag3 == 0
           I_flag = 0;
            if R2 <= R3
              BET = BET2;
                 Yh = Ymh2 + TB;  
                 R = R2;
            else
              BET = BET3;
               Yh = Ymh3 + TB;  
                R = R3;
            end
            return;
       else     
           I_flag = 0;
            if R2 == min([R1,R2,R3])
               BET = BET2;
                Yh = Ymh2 + TB;  
                 R = R2; 
            elseif R1 == min([R1,R2,R3]) 
               BET = BET1;
                Yh = Ymh1 + TB;  
                 R = R1;
            else
                BET = BET3;
                 Yh = Ymh3 + TB;  
                 R = R3;
            end
            return;
       end
       
       
       
       return;
