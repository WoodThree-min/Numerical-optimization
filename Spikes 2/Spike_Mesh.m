
load ys.mat;
%load nys.mat;    ys = YY;


m  =   90;
n  =   85;

ZB0 = 61; 
ZB1 = 75.9;

TB0 = 390;
TB1 = 450;

ZBS  =  ZB0 + (1:n)/(n+1) * (ZB1 - ZB0);
TBS  =  TB0 + (1:m)/(m+1) * (TB1 - TB0);

RS   = zeros(m,n);

for ii = 1:m
   TB = TBS(ii); 
    for jj = 1:n,
        ZB = ZBS(jj); 
        
         I_start = 0;                   %  Backup ;
         Rb      = 0;     
         
         for i = 1:5,
            par   =   [i, ys];           
             [I_flag, BET, Yh, R] = QS(ZB, TB, par); 
                       % evaluation over ordered layers; 
               if i == 1 && I_flag == 0 
                     Rb   = R;     break; 
               elseif  I_flag == 0 && I_start == 0
                 I_start = 1;
                    Rb   = R; 
               elseif I_flag == 0 && R < Rb 
                    Rb   = R; 
               end
         end
         
      RS(ii,jj) = Rb;    
  end
end

%%

figure (9);  
            subplot(1,2,2);surfc(ZBS,TBS,RS); shading flat; alpha(0.6);
          subplot(1,2,1); contour(ZBS,TBS,RS,250);