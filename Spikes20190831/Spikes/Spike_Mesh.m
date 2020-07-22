
load ys.mat;


m  =   150;
n  =   130;

ZB0 = 68; 
ZB1 = 68.5;

TB0 = 421;
TB1 = 428;

ZBS  =  ZB0 + (1:n)/(n+1) * (ZB1 - ZB0);
TBS  =  TB0 + (1:m)/(m+1) * (TB1 - TB0);

RS   = zeros(m,n);

for ii = 1:m
   TB = TBS(ii); 
    for jj = 1:n,
        ZB = ZBS(jj); 
         for i = 1:5,
             par   =   [i, ys];           
            [I_flag, BET, Yh, R] = QS(ZB, TB, par); 
                                 % evaluation over ordered layers;  
              if I_flag == 0,
                     break; 
              end
         end
      RS(ii,jj) = R;    
  end
end

%%

figure (9);  
            subplot(1,2,2);surfc(ZBS,TBS,RS); shading flat; alpha(0.6);
          subplot(1,2,1); contour(ZBS,TBS,RS,250);