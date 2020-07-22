
load ys.mat;


m  =   150;
n  =   130;

ZB0 = 67; 
ZB1 = 70;

TB0 = 400;
TB1 = 450;

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

figure (9); meshc(ZBS,TBS,RS);
figure (10); contour(ZBS,TBS,RS,250)