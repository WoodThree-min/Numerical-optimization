 
function hhh = BTP_Mesh_Polish(ys,ZB0,ZB1,TB0,TB1,m,n, I_figure)


% load ys.mat;
% %load nys.mat;    ys = YY;
% 
% %%---------------- Setting up the range. Depict the function surface
% 
% ZB0 = 67;  %[55, 76];    %---- Range of candidate ZB  
% ZB1 = 72; 
% 
% TB0 = 380;              % ---- Range of candidate TB
% TB1 = 450;
% 
% m  =   35;    % -------- section number for rows [ TB ];
% n  =   50;    % -------- section number for columns [ ZB ];
% I_figure = 10;

 ZBS  =  ZB0 + (1:n)/(n+1) * (ZB1 - ZB0);
 TBS  =  TB0 + (1:m)/(m+1) * (TB1 - TB0);

RS   = zeros(m,n);

for ii = 1:m
   TB = TBS(ii); 
    for jj = 1:n,
        ZB = ZBS(jj); 
        
%          I_start = 0;                   %  Backup ;
%          Rb      = 0;     
%          
%          for i = 1:5,
%             par   =   [i, ys];           
%              [I_flag, BET, Yh, R] = QS(ZB, TB, par); 
%                        % evaluation over ordered layers; 
%                if i == 1 && I_flag == 0 
%                      Rb   = R;     break; 
%                elseif  I_flag == 0 && I_start == 0
%                  I_start = 1;
%                     Rb   = R; 
%                elseif I_flag == 0 && R < Rb 
%                     Rb   = R; 
%                end
%          end
%          
%       RS(ii,jj) = Rb;    

  [BET, Yh, R] = Q(ZB, TB, ys);
  RS(ii,jj)    = R;
  end
end

%%

figure (I_figure);  
            subplot(1,2,2);surfc(ZBS,TBS,RS); shading flat; alpha(0.6);
          subplot(1,2,1); contour(ZBS,TBS,RS,450);
          
          drawnow;
          
   hhh = 0;
          
    return
          