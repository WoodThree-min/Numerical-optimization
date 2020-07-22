
%% Data uploading


load ys.mat;
%load nys.mat;    ys = YY;
 
 alfs = [0,1,3,5,7,10,13,16,19,22,25,28,31,34,37,40,43,...
        46,49,52,55,58,61,64,67,70,72,74];
 bets = [alfs(2:28), 76];
 
z = alfs + .5 *( bets -alfs);


%% ----------- Parameter Setting ------------------
N    = 2400;                       % Number of iterations;
 
%TB   =   425 ;
ZB   =   68.5 ;

 TB0 = 200;  %[150,500]; [370,400]; [380,390] 
 TB1 = 450;     
 
 
 
 
 
 
 
 
 
 %%  --------------  Real Time Animations
 TBS = TB0 + (TB1-TB0)* (0:(N-1))/(N-1); 
 
 RS   =   (1:N);                 % Records of SSR;
 LS   =   RS;                    % Records of LSE levels;

 JS   = RS;                      % Records of Spikes and Kinks;
 NJ   = 0;
 
 m = 200;
 Z   =  55 + (0:(m-1))/(m-1) * (76 - 55);
  
 
 
 
figure (1); subplot(1,3,[2,3]);
 plot(z(2:28),ys(2:28),'.-','MarkerSize',20);
 h_str = text(10,350,'counting');
for ii = 1:N,
%TB   =   425 + 15* randn(1);
% ZB   =   68.5 + 0.5 *sin((ii/N)* 8 * pi) ;%+0.025 * randn(1);

  TB = TBS(ii);
  x = [z(21:24),ZB,z(26:28)];
  
for i = 1:5,
             par   =   [i, ys];           
  [I_flag, BET, Yh, R] = QS(ZB, TB, par);  % evaluation over ordered layers; 
    
      if I_flag == 0,
        break; 
     end
end
  str1 = ['  Level of LSEs :          [   ', num2str(i), '  ]'];
  str2 = ['  SSE   =    ', num2str(R) ];
  str3 = ['Total No. of iterations: [ ',num2str(N),'  ] '];
  str4 = ['The iteration No.:         [ ',  num2str(ii), ' ]'];    
    
  %   fittings:
  
     a1 = BET(1); b1 = BET(2); a2 = BET(3); b2 = BET(4);
      c1 = TB + a1 * ( - ZB^2) + b1 *(-ZB);   % Ignore a * (delta ^2 /3);
       c2 = TB + a2 * ( - ZB^2) + b2 *(-ZB); 
     y1 = a1 * Z.^2 + b1 * Z + c1;
     y2 = a2 * Z.^2 + b2 * Z + c2;
     
    
  RS(ii) = R;    % SSEs
  LS(ii) = i;    % LSE levels; 
  
  
  % Jumps (possible);
  
  if ii >= 3 
      if  abs(RS(ii)-RS(ii-1)) > 3 * abs(RS(ii-1)-RS(ii-2)) 
         NJ = NJ + 1;
         JS(NJ) = ii;
      end
  end
 
 subplot(1,3,[2,3]); 
    
 plot(z(2:28),ys(2:28),'.-','MarkerSize',20); hold on;
  plot(x,Yh,'o-r','MarkerSize',10);hold off;
       line([40,76],[TB0,TB0],'Color','r');
       line([40,76],[TB1,TB1],'Color','r');
  text(5,400,str1,'FontSize',20);
                   text(5,350,str2,'FontSize',20);
                   text(5,300,str3,'FontSize',20);
                   text(5,250,str4, 'FontSize',20);
       line(Z,y1,'Color','k');
       line(Z,y2,'Color','k');
       
       %--- Jie Qing lines ------------------
       if a1 < - 0.01
           M1 = - 0.5 * b1 /a1; 
           line([M1,M1],[50,500],'Color','b');
       end
       if a2 < - 0.01
           M2 = - 0.5 * b2 /a2; 
           line([M2,M2],[50,500],'Color','r');
       end

                   
 xlim([1,77]);
 ylim([50,500]);
 grid on; grid minor;
 drawnow;
 %pause(0.1);
 
 subplot(1,3,1);
  if ii==1
     plot(TB,R,'.'); xlim([TB0,TB1]); ylim([0,3000]); 
     title('SSE');
     grid on; grid minor;
     hold on;
  else
      plot(TB,R,'.');  hold on;
  end

end

  JS =JS(1:NJ);
  
subplot(1,3,1); 
   str5 = ['No. of SpikesKinks:      ',num2str(NJ)];
   text(TB0 + 0.1 * (TB1-TB0), 2500, str5, 'FontSize',20);
line(TBS,RS); hold off;
%%
figure (8); 
 subplot(2,1,1);
      plot(TBS,RS,'.-');  title('SSE');
 subplot(2,1,2); 
      plot(TBS,LS,'.-'); title('LSE Levels');
        xlim([TB0,TB1]);
        ylim([0,6]);
        grid on; grid minor;
