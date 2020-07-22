
%% Data uploading

%load ys.mat;
load nys.mat;    ys = YY;
 
 alfs = [0,1,3,5,7,10,13,16,19,22,25,28,31,34,37,40,43,...
        46,49,52,55,58,61,64,67,70,72,74];
 bets = [alfs(2:28), 76];
 
 z = alfs + .5 *( bets -alfs);

%% ----------- Parameter Setting ------------------
N    = 1200;                       % Number of iterations;
 
ZB   =   68.51 ;

 TB0 = 370;  %[150,500]; [370,400]; [380,390] 
 TB1 = 450;     
  
 
 %%  --------------  Real Time Animations
 TBS = TB0 + (TB1-TB0)* (0:(N-1))/(N-1); 
 
 RS   =   (1:N);                 % Records of SSR;
 LS   =   RS;                    % Records of LSE levels;

 
 m  =   200;
 Z   =  55 + (0:(m-1))/(m-1) * (76 - 55);
  
 Rm = 1E+5;
 ZBm = 68;
 TBm = 450;
 
 
figure (3); subplot(1,3,[2,3]);
 plot(z(2:28),ys(2:28),'.-','MarkerSize',20);
 h_str = text(10,350,'counting');
for ii = 1:N, 
  TB = TBS(ii);  
  ZB   =   68.5 + 4 *sin((ii/N)* 16 * pi) ; %+0.0015 * randn(1);
 
  I_start = 0;                   %  Backup ;
  BETb    = ones(4,1);   
  Yhb     = ys;
  Rb      = 0;              
  ib      = 0;
  
for i = 1:5,
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
      
         if R < Rm
               Rm = R;
              ZBm = ZB;
              TBm = TB;
         end
         
  str1 = ['Total No. of iterations:   [ ',num2str(N),'  ] '];
  str2 = ['The iteration No. :         [ ',  num2str(ii), ' ]'];             
  str3 = ['Level of LSEs :      [   ', num2str(ib), '   ]'];
  str4 = ['SSE   =                  [  ', num2str(R),' ]' ];
  str5 = ['The Min. SSE =     [ ',num2str(Rm),' ]'];
  str6 = ['Est. [ ZB, TB ] =     [ ',num2str(ZBm),', ',...
            num2str(TBm),' ]'];
  
 %  -------- Fittings:
  
     a1 = BET(1); b1 = BET(2); a2 = BET(3); b2 = BET(4);
      c1 = TB + a1 * ( - ZB^2) + b1 *(-ZB);   % Ignore a * (delta ^2 /3);
      c2 = TB + a2 * ( - ZB^2) + b2 *(-ZB); 
     y1 = a1 * Z.^2 + b1 * Z + c1;
     y2 = a2 * Z.^2 + b2 * Z + c2;
    
  RS(ii) = R;    % SSEs
  LS(ii) = ib;    % LSE levels; 
  
 
  
  
 subplot(1,3,[2,3]); 
    
 plot(z(2:28),ys(2:28),'.-','MarkerSize',30); hold on; 
 plot(z(21:28),Yh,'o-r','MarkerSize',10);
 plot(ZB,TB,'.','MarkerSize',30,'Color','r');hold off;

     line([40,76],[TB0,TB0],'Color','r');   % Range of temperature;
     line([40,76],[TB1,TB1],'Color','r');
     
   text(5,450,str1,'FontSize',20);
   text(5,400,str2,'FontSize',20);
   text(5,350,str3,'FontSize',20);
   text(5,300,str4,'FontSize',20);
   text(5,250,str5,'FontSize',20);
   text(5,200,str6,'FontSize',20); 
       
     %--- Jie Qing lines ------------------
       if a1 < - 0.01
           M1 = - 0.5 * b1 /a1; 
           line([M1,M1],[50,500],'Color','b');
       end
       if a2 < - 0.01
           M2 = - 0.5 * b2 /a2; 
           line([M2,M2],[50,500],'Color','b'); 
       end
       
    line(Z,y1,'Color','k');          % Background fitting curves;
         line(Z,y2,'Color','k');
         line([ZB,ZB],[50,TB],'Color','r');   % Current [ZB,TB];
       for jj = 20:28,                     % Marking flue gas boxes;
          line([bets(jj),bets(jj)],[50,60],'LineWidth',5,'Color','r'); 
       end          
       
       xlim([1,77]);
       ylim([50,500]);
       grid on; grid minor;
    drawnow; 
 
 subplot(1,3,1);
  if ii==1
     plot(TB,R,'.','MarkerSize',20); xlim([TB0,TB1]); ylim([0,25000]); 
     ht = title('SSE');
     grid on; grid minor;
     hold on;
  else
      plot(TB,R,'.','MarkerSize',20);  hold on;
  end

end

 
  
subplot(1,3,1);  
     line(TBS,RS);   hold off;

     %%
figure (4); 
 subplot(2,1,1);
      plot(TBS,RS,'.-');  title('SSE');
       xlim([TB0,TB1]); ylim([0,max(RS)]); 
         grid on; grid minor;
 subplot(2,1,2); 
      plot(TBS,LS,'.-'); title('LSE Levels');
        xlim([TB0,TB1]);
        ylim([0,6]);
        grid on; grid minor;
