

%%------------------------------  Loading Data -----------

load ys.mat;
% load nys.mat;    ys = YY;
% load YY.mat;     ys = YY;
% load YY1.mat;    ys = YY;


figure (100); plot(ys,'o-'); drawnow;

%%---------------- Setting up the range. Depict the function surface


 ZB_00 =  68.5 ;  %71; %68.5;  %[61, 75.9]; 
 TB_00 = max(ys);    
 
 
ZB0 = 67;  %[55, 76];    %---- Range of candidate ZB  
ZB1 = 73; 

TB0 =    TB_00 - 20; %380;              % ---- Range of candidate TB
TB1 =    TB_00 + 20; %450;

m  =   35;    % -------- section number for rows [ TB ];
n  =   33;    % -------- section number for columns [ ZB ];

I_figure = 10;

% hhh = BTP_Mesh(ys,ZB0,ZB1,TB0,TB1,m,n,I_figure);
  hhh = BTP_Mesh_Polish(ys,ZB0,ZB1,TB0,TB1,m,n,I_figure);

 figure (I_figure); subplot(1,2,1); hold on; 
  
 hcount = title(' counting');


%% ----------- Evaluating the virtual BTP: [ZB, TB], R = SSE

 tic;
 

NN = 100;
 for III = 1:NN
     ZB0  =  ZB_00 +  0.8 * (rand(1)-0.5) * 2;
     TB0  =  TB_00 +   10 * (rand(1)-0.5) * 2;
     
 [ZB, TB, BET, Yh, R, xs] = BTP_Search(ZB0,TB0,ys);  % Evaluating Q-Function;
  
    plot(ZB0,TB0,'o','MarkerSize',5);
    plot(ZB,TB,'.','MarkerSize',40); 
      %line([ZB0,ZB], [TB0, TB]);
          plot(xs(1,:),xs(2,:),'.-');
       
   str_t = ['Counting ',num2str(III), ' of the total ',num2str(NN)...
       '.            Time used: ', num2str(floor(toc)),' seconds'];
   set(hcount,'String',str_t);
     %%%%    disp([' ZB = ',num2str(ZB),',   TB= ',num2str(TB)]);
    drawnow;
 end
     t_end = toc; t_min = floor(t_end / 60); 
                  t_sec = floor (t_end - t_min * 60);
    t_gap =['  Total searching time: ',...
           num2str(t_min ),' minutes  ',...
           num2str(t_sec), ' seconds.'];
    disp(t_gap);
    
    hold off;