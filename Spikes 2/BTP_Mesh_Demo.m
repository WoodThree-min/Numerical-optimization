%%------------------------------  Loading Data -----------

%load ys.mat;
 load nys.mat;    ys = YY;
%load YY.mat;     ys = YY;
%load YY1.mat;    ys = YY;


%%---------------- Setting up the range. Depict the function surface
ZB_00 = 68.5;  %[61, 75.9]; 
 TB_00 = max(ys);    
 
 
ZB0 = 66;  %[55, 76];    %---- Range of candidate ZB  
ZB1 = 73; 

TB0 =    TB_00 - 20; %380;              % ---- Range of candidate TB
TB1 =    TB_00 + 20; %450;

m  =   85;    % -------- section number for rows [ TB ];
n  =   90;    % -------- section number for columns [ ZB ];

I_figure = 9;
    hhh = BTP_Mesh(ys,ZB0,ZB1,TB0,TB1,m,n,I_figure);
    
I_figure = 10;
    hhh = BTP_Mesh_Polish(ys,ZB0,ZB1,TB0,TB1,m,n,I_figure);

 