%=========================================================================
%|                                                                       |
%|           NR2D_Demo.m                                                 |
%|                                                                       |
%|  This program performs the Newton-Raphson minimization method upon    |
%| the target function Q in a neighbourhood of [x].                      |
%|                                                                       |
%|    x            :   point in Euclidean R^2;                           |
%|    par          :   parameters to be transmitted;                     |
%|    Q(x,par)     :   target function.                                  |
%|    NR2D(x,par)  :   Newton-Raphson optimization for Q-function        |
%|                                                                       |
%=========================================================================


Q_Function_Demo;   % Resulting : {x, y, Y} mesh of [Q] function.
 
  x0     = [2,2]';      % Initial point;
  par   = 0; 
  
xs = NR2D(x0,par); % Resulting trajectories of Newton-Raphson optimization.
 
figure (2); subplot(1,2,2);
   contour(x,y,Y,50);
   hold on;
    [~,m]  = size(xs);
    plot(xs(1,1),xs(2,1),'o');
    plot(xs(1,m),xs(2,m),'.','MarkerSize',30);
  xlim([-8,8]); 
  ylim([-8,8]); 
   axis square;
   str = num2str(0);
   hhh= title(str);
   
  subplot(1,2,1);
  surfc(x,y,Y); shading flat; alpha(0.8); 
  view([-8,-8,4]);
     colormap  hot; 
  hold on;
  
 for j = 1:250,
   xx = x0 + 1 *randn(2,1);
   xs = NR2D(xx,par);  
   [~,m]  = size(xs);
   subplot(1,2,2);
  %  plot(xs(1,:),xs(2,:),'.-');    % plot the optimization trajectory;
      plot(xs(1,1),xs(2,1),'o');      % Plot only the start-end points.
      plot(xs(1,m),xs(2,m),'.','MarkerSize',30);
      str = num2str(j);
      set(hhh,'String',str);
   drawnow;
   subplot(1,2,1);;
    plot3(xs(1,1),xs(2,1),-100,'o');       % Plot only the start-end points.
    plot3(xs(1,m),xs(2,m),-100,'.','MarkerSize',10);
    xlim([-8,8]);ylim([-8,8]);zlim([-120,120]);
   drawnow;
end 
  
subplot(1,2,1);  hold off;
subplot(1,2,2); hold off;