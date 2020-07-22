
%=========================================================================
%|                                                                       |
%|           NR2D(x,par).m                                               |
%|                                                                       |
%|  This program performs the Newton-Raphson optimization method upon    |
%| the target function Q in a neighbourhood of [x].                      |
%|                                                                       |
%|    x        :   point in Euclidean R^2;                               |
%|    par      :   parameters to be transmitted;                         |
%|    Q(x,par) :   target function.                                      |
%|                                                                       |
%=========================================================================



  x     = [3,3]';      % Initial point;
  par   = 0; 
  
  xs = NR2D(x,par);