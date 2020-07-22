


%=========================================================================
%|                                                                       |
%|           Q(x,par).m                                                  |
%|                                                                       |
%|    Target function to test the Newton-Raphson optimization method.    |
%|                                                                       |
%|    x    :   point in Euclidean R^n;                                   |
%|    par  :   parameters to be transmitted.                             |
%|                                                                       |
%=========================================================================

function y = Q(x,par)

  y = fun_alt(x,par);
  
  return