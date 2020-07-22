


function z = fun_alt(x,par);

     
     z = -80 * sin(x(1))*cos(x(2))*exp( - (x-4)' * [1/40,0;0,1/60] *(x-4) );
         
     %z=   -  100 *exp( - (x-4)' * [1/40,0;0,1/60] *(x-4) );
     
     
     
     return
     