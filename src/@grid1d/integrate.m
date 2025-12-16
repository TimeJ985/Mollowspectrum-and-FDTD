function int = integrate( obj, y )
%  INTEGRATE - Integrate function over grid domain.
%
%  Usage for obj = grid1d :
%    int = integrate( obj, y )
%  Input
%    y      :  integrand
%  Output
%    int    :  integrated function

int = trapz( obj.pos, y, 1 );