function w = weight( obj )
%  WEIGHT - Integration weights for Crank-Nicolson steps.
%
%  Usage for obj = grid1d :
%    w = weight( obj )
%  Output
%    w      :  integration weights such that integral is sum( w .* fun )

w = obj.pos( 2 ) - obj.pos( 1 );
