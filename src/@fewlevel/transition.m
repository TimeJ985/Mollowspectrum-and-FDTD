function t = transition( obj, f, i, val )
%  TRANSITION - Transition matrix for Lindblad operators.
%  nanoopt toolbox Hohenester, Ulrich. Nano and quantum optics: 
%  an introduction to basic principles and theory. Springer Nature, 2019.

%  Usage for obj = fewlevel :
%    t = transition( obj, f, i, val )
%  Input
%    f    :  final state
%    i    :  initial state
%    val  :  value
%  Output
%    t    :  transition operator

if ~exist( 'val', 'var' ),  val = 1;  end
%  sparse transition matrix
t = sparse( f, i, val, obj.n, obj.n );
