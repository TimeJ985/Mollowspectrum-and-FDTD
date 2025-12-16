function varargout = subarray( a, s )
%  SUBARRAY - Pass arguments to subsref.
%  nanoopt toolbox Hohenester, Ulrich. Nano and quantum optics: 
%  an introduction to basic principles and theory. Springer Nature, 2019.

if length( s ) > 1
  [ varargout{ 1 : nargout } ] = subsref( a, s( 2 : end ) ); 
else
  varargout{ 1 } = a;
end