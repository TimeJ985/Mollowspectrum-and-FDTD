function varargout = subsref( obj, s )
%  Derived properties for objects of class fewlevel.
%  nanoopt toolbox Hohenester, Ulrich. Nano and quantum optics: 
%  an introduction to basic principles and theory. Springer Nature, 2019.

switch s( 1 ).type
case '.'
  switch s( 1 ).subs
    case 'trans'
      varargout{ 1 } = transition( obj, s( 2 ).subs{ : } );
    otherwise
      [ varargout{ 1 : nargout } ] = builtin( 'subsref', obj, s );  
  end
end
