function fun = eval( obj, fun )
%  EVAL - Evaluate function .
%
%  Usage for obj = grid1d :
%    n = norm( obj, fun )
%  Input
%    fun    :  function f( x )
%  Output
%    fun    :  evaluated function

fun = fun( obj.pos );
