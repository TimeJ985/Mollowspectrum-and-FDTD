function  obj = init( obj )
%  INIT - Initialize 1D grid.

%  spatial increment
h = obj.pos( 2 ) - obj.pos( 1 );
%  number of space points
obj.n = length( obj.pos );

%  real space derivatives (2nd order)
obj.grad = gallery( 'tridiag', obj.n, - 1,   0, 1 ) / ( h * 2 );
obj.lap  = gallery( 'tridiag', obj.n,   1, - 2, 1 ) / ( h ^ 2 );
%  real space derivatives (4th order)
obj.grad4 =  ...
  gallery( 'toeppen', obj.n,   1, - 8,    0,  8, - 1 ) / ( 12 * h     );
obj.lap4   =  ...
  gallery( 'toeppen', obj.n, - 1,  16, - 30, 16, - 1 ) / ( 12 * h ^ 2 );
%  wavenumber space
obj.wav = 2 * pi * ( 0 : obj.n - 1 )' / obj.n;
obj.igrad = 1i *       sin( obj.wav )   / h;
obj.ilap = - 2 * ( 1 - cos( obj.wav ) ) / h ^ 2;
%  wavenumber space (4th order)
obj.igrad4 = 1i * ( - sin( 2 * obj.wav ) + 8 * sin( obj.wav ) ) / ( 6 * h );
obj.ilap4 = ( - cos( 2 * obj.wav ) + 16 * cos( obj.wav ) - 15 ) / ( 6 * h ^ 2 ); 
