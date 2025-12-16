function  L = liouville( obj, varargin )
%  LIOUVILLE - Compute Liouville operator.
%  nanoopt toolbox Hohenester, Ulrich. Nano and quantum optics: 
%  an introduction to basic principles and theory. Springer Nature, 2019.

%  Usage for obj = fewlevel :
%    L = liouville( obj )
%    L = liouville( obj, t )
%  Input
%    t    :  time argument
%  Output
%    L    :  Liouville matrix

%  effective Hamiltonian
h = heff( obj, varargin{ : } );

%  number of atomic states
n = obj.n;
%  index for atomic states
ind = reshape( 1 : n ^ 2, [ n, n ] );
%  allocate matrix
L = zeros( n ^ 2 );

%  heff * rho - rho * heff' - ...
for i = 1 : n
for j = 1 : n
for k = 1 : n
  L( ind( i, j ), ind( k, j ) ) =   ...
  L( ind( i, j ), ind( k, j ) ) + h( i, k );
  L( ind( i, j ), ind( i, k ) ) =   ...
  L( ind( i, j ), ind( i, k ) ) - h( j, k )'; 
end
end
end
%  ... + 1i * L_mu * rho * L_mu'
for mu = 1 : numel( obj.lin )
  lin = obj.lin{ mu };
  [ first, second ] = find( lin );
  for it = 1 : numel( first );  i = first( it );  k = second( it );
  for jt = 1 : numel( first );  j = first( jt );  l = second( jt );
    L( ind( i, j ), ind( k, l ) ) =  ...
    L( ind( i, j ), ind( k, l ) ) + 1i * lin( i, k ) * lin( j, l )';
  end
  end
end
