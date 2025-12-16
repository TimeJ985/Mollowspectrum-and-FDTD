function h = heff( obj, varargin )
%  HEFF - Compute effective Hamiltonian.
%  nanoopt toolbox Hohenester, Ulrich. Nano and quantum optics: 
%  an introduction to basic principles and theory. Springer Nature, 2019.

%  Usage for obj = fewlevel :
%    h = heff( obj )
%    h = heff( obj, t )
%  Input
%    t    :  time argument
%  Output
%    h    :  effective Hamiltonian

h = obj.ham;
%  evaluate function ?
if isa( h, 'function_handle' ),  h = h( varargin{ : } );  end

for mu = 1 : numel( obj.lin )
  lin = obj.lin{ mu };
  h = h - 0.5i * ( lin' ) * lin;
end
