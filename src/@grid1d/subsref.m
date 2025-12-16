function varargout = subsref( obj, s )
%  Access to class properties and methods.
%
%  Usage for obj = grid1d :
%    obj.x
%    obj.y
%    obj.z                      :  positions of grid
%    obj.kx
%    obj.ky 
%    obj.kz                     :  wavenumbers of grid
%    obj.fft(  fun )            :  fast Fourier transform on grid
%    obj.ifft( fun )            :  inverse fast Fourier transform on grid
%    obj.eval( fun )            :  evaluate function
%    obj.inner( lhs, rhs )      :  inner product
%    obj.integrate( fun )       :  integrate function over grid domain
%    obj.norm( psi )            :  norm of wavefunction
%    obj.normalize( psi )       :  normalize wavefunction
%    obj.weight                 :  integration weight

switch s( 1 ).type
  
case '.'
  switch s( 1 ).subs
    %  positions
    case { 'x', 'y', 'z' }
      varargout{ 1 } = subarray( obj.pos, s );
    %  wavenumbers
    case { 'kx', 'ky', 'kz' }
      varargout{ 1 } = subarray( obj.wav, s );
    %  class methods
    case { 'fft', 'ifft', 'eval', 'inner', 'integrate', 'norm', 'normalize' }
      [ varargout{ 1 : nargout } ] =  ...
        feval( s( 1 ).subs, obj, s( 2 ).subs{ : } );
    %  integration weight
    case 'weight'
      varargout{ 1 } = subarray( weight( obj ), s );
    %  default subsref operator
    otherwise
      [ varargout{ 1 : nargout } ] = builtin( 'subsref', obj, s );
  end
  
end
    
