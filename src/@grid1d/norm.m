function n = norm( obj, psi )
%  NORM - Norm of wavefunction.
%
%  Usage for obj = grid1d :
%    n = norm( obj, psi )
%  Input
%    psi    :  wavefunction
%  Output
%    n      :  norm of wavefunction

n = sqrt( integrate( obj, abs( psi ) .^ 2 ) );
