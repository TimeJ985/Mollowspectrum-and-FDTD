function y0 = rho_steady( D, V )
% rho steady state
%  nanoopt toolbox Hohenester, Ulrich. Nano and quantum optics: 
%  an introduction to basic principles and theory. Springer Nature, 2019.
    y0 = [ 0, 0; 0, 1 ];
    y0 = V * ( exp( - 1i * diag( D ) * 1e10 ) .* ( V \ y0( : ) ) );
    y0 = reshape( y0, 2, 2 );
end