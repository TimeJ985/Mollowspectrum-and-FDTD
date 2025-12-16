%% demofdtd_01
% stochastic trajectories coupled to 1D-FDTD 1-st order Mur-boundary condition 


%% parameters TLS
det = 0.9; % detuning
w = 8; % Rabi
N = 1000;  % walkers
T1 = 1; % relaxation
T2 = 2; % dephasing

%% general parameters
om_c = 100;  % central frequency
T = 2*pi/om_c; % period
dt = T/100; % time increments
nt = 10000; % length time vector
tout = ( 0:( nt -1 ) )*dt; % time vector
%% parameters FDTD
dx = dt; % for 1D-FDTD dt = dx/c, c=1;
eps = 1; % dielectric
xmin = -om_c * dx;
xmax = om_c * dx;
m = round( (xmax - xmin) / dx ) + 1;

grid = grid1d( xmin, xmax, m );
[ n, dx1 ] = deal( grid.n, grid.x( 2 ) - grid.x( 1 ) );
grad = grid.grad;

%% dynamics

 % drift matrix 
 A = [ -1i * det - 1/(2*T1) - 1/(T2) , 0, -0.5i * w; 
           0, 1i * det - 1/(2*T1) - 1/(T2) ,  0.5i * w; 
                -1i * w, 1i * w, -1/T1 ];
 % inhomogenity
 inh = [ 0; 0; -1 ] ;

 % steady state vector
 u0 = - A \ inh;
        
 % Set up two-level system
 few = fewlevel(2); 
 few.ham = [ det, - 0.5 * w; - 0.5 * w, 0 ]; 
 few.lin = { few.trans( 2, 1, sqrt( 1 / T1 ) ), few.trans( 1, 1, sqrt( 2 / T2)) };
            
% diagonalize Liouville operator
[ V, D ] = eig( liouville( few ) );
y0 = rho_steady( D, V );

% operators
sig{1} = few.trans( 2, 1, 1 );   % -
sig{2} = few.trans( 1, 2 ,1 );   % + 
sig{3} = few.trans( 1, 1, 1 ) - few.trans( 2, 2, 1 );  % z

% initial Bloch vector from the Lindblad approach
u1 = cellfun(@(s) trace(y0 * s), sig).';

% initial covariance
covm = covariance(y0, sig, u0);

sigmaa= 0.0014; % variance Gaussian
g = normpdf(grid.x,0,sigmaa); % normpdf (Gaussian)
g_n = g / sum(g*dx1) ; % normalized
%% SDE, FDTD solver
       
% allocation 
[ s_1, s_2 ] = deal(zeros(3, N, numel(tout)));
[ e_endm, e_endp ] = deal( zeros( N, numel(tout) ) );
sde = sdesolver_T1T2(A, covm, dt,inh, N, u0 );
fdtd1 = fdtdsolver( n, dx1, dt, eps, N , grad);
fdtd2 = fdtdsolver( n, dx1, dt, eps, N , grad);

for it = - 1000 : numel(tout)

   sde = sde.step();
   curr_p = g_n * sde.s1(1,:) * exp(-1i * om_c * it * dt ); % source current J-
   curr_m = g_n * sde.s2(2,:) * exp(1i * om_c * it * dt );% source current J+

   fdtd1 = fdtd1.step( curr_p );
   fdtd2 = fdtd2.step( curr_m );

        if it > 0,  s_1(:, :, it) = sde.s1; 
                    s_2(:, :, it) = sde.s2;  
                    e_endp( :, it ) = fdtd1.e( n-5, : ); % E-fields +
                    e_endm( :, it ) = fdtd2.e( n-5, : ); % E-fields -
        end
end


%% correlation E-field

% rescaling E-field
Em = e_endm .* exp( -1i*om_c*tout );
Ep = e_endp .* exp( +1i*om_c*tout );

% steady-state mean E-fields
mu_m = mean(Em(:,end-10:end), 'all');  
mu_p = mean(Ep(:,end-10:end), 'all');  

% correlation E-fields
C_E = zeros( numel(tout),1  );

for it = 1:numel(tout)
        C_E(it) = mean(Em(:,it).*Ep(:,1),1) ;
end

% steady-state fluctiations E-fields
C_E_fin = C_E - mu_p *mu_m;

%scaled steady-state fluctiations E-field 
CE_s = 1/0.25*C_E_fin;
%% correlation atom

% correlation stochastic solutions
corr_sto = C_sto( s_1, s_2, tout, u0 );

corr_stdy_pm = squeeze(corr_sto(2,1,:));

%% Mollow

% frequency range spectrum integral
om_int = linspace( -0.2*om_c, 0.2*om_c, 1000); 

E = exp(1i * om_int.' * tout);  

% compute the Fourier integrals
S_omega_atom = trapz(tout, E .* corr_stdy_pm', 2).' ;
S_omega_fdtd = trapz(tout, E .* CE_s', 2).';

 %% final plots

figure_fdtdplot = figure('Name','fdtdplot');

tiledlayout(1,1)

nexttile
plot( om_int, real( S_omega_atom ) ); hold on
plot( om_int, real( S_omega_fdtd ), 'k' ) 
legend('atom','fields');
xlabel('$\omega \cdot T_{1}$');
ylabel('$S(\omega)$');
pbaspect([1 1 1]);


set(findall(figure_fdtdplot, '-property', 'Linewidth'), 'Linewidth', 1.1);
set(findall(figure_fdtdplot, '-property', 'FontSize'), 'FontSize', 15);
set(findall(figure_fdtdplot, '-property', 'Interpreter'), 'Interpreter', 'latex');
set(findall(figure_fdtdplot, '-property', 'TickLabelInterpreter'), 'TickLabelInterpreter', 'latex');  

