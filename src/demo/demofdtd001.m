%% demofdtd001
%% Parameters
% two-level system
det = -0.08;        % detuning
w   = 1.0;          % Rabi freq
N   = 10;         % walkers
T1  = 12;           % relaxation time
T2  = 24;           % dephasing time

% propagation
om_c          = 80; % carrier frequency
eps           = 1;  % dielectric function
pts_per_cycle = 30; 
lambda        = 2*pi / om_c;
dt            = 2*pi / (om_c * pts_per_cycle);

% spatial grid
n_lambda = 3;
L        = n_lambda * lambda;
sigma    = 0.09 * lambda;
grid     = grid1d(-L/2, L/2, round(L/dt) + 1);
[n, dx, grad] = deal(grid.n, grid.x(2) - grid.x(1), grid.grad);

% time grid
T_total  = 200;
nt       = round(T_total / dt);
tout     = (0:nt-1) * dt;

% observation
n_warmup = round(4 * T1 / dt);
i_obs    = round(0.5 * grid.n);

% gaussian launch
g   = normpdf(grid.x, 0, sigma);
g_n = g / sum(g * dx);


%% dynamics

 % drift matrix 
 A = [ -1i * det - 1/(2*T1) - 1/(T2) , 0, -0.5i * w; 
           0, 1i * det - 1/(2*T1) - 1/(T2) ,  0.5i * w; 
                -1i * w, 1i * w, -1/T1 ];

 inh = [ 0; 0; -1/T1 ] ;

 u0 = - A \ inh;
        
 % Set up two-level system
 few = fewlevel(2); 
 few.ham = [ det, - 0.5 * w; - 0.5 * w, 0 ]; 
 few.lin = { few.trans( 2, 1, sqrt( 1 / T1 ) ), few.trans( 1, 1, sqrt( 2 / T2)) };
            

[ V, D ] = eig( liouville( few ) );
y0 = rho_steady( D, V );


sig{1} = few.trans( 2, 1, 1 );   % -
sig{2} = few.trans( 1, 2 ,1 );   % + 
sig{3} = few.trans( 1, 1, 1 ) - few.trans( 2, 2, 1 );  % z

covm = covariance(y0, sig, u0);

%% SDE, FDTD solver
       
% allocation
om_dt = om_c * dt;
curr_m = zeros(n, N);
curr_p = zeros(n, N);
[ s_1, s_2 ] = deal(zeros(3, N, numel(tout)));
[ e_endm, e_endp ] = deal( zeros( N, numel(tout) ) );

% intialization
sde = sdesolver_T1T2(A, covm, dt,inh, N, u0 );
fdtd1 = fdtdsolver( n, dx, dt, eps, N , grad);
fdtd2 = fdtdsolver( n, dx, dt, eps, N , grad);


%% main loop
for it = -n_warmup : numel(tout)
   sde  = sde.step();
   curr_m = g_n * sde.s1(1,:) * exp(-1i * om_dt * it); 
   curr_p = g_n * sde.s2(2,:) * exp( 1i * om_dt * it);
   fdtd1  = fdtd1.step( curr_m );
   fdtd2  = fdtd2.step( curr_p );
   if it > 0
       s_1(:, :, it) = sde.s1; 
       s_2(:, :, it) = sde.s2;  
       e_endm(:, it) = fdtd1.e( i_obs, : );  
       e_endp(:, it) = fdtd2.e( i_obs, : );
   end
end

%% correlation E-field

Ep = e_endp .* exp(-1i*om_c*tout);  
Em = e_endm .* exp( 1i*om_c*tout);

C_E = C_fields(Ep, Em, tout) ; 
%% correlation SDE

corr_sto = C_sto_full( s_1, s_2, tout );

corr_stdy_pm = squeeze( corr_sto(2,1,:) );

diff1 = corr_stdy_pm(1)/C_E(1);

C_E_s= real(diff1)*C_E;                 

%% Mollow

om_int = linspace( -2*w, 2*w, 1000); 

E = exp(1i * om_int.' * tout);        

S_omega_atom = trapz(tout, E .* corr_stdy_pm', 2).' ;
S_omega_fdtd = trapz(tout, E .* C_E_s', 2).';


%% final plot

figure_fdtdplot = figure('Name','fdtdplot');

tiledlayout(1,1)

nexttile
plot( om_int, real( S_omega_atom ) ); hold on
plot( om_int, real( S_omega_fdtd ), 'k' ) 
legend('atom','fields');
xlabel('$\omega / \Omega_{R}$');
ylabel('$S(\omega)$');
pbaspect([1 1 1]);


set(findall(figure_fdtdplot, '-property', 'Linewidth'), 'Linewidth', 1.1);
set(findall(figure_fdtdplot, '-property', 'FontSize'), 'FontSize', 15);
set(findall(figure_fdtdplot, '-property', 'Interpreter'), 'Interpreter', 'latex');
set(findall(figure_fdtdplot, '-property', 'TickLabelInterpreter'), 'TickLabelInterpreter', 'latex');  