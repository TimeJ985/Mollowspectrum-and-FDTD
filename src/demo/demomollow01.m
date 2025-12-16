%% demomollow01
% Mollow spectra obtained from qrt, grn and sto approach
% w/o phonon coupling
% Fix rabi and vary detuning or fix detuning and vary Rabi


%% parameters

detuning = linspace(-0.030,0.030,40); % meV
rabi = 0.03; % meV

N = 2000; % walkers
tout = linspace(0, 2500, 1000); % (ps)
dt = tout(2) - tout(1);

 
T1 = 400; % relaxation (ps)
T2 = 800; % dephasing  (ps)

%% steady-state fluctuation correlations

% preallocate
[corr_qrt_stdypm,corr_grn_stdypm, corr_sto_stdypm] = deal( cell(numel(detuning), numel(rabi)) );

for idx = 1 : numel(detuning)
    for idy = 1 : numel(rabi)
        
        det = detuning(idx);
        w = rabi(idy);

        A = [ -1i * det - 1/(2*T1) - 1/(T2) , 0, -0.5i * w; 
                0, 1i * det - 1/(2*T1) - 1/(T2) ,  0.5i * w; 
                -1i * w, 1i * w, -1/T1 ];
        inh = [ 0; 0; -1/T1 ] ;
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

        % initialize sdesolver
        sde = sdesolver_T1T2(A, covm, dt,inh, N, u0 );

        % allocation for the stochastic simulations
        [ s_1, s_2 ] = deal(zeros(3, N, numel(tout)));

        for it = -600 : numel(tout)
            sde = sde.step();  % step in stochastic solver
            if it > 0
                s_1(:, :, it) = sde.s1; 
                s_2(:, :, it) = sde.s2;  
            end
        end

        % correlation functions
        corr_qrt = C_qrt( sig, D, tout, V, y0, u0 );
        corr_grn = C_green( covm, A, tout);
        corr_sto = C_sto( s_1, s_2, tout, u0 );

        % steady-state fluctuation correlation
        corr_qrt_stdy = corr_qrt;
        corr_grn_stdy = corr_grn;
        corr_sto_stdy = corr_sto;

        corr_qrt_stdypm{idx, idy} = squeeze( corr_qrt_stdy(2,1,:));
        corr_grn_stdypm{idx, idy} = squeeze( corr_grn_stdy(2,1,:));
        corr_sto_stdypm{idx, idy} = squeeze( corr_sto_stdy(2,1,:));

    end
end

%% spectrum

% frequency range for Mollow spectrum
om = linspace(-2*rabi(end), 2*rabi(end), 1000);  

% coherent spectrum
coh_amp = abs(u0(1))^2;

% coherent spectrum as a spike at omega = 0 (approximate delta)
[~, idx0] = min(abs(om));   
S_coh = zeros(size(om));
S_coh(idx0) = coh_amp; 

[ S_om_qrt_f, S_om_grn_f, S_om_sto_f ] = deal( cell(numel(detuning), numel(rabi)) );


% incoherent spectrum
for idx = 1 : numel(detuning)
    for idy = 1 : numel(rabi)
        
        det = detuning(idx);
        w = rabi(idy);

        corr_qrt_stypm_int = corr_qrt_stdypm{idx, idy};  % correlation for qrt
        corr_grn_stypm_int = corr_grn_stdypm{idx, idy};  % correlation for green
        corr_sto_stypm_int = corr_sto_stdypm{idx, idy};  % correlation for sto

        % exponential matrix for fourier integral
        E = exp(1i * om.' * tout);  

        % fourier integrals
        S_om_qrt = trapz(tout, E .* corr_qrt_stypm_int', 2).' + S_coh*pi;
        S_om_grn = trapz(tout, E .* corr_grn_stypm_int', 2).' + S_coh*pi;
        S_om_sto = trapz(tout, E .* corr_sto_stypm_int', 2).' + S_coh*pi;

        % store the Mollow spectra 
        S_om_qrt_f{idx, idy} = S_om_qrt;
        S_om_grn_f{idx, idy} = S_om_grn;
        S_om_sto_f{idx, idy} = S_om_sto;
    end
end


%% Mollow spectrum matrix

[S_mat_green, S_mat_qrt, S_mat_sto] = deal( zeros(numel(detuning), numel(rabi), numel(om)) );
    
% populate matrix
for i = 1 : numel(detuning)
   for j = 1 : numel(rabi)
      S_mat_qrt(i, j, :) = S_om_qrt_f{i, j}; 
      S_mat_green(i, j, :) = S_om_grn_f{i, j};
      S_mat_sto(i, j, :) = S_om_sto_f{i, j};
   end
end

% extract for fixed Rabi and varying detuning
s_qrt = real( squeeze( S_mat_qrt(:,1,:)));
s_grn = real( squeeze( S_mat_green(:,1,:)));
s_sto = real( squeeze( S_mat_sto(:,1,:)));

% normalize to own maximum
s_qrt_norm = s_qrt/ max( s_qrt(:));
s_grn_norm = s_grn/max( s_grn(:));
s_sto_norm = s_sto/max( s_sto(:));


%% plot

figure_plot = figure('Name','plot');

tiledlayout(2,3,TileSpacing="compact");

nexttile
plot( om*1e3, s_qrt_norm(1,:) ); hold on;
plot( om*1e3, s_grn_norm(1,:),"LineStyle","--" );
plot( om*1e3, s_sto_norm(1,:) );
legend('qrt', 'grn', 'sto');
xlabel('$\hbar \omega$');
ylabel('$S(\omega)$');
title('all plots in $\mu eV$')
nexttile
plot( om*1e3, s_qrt_norm( numel(detuning)/2,:) ); hold on;
plot( om*1e3, s_grn_norm( numel(detuning)/2,:),"LineStyle","--" );
plot( om*1e3, s_sto_norm( numel(detuning)/2,:) );
legend('qrt', 'grn', 'sto');
xlabel('$\hbar \omega$');
ylabel('$S(\omega)$');

nexttile
plot( om*1e3, s_qrt_norm( numel(detuning),:) ); hold on;
plot( om*1e3, s_grn_norm( numel(detuning),:),"LineStyle","--" );
plot( om*1e3, s_sto_norm( numel(detuning),:) );
legend('qrt', 'grn', 'sto');
xlabel('$\hbar \omega$');
ylabel('$S(\omega)$');

nexttile
imagesc(om*1e3, detuning*1e3, s_qrt_norm );
xlabel('$\hbar \omega$');
ylabel('$\hbar \Delta$');
title('qrt');
axis xy;

nexttile
imagesc(om*1e3, detuning*1e3, s_grn_norm );
xlabel('$\hbar \omega$');
ylabel('$\hbar \Delta$');
title('grn');
axis xy;

nexttile
imagesc(om*1e3, detuning*1e3, s_sto_norm );
xlabel('$\hbar \omega$');
ylabel('$\hbar \Delta$');
title('sto');
axis xy;

set(findall(figure_plot, '-property', 'Linewidth'), 'Linewidth', 1.2);
set(findall(figure_plot, '-property', 'FontSize'), 'FontSize', 12);
set(findall(figure_plot, '-property', 'Interpreter'), 'Interpreter', 'latex');
set(findall(figure_plot, '-property', 'TickLabelInterpreter'), 'TickLabelInterpreter', 'latex');


