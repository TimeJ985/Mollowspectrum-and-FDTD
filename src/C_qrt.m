function C = C_qrt(sig, D, tout, V, rho_ss, u0)
    % QRT: computes two-time correlations
   

    num_ops = numel(sig);
    C = zeros(num_ops, num_ops, numel(tout));

    exp_factors = exp(-1i * diag(D) * tout);

    
    for i = 1:num_ops
        for j = 1:num_ops
            
            pert0 = sig{j} * rho_ss;       % matrix form
            pert0_vec = pert0(:);         % vectorize

           
            for tt = 1:numel(tout)
                y_t_vec = V * (exp_factors(:,tt) .* (V \ pert0_vec));
                y_t = reshape(y_t_vec, size(rho_ss));
                C(i, j, tt) = trace(sig{i} * y_t);
            end
        end
    end
    C = C- u0 * u0.';
end
