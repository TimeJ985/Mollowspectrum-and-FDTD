function C1 = C_green(Cov, A, tout)
    % greens functions propagator 
    Ccell = arrayfun(@(t) expm(A * t) * Cov, tout, 'UniformOutput', false);

    C1 = cat(3, Ccell{:});
end
