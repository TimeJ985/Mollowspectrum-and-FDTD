function C3 = C_sto_full(s_1, s_2, tout)
% s_1, s_2: [3 x walkers x time]
% tout:     [1 x nt]
% C3: [3 x 3 x nt]

n_c = 3;
nt     = numel(tout);
nfft   = 2^nextpow2(2 * nt);      
norm_v = nt : -1 : 1;              

C3 = zeros(n_c, n_c, nt);

for i = 1:n_c
    for j = 1:n_c
        
        a = squeeze(s_1(i, :, :));
        b = squeeze(s_2(j, :, :));
        
        A  = fft(a,       nfft, 2);    % FFT along time axis
        B  = fft(conj(b), nfft, 2);    
        cc = ifft(A .* conj(B), [], 2);
        
    
        C3(i, j, :) = sum(cc(:, 1:nt), 1) ./ (size(a,1) * norm_v);
    end
end

% subtract steady-state offset
half   = round(nt / 2);
C_mean = mean(C3(:, :, half:end), 3);
C3     = C3 - C_mean;
end