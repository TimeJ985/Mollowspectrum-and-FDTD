function C_E = C_fields(Em, Ep, tout)
% Em, Ep : [walkers x nt]
% tout   : [1 x nt]
% C_E    : [nt x 1]

nt   = numel(tout);
nfft = 2^nextpow2(2 * nt);
norm_v = nt : -1 : 1;

A  = fft(Em,       nfft, 2);
B  = fft(conj(Ep), nfft, 2);
cc = ifft(A .* conj(B), [], 2);

% sum over walkers
C_E = sum(cc(:, 1:nt), 1).' ./ (size(Em, 1) * norm_v.');

% subtract steady-state offset
half = round(nt / 2);
C_E  = C_E - mean(C_E(half:end));
end