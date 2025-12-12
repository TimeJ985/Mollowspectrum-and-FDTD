function C0 = covariance(y0, sig, u0)
% output: second-order cumulants
C0 = zeros(3, 3);  % Preallocate the matrix

for i = 1:3
    for j = 1:3
        C0(i, j) = trace(y0 * sig{i} * sig{j});
    end
end

C0 = C0 - u0 * u0.';

end