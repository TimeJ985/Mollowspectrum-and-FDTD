function C3 = C_sto( s_1, s_2, tout, u0)

% allocation
C3 = zeros( 3 ,3, numel( tout )  );

% correlation stochastic solutions

for i = 1:3
    for j = 1:3
        for it = 1:numel(tout)
        C3(i,j,it) = mean(s_1(i,:,it).*s_2(j,:,1),2);
        end
    end
end
C3 = C3 - u0 * u0.';
end