classdef sdesolver_T1T2
    % Solves dynamics (one time step) of two level system with stochastic 
    % differential equations. 
    % Input: detuning, Rabi freqeuncy, # of walkers, time increment,
    % relaxation time;
    % Output: stochastic solution for sigma_-, sigma_+ and sigma_z for one
    % time step ✓

    properties
        A       % drift matrix
        covm     % covariance matrix
        N       % # of walkers
        dt      % time increment
        P       % propagator
        inh     % inhomogenity
        u0      % inital Bloch vector
        B1      % noise matrix 1
        B2      % noise matrix 2
        s1      % stochastic solution 1
        s2      % stochastic solution 2
        D       % diffusion matrix
    end

    methods

        function obj = sdesolver_T1T2(A, covm, dt,inh, N, u0 )

            
            obj.A = A;
            obj.covm = covm;
            obj.u0 = u0;
            obj.dt = dt;
            obj.N = N;
            obj.inh = inh;
            obj = obj.initialize();

        end

       function obj = initialize(obj)

            % propagator
            obj.P = expm( obj.A * obj.dt );

            %  diffusion matrix
            obj.D = - obj.A * obj.covm - obj.covm * obj.A .';

            % singular value decomposition
            [U,S,V] = svd(obj.D);

            % noise matrix
            obj.B1 = U*sqrt(S);
            obj.B2 = conj(V)*sqrt(S);

              % allocation stochastic solution
            obj.s1 = repmat(obj.u0, 1, obj.N);  
            obj.s2 = repmat(obj.u0, 1, obj.N);  

       end
        

       function obj = step(obj)

           % rnd vector
           xi = randn(3, obj.N);

           % one time step stochastic solution
           obj.s1 = obj.P * obj.s1 + obj.inh * obj.dt + sqrt(obj.dt) * obj.B1 * xi;
           obj.s2 = obj.P * obj.s2 + obj.inh * obj.dt + sqrt(obj.dt) * obj.B2 * xi; 

       end
    end
end