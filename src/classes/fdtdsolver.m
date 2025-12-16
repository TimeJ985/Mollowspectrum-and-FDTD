classdef fdtdsolver
    % 1-D FDTD simulation 

    properties
       n        % # grid points
       dx       % grid step size
       dt       % timestep size
       eps      % dielectric
       grid     % grid object
       N        % Walkers
       grad     % gradient
       e           % electric field
       h           % magnetic field
      
    end

    methods
        function obj = fdtdsolver( n, dx, dt, eps, N, grad)

          obj.n = n;
          obj.dx = dx;
          obj.dt = dt;
          obj.eps = eps;
          obj.grad = grad;
          obj.N = N;
          [obj.e, obj.h] = deal(zeros(obj.n, N, 1));
          

        end

        function obj = step( obj, curr )
                eold = obj.e;
                hold = obj.h;
           
                % Propagate electric field
                obj.e = obj.e + obj.dt * (curr -  obj.grad  * obj.h) ./ obj.eps;
                
                % 1st order Mur boundary condition (electric field)
                obj.e(1, :) = eold(2, :) + (obj.dt - obj.dx) / (obj.dt + obj.dx) * (obj.e(2, :) - eold(1, :));
                obj.e(obj.n, :) = eold(obj.n - 1, :) + (obj.dt - obj.dx) / (obj.dt + obj.dx) * (obj.e(obj.n - 1, :) - eold(obj.n, :));
                
                % Propagate magnetic field
                obj.h = obj.h - obj.dt * ( obj.grad  * obj.e);
                
                % 1st order Mur boundary condition (magnetic field)
                obj.h(1, :) = hold(2, :) + (obj.dt - obj.dx) / (obj.dt + obj.dx) * (obj.h(2, :) - hold(1, :));
                obj.h(obj.n, :) = hold(obj.n - 1, :) + (obj.dt - obj.dx) / (obj.dt + obj.dx) * (obj.h(obj.n - 1, :) - hold(obj.n, :));
            
        end 
    end 
end
 