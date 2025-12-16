classdef fewlevel
  %  Define few level system.
  %  nanoopt toolbox Hohenester, Ulrich. Nano and quantum optics: 
  %  an introduction to basic principles and theory. Springer Nature, 2019.

  properties
    n           %  number of states
    ham         %  Hamiltonian
    lin = {};   %  Lindblad operators
  end
  
  methods  
    function obj = fewlevel( n )
      %  Initialize few-level system.
      % 
      %  Usage :
      %    obj = fewlevel( n )
      %  Input
      %    n    :  number of states
      obj.n = n;
    end
  end
end
