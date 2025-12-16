classdef grid1d
  %  One-dimensional grid for solution of FDTD equations.
  
%%  Properties
  properties
    n           %  number of positions
    pos         %  positions of grid
    wav         %  wavenumbers for grid
    grad        %  derivative operator
    grad4       %  derivative operator (4th order accuracy)    
    lap         %  Laplace operator
    lap4        %  Laplace operator (4th order accuracy)
    igrad       %  Fourier transform of derivative operator
    igrad4      %  Fourier transform of derivative operator (4th order accuracy)
    ilap        %  Fourier transform of Laplace operator
    ilap4       %  Fourier transform of Laplace operator (4th order accuracy)
  end
  
%%  Methods
  methods
    
    function obj = grid1d( varargin )
      %  Initialize 1D grid.
      %
      %  Usage :
      %    obj = grid1d( x )
      %    obj = grid1d( xmin, xmax, n )
      %  Input
      %    x        :  array with grid positions
      %    xmin     :  smallest grid point
      %    xmax     :  largest  grid point
      %    n        :  number of grid points
      switch nargin
        case 1
          obj.pos = varargin{ 1 };
        case 3
          obj.pos = reshape( linspace( varargin{ : } ), [], 1 );
      end
      %  initialization
      obj = init( obj );
    end
  
    function display( obj )
      %  Command window display.
      disp( 'grid1d :' );
      disp( struct( 'pos', [ min( obj.pos ), max( obj.pos ), obj.n ] ) );
    end
    
    function [] = plot( obj, y, varargin )
      %  Plot function on grid.
      plot( obj.pos, y, varargin{ : } );
    end
    
  end
  
  methods (Access = private)
    obj = init( obj );
  end
  
end
