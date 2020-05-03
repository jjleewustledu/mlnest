classdef GammaDistributions2 < mlnest.GammaDistributions
	%% GAMMADISTRIBUTIONS2 has ~19% speed up compared to GammaDistributions but may have degraded results
    %  if native |params| >> 1

	%  $Revision$
 	%  was created 28-Jan-2020 20:52:09 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlnest/src/+mlnest.
 	%% It was developed on Matlab 9.7.0.1261785 (R2019b) Update 3 for MACI64.  Copyright 2020 John Joowon Lee.meInterpolants

	methods 
        function est  = Estimation(this, Obj)
            est = this.estimatorGamma_(Obj);
        end
		  
 		function this = GammaDistributions2(varargin)
 			%% GAMMADISTRIBUTIONS2
            
            this = this@mlnest.GammaDistributions(varargin{:});
 		end
    end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

