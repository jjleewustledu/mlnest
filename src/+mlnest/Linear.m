classdef Linear < mlnest.AbstractApply 
	%% LINEAR implements application ideas from "Data Analysis:  A Bayesian Tutorial, Second Edition"
    %  by D.S. Sivia and J. Skilling, section 9.3.2

	%  $Revision$ 
 	%  was created $Date$ 
 	%  by $Author$,  
 	%  last modified $LastChangedDate$ 
 	%  and checked into repository $URL$,  
 	%  developed on Matlab 8.4.0.150421 (R2014b) 
 	%  $Id$  	 
    
	properties 
 		 n   = 100
         MAX = 3000
         
         Measurement % tracer concentration of single voxel over time, using dimensionless, scaled parameters
         timeFinal
         dt
         timeInterpolants
         timeLength         
         map
 	end 

	methods
        function est  = Estimation(this, Obj)
            Obj = this.Obj2native(Obj);            
            est = Obj.slope*this.timeInterpolants + Obj.intercept;
        end
        
        function this = Linear
            SLOPE = 2;
            INTERCEPT = -20;
            this.Measurement = SLOPE*(0:119) + INTERCEPT;
            this.timeFinal = 119;
            
            this.dt = this.timeFinal/length(this.Measurement);
            this.timeInterpolants = linspace(0, this.timeFinal, length(this.Measurement));
            this.timeLength = length(this.Measurement);
            
            this.map = containers.Map;
            this.map('slope') = struct('min', -10, 'max', 10, 'init', 1);
            this.map('intercept') = struct('min', -100, 'max', 100, 'init', 0);
            
            this.sigma0 = 0.01;
        end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

