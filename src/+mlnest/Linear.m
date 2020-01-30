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
         MAX = 1000
         
         Measurement % tracer concentration of single voxel over time, using dimensionless, scaled parameters
         timeFinal
         dt
         timeInterpolants
         timeLength         
         map
 	end 

	methods
        function est  = Estimation(this, Obj)
            slope     = this.uniform2limits(Obj.slope,     this.limits('slope'));
            intercept = this.uniform2limits(Obj.intercept, this.limits('intercept'));
            
            est = slope*this.timeInterpolants + intercept;
        end
        
        function this = Linear
            this.Measurement = 0:119;
            this.timeFinal = 119;
            
            this.dt = this.timeFinal/length(this.Measurement);
            this.timeInterpolants = linspace(0, this.timeFinal, length(this.Measurement));
            this.timeLength = length(this.Measurement);
            
            this.map = containers.Map;
            this.map('slope') = struct('min', 0.1, 'max', 10);
            this.map('intercept') = struct('min', -1, 'max', 1);
            
            this.sigma0 = 0.01;
        end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

