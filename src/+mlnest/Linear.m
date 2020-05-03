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
        ignoredObjFields = {'logL' 'logWt'}
        MAX = 500          % # of nested sampling loops, similar to temperature for s.a.
        MCMC_Counter = 50  % MCMC counter (pre-judged # steps)
        n = 25             % # of sampling particles \sim (log width of outer prior mass)^{-1}; reduces sampling space
        STEP_Initial = 0.2 % Initial guess suitable step-size in (0,1)
        
        dt
        Measurement % tracer concentration of single voxel over time, using dimensionless, scaled parameters
        timeFinal
        timeInterpolants
        timeLength
        map
        visualize = true
 	end 

	methods
        function est  = Estimation(this, Obj)
            Obj = this.Obj2native(Obj);            
            est = Obj.slope*this.timeInterpolants + Obj.intercept;
        end
        
        function this = Linear(varargin)
            this = this@mlnest.AbstractApply(varargin{:});
            
            SLOPE = 2;
            INTERCEPT = -20;
            this.Measurement = SLOPE*(0:119) + INTERCEPT;
            this.timeFinal = 119;
            
            this.dt = this.timeFinal/length(this.Measurement);
            this.timeInterpolants = linspace(0, this.timeFinal, length(this.Measurement));
            this.timeLength = length(this.Measurement);
            
            this.map = containers.Map;
            this.map('slope')     = struct('min',  -10, 'max',  10, 'init', 0);
            this.map('intercept') = struct('min', -100, 'max', 100, 'init', 0);
            
            this.sigma0 = 0.01;
        end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

