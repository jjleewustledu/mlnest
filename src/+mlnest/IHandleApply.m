classdef IHandleApply < handle
	%% IHANDLEAPPLY is an interface for application ideas from "Data Analysis:  A Bayesian Tutorial, Second Edition"
    %  by D.S. Sivia and J. Skilling, section 9.3.2  

	%  $Revision$ 
 	%  was created $Date$ 
 	%  by $Author$,  
 	%  last modified $LastChangedDate$ 
 	%  and checked into repository $URL$,  
 	%  developed on Matlab 8.4.0.150421 (R2014b) 
 	%  $Id$ 
 	 

	properties (Abstract)
        ignoredObjFields
        MAX          % # of nested sampling loops, similar to temperature for s.a.
                     % ~O(Ntime)
        MCMC_Counter % counter for explorting single particle (pre-judged # steps); 
                     % nonlinearly affects precision;
                     % ~O(Ntime)
        Measurement  % external data
        Object       % struct for sampling particles
        n            % # of sampling particles \sim (log width of outer prior mass)^{-1}; 
                     % values >\approx 100 can significantly reduce sampling space & broaden histograms;
                     % ~O(1)
        results
        sigma0       % of model estimation used by logLhood()
        STEP_Initial % Initial guess suitable step-size in (0,1); 
                     % 0.01*MCMC_Counter^{-1} < STEP < 10*MCMC_Counter^{-1} improve precision;
                     % ~O(1)
 	end 

	methods (Abstract)
        e = Estimation(this)
        logL = logLhood(this)
        Obj = Prior(this)
        Explore(this, Obj, logLstar)
        Results(this, Samples, nest, logZ)
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

