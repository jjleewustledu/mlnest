classdef IApply  
	%% IAPPLY is an interface for application ideas from "Data Analysis:  A Bayesian Tutorial, Second Edition"
    %  by D.S. Sivia and J. Skilling, section 9.3.2  

	%  $Revision$ 
 	%  was created $Date$ 
 	%  by $Author$,  
 	%  last modified $LastChangedDate$ 
 	%  and checked into repository $URL$,  
 	%  developed on Matlab 8.4.0.150421 (R2014b) 
 	%  $Id$ 
 	 

	properties (Abstract)
        MAX    % int
        Measurement
        Object % struct
        n      % int
        sigma0
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

