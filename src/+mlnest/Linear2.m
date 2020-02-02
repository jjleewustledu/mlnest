classdef Linear2 < mlnest.Linear
	%% LINEAR2 implements application ideas from "Data Analysis:  A Bayesian Tutorial, Second Edition"
    %  by D.S. Sivia and J. Skilling, section 9.3.2

	%  $Revision$ 
 	%  was created $Date$ 
 	%  by $Author$,  
 	%  last modified $LastChangedDate$ 
 	%  and checked into repository $URL$,  
 	%  developed on Matlab 8.4.0.150421 (R2014b) 
 	%  $Id$  	

	methods
        function est  = Estimation(this, Obj)          
            est = Obj.slope*this.timeInterpolants + Obj.intercept;
        end
        
        function this = Linear2(varargin)
            this = this@mlnest.Linear(varargin{:});
        end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

