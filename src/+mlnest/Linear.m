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
         k = 1
         
         limits
 	end 

	methods
        function this = Linear
            this.Measurement = 0:119;
            this.timeFinal = 119;
            
            this.dt = this.timeFinal/length(this.Measurement);
            this.timeInterpolants = linspace(0, this.timeFinal, length(this.Measurement));
            this.timeLength = length(this.Measurement);
            
            this.limits.slope     = [0.1     10];
            this.limits.intercept = [-1      1];
        end
        function est  = Estimation(this, Obj)
            slope     = this.uniform2limits(Obj.slope,     this.limits.slope);
            intercept = this.uniform2limits(Obj.intercept, this.limits.intercept);
            
            est = slope*this.timeInterpolants + intercept;
        end
        function logL = logLhood(this, Obj)
            Estimation = this.Estimation(Obj);       
            DiffSquare = (this.Measurement - Estimation).^2;
            S          = sqrt(sum(DiffSquare)/this.timeLength);
            logL       = -this.timeLength*(log(S) + this.logSqrt2pi) - ...
                          (1/(2*S^2))*sum(DiffSquare);
        end
        function logL = logLhood_Cauchy(this, Obj)
            Y0         = this.timeFinal/2;
            Estimation = this.Estimation(Obj);
            logL       = sum(log((Y0/pi)./((this.Measurement - Estimation).^2 + Y0^2)));
        end
        function Obj  = Prior(this, ~)
            Obj = struct( ...
                'slope',     this.UNIFORM, ...
                'intercept', this.UNIFORM, ...
                'logL',  [], ...
                'logWt', []);
            Obj.logL = this.logLhood(Obj);
        end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

