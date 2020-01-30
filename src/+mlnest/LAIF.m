classdef LAIF < mlnest.AbstractApply 
	%% LAIF implements application ideas from "Data Analysis:  A Bayesian Tutorial, Second Edition"
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
        function this = LAIF(magn)
            this.Measurement = magn;
            this.timeFinal = 119;
            
            this.dt = this.timeFinal/length(this.Measurement);
            this.timeInterpolants = linspace(0, this.timeFinal, length(this.Measurement));
            this.timeLength = length(this.Measurement);
            
            Tf = this.timeFinal;
            M0 = this.Measurement(1);
            this.limits.alpha = [0.1     10];
            this.limits.bet   = [1       2];
            this.limits.beta  = [0.1     10];
            this.limits.Delta = [1/Tf    1];
            this.limits.eps   = [0.01    0.5];
            this.limits.F     = [0.01    100];
            this.limits.gamma = [1/Tf    1];
            this.limits.M0    = [0.01*M0 100*M0];
            this.limits.t0    = [1       Tf/2];
        end
        function est  = Estimation(this, Obj)
            F   = this.uniform2limits(Obj.F,  this.limits.F);
            M0  = this.uniform2limits(Obj.M0, this.limits.M0);
            
            est = F*conv(this.bolus(Obj), this.residue(Obj));
            est = est(1:this.timeLength);            
            est = M0*exp(-this.k.*est);
        end
        function r    = residue(this, Obj)
            D   = this.uniform2limits(Obj.D, this.limits.D);
            b   = this.uniform2limits(Obj.b, this.limits.b);
            r = exp((-D*this.timeInterpolants).^b);
        end
        function conc = bolus(this, Obj)
            e    = this.uniform2limits(Obj.eps, this.limits.eps);
            conc = (1 - e)*this.gammaVariate(Obj) + e*this.steadyState(Obj);
        end 
        function conc = gammaVariate(this, Obj)
            a      = this.uniform2limits(Obj.alpha, this.limits.alpha);
            b      = this.uniform2limits(Obj.beta,  this.limits.beta);
            t0     = this.uniform2limits(Obj.t0,    this.limits.t0);
            t0_idx = ceil(t0);
            g_norm = gamma(a + 1)/b^(a + 1);
            conc0  = this.timeInterpolants.^a.*exp(-b*this.timeInterpolants)/g_norm;      
            conc   = zeros(1, this.timeLength);    
            conc(t0_idx:end) = conc0(1:end+1-t0_idx);
        end
        function conc = steadyState(this, Obj)
            g       = uniform2limits(Obj.gamma, this.limits.gamma);
            t0      = uniform2limits(Obj.t0,    this.limits.t0);
            t0_idx  = ceil(t0);
            ss_norm = this.timeFinal + (exp(-g*this.timeFinal) - 1)/g;
            conc0   = (1 - exp(-g*this.timeInterpolants))/ss_norm;    
            conc    = zeros(1, this.timeLength);
            conc(t0_idx:end) = conc0(1:end+1-t0_idx);
        end        
        function logL = logLhood(this, Obj)
            Estimation = this.Estimation(Obj);       
            DiffSquare = (this.Measurement - Estimation).^2;
            S          = sqrt(sum(DiffSquare)/this.timeLength);
            logL       = -this.timeLength*(log(S) + this.logSqrt2pi) - ...
                          (1/(2*S^2))*sum(DiffSquare);
        end
        function logL = logLhood_Cauchy(this, Obj)
            M0   = this.uniform2limits(Obj.M0, this.limits.M0);
            Est  = this.Estimation(Obj);
            logL = sum(log((M0/pi)./((this.Measurement - Est).^2 + M0^2)));
        end
        function Obj  = Prior(this, ~)
            Obj = struct( ...
                'alpha', rand(), ...
                'bet',   rand(), ...
                'beta',  rand(), ...
                'Delta', rand(), ...
                'eps',   rand(), ...
                'F',     rand(), ...
                'gamma', rand(), ...
                'M0',    rand(), ...
                't0',    rand(), ...
                'logL',  [], ...
                'logWt', []);
            Obj.logL = this.logLhood(Obj);
        end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

