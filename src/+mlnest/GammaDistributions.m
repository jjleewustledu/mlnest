classdef GammaDistributions < handle & mlnest.AbstractApply
	%% GAMMADISTRIBUTIONS  

	%  $Revision$
 	%  was created 28-Jan-2020 20:52:09 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlnest/src/+mlnest.
 	%% It was developed on Matlab 9.7.0.1261785 (R2019b) Update 3 for MACI64.  Copyright 2020 John Joowon Lee.
 	
	properties
        ignoredObjFields = {'logL' 'logWt'}
        MAX = 500          % # of nested sampling loops, similar to temperature for s.a.
        MCMC_Counter = 50  % MCMC counter (pre-judged # steps)
        n = 25             % # of sampling particles \sim (log width of outer prior mass)^{-1}; reduces sampling space
        STEP_Initial = 0.2 % Initial guess suitable step-size in (0,1)
        
        map % containers.Map containing model params as structs with fields:  min, max ,init
        Measurement % tracer concentration of single voxel over time, using dimensionless, scaled parameters
        timeInterpolants
 	end

	methods 
        function est  = Estimation(this, Obj)
            est = this.estimatorGamma_(this.Obj2native(Obj));
        end
        function k = estimatorGamma(this, Obj)
            a = Obj.a;
            b = Obj.b;
            t0 = Obj.t0;
            t = this.timeInterpolants;
            
            if (t(1) >= t0)
                t_   = t - t0;
                k = t_.^a .* exp(-b*t_);
                k = abs(k);
            else 
                t_   = t - t(1);
                k = t_.^a .* exp(-b*t_);
                k = mlnest.GammaDistributions.slide(abs(k), t, t0 - t(1));
            end
            %k = k*b^(a+1)/gamma(a+1); % lossy
            sumk = sum(k);          
            if sumk > eps
                k = k/sumk;
            end
        end
        function k = estimatorGammaP(this, Obj)
            a = Obj.a;
            b = Obj.b;
            t0 = Obj.t0;
            t = this.timeInterpolants;
            
            if (t(1) >= t0)
                t_   = t - t0;
                k = t_.^a .* exp(-b*t_);
                k = abs(k);
            else 
                t_   = t - t(1);
                k = t_.^a .* exp(-b*t_);
                k = mlnest.GammaDistributions.slide(abs(k), t, t0 - t(1));
            end 
            
            k = k .* (1 + Obj.w*this.timeInterpolants);
            sumk = sum(k);
            if sumk > eps
                k = k/sumk;
            end
        end
        function k = estimatorGenGamma(this, Obj)
            a = Obj.a;
            b = Obj.b;
            p = Obj.p;
            t0 = Obj.t0;
            t = this.timeInterpolants;
            
            if (t(1) >= t0) % saves extra flops from slide()
                t_   = t - t0;
                k = t_.^a .* exp(-(b*t_).^p);
                k = abs(k);
            else
                t_   = t - t(1);
                k = t_.^a .* exp(-(b*t_).^p);
                k = mlnest.GammaDistributions.slide(abs(k), t, t0 - t(1));
            end
            %k = abs(k*p*b^(a+1)/gamma((a+1)/p)); % lossy
            sumk = sum(k);         
            if sumk > eps
                k = k/sumk;
            end
        end
        function k = estimatorGenGammaP(this, Obj)
            a = Obj.a;
            b = Obj.b;
            p = Obj.p;
            t0 = Obj.t0;
            t = this.timeInterpolants;
            
            if (t(1) >= t0) % saves extra flops from slide()
                t_   = t - t0;
                k = t_.^a .* exp(-(b*t_).^p);
                k = abs(k);
            else
                t_   = t - t(1);
                k = t_.^a .* exp(-(b*t_).^p);
                k = mlnest.GammaDistributions.slide(abs(k), t, t0 - t(1));
            end
            
            k = k .* (1 + Obj.w*this.timeInterpolants);            
            sumk = sum(k);
            if sumk > eps
                k = k/sumk;
            end
        end
		  
 		function this = GammaDistributions(varargin)
 			%% GAMMADISTRIBUTIONS
 			%  @param measurement
            %  @param timeInterpolants
            %  @param paramMap
            %  @param MAX
            %  @param MCMC_Counter
            %  @param n
            %  @param STEP_Initial
            %  @param sigma0
            %  @param modelName
            
            ip = inputParser;            
            ip.KeepUnmatched = true;
            addParameter(ip, 'measurement', [], @isnumeric)
            addParameter(ip, 'timeInterpolants', [], @isnumeric)
            addParameter(ip, 'paramMap', containers.Map, @(x) isa(x, 'containers.Map'))
            addParameter(ip, 'MAX', [], @isnumeric)
            addParameter(ip, 'MCMC_Counter', [], @isnumeric)
            addParameter(ip, 'n', [], @isnumeric)
            addParameter(ip, 'STEP_Initial', [], @isnumeric)
            addParameter(ip, 'sigma0', [], @isnumeric)
            addParameter(ip, 'modelName', 'GeneralizedGammaDistributionP', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;

            this.Measurement = ipr.measurement;
            this.timeInterpolants = ipr.timeInterpolants; 
            this.map = ipr.paramMap;
            this.MAX = ipr.MAX;
            this.MCMC_Counter = ipr.MCMC_Counter;
            this.n = ipr.n;
            this.STEP_Initial = ipr.STEP_Initial;
            this.sigma0 = ipr.sigma0;
            this.modelName_ = ipr.modelName;
            
            switch this.modelName_
                case 'GammaDistribution'
                    this.ignored = [this.ignored {'p' 'w'}];
                    this.estimatorGamma_ = @(obj_) this.estimatorGamma(obj_);
                case 'GammaDistributionP'
                    this.ignored = [this.ignored 'p'];
                    this.estimatorGamma_ = @(obj_) this.estimatorGammaP(obj_);
                case 'GeneralizedGammaDistribution'
                    this.ignored = [this.ignored 'w'];
                    this.estimatorGamma_ = @(obj_) this.estimatorGenGamma(obj_);
                case 'GeneralizedGammaDistributionP'
                    this.estimatorGamma_ = @(obj_) this.estimatorGenGammaP(obj_);
                otherwise
                    error('mlnest:NotImplementedError', 'GammaDistributions.ctor.modelName_->%s', this.modelName_)                
            end
 		end
    end 
    
    %% PROTECTED
    
    properties (Access = protected)
        estimatorGamma_
        modelName_
    end
    
    methods (Static, Access = protected)  
        function [vec,T] = ensureRow(vec)
            if (~isrow(vec))
                vec = vec';
                T = true;
                return
            end
            T = false; 
        end      
        function conc    = slide(conc, t, Dt)
            %% SLIDE slides discretized function conc(t) to conc(t - Dt);
            %  Dt > 0 will slide conc(t) towards later times t.
            %  Dt < 0 will slide conc(t) towards earlier times t.
            %  It works for inhomogeneous t according to the ability of pchip to interpolate.
            %  It may not preserve information according to the Nyquist-Shannon theorem.  
            
            import mlnest.GammaDistributions;
            [conc,trans] = GammaDistributions.ensureRow(conc);
            t            = GammaDistributions.ensureRow(t);
            
            tspan = t(end) - t(1);
            tinc  = t(2) - t(1);
            t_    = [(t - tspan - tinc) t];   % prepend times
            conc_ = [zeros(size(conc)) conc]; % prepend zeros
            conc_(isnan(conc_)) = 0;
            conc  = pchip(t_, conc_, t - Dt); % interpolate onto t shifted by Dt; Dt > 0 shifts to right
            
            if (trans)
                conc = conc';
            end
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

