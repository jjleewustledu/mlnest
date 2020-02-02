classdef (Abstract) AbstractApply < handle & matlab.mixin.Copyable & mlnest.IApply 
	%% ABSTRACTAPPLY implements application ideas from "Data Analysis:  A Bayesian Tutorial, Second Edition"
    %  by D.S. Sivia and J. Skilling, section 9.3.2  

	%  $Revision$ 
 	%  was created $Date$ 
 	%  by $Author$,  
 	%  last modified $LastChangedDate$ 
 	%  and checked into repository $URL$,  
 	%  developed on Matlab 8.4.0.150421 (R2014b) 
 	%  $Id$     
    
    properties (Constant)
        logSqrt2pi   = 0.9189385332046727;
    end
    
    properties
        Object
        sigma0
    end
    
    properties (Dependent)
        results
    end
    
    methods (Static)
        function main = run(application)
            assert(isa(application, 'mlnest.IApply'))
            main = mlnest.NestedSamplingMain(application);
%            save(sprintf('%s_run_%s.mat', ...
%                strrep(class(application), '.', '_'), datestr(now, 'yyyymmddHHMMSS')));
        end
    end
    
    methods
        
        %% GET
        
        function g = get.results(this)
            g = this.results_;
        end
        
        %%
        
        function [Obj,acceptRejectRatio] = Explore(this, Obj, logLstar)
            %% EXPLORE evolves sampling object within likelihood constraint
            %  Usage:  obj = this.Explore(Obj, log_likelihood_star)
            %                             ^ objects being evolved
            %                                  ^ likelihood constraint L > Lstar
            
            step   = this.STEP_Initial;              
            accept = 0;                 % # MCMC acceptances
            reject = 0;                 % # MCMC rejections
            Try    = Obj;
            for m = this.MCMC_Counter:-1:1

                %% Trial object
                flds = fields(Try);
                flds = this.ignoreFields(flds);
                for f = 1:length(flds)                    
                    Try.(flds{f}) = Obj.(flds{f}) + step*(2*rand() - 1.);  % |move| < step
                    Try.(flds{f}) = Try.(flds{f}) - floor(Try.(flds{f}));  % wraparound to stay within (0,1)
                end
                Try.logL = this.logLhood(Try); % trial likelihood value
                
                %% Accept if trial likelihood > previous likelihood * UNIFORM(0,1); evaluate logs.
                %  Cf. Sivia sec. 9.4.4.
                if (Try.logL > logLstar + log(rand()))
                    Obj = Try;
                    accept = accept + 1;  
                else
                    reject = reject + 1;
                end
       
                %% Refine step-size to let acceptance ratio converge around 50%
                if (accept > reject); step = step*exp(1/accept); end
                if (accept < reject); step = step/exp(1/reject); end
            end            
            acceptRejectRatio = accept/reject;
        end
        function logL = logLhood(this, Obj)
            Estimation = this.Estimation(Obj);
            positive   = this.Measurement > 0;
            EoverM     = Estimation(positive)./this.Measurement(positive);
            Q          = sum((1 - EoverM).^2);
            logL       = -0.5*Q/this.sigma0^2;
        end
        function Obj  = Prior(this, ~)
            Obj = struct('logL',  [], 'logWt', []);
            keys = this.map.keys;
            for k = 1:length(keys)
                Obj.(keys{k}) = this.priorValue(this.map(keys{k}), keys{k});
            end
            Obj.logL = this.logLhood(Obj);
        end
        function val  = priorValue(this, mapStruct, key)
            if lstrfind(this.ignoredObjFields, key)
                val = mapStruct.init;
                return
            end
            
            init_ = this.vec2uniform(mapStruct.init, this.limits(key));
            min_  = this.vec2uniform(mapStruct.min,  this.limits(key));
            max_  = this.vec2uniform(mapStruct.max,  this.limits(key));
            std_  = 0.25*abs(max_ - min_);
            val   = init_ + randn()*std_;
            while val < min_ || max_ < val % from Josh and Larry
                val = init_ + randn()*std_;
            end
            
            % surprisingly inferior:
            % val = val - floor(val);
            % val = val - (val - mod(val - min_, max_ - min_)) + min_;
        end
        function [r,this] = Results(this, Samples, nest, logZ)
            %% RESULTS prints the posterior properties; here mean and stddev of x, y
            %  @return struct with fields:  
            %          chains, ObjMoment1, ObjMoment2:  in reduced uniform variables;
            %          moment1, moment2:  in native variable units;
            %          flds:  fields of Obj.
            %  @return this with updated this.results.
            %
            %  Usage:  this.Results(Samples, nest, logZ)
            %                       ^ objects defining posterior
            %                                ^ # samples 
            %                                      ^ evidence (= total weight = SUM[samples] weight)
            
            flds = fields(Samples{1});
            flds = this.ignoreFields(flds);
            Obj1 = structfun(@(x) 0, Samples{1}, 'UniformOutput', false); 
            Obj2 = Obj1;
            chains = zeros(nest, length(flds));
            
            for ni = 1:nest
                
                w  = exp(Samples{ni}.logWt - logZ); % proportional weight
                
                for f = 1:length(flds)                    
                    chains(ni, f) = Samples{ni}.(flds{f});                                          
                    Obj1.(flds{f}) = Obj1.(flds{f}) + w*Samples{ni}.(flds{f});
                    Obj2.(flds{f}) = Obj2.(flds{f}) + w*(Samples{ni}.(flds{f}))^2;
                end
            end
            
            %% gather
            
            r.flds = flds;            
            r.Obj1 = Obj1;
            r.Obj2 = Obj2;
            r.chains = chains;
            r.moment1 = this.Obj2vec(this.Obj2native(Obj1));
            r.moment2 = this.Obj2vec(this.Obj2native(Obj2));
            this.results_ = r;
        end
        
        %% UTILITY
        
        function flds = ignoreFields(this, flds)
            for ig = this.ignoredObjFields
                flds = flds(~strcmp(flds, ig{1}));
            end
        end
        function vec = limits(this, key)
            vec = [this.map(key).min this.map(key).max];
        end
        function o = Obj2native(this, o)
            %% rescale to native units
            
            flds = fields(o);
            flds = this.ignoreFields(flds);
            for f = flds'             
                lims = this.limits(f{1});
                o.(f{1}) = lims(2)*o.(f{1}) + lims(1)*(1 - o.(f{1}));
            end
        end
        function o = Obj2uniform(this, o)
            %% rescale to (0,1)
            
            flds = fields(o);
            flds = this.ignoreFields(flds);
            for f = flds'
                lims = this.limits(f{1});
                u = (o.(f{1}) - lims(1))/(lims(2) - lims(1));
                o.(f{1}) = u;
            end
        end
        function vec = Obj2vec(this, obj)
            flds = fields(obj);
            flds = this.ignoreFields(flds);
            vec = zeros(1, length(flds));
            for idx = 1:length(vec)
                vec(idx) = obj.(flds{idx});
            end
        end
        function h = plotResults(this)
            figure;
            est = this.Estimation(this.results.Obj1);
            h = plot(1:length(this.Measurement), this.Measurement, 'o', ...
                     1:length(est), est, '-+');
            title([class(this) '.plotResults()'])
            legend('measurement', 'estimation')
        end
    end
    
    %% PROTECTED
    
    properties (Access = protected)
        results_
    end
    
    methods (Static, Access = protected)
        function z = PLUS(x, y)
            %% logarithmic addition log(exp(x) + exp(y))
            %  protects against over/underflow of exponential quantities such as likelihood \mathcal{L}^*
            %  by storing logarithms
            
            if (x > y)
                z = x + log(1 + exp(y - x));
            else
                z = y + log(1 + exp(x - y));
            end
        end
        function u = UNIFORM
            %% uniform inside (0,1)
           
            u = rand;
        end
        function y = vec2native(u, lims)
            y = lims(2)*u + lims(1)*(1 - u);
        end
        function u = vec2uniform(y, lims)
            u = (y - lims(1))/(lims(2) - lims(1));
        end
    end 
    
    %% HIDDEN & DEPRECATED
    
    methods (Hidden)
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
    
end

