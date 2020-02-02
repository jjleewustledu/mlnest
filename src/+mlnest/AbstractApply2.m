classdef (Abstract) AbstractApply2 < handle & mlnest.AbstractApply
	%% ABSTRACTAPPLY2 implements application ideas from "Data Analysis:  A Bayesian Tutorial, Second Edition"
    %  by D.S. Sivia and J. Skilling, section 9.3.2  

	%  $Revision$ 
 	%  was created $Date$ 
 	%  by $Author$,  
 	%  last modified $LastChangedDate$ 
 	%  and checked into repository $URL$,  
 	%  developed on Matlab 8.4.0.150421 (R2014b) 
 	%  $Id$     
    
    methods (Static)
        function main = run(application)
            assert(isa(application, 'mlnest.IApply'))
            main = mlnest.NestedSamplingMain2(application);
%            save(sprintf('%s_run_%s.mat', ...
%                strrep(class(application), '.', '_'), datestr(now, 'yyyymmddHHMMSS')));
        end
    end
    
    methods
        function [Obj,acceptRejectRatio] = Explore(this, Obj, logLstar)
            %% EXPLORE evolves sampling object within likelihood constraint
            %  Usage:  obj = this.Explore(Obj, log_likelihood_star)
            %                             ^ object to evolve
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
                    min_ = this.map(flds{f}).min;
                    max_ = this.map(flds{f}).max;
                    Try.(flds{f}) = Obj.(flds{f}) + step*(max_ - min_)*(2*rand() - 1.);  % |move| < step
                    val = Try.(flds{f});
                    Try.(flds{f}) = val - (val - mod(val - min_, max_ - min_)) + min_; % wraparound to stay within (min_, max_)
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
        function val  = priorValue(this, mapStruct, key)            
            if lstrfind(this.ignoredObjFields, key)
                val = mapStruct.init;
                return
            end
            
            std_  = 0.25*abs(mapStruct.max - mapStruct.min);
            val   = mapStruct.init + randn()*std_;
            while val < mapStruct.min || mapStruct.max < val
                val = mapStruct.init + randn()*std_;
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
            r.moment1 = this.Obj2vec(Obj1);
            r.moment2 = this.Obj2vec(Obj2);
            this.results_ = r;
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
    
end

