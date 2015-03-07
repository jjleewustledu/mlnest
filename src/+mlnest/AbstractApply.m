classdef AbstractApply < mlnest.IApply 
	%% ABSTRACTAPPLY   

	%  $Revision$ 
 	%  was created $Date$ 
 	%  by $Author$,  
 	%  last modified $LastChangedDate$ 
 	%  and checked into repository $URL$,  
 	%  developed on Matlab 8.4.0.150421 (R2014b) 
 	%  $Id$     

    properties (Constant)
        logSqrt2pi = 0.9189385332046727;
    end
    
    methods        
        function u    = uniform2limits(~, u, lims)
            u = lims(2)*u + lims(1)*(1 - u);
        end        
        function Obj  = Explore(this, Obj, logLstar)
            %% EXPLORE evolves object within likelihood constraint
            %  Usage:  obj = this.Explore(obj, log_likelihood_star)
            
            step   = 0.1; % Initial guess suitable step-size in (0,1)
            m      = 20;  % MCMC counter (pre-judged # steps)
            accept = 0;   % # MCMC acceptances
            reject = 0;   % # MCMC rejections
            Try = Obj;
            while (m > 0)
                m = m - 1;

                %% Trial object
                flds = fields(Try);
                for f = 1:length(flds)
                    Try.(flds{f}) = Obj.(flds{f}) + step*(2*this.UNIFORM - 1.);  % |move| < step
                    Try.(flds{f}) = Try.(flds{f}) - floor(Try.(flds{f}));        % wraparound to stay within (0,1)
                end                
                Try.logL = this.logLhood(Try); % trial likelihood value
                
                %% Accept if and only if within hard likelihood constraint
                if (Try.logL > logLstar)
                    Obj = Try;
                    accept = accept + 1;  
                else
                    reject = reject + 1;
                end
       
                %% Refine step-size to let acceptance ratio converge around 50%
                if (accept > reject); step = step*exp(1/accept); end
                if (accept < reject); step = step/exp(1/reject); end
            end
        end
        function        Results(this, Samples, nest, logZ)
            %% RESULTS prints the posterior properties; here mean and stddev of x, y
            %  Usage:  this.Results(Samples, nest, logZ)
            %                       ^ objects defining posterior
            %                                ^ # samples 
            %                                      ^ evidence (= total weight = SUM[samples] weight)

            flds    = fields(Samples{1});
            moment1 = zeros(1, length(flds));
            moment2 = zeros(1, length(flds));
            for i = 1:nest
                w  = exp(Samples{i}.logWt - logZ); % proportional weight
                for f = 1:length(flds)
                    if (~strcmp('logL', flds{f}) && ~strcmp('logWt', flds{f}))
                        fvalue     = this.uniform2limits(Samples{i}.(flds{f}), this.limits.(flds{f}));
                        moment1(f) = moment1(f) + w*fvalue;
                        moment2(f) = moment2(f) + w*fvalue^2;
                    end
                end
            end
                fprintf('Estimated parameter = mean +/- stddev:\n');
                for f = 1:length(flds)
                    fprintf('%s = %g +/- %g\n', flds{f}, moment1(f), sqrt(moment2(f) - moment1(f)^2));
                end
        end
    end
    
    methods (Static, Access = 'protected')
        function u = UNIFORM
            %% uniform inside (0,1)
           
            u = rand;
        end
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
    end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

