classdef HandleNestedSamplingMain < handle & matlab.mixin.Copyable
	%% HANDLENESTEDSAMPLINGMAIN implements main.c from "Data Analysis:  A Bayesian Tutorial, Second Edition"
    %  by D.S. Sivia and J. Skilling, section 9.2.4.  
    %  (GNU General Public License software, (C) Sivia and Skilling 2006)

	%  $Revision$ 
 	%  was created $Date$ 
 	%  by $Author$,  
 	%  last modified $LastChangedDate$ 
 	%  and checked into repository $URL$,  
 	%  developed on Matlab 8.4.0.150421 (R2014b) 
 	%  $Id$ 

	properties 
        nCheckStoppingCondition = 5
        nReports = 10
    end 
    
    properties (Dependent)
        apply
        map
        MAX          % # of nested sampling loops, similar to temperature for s.a.
        MCMC_Counter % counter for explorting single particle (pre-judged # steps)
        Measurement  % external data
        n            % # of sampling particles \sim (log width of outer prior mass)^{-1}
        results
        sigma0       % of model estimation
        STEP_Initial % Initial guess suitable step-size in (0,1)
    end
    
    methods %% GET/SET
        function g = get.apply(this)
            g = this.apply_;
        end
        function m = get.map(this)
            m = this.apply.map;
        end
        function m = get.MAX(this)
            m = this.apply.MAX;
        end
        function m = get.MCMC_Counter(this)
            m = this.apply.MCMC_Counter;
        end
        function m = get.Measurement(this)
            m = this.apply.Measurement;
        end
        function m = get.n(this)
            m = this.apply.n;
        end
        function m = get.results(this)
            m = this.apply.results;
        end
        function g = get.STEP_Initial(this)
            g = this.apply.STEP_Initial;
        end
        function g = get.sigma0(this)
            g = this.apply.sigma0;
        end
    end

	methods 
        function est  = Estimation(this, varargin)
            est = this.apply.Estimation(varargin{:});
        end
        function logL = logLhood(this, varargin)
            logL = this.apply.logLhood(varargin{:});
        end
        function Obj  = Prior(this)
            Obj = this.apply.Prior;
        end
        function [obj,acceptRejectRatio] = Explore(this, obj, logLstar)
            [obj,acceptRejectRatio] = this.apply.Explore(obj, logLstar);
        end   
        
        %% UTILITIES
        
        function        fprintfModel(this)
            this.apply.fprintfModel();
        end        
        function        printResults(this, nest, logZ, H, Samples)
            this.apply.Results(Samples, nest, logZ);
            fprintf('-----------------------------------------------------------\n');
            fprintf('# iterates ~ nest = %i <= MAX = %i\n', nest, this.MAX);
            fprintf('# sampling particles ~ n = %f\n', this.n);
            fprintf('MCMC_Counter = %f\n', this.MCMC_Counter);
            fprintf('STEP_Initial = %f\n', this.STEP_Initial);
            fprintf('Stopping criteria = %f\n', nest/(H * this.n))
            fprintf('Evidence:  ln(Z) = %f +/- %f\n', logZ, sqrt(H/this.n));
            fprintf('Information:  H = %f nats = %f bits\n', H, H/log(2));
            this.fprintfModel()
            fprintf('\t%s(k = %i) = %f\n', 'sampled logL ', nest, Samples{nest}.logL);
            fprintf('\t%s(k = %i) = %f\n', 'sampled logWt', nest, Samples{nest}.logWt);
            fprintf('\n')
        end
        function h    = plotMap(this)
            h = this.apply.plotMap();
        end
        function h    = plotMatrix(this)
            h = this.apply.plotMatrix();
        end
        function h    = plotObjs(this, varargin)
            h = this.apply.plotObjs(varargin{:});
        end
        function h    = plotResults(this)
            h = this.apply.plotResults();
        end
        function s =    sprintfModel(this)
            s = this.apply.sprintfModel();
        end
        
        %%
        
  		function this = HandleNestedSamplingMain(app) 
 			%% HANDLENESTEDSAMPLINGMAIN 
            %  @param app is an application object implementing mlnest.IHandleApply
            %  @return this
            %
 			%  Usage:  this = HandleNestedSamplingMain(mlnest.IHandleApply object)
            %  Internally:
            %     Object Obj[n];        // Collection of n objects
            %     Object Samples[MAX];  // Objects stored for posterior results
            %     double logwidth;      // ln(width in prior mass)
            %     double logLstar;      // ln(Likelihood constraint)
            %     double H    = 0.0;    // Information, initially 0
            %     double logZ =-DBL_MAX;// ln(Evidence Z, initially 0)
            %     double logZnew;       // Updated logZ
            %     int    i;             // Object counter
            %     int    copy;          // Duplicated object
            %     int    worst;         // Worst object
            %     int    nest;          // Nested sampling iteration count
            
            this.apply_ = app;
            Obj         = cell(1, this.n);
            Samples     = cell(1, this.MAX);
            H           = 0;
            logZ        = -realmax;
            
            %% set prior objects
            for i = 1:this.n
                Obj{i} = this.Prior; 
            end
            
            %% outermost interval of prior mass
            logwidth = log(1 - exp(-1/this.n)); 
            try
                
                %% NESTED SAMPLING LOOP
                for nest = 1:this.MAX

                    %% worst object in collection with weight = width * likelihood
                    worst = 1;
                    for i = 2:this.n
                        if (Obj{i}.logL < Obj{worst}.logL); worst = i; end
                    end
                    Obj{worst}.logWt = logwidth + Obj{worst}.logL;

                    %% update evidence Z and information H
                    logZnew = this.PLUS(logZ, Obj{worst}.logWt);
                    H = exp(Obj{worst}.logWt - logZnew) * Obj{worst}.logL + ...
                        exp(logZ - logZnew) * (H + logZ) - ...
                        logZnew;
                    logZ = logZnew;

                    %% posterior samples (optional, for reporting)
                    Samples{nest} = Obj{worst};

                    %% kill worst object in favour of copy of different survivor
                    copy = mod(floor(this.n * rand()), this.n) + 1;
                    while (copy == worst && this.n > 1)
                        copy = mod(floor(this.n * rand()), this.n) + 1; 
                    end
                    logLstar = Obj{worst}.logL;
                    Obj{worst} = Obj{copy};

                    %% evolve copied object within constraints
                    [Obj{worst},acceptRejectRatio] = this.Explore(Obj{worst}, logLstar); %#ok<ASGLU>

                    %% shrink interval                
                    logwidth = logwidth - 1/this.n;
                    
                end %% NESTED SAMPLING LOOP (might be ok to terminate early)
            catch ME
                handexcept(ME);
            end
            
            %% finalize with evidence Z, information H, and optional posterior Samples
            this.finalize(nest, logZ, H, Samples);            
        end 
    end 
    
    %% PROTECTED  
    
    properties (Access = protected)        
        apply_
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
    end 
    
    methods (Access = 'protected')    
        function tf = stoppingConditionMet(this, nest, H)
            tf = nest/(H * this.n) > 4;
        end
        function      finalize(this, nest, logZ, H, Samples)
            this.apply_.Results(Samples, nest, logZ);
            this.printResults(nest, logZ, H, Samples);
            this.plotResults();
            this.plotMatrix();
        end
        function printAcceptsRejects(this, nest, acceptRejectRatio)
            if (1        == nest || ...
                this.MAX == nest || ...    
                0        == mod(nest, floor(this.MAX/this.nReports^2)))
                fprintf('# accepts/# rejects = %f\n', acceptRejectRatio);
            end
        end
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

