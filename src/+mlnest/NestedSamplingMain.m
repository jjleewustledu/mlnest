classdef NestedSamplingMain < mlnest.AbstractApply
	%% NESTEDSAMPLINGMAIN implements main.c from "Data Analysis:  A Bayesian Tutorial, Second Edition"
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
        apply
        nReports = 10
        results
    end 
    
    properties (Dependent)
        map
        MAX
        n
    end
    
    methods %% GET/SET
        function m = get.map(this)
            m = this.apply.map;
        end
        function m = get.MAX(this)
            m = this.apply.MAX;
        end
        function m = get.n(this)
            m = this.apply.n;
        end
    end

	methods 
  		function this = NestedSamplingMain(app) 
 			%% NESTEDSAMPLINGMAIN 
            %  @param app is an application object implementing mlnest.IApply
            %  @return this
 			%  Usage:  this = NestedSamplingMain(mlnest.IApply object)
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
            
            this.apply = app;
            Obj        = cell(1, this.n);
            Samples    = cell(1, this.MAX);
            H          = 0;
            logZ       = -realmax;
            
            %% set prior objects
            for i = 1:this.n
                Obj{i} = this.Prior; end
            
            %% outermost interval of prior mass
            logwidth = log(1 - exp(-1/this.n)); 
            try
                
                %% NESTED SAMPLING LOOP
                for nest = 1:this.MAX

                    %% worst object in collection with weight = width * likelihood
                    worst = 1;
                    for i = 2:this.n
                        if (Obj{i}.logL < Obj{worst}.logL); worst = i; end; end
                    Obj{worst}.logWt = logwidth + Obj{worst}.logL;

                    %% update evidence Z and information H
                    logZnew = this.PLUS(logZ, Obj{worst}.logWt);
                    H = exp(Obj{worst}.logWt - logZnew) * Obj{worst}.logL + ...
                        exp(logZ - logZnew) * (H + logZ) - ...
                        logZnew;
                    logZ = logZnew;

                    %% posterior samples (optional)
                    Samples{nest} = Obj{worst};
                    this = this.printResults(nest, logZ, H, Samples);

                    %% kill worst object in favour of copy of different survivor
                    copy = mod(floor(this.n * this.UNIFORM), this.n) + 1;
                    while (copy == worst && this.n > 1)
                        copy = mod(floor(this.n * this.UNIFORM), this.n) + 1; end
                    logLstar = Obj{worst}.logL;
                    Obj{worst} = Obj{copy};

                    %% evolve copied object within constraints
                    [Obj{worst},acceptRejectRatio] = this.Explore(Obj{worst}, logLstar);
                    this.printAcceptsRejects(nest, acceptRejectRatio);

                    %% shrink interval                
                    logwidth = logwidth - 1/this.n;
                end %% NESTED SAMPLING LOOP (might be ok to terminate early)
            catch ME
                handexcept(ME);
            end
            
            %% exit with evidence Z, information H, and optional posterior Samples
            this = this.exit(nest, logZ, H, Samples);            
        end 
        function Obj = Prior(this)
            Obj = this.apply.Prior;
        end
        function [obj,acceptRejectRatio] = Explore(this, obj, logLstar)
            [obj,acceptRejectRatio] = this.apply.Explore(obj, logLstar);
        end
        function dat = Results(this, smpls, nest, logZ)
            dat = this.apply.Results(smpls, nest, logZ);
        end
    end 
    
    %% PRIVATE
    
    methods (Access = 'private')        
        function this = printResults(this, nest, logZ, H, Samples)
            if (1        == nest || ...
                this.MAX == nest || ...    
                0        == mod(nest, floor(this.MAX/this.nReports)))
                fprintf('-----------------------------------------------------------\n');
                fprintf('# iterates = %i\n', nest);
                fprintf('Evidence:  ln(Z) = %g +- %g\n', logZ, sqrt(H/this.n));
                fprintf('Information:  H = %g nats = %g bits\n', H, H/log(2));
                this.results = this.Results(Samples, nest, logZ);
                
                for f = 1:length(this.results.flds)
                    if (~strcmp('logL', this.results.flds{f}) && ~strcmp('logWt', this.results.flds{f}))
                        fprintf('\t%s = %g +/- %g\n', ...
                            this.results.flds{f}, this.results.moment1(f), ...
                            sqrt(this.results.moment2(f) - this.results.moment1(f)^2));
                    end
                end
                fprintf('%s(k = %i) = %g\n', 'logL',  nest, Samples{nest}.logL);
                fprintf('%s(k = %i) = %g\n', 'logWt', nest, Samples{nest}.logWt);                
                this.checkStoppingCondition(nest,H);
            end
        end 
        function printAcceptsRejects(this, nest, acceptRejectRatio)        
            if (1        == nest || ...
                this.MAX == nest || ...    
                0        == mod(nest, floor(this.MAX/this.nReports^2)))
                fprintf('# accepts/# rejects = %g\n', acceptRejectRatio);
            end
        end
        function tf = checkStoppingCondition(this, nest, H)
            tf = nest/(H * this.n) > 4;
            if (tf)                
                error('mlnest:stoppingCondition', ...
                      'NestedSamplingMain.checkStoppingCondition: stopping at nest->%g, H->%g\n', nest, H);
            end
        end
        function this = exit(this, nest, logZ, H, Samples)
            fprintf('-----------------------------------------------------------\n');
            fprintf('# iterates = %i\n', nest);
            fprintf('Evidence:  ln(Z) = %g +- %g\n', logZ, sqrt(H/this.n));
            fprintf('Information:  H = %g nats = %g bits\n', H, H/log(2));
            this.results = this.Results(Samples, nest, logZ);
        end
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

