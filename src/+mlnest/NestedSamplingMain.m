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
  		function [this,Samples] = NestedSamplingMain(app) 
 			%% NESTEDSAMPLINGMAIN 
 			%  Usage:  [this,Samples] = NestedSamplingMain(application_object)
            %                                              ^ mlnest.IAppy
            %                ^ large, working space of data

            this.apply = app;
            Obj        = cell(1, this.n);
            Samples    = cell(1, this.MAX);
            for i = 1:this.n
                Obj{i} = this.Prior; end
            
            H        = 0;
            logZ     = -realmax;
            logwidth = log(1 - exp(-1/this.n)); % outermost interval of prior mass
            for nest = 1:this.MAX
               
                %% worst object in collection with weight = width * likelihood
                worst = 1;
                for i = 2:this.n
                    if (Obj{i}.logL < Obj{worst}.logL)
                        worst = i; end
                end
                Obj{worst}.logWt = logwidth + Obj{worst}.logL;
                
                %% update evidence Z and information H
                logZnew = this.PLUS(logZ, Obj{worst}.logWt);
                H = exp(Obj{worst}.logWt - logZnew) * Obj{worst}.logL + ...
                         exp(logZ - logZnew) * (H + logZ) - ...
                         logZnew;
                logZ = logZnew;
                
                %% posterior samples
                Samples{nest} = Obj{worst};                
                this = this.printResults(Samples, nest, logZ, H, Samples{nest}.logL, Samples{nest}.logWt);
                if (this.stopCondition(nest,H)); break; end

                %% kill worst object in favour of copy of different survivor
                copy = mod(floor(this.n * this.UNIFORM), this.n) + 1;
                while (copy == worst && this.n > 1)
                    copy = mod(floor(this.n * this.UNIFORM), this.n) + 1; end
                logLstar = Obj{worst}.logL;
                Obj{worst} = Obj{copy};
                
                %% evolve copied object within constraints
                Obj{worst} = this.Explore(Obj{worst}, logLstar);
                
                %% shrink interval                
                logwidth = logwidth - 1/this.n;
            end
        end 
        function obj = Prior(this)
            obj = this.apply.Prior;
        end
        function obj = Explore(this, obj, logLstar)
            obj = this.apply.Explore(obj, logLstar);
        end
        function dat = Results(this, smpls, nest, logZ)
            dat = this.apply.Results(smpls, nest, logZ);
        end
    end 
    
    %% PRIVATE
    
    methods (Access = 'private')        
        function this = printResults(this, smpls, nest, logZ, H, logL, logWt)        
            if (1        == nest || ...
                this.MAX == nest || ...    
                0        == mod(nest, floor(this.MAX/this.nReports)))
                fprintf('-----------------------------------------------------------\n');
                fprintf('# iterates = %i\n', nest);
                fprintf('Evidence: ln(Z) = %g +- %g\n', logZ, sqrt(H/this.n));
                fprintf('Information:  H = %g nats = %g bits\n', H, H/log(2));                
                fprintf('Estimated parameter = mean +/- stddev:\n');
                
                this.results = this.Results(smpls, nest, logZ);
                
                for f = 1:length(this.results.flds)
                    if (~strcmp('logL', this.results.flds{f}) && ~strcmp('logWt', this.results.flds{f}))
                        fprintf('%s = %g +/- %g\n', ...
                            this.results.flds{f}, this.results.moment1(f), ...
                            sqrt(this.results.moment2(f) - this.results.moment1(f)^2));
                    end
                end
                fprintf('%s = %g\n', 'logL',  logL);
                fprintf('%s = %g\n', 'logWt', logWt);
            end
        end
        function tf = stopCondition(this, nest, H)
            tf = nest/(H * this.n) > 10;
            fprintf('NestedSamplingMain.stopCondition:  stopping at nest->%g, H->%g\n', nest, H);
        end
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

