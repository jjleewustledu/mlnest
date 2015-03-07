classdef NestedSamplingMain < mlnest.AbstractApply
	%% NESTEDSAMPLINGMAIN implements main.c from "Data Analysis:  A BAyesian Tutorial, Second Edition"
    %  by D.S. Sivia and J. Skilling, section 9.2.4
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
    end 
    
    properties (Dependent)
        n
        MAX
    end
    
    methods %% GET/SET
        function g = get.n(this)
            g = this.apply.n;
        end
        function g = get.MAX(this)
            g = this.apply.MAX;
        end
    end

	methods 
  		function this = NestedSamplingMain(app) 
 			%% NESTEDSAMPLINGMAIN 
 			%  Usage:  this = NestedSamplingMain(application_object)
            %                                    ^ mlnest.IAppy

            this.apply = app;
            Obj        = cell(1, this.n);
            Samples    = cell(1, this.MAX);
            for i = 1:this.n
                Obj{i} = this.Prior(Obj{i}); end
            
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
                
                %% kill worst object in favour of copy of different survivor
                copy = ceil(mod(this.n * this.UNIFORM, this.n));
                while (copy == worst && this.n > 1)
                    copy = ceil(mod(this.n * this.UNIFORM, this.n)); end
                logLstar = Obj{worst}.logL;
                Obj{worst} = Obj{copy};
                
                %% evolve copied object within constraints
                Obj{worst} = this.Explore(Obj{worst}, logLstar);
                
                %% shrink interval                
                logwidth = logwidth - 1/this.n;
            end
            this.printResults(Samples, nest, logZ, H);
        end 
        function obj = Prior(this, obj)
            obj = this.apply.Prior(obj);
        end
        function obj = Explore(this, obj, logLstar)
            obj = this.apply.Explore(obj, logLstar);
        end
        function       Results(this, smpls, nest, logZ)
            this.apply.Results(smpls, nest, logZ);
        end
        function       printResults(this, smpls, nest, logZ, H)
            fprintf('# iterates = %i\n', nest);
            fprintf('Evidence: ln(Z) = %g +- %g\n', logZ, sqrt(H/this.n));
            fprintf('Information:  H = %g nats = %g bits\n', H, H/log(2));
            this.Results(smpls, nest, logZ);
        end
    end 
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

