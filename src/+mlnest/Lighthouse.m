classdef Lighthouse < mlnest.AbstractApply 
	%% LIGHTHOUSE implements lighthouse from "Data Analysis:  A BAyesian Tutorial, Second Edition"
    %  by D.S. Sivia and J. Skilling, section 9.3.2

    %  apply.c     "LIGHTHOUSE" NESTED SAMPLING APPLICATION
    %  (GNU General Public License software, (C) Sivia and Skilling 2006)
    %               u=0                                 u=1
    %                -------------------------------------
    %           y=2 |:::::::::::::::::::::::::::::::::::::| v=1
    %               |::::::::::::::::::::::LIGHT::::::::::|
    %          north|::::::::::::::::::::::HOUSE::::::::::|
    %               |:::::::::::::::::::::::::::::::::::::|
    %               |:::::::::::::::::::::::::::::::::::::|
    %           y=0 |:::::::::::::::::::::::::::::::::::::| v=0
    %  --*--------------*----*--------*-**--**--*-*-------------*--------
    %              x=-2          coastline -.east      x=2
    %  Problem:
    %   Lighthouse at (x,y) emitted n flashes observed at D[.] on coast.
    %  Inputs:
    %   Prior(u)    is uniform (=1) over (0,1), mapped to x = 4*u - 2; and
    %   Prior(v)    is uniform (=1) over (0,1), mapped to y = 2*v; so that
    %   Position    is 2-dimensional -2 < x < 2, 0 < y < 2 with flat prior
    %   Likelihood  is L(x,y) = PRODUCT[k] (y/pi) / ((D[k] - x)^2 + y^2)
    %  Outputs:
    %   Evidence    is Z = INTEGRAL L(x,y) Prior(x,y) dxdy
    %   Posterior   is P(x,y) = L(x,y) / Z estimating lighthouse position
    %   Information is H = INTEGRAL P(x,y) log(P(x,y)/Prior(x,y)) dxdy

    %  Results from original C source codes:
    %  # iterates = 1000
    %  Evidence: ln(Z) = -160.557 +- 0.176466
    %  Information: H = 3.11403 nats = 4.49259 bits
    %  mean(x) = 1.27463, stddev(x) = 0.173153
    %  mean(y) = 0.990962, stddev(y) = 0.178608
    
    %  For D = D + 1:    
    %  # iterates = 1000
    %  Evidence: ln(Z) = -162.927 +- 0.195109
    %  Information:  H = 3.80677 nats = 5.492 bits
    %  mean(x) = 1.91361, stddev(x) = 0.0753362
    %  mean(y) = 1.09591, stddev(y) = 0.21325

	%  $Revision$ 
 	%  was created $Date$ 
 	%  by $Author$,  
 	%  last modified $LastChangedDate$ 
 	%  and checked into repository $URL$,  
 	%  developed on Matlab 8.4.0.150421 (R2014b) 
 	%  $Id$  
    
    properties (Dependent)
        map
    end

	properties 
        ignoredObjFields
         MAX = 1000
         MCMC_Counter
         Measurement
 		 n = 100
         STEP_Initial
    end 
    
    methods (Static)
        function this = run
            import mlnest.*;
            this = NestedSamplingMain(Lighthouse);
        end
    end

	methods 
        
        %% GET
        
        function m   = get.map(~)
            m = containers.Map;
            m('u') = struct('min', 0, 'init', 0.5, 'max', 1);
            m('v') = struct('min', 0, 'init', 0.5, 'max', 1);
            m('x') = struct('min',-2, 'init', 0,   'max', 2);
            m('y') = struct('min', 0, 'init', 1,   'max', 2);
        end
        
        %%
        
        function logL = logLhood(this, x, y)
            %% LOGLHOOD
            %  Usage:  log_likelihood = this.logLhood(x, y)
            %                                         ^ easterly position
            %                                            ^ northerly position
            
            D = this.Measurement;
            logL = 0; % logLikelihood accumulator
            for k = 1:length(D)
                logL = logL + log((y/3.1416) / ((D(k)-x)*(D(k)-x) + y*y)); end
        end
        function Obj  = Prior(this, ~)
            Obj = struct( ...
                'u', rand(), ...
                'v', rand(), ...
                'x', [], ...
                'y', [], ...
                'logL', [], ...
                'logWt', []);
            Obj.x    = 4*Obj.u - 2;
            Obj.y    = 2*Obj.v;
            Obj.logL = this.logLhood(Obj.x, Obj.y);
        end
        function [Obj,acceptRejectRatio] = Explore(this, Obj, logLstar)
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
                Try.u = Obj.u + step * (2*rand() - 1.);  % |move| < step
                Try.v = Obj.v + step * (2*rand() - 1.);  % |move| < step
                Try.u = Try.u - floor(Try.u);      % wraparound to stay within (0,1)
                Try.v = Try.v - floor(Try.v);      % wraparound to stay within (0,1)
                Try.x = 4.0 * Try.u - 2.0;         % map to x
                Try.y = 2.0 * Try.v;               % map to y
                Try.logL = this.logLhood(Try.x, Try.y); % trial likelihood value
                
                %% Accept if and only if within hard likelihood constraint
                if( Try.logL > logLstar )
                    Obj = Try;
                    accept = accept + 1;  
                else
                    reject = reject + 1;
                end
       
                %% Refine step-size to let acceptance ratio converge around 50%
                if ( accept > reject ); step = step * exp(1.0 / accept); end
                if ( accept < reject ); step = step / exp(1.0 / reject); end
                
                acceptRejectRatio = accept/reject;
            end
        end
        function [r,this] = Results(this, Samples, nest, logZ)
            %% RESULTS prints the posterior properties; here mean and stddev of x, y
            %  Usage:  this.Results(Samples, nest, logZ)
            %                       ^ objects defining posterior
            %                                ^ # samples 
            %                                      ^ evidence (= total weight = SUM[samples] weight)

            x = 0.0; xx = 0.0;   % 1st and 2nd moments of x
            y = 0.0; yy = 0.0;   % 1st and 2nd moments of y
            for i = 1:nest
                w  = exp(Samples{i}.logWt - logZ); % proportional weight
                x  = x  + w * Samples{i}.x;
                xx = xx + w * Samples{i}.x * Samples{i}.x;
                y  = y  + w * Samples{i}.y;
                yy = yy + w * Samples{i}.y * Samples{i}.y;
            end
            fprintf('mean(x) = %g, stddev(x) = %g\n', x, sqrt(xx-x*x));
            fprintf('mean(y) = %g, stddev(y) = %g\n', y, sqrt(yy-y*y));
            
            r.flds = fields(Samples{1});
            r.moment1 = [0 0 x y];
            r.moment2 = [0 0 xx yy];
            this.results_ = r;
        end
        
        function this = Lighthouse(varargin)            
            this.Measurement =       [ 4.73,  0.45, -1.73,  1.09,  2.19,  0.12, ...
                  1.31,  1.00,  1.32,  1.07,  0.86, -0.49, -2.59,  1.73,  2.11, ...
                  1.61,  4.98,  1.71,  2.23,-57.20,  0.96,  1.25, -1.56,  2.45, ...
                  1.19,  2.17,-10.66,  1.91, -4.16,  1.92,  0.10,  1.98, -2.51, ...
                  5.55, -0.47,  1.91,  0.95, -0.78, -0.84,  1.72, -0.01,  1.48, ...
                  2.70,  1.21,  4.41, -4.79,  1.33,  0.81,  0.20,  1.58,  1.29, ...
                 16.19,  2.75, -2.38, -1.79,  6.50,-18.53,  0.72,  0.94,  3.64, ...
                  1.94, -0.11,  1.57,  0.57]; % arrival positions
        end
    end 
    
    %% HIDDEN
    
    methods (Hidden)        
        function obj  = Estimation(~, obj)
        end
        function plotResults(~)
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

