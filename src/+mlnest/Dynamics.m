classdef Dynamics < mlnest.AbstractApply
	%% DYNAMICS  

	%  $Revision$
 	%  was created 11-Apr-2020 20:20:51 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlnest/src/+mlnest.
 	%% It was developed on Matlab 9.7.0.1319299 (R2019b) Update 5 for MACI64.  Copyright 2020 John Joowon Lee.
 	
	properties
        background = 5       % Bq/mL
        context              % client of mlnest.Dynamics
        ignoredObjFields = {'logL' 'logWt'}
        map                  % containers.Map containing model params as structs with fields:  min, max ,init
        MAX = 1000           % # of nested sampling loops, similar to temperature for s.a.
        MCMC_Counter = 50    % counter for explorting single particle (pre-judged # steps); nonlinearly affects precision
        model                % 
        Measurement          % external data
        n = 10               % # of sampling particles \sim (log width of outer prior mass)^{-1}; reduces sampling space
        STEP_Initial         % Initial guess suitable step-size in (0,1); 0.01*MCMC_Counter^{-1} < STEP < 10*MCMC_Counter^{-1} improve precision        
        times_sampled        % numeric times for Measurement; midpoints of frames for PET        
        visualize = true
 	end
    
    methods (Static)
        function conc = slide_fast(conc, Dt)
            %% SLIDE_FAST slides discretized function conc(t) to conc(t - Dt);
            %  @param conc is row vector without NaN.
            %  @param t is row vector with same size as conc.
            %  @param Dt is scalar rounded to integer.
            %
            %  Dt > 0 will slide conc(t) towards later times t.
            %  Dt < 0 will slide conc(t) towards earlier times t.
            
            Dt = round(Dt);
            if Dt == 0
                return
            end
            if Dt < 0
                T = length(conc);
               conc_ = conc(end)*ones(1, length(conc));
               conc_(1:T+Dt) = conc(1-Dt:end);
               conc = conc_;
               return
            end
            conc_ = zeros(size(conc));
            conc_(1+Dt:end) = conc(1:end-Dt);
            conc = conc_;
        end
    end

	methods
 		function this = Dynamics(varargin)
 			%% DYNAMICS; see also mlnest.IApply for parameter specifications.
            %  @param context is mlglucose.Huang1980.
            %  @param MAX.
            %  @param MCMC_Counter.
            %  @param n.
            %  @param STEP_Initial.
            %  @param sigma0.
            
 			this = this@mlnest.AbstractApply(varargin{:});
            
            ip = inputParser;            
            ip.KeepUnmatched = true;
            addParameter(ip, 'context', [], @(x) isa(x, 'mlglucose.Huang1980'))
            addParameter(ip, 'MAX', this.MAX, @isnumeric)
            addParameter(ip, 'MCMC_Counter', this.MCMC_Counter, @isnumeric)
            addParameter(ip, 'n', this.n, @isnumeric)
            addParameter(ip, 'STEP_Initial', [], @isnumeric)
            addParameter(ip, 'sigma0', this.sigma0, @(x) isscalar(x) && x <= 1)
            parse(ip, varargin{:})
            ipr = ip.Results;

            this.context = ipr.context;            
            this.model = this.context.model;               % copy objects for speed
            this.map = this.model.map;                     %
            this.Measurement = this.context.measurement;   %
            this.times_sampled = this.model.times_sampled; %           
            this.MAX = ipr.MAX;
            this.MCMC_Counter = ipr.MCMC_Counter;
            this.n = ipr.n;
            if isempty(ipr.STEP_Initial)
                this.STEP_Initial = 0.1/this.MCMC_Counter;
            end
            this.sigma0 = ipr.sigma0;
        end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

