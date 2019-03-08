classdef C11Glucose < mlnest.AbstractApply
	%% C11GLUCOSE implements application ideas from "Data Analysis:  A Bayesian Tutorial, Second Edition"
    %  by D.S. Sivia and J. Skilling, section 9.3.2

	%  $Revision$
 	%  was created 20-Oct-2016 21:18:43
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlnest/src/+mlnest.
 	%% It was developed on Matlab 9.0.0.341360 (R2016a) for MACI64.
 	
	 

	properties         
        n   = 100
        MAX = 10000
         
        k12 = 0.00525
        k21 = 0.0860
        k32 = 0.00344
        k43 = 0.000302
        t0  = 44.1
        sigma = nan
        
        pnumber
        scanIndex
        xLabel = 'times/s'
        yLabel = 'concentration/(wellcounts/mL)'
        
        dta
        Measurements % tracer concentration of single voxel over time, using dimensionless, scaled parameters
        region 
        times
    end 
    
    properties (Dependent)        
        baseTitle
        detailedTitle
        gluTxlsxInfo
        map
        mode
        VB
        FB
    end
    
    methods %% GET
        function bt  = get.baseTitle(this)
            bt = sprintf('%s %s', class(this),str2pnum(pwd));
        end
        function dt  = get.detailedTitle(this)
            dt = sprintf('%s:\nk21 %g, k12 %g, k32 %g, k43 %g, t0 %g, VB %g, FB %g', ...
                         this.baseTitle, ...
                         this.k21, this.k12, this.k32, this.k43, this.t0, this.VB, this.FB);
        end
        function inf = get.gluTxlsxInfo(this)
            inf = this.gluTxlsx_.pid_map(this.pnumber).(sprintf('scan%i', this.scanIndex));
        end       
        function m   = get.map(this)
            m = containers.Map;
            fL = 1; fH = 1;
            meanMeas = mean(this.Measurements);
            m('k12') = struct('fixed', 0, 'min', 0.00192*fL, 'mean', this.k12, 'max', 0.0204*fH);  % Powers' monkey paper
            m('k21') = struct('fixed', 0, 'min', 0.0435*fL,  'mean', this.k21, 'max', 0.0942*fH);  % "
            m('k32') = struct('fixed', 0, 'min', 0.0015*fL,  'mean', this.k32, 'max', 0.5589*fH);  % " excluding last 2 entries
            m('k43') = struct('fixed', 0, 'min', 2.03e-5*fL, 'mean', this.k43, 'max', 3.85e-4*fH); % "
            m('t0' ) = struct('fixed', 0, 'min', 0*fL,       'mean', this.t0,  'max', 5e2*fH);   
            m('sigma') = struct('fixed', 0, 'min', 0,        'mean', 0.5*meanMeas, 'max', meanMeas);  
        end
        function m   = get.mode(this)
            m = this.gluTxlsx_.mode;
        end
        function v   = get.VB(this)
            % fraction
            v = this.gluTxlsxInfo.cbv;
            if (v > 1)
                v = v/100; end
        end
        function f   = get.FB(this)
            % fraction/s
            f = this.gluTxlsxInfo.cbf;
            assert(~isnumeric(f) || ~isnan(f), 'mlnext:nan', 'C11Glucose.get.FB');
            f = 1.05 * f / 6000; % mL/min/100g to 1/s
        end
    end

    methods (Static)
        function [these,dt] = parRun

            t0 = tic %#ok<*NOPRT>
            subjectsPth = fullfile('/scratch/jjlee/powers', 'pet6_c11monkey', '');
            cd(subjectsPth);
            runLabel = fullfile(subjectsPth, sprintf('parRun_C11Glucose_%s', mydatetimestr(now)));
            %diary([runLabel '.log']);

            import mlnext.* mlsystem.*;   

            dt    = DirTool('M*'); dns = dt.dns;
            assert(~isempty(dt.dns));
            these = cell(length(dt.dns));
            
            parfor d = 1:length(dns)
                try
                    pth = fullfile(subjectsPth, dns{d}, '');
                    cd(pth);
                    fprintf('-------------------------------------------------------------------------------------------------------------------------------\n');
                    fprintf('parRun:  working in %s\n', pth);
                    these{d} = C11Glucose.run;
                catch ME
                    handwarning(ME)
                end
            end

            cd(subjectsPth);
            save([runLabel '.mat']);            
            tf = toc(t0) %#ok<NASGU>
            %diary off
        end
        function this = run
            import mlnest.*;
            this = NestedSamplingMain(C11Glucose);
        end
        function Q_sampl = concentrationQ(k04, k12, k21, k32, k43, t0, dta, VB, t_sampl)
            t                    = dta.timeInterpolants; % use interpolants internally            
            t0_idx               = floor(t0/dta.dt) + 1;
            cart                 = dta.wellCountInterpolants(end) * ones(1, length(t));
            cart(1:end-t0_idx+1) = dta.wellCountInterpolants(t0_idx:end); % shift cart earlier in time
            
            k22 = k12 + k32;
            q1_ = VB * cart;
            q2_ = VB * k21 * exp(-k22*t);
            q3_ = VB * k21 * k32 * (k22 - k43)^-1 * (exp(-k43*t) - exp(-k22*t));
            q4_ = VB * k21 * k32 * k43 * ( ...
                     exp(-k22*t)/((k04 - k22)*(k43 - k22)) + ...
                     exp(-k43*t)/((k22 - k43)*(k04 - k43)) + ...
                     exp(-k04*t)/((k22 - k04)*(k43 - k04)));
                 
            q234    = conv(q2_ + q3_ + q4_, cart);
            Q       = q1_ + q234(1:length(t)); % truncate convolution         
            Q_sampl = pchip(t, Q, t_sampl); % resample interpolants
        end
    end
    
	methods
        function this = C11Glucose(aTsc, aDta, varargin)
            
            ip = inputParser;
            addRequired(ip, 'aTsc', @(x) isa(x, 'mlpet.TSC'));
            addRequired(ip, 'aDta', @(x) isa(x, 'mlpet.DTA'));
            parse(ip, aTsc, aDta, varargin{:});
            
            this.pnumber      = ip.Results.aTsc.pnumber;
            this.scanIndex    = ip.Results.aDta.scanIndex;
            this.dta          = ip.Results.aDta;
            this.Measurements = ip.Results.aTsc.activity;
            this.region       = ip.Results.region;
            this.times        = ip.Results.aTsc.times;            
            this.gluTxlsx_    = mlpowers.GluTxlsx('Mode', 'WholeBrain');
            this.k04_         = this.FB / this.VB;
        end
        function Q    = itsConcentrationQ(this)
            Q = this.concentrationQ(this.k04_, this.k12, this.k21, this.k32, this.k43, this.t0, this.dta, this.VB, this.times);
        end    
        function logL = logLhood(this, Obj)
            this.sigma = this.uniform2limits(Obj.sigma, this.limits('sigma'));
            logL       = this.logLhood_Cauchy(Obj);
        end
        function logL = logLhood_Chi2(this, Obj)
            DiffSquare = (this.Measurements - this.Estimation(Obj)).^2;
            len        = length(this.Measurements);
            S          = sqrt(sum(DiffSquare)/len);
            logL       = -len*(log(S) + this.logSqrt2pi) - (1/(2*S^2))*sum(DiffSquare);
        end
        function logL = logLhood_Cauchy(this, Obj)
            logL     = sum(log((this.sigma/pi)./((this.Measurements - this.Estimation(Obj)).^2 + this.sigma^2)));
        end
        function est  = Estimation(this, Obj)
            this.k12 = this.uniform2limits(Obj.k12, this.limits('k12'));
            this.k21 = this.uniform2limits(Obj.k21, this.limits('k21'));
            this.k32 = this.uniform2limits(Obj.k32, this.limits('k32'));
            this.k43 = this.uniform2limits(Obj.k43, this.limits('k43')); 
            this.t0  = this.uniform2limits(Obj.t0,  this.limits('t0'));            
            est      = this.itsConcentrationQ;
        end  
        function Obj  = Prior(this)
            Obj = struct( ...
                'k12', this.UNIFORM, ...
                'k21', this.UNIFORM, ...
                'k32', this.UNIFORM, ...
                'k43', this.UNIFORM, ...
                't0',  this.UNIFORM, ...
                'sigma', this.UNIFORM, ...
                'logL', [], ...
                'logWt', []);
            Obj.logL = this.logLhood(Obj);
        end
 	end 

    %% PRIVATE

    properties (Access = 'private')
        gluTxlsx_
        k04_
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

