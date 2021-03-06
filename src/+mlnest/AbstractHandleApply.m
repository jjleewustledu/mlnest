classdef (Abstract) AbstractHandleApply < handle & matlab.mixin.Copyable & mlnest.IHandleApply 
	%% ABSTRACTHANDLEAPPLY implements application ideas from "Data Analysis:  A Bayesian Tutorial, Second Edition"
    %  by D.S. Sivia and J. Skilling, section 9.3.2  

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
    
    properties
        fileprefix
        Object
        sigma0 % fraction of Measurement < 1
    end
    
    properties (Dependent)
        results
    end
    
    methods (Static)
        function main = run(application)
            assert(isa(application, 'mlnest.IHandleApply'))
            main = mlnest.HandleNestedSamplingMain(application);
        end
        function        run_varying(application, par, rng)
            assert(ischar(par))
            assert(isnumeric(rng))
            
            fp = sprintf('%s_run_varying_%s%g-%g', strrep(class(application), '.', '_'), par, rng(1), rng(end));            
            application.fileprefix = fp;
            ensuredir(fp)
            diary(fullfile(fp, [fp '.log']))
            application.varying_ = true;
            par0 = application.(par);
            for r = rng
                fprintf('\n%s -> %g\n', par, r)
                application.(par) = r;
                tic
                mlnest.NestedSamplingMain(application);
                toc
            end
            application.(par) = par0;
            saveFigures(fp, 'closeFigure', false)
            diary('off')
        end        
        
        function y    = vec2native(u, lims)
            y = lims(2)*u + lims(1)*(1 - u);
        end
        function u    = vec2uniform(y, lims)
            u = (y - lims(1))/(lims(2) - lims(1));
        end
    end
    
    methods
        
        %% GET
        
        function g = get.results(this)
            g = this.results_;
        end
        
        %%
        
        function this = AbstractHandleApply(varargin)
            this.fileprefix = ...
                sprintf('%s_ctor_%s', strrep(class(this), '.', '_'), datestr(now, 'yyyymmddHHMMSS'));
            
            ip = inputParser;            
            ip.KeepUnmatched = true;
            addParameter(ip, 'fileprefix', this.fileprefix, @ischar)
            parse(ip, varargin{:})
            this.fileprefix = ip.Results.fileprefix;
        end
        
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
                flds = ignoreFields(this, flds);
                for f = 1:length(flds)                    
                    Try.(flds{f}) = Obj.(flds{f}) + step*(2*rand() - 1.);  % |move| < step
                    Try.(flds{f}) = Try.(flds{f}) - floor(Try.(flds{f}));  % wraparound to stay within (0,1)
                end
                Try.logL = logLhood(this, Try); % trial likelihood value
                
                %% Accept if and only if within hard likelihood constraint
                %  Cf. Sivia sec. 9.4.4.
                if Try.logL > logLstar % + log(rand())
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
            estimation  = Estimation(this, Obj);
            positive    = this.Measurement > 0;
            measurement = this.Measurement(positive);
            eoverm      = estimation(positive)./measurement;
            Q           = sum((1 - eoverm).^2);
            logL        = -0.5*Q/this.sigma0^2 - sum(log(this.sigma0*measurement));
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
            assert(max_ > min_)
            if init_ < min_ || max_ > init_
                init_ = 0.5*(min_ + max_);
            end
            std_  = 0.25*abs(max_ - min_);
            val   = init_ + randn()*std_;
            while val < min_ || max_ < val % from Josh and Larry
                val = init_ + randn()*std_;
            end
            assert(~isnan(val))
            assert(isfinite(val))
            
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
            ws = zeros(nest, 1);
            chains = zeros(nest, length(flds));
            %wchains = zeros(nest, length(flds));
            Obj1 = structfun(@(x) 0, Samples{1}, 'UniformOutput', false); 
            Obj2 = Obj1;
            moment1 = zeros(1, length(flds));
            moment2 = zeros(1, length(flds));
            
            for ni = 1:nest
                
                w = exp(Samples{ni}.logWt - logZ); % proportional weight
                ws(ni) = w;
                
                for f = 1:length(flds)                    
                    chains(ni, f) = Samples{ni}.(flds{f});
                    %wchains(ni, f) = w*Samples{ni}.(flds{f});
                    
                    Obj1.(flds{f}) = Obj1.(flds{f}) + w* Samples{ni}.(flds{f});
                    Obj2.(flds{f}) = Obj2.(flds{f}) + w*(Samples{ni}.(flds{f}))^2;
                    
                    Obn = this.Obj2native(Samples{ni}); % needed by mean +/- std reports
                    moment1(f) = moment1(f) + w* Obn.(flds{f});
                    moment2(f) = moment2(f) + w*(Obn.(flds{f}))^2;
                end
            end
            
            %% gather
            
            r.flds = flds;  
            r.ws = ws;
            r.chains = chains;  
            %r.wchains = wchains;
            r.Obj1 = Obj1;
            r.Obj2 = Obj2;
            r.moment1 = moment1;
            r.moment2 = moment2;
            this.results_ = r;
        end
        
        %% UTILITIES
        
        function        fprintfModel(this)
            fprintf('Model:\n');
            for f = 1:length(this.results.flds)
                fprintf('\t%s = %f +/- %f\n', ...
                    this.results.flds{f}, ...
                    this.results.moment1(f), ...
                    sqrt(this.results.moment2(f) - this.results.moment1(f)^2));
            end
            fprintf('\tsigma0 = %f\n', this.sigma0);
            for f = asrow(this.results.flds)
                fprintf('\tmap(''%s'') => %s\n', f{1}, struct2str(this.map(f{1})));
            end
        end
        function flds = ignoreFields(this, flds)
            for ig = this.ignoredObjFields
                flds = flds(~strcmp(flds, ig{1}));
            end
        end  % OPTIMIZE?
        function vec  = limits(this, key)
            vec = [this.map(key).min this.map(key).max];
        end
        function o    = Obj2native(this, o)
            %% rescale to native units
            
            import mlnest.AbstractHandleApply.vec2native
            flds = fields(o);
            flds = ignoreFields(this, flds);
            for f = flds'             
                o.(f{1}) = vec2native(o.(f{1}), limits(this, f{1}));
            end
        end
        function o    = Obj2uniform(this, o)
            %% rescale to (0,1)
            
            import mlnest.AbstractHandleApply.vec2uniform
            flds = fields(o);
            flds = ignoreFields(this, flds);
            for f = flds'
                o.(f{1}) = vec2uniform(o.(f{1}), limits(this, f{1}));
            end
        end
        function vec  = Obj2vec(this, obj)
            flds = fields(obj);
            flds = ignoreFields(this, flds);
            vec = zeros(1, length(flds));
            for idx = 1:length(vec)
                vec(idx) = obj.(flds{idx});
            end
        end
        function h    = plotMap(this)
            m = this.map;
            for k = m.keys
                objs(1).(k{1}) = m(k{1}).min;
                objs(2).(k{1}) = m(k{1}).init;
                objs(3).(k{1}) = m(k{1}).max;
            end
            for j = 1:3
                objs(j) = this.Obj2uniform(objs(j));
            end
            h = this.plotObjs(objs);
        end
        function h    = plotMatrix(this)
            if isfield(this.results, 'chains')
                figure
                h = plotmatrix(this.results.chains);
                flds = this.results.flds;
                title({ this.titlePlot('plotMatrix()') cell2str(flds) })
            end
        end
        function h    = plotObjs(this, objs)
            figure;
            hold on
            h = plot(1:length(this.Measurement), this.Measurement, 'o');
            
            lbls = cell(1, length(objs));
            for i = 1:length(objs)
                lbls{i} = struct2str(this.Obj2native(objs(i)));
                est = this.Estimation(objs(i));
                plot(1:length(est), est, '-');
            end           
            title(this.titlePlot('plotObjs()'))
            legend(['measurement' lbls])
            hold off
        end
        function h    = plotResults(this)
            figure;
            est = this.Estimation(this.results.Obj1);
            h = plot(1:length(this.Measurement), this.Measurement, 'o', ...
                     1:length(est), est, '-');
            title(this.titlePlot('plotResults()'))
            legend('measurement', 'estimation')
            annotation('textbox', [.175 .25 .3 .3], 'String', this.sprintfModel(), 'FitBoxToText', 'on')
        end
        function        save(this)
            save([this.fileprefix '.mat'], this);
        end
        function        saveas(this, fn)
            save(fn, this);
        end
        function s    = sprintfModel(this)
            s = sprintf('Model:\n');
            for f = 1:length(this.results.flds)
                s = [s sprintf('\t%s = %f +/- %f\n', ...
                          this.results.flds{f}, ...
                          this.results.moment1(f), ...
                          sqrt(this.results.moment2(f) - this.results.moment1(f)^2))]; %#ok<AGROW>
            end
            s = [s sprintf('\tsigma0 = %f\n', this.sigma0)];
            for f = asrow(this.results.flds)
                s = [s sprintf('\tmap(''%s'') => %s\n', f{1}, struct2str(this.map(f{1})))];
            end
        end
        function t    = titlePlot(this, client)
            assert(ischar(client))
            t = [client ': ' strrep(this.fileprefix, '_', ' ')];
        end
    end
    
    %% PROTECTED
    
    properties (Access = protected)
        results_
        varying_ = false
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
    
    %% HIDDEN & DEPRECATED
    
    methods (Hidden)
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
    
end

