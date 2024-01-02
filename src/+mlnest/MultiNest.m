classdef MultiNest < mlio.AbstractIO
    %% is a value class for performance, using nested_sampler(): 
    %
    % function [logZ, nest_samples, post_samples] = nested_sampler(data, ...
    %           Nlive, Nmcmc, tolerance, likelihood, model, prior, extraparams)
    %
    % This function performs nested sampling of the likelihood function from
    % the given prior (given a set of data, a model, and a set of extra model
    % parameters).
    %
    % By default the algorithm will draw new samples from a set of bounding
    % ellipsoids constructed using the MultiNest algorithm for partitioning
    % live points. However, if the optional 'Nmcmc' argument is set and
    % Nmcmc > 0, new samples will be drawn from a proposal using an MCMC. This
    % method is based on that of Veitch & Vecchio. For both methods the
    % sampling will stop once the tolerance critereon has been reached.
    %
    % The likelihood should be the function handle of a likelihood function to
    % use. This should return the log likelihood of the model parameters given
    % the data.
    %
    % The model should be the function handle of the model function to be
    % passed to the likelihood function.
    %
    % The prior should be a cell array with each cell containing five values:
    %   parameter name (string)
    %   prior type (string) e.g. 'uniform', 'gaussian' of 'jeffreys'
    %   minimum value (for uniform prior), or mean value (for Gaussian prior)
    %   maximum value (for uniform prior), or width (for Gaussian prior)
    %   parameter behaviour (string):
    %       'reflect' - if the parameters reflect off the boundaries
    %       'cyclic'  - if the parameter space is cyclic
    %       'fixed'   - if the parameters have fixe boundaries
    %       ''        - for gaussian priors
    %   e.g., prior = {'h0', 'uniform', 0, 1, 'reflect';
    %                  'r', 'gaussian', 0, 5, '';
    %                  'phi', 'uniform', 0, 2*pi, 'cyclic'};
    %
    % extraparams is a cell array of fixed extra parameters (in addition
    % to those specified by prior) used by the model
    % e.g.  extraparams = {'phi', 2;
    %                      'x', 4};
    %
    % Optional arguments:
    %  Set these via e.g. 'Nmcmc', 100
    %   Nmcmc - if this is set then MultiNest will not be used as the sampling
    %           algorithm. Instead an MCMC chain with this number of iterations
    %           will be used to draw the number nested sample point.
    %   Nsloppy - if this is set then during the MCMC the likelihood will only
    %             be evaluted once every Nsloppy points rather than at every
    %             iteration of the chain.
    %   covfrac - the relative fraction of the iterations for which the MCMC
    %             proposal distribution will be based on a Students-t
    %             distribution defined by the covariance of the current live
    %             points.
    %   diffevfrac - the relative fraction of the iterations that will use
    %                differential evolution to draw the new sample.
    %   stretchfrac - the relative fraction of the iterations that will use the
    %                 affine invariant ensemble stretch method for drawing a
    %                 new sample
    %   walkfrac - the relative fraction of the iterations that will use the
    %              affine invariant ensemble walk method for drawing a new
    %              sample
    %   propscale - the scaling factor for the covariance matrix used by the
    %               'covfrac' Students-t distribution proposal. This defaults
    %               to 0.1.
    %
    % E.g. if covfrac = 10 then diffevfrac = 5 the Students-t proposal will be
    % used 2/3s of the time and differential evolution 1/3. The default is to
    % use the affine invariant samplers with the stretch move 75% of the time
    % and the walk move 25% of the time.
    %  
    % Created 18-Dec-2023 23:59:32 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlnest/src/+mlnest.
    % Developed on Matlab 23.2.0.2459199 (R2023b) Update 5 for MACA64.  Copyright 2023 John J. Lee.
    
	properties
        context         % client, must implement:
                        % prior, simulate(), Data, measurement, measurement_sigma, times_sampled
        
        Data            % extra data (metadata) as struct
        Measurement     % sampled measurements to be inferred
        Prior           % cell array
        Sigma           % estimated sigma for Measurement
        TimesSampled    % sampled times for Measurement
    end

    properties (Dependent)
        Data_cell       % primary data as cell array structured as expected by nested_sampler()
        Data_cell_extra % extra data as cell array structured as expected by nested_sampler()
        preferred_names
        product
    end

    methods % GET
        function g = get.Data_cell(this)
            g{1} = ascol(this.TimesSampled);
            g{2} = ascol(this.Measurement);
            Sigma_ = ascol(this.Sigma);
            g{3} = diag((Sigma_.^2)); % convert sigmas to a covariance matrix for logL_gaussian
        end
        function g = get.Data_cell_extra(this)
            g = [fields(this.Data), struct2cell(this.Data)];
        end
        function g = get.preferred_names(this)
            if isemptytext(this.context.preferred_names)
                g = this.Prior{:,1};
                return
            end
            g = this.context.preferred_names;
        end
        function g = get.product(this)
            g = this.product_;
        end
    end

    methods
        function k = ks(this)
            %% returns row vector

            k = zeros(1, size(this.Prior, 1));
            for ki = 1:size(this.Prior, 1)
                k(ki) = this.posteriors(this.post_samples_, ki, {this.Prior{ki,1}}, show_plot=false);
            end
        end
        function loss_ = loss(this)
            %% smaller positive is better ~ sum(abs(fitted/measured - 1))

            arguments
                this mlnest.MultiNest
            end

            % params_names = [this.Prior(:,1); this.Data_cell_extra(:,1)];
            % params_values = [num2cell(ascol(this.ks())); this.Data_cell_extra(:,2)];
            % loss_ = this.logL_gaussian( ...
            %     this.Data_cell, opts.signal_model, params_names, params_values);

            try
                product0 = struct( ...
                    'ks', this.ks(), ...
                    'logZ', this.logZ_, ...
                    'loss', NaN);
                fitted = asrow(double(this.simulate(product0))); % simulate() returns [signal, ideal]
                measured = asrow(this.Measurement);
                positive = fitted + measured > 0;
                fp = fitted(positive);
                mp = measured(positive);
                loss_ = mean(abs(2*(fp - mp)./(fp + mp))); % unbiased mean abs residual
            catch ME
                handwarning(ME)
                loss_ = [];
            end
        end
        function plot_posteriors(this, opts)
            arguments
                this mlnest.MultiNest
                opts.singles logical = true
                opts.montage logical = true
                opts.do_save logical = false
            end

            pwd0 = pushd(this.filepath);
            
            if opts.singles
                for ki = 1:size(this.Prior, 1)
                    pname = this.preferred_names(ki);
                    this.posteriors(this.product.post_samples, ki, pname);
                end
            end
            if opts.montage
                wp = 1:size(this.Prior, 1);
                pnames = asrow(this.preferred_names);
                this.posteriors(this.product.post_samples, wp, pnames);
            end
            if opts.do_save
                saveFigures(this.filepath)
            end

            popd(pwd0);
        end
        function est = rescaleModelEstimate(this, est, opts)
            arguments
                this mloptimization.SimulatedAnnealing
                est {mustBeNumeric}
                opts.norm {mustBeTextScalar} = "max"
            end

            if strcmpi(opts.norm, "max")
                M0 = max(this.Measurement);
                est = M0*est;
                return
            end

            % rescale by this.int_dt_M
            int_dt_M = trapz(this.TimesSampled, this.Measurement);
            int_dt_est = trapz(this.TimesSampled, est);
            est = est*int_dt_M/int_dt_est;
        end
        function save(this)
            try
                fp = mlpipeline.Bids.adjust_fileprefix( ...
                    this.context.artery.fileprefix, ...
                    post_proc=stackstr(3, use_dashes=true)+"-"+this.context.datestr());
            catch ME
                handwarning(ME)
                fp = this.fileprefix+"-"+this.context.datestr();
            end
            fqfp = fullfile(this.filepath, fp);
            save(strcat(fqfp, ".mat"), "this");
            
            warning("off", "MATLAB:structOnObject");
            struct_this = struct(this);
            save(strcat(fqfp, "_struct.mat"), "struct_this"); % anticipating changing MultiNest src.
            warning("on", "MATLAB:structOnObject");
        end
        function saveas(this, fn)
            save(fn, "this");
        end
        function ic = simulate(this, product)
            arguments
                this mlnest.MultiNest
                product struct = this.product
            end

            ic = this.context.simulate(product);
        end
        function this = solve(this, opts)
            % example: (from hogg et al., 1008.4686)
            % fit a signal model to data keeping outlier points

            % provide sample data for fitting a line
            % NOTE: data{1} = x_i, data{2} = y_i, data{3} = sigma_yi

            arguments
                this mlnest.MultiNest
                opts.DEBUG double = 0
                opts.likelihood function_handle = @mlnest.MultiNest.logL_gaussian
                opts.Nlive double = 50
                opts.Nmcmc double = 0 % number of iterations for MCMC sampling: (enter 0 for multinest sampling)
                opts.signal_model function_handle
                opts.tol double = 0.1
                opts.verbose double = 0
            end
            assert(opts.Nlive > size(this.Prior, 1) + size(this.Data_cell_extra, 1))
            
            global verbose; %#ok<*GVMIS>
            verbose = opts.verbose;
            global DEBUG;
            DEBUG = opts.DEBUG;
            
            % called nested sampling routine 
            [this.logZ_, this.nest_samples_, this.post_samples_] = ...
                mlnest.MultiNest.nested_sampler( ...
                this.Data_cell, opts.Nlive, ...
                opts.tol, opts.likelihood, opts.signal_model, this.Prior, this.Data_cell_extra, Nmcmc=opts.Nmcmc);

            this.product_ = struct( ...
                'ks', this.ks(), ...
                'likelihood', func2str(opts.likelihood), ...
                'logZ', this.logZ_, ...
                'loss', this.loss(), ...
                'nest_samples', this.nest_samples_, ...
                'Nlive', opts.Nlive, ...
                'Nmcmc', opts.Nmcmc, ...
                'post_samples', this.post_samples_, ...
                'signal_model', func2str(opts.signal_model), ...
                'tol', opts.tol);
        end

        function this = MultiNest(opts)
            arguments
                opts.context mlsystem.IHandle
                opts.filepath {mustBeFolder} = pwd
                opts.fileprefix {mustBeTextScalar} = stackstr(3)
            end

            % first-order vars
            this.context = opts.context;
            this.filepath = opts.filepath;
            this.fileprefix = opts.fileprefix;

            % second-order vars; copy objects for speed
            this.Data = this.context.Data;
            this.Measurement = this.context.measurement;
            this.Prior = this.context.prior(this.Data);
            this.Sigma = this.context.measurement_sigma;
            this.TimesSampled = this.context.times_sampled;
        end
    end

    methods (Static)
        function logL = logL_gaussian(data, model, parnames, parvals)
            %% logL = logL_gaussian(data, model, parnames, parvals)
            %
            % This function will compute the log likelihood of a multivariate
            % gaussian:
            %
            %     L = 1/sqrt((2 pi)^N det C)
            %         exp[-0.5*(y - model(x,params))^T * inv(C) * (y - model(x,params))]
            %
            % The input parameters are:
            %     data - a cell array with three columns
            %            { x values, y values, C: covariance matrix }
            %     NOTE: if C is a single number, convert to a diag covariance matrix
            %     model - the function handle for the signal model.
            %     parnames - a cell array listing the names of the model parameters
            %     parvals - a cell array containing the values of the parameters given
            %         in parnames. These must be in the same order as in parnames.
            %         If parvals is an empty vector the noise-only likelihood will be
            %         calculated.
            %
            % -------------------------------------------------------------------------
            %           This is the format required by nested_sampler.m.
            % -------------------------------------------------------------------------
            
            % check whether model is a string or function handle
            if ischar(model)
                fmodel = str2func(model);
            elseif isa(model, 'function_handle')
                fmodel = model;
            else
                error('Error... Expecting a model function!');
            end
            
            % get data values from cell array
            x = data{1};
            y = data{2};
            C = data{3};
            N = length(x);
            
            % evaluate the model
            if isempty(parvals)
                % if parvals is not defined get the null likelihood (noise model
                % likelihood)
                md = 0;
            else
                md = feval(fmodel, x, parnames, parvals);
            
                % if the model returns a NaN then set the likelihood to be zero (e.g.
                % loglikelihood to be -inf
                if isnan(md)
                    logL = -inf;
                    return;
                end
            end
            
            % get inverse of covariance matrix and the log of the determinant
            if size(C,1)==1 || N == 1
                % variance is just a single value
                invC = 1/C;
                lDetC = N*log(C);
                % calculate the log likelihood
                logL = -0.5*(y - md)' * invC * (y - md);
            elseif size(C,1) == N && size(C,2) == N && N > 1
                % use trick from http://xcorr.net/2008/06/11/log-determinant-of-positive-definite-matrices-in-matlab/
                % to calculate log of determinant and avoid infinities
                Cchol = chol(C); % Cholesky decomposition
                lDetC = 2*sum(log(diag(Cchol)));
            
                % calculate the log likelihood (use / for matrix inverse)
                logL = -0.5*((y - md)'/C)*(y - md);
            end
            
            % calculate the log likelihood
            logL = logL - 0.5*N*log(2*pi) - 0.5*lDetC;
            
            if isnan(logL)
                error('Error: log likelihood is NaN!');
            end
        end
        function [logZ, nest_samples, post_samples] = nested_sampler(data, ...
                Nlive, tolerance, likelihood, model, prior, extraparams, ...
                varargin)

            % function [logZ, nest_samples, post_samples] = nested_sampler(data, ...
            %           Nlive, Nmcmc, tolerance, likelihood, model, prior, extraparams)
            %
            % This function performs nested sampling of the likelihood function from
            % the given prior (given a set of data, a model, and a set of extra model
            % parameters).
            %
            % By default the algorithm will draw new samples from a set of bounding
            % ellipsoids constructed using the MultiNest algorithm for partitioning
            % live points. However, if the optional 'Nmcmc' argument is set and
            % Nmcmc > 0, new samples will be drawn from a proposal using an MCMC. This
            % method is based on that of Veitch & Vecchio. For both methods the
            % sampling will stop once the tolerance critereon has been reached.
            %
            % The likelihood should be the function handle of a likelihood function to
            % use. This should return the log likelihood of the model parameters given
            % the data.
            %
            % The model should be the function handle of the model function to be
            % passed to the likelihood function.
            %
            % The prior should be a cell array with each cell containing five values:
            %   parameter name (string)
            %   prior type (string) e.g. 'uniform', 'gaussian' of 'jeffreys'
            %   minimum value (for uniform prior), or mean value (for Gaussian prior)
            %   maximum value (for uniform prior), or width (for Gaussian prior)
            %   parameter behaviour (string):
            %       'reflect' - if the parameters reflect off the boundaries
            %       'cyclic'  - if the parameter space is cyclic
            %       'fixed'   - if the parameters have fixe boundaries
            %       ''        - for gaussian priors
            %   e.g., prior = {'h0', 'uniform', 0, 1, 'reflect';
            %                  'r', 'gaussian', 0, 5, '';
            %                  'phi', 'uniform', 0, 2*pi, 'cyclic'};
            %
            % extraparams is a cell array of fixed extra parameters (in addition
            % to those specified by prior) used by the model
            % e.g.  extraparams = {'phi', 2;
            %                      'x', 4};
            %
            % Optional arguments:
            %  Set these via e.g. 'Nmcmc', 100
            %   Nmcmc - if this is set then MultiNest will not be used as the sampling
            %           algorithm. Instead an MCMC chain with this number of iterations
            %           will be used to draw the number nested sample point.
            %   Nsloppy - if this is set then during the MCMC the likelihood will only
            %             be evaluted once every Nsloppy points rather than at every
            %             iteration of the chain.
            %   covfrac - the relative fraction of the iterations for which the MCMC
            %             proposal distribution will be based on a Students-t
            %             distribution defined by the covariance of the current live
            %             points.
            %   diffevfrac - the relative fraction of the iterations that will use
            %                differential evolution to draw the new sample.
            %   stretchfrac - the relative fraction of the iterations that will use the
            %                 affine invariant ensemble stretch method for drawing a
            %                 new sample
            %   walkfrac - the relative fraction of the iterations that will use the
            %              affine invariant ensemble walk method for drawing a new
            %              sample
            %   propscale - the scaling factor for the covariance matrix used by the
            %               'covfrac' Students-t distribution proposal. This defaults
            %               to 0.1.
            %
            % E.g. if covfrac = 10 then diffevfrac = 5 the Students-t proposal will be
            % used 2/3s of the time and differential evolution 1/3. The default is to
            % use the affine invariant samplers with the stretch move 75% of the time
            % and the walk move 25% of the time.
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            global verbose;
            global DEBUG;

            Nmcmc = 0; % default with this set to zero is to use MultiNest
            Nsloppy = 0;
            covfrac = 0;
            diffevfrac = 0;
            walkfrac = 25;
            stretchfrac = 75;
            propscale = 0.1;

            % get optional input arguments
            optargin = size(varargin,2);

            % get optional arguments
            if optargin > 1
                for i = 1:2:optargin
                    if strcmpi(varargin{i}, 'Nmcmc') % number of MCMC samples
                        if ~isempty(varargin{i+1})
                            if varargin{i+1} < 1
                                fprintf(1, 'Using MultiNest algorithm\n');
                            else
                                Nmcmc = varargin{i+1};
                            end
                        end
                    elseif strcmpi(varargin{i}, 'Nsloppy') % number of burn in samples
                        if ~isempty(varargin{i+1})
                            if varargin{i+1} < 0
                                fprintf(1, 'Number of \"sloppy\" samples is silly. Setting to zero\n');
                            else
                                Nsloppy = varargin{i+1};
                            end
                        end
                    elseif strcmpi(varargin{i}, 'covfrac') % fraction of MCMC iterations using Students't proposal
                        if varargin{i+1} > 0
                            covfrac = varargin{i+1};
                        end
                    elseif strcmpi(varargin{i}, 'diffevfrac') % fraction of MCMC iterations using differential evolution
                        if varargin{i+1} > 0
                            diffevfrac = varargin{i+1};
                        end
                    elseif strcmpi(varargin{i}, 'walkfrac') % fraction of MCMC iterations using walk move
                        if varargin{i+1} >= 0
                            walkfrac = varargin{i+1};
                        end
                    elseif strcmpi(varargin{i}, 'stretchfrac') % fraction of MCMC iterations using stretch move
                        if varargin{i+1} >= 0
                            stretchfrac = varargin{i+1};
                        end
                    elseif strcmpi(varargin{i}, 'propscale') % the scaling factor for the covariance matrix
                        if varargin{i+1} > 0
                            propscale = varargin{i+1};
                        end
                    end
                end
            end

            % get the number of parameters from the prior array
            D = size(prior,1);

            % get all parameter names
            parnames = prior(:,1);

            if ~isempty(extraparams)
                extraparnames = extraparams(:,1);
                extraparvals = extraparams(:,2);
                parnames = cat(1, parnames, extraparnames);
            else
                extraparvals = [];
            end

            % draw the set of initial live points from the prior
            livepoints = zeros(Nlive, D);

            for i=1:D
                priortype = prior{i,2};
                p3 = prior{i,3};
                p4 = prior{i,4};

                % currently only handles uniform or Gaussian priors
                if strcmp(priortype, 'uniform')
                    livepoints(:,i) = p3 + (p4-p3)*rand(Nlive,1);
                elseif strcmp(priortype, 'gaussian')
                    livepoints(:,i) = p3 + p4*randn(Nlive,1);
                elseif strcmp(priortype, 'jeffreys')
                    % uniform in log space
                    livepoints(:,i) = 10.^(log10(p3) + (log10(p4)-log10(p3))*rand(Nlive,1));
                end
            end

            % check whether likelihood is a function handle, or a string that is a
            % function name
            if ischar(likelihood)
                flike = str2func(likelihood);
            elseif isa(likelihood, 'function_handle')
                flike = likelihood;
            else
                error('Error... Expecting a model function!');
            end

            % calculate the log likelihood of all the live points
            logL = zeros(Nlive,1);

            for i=1:Nlive
                parvals = cat(1, num2cell(livepoints(i,:)'), extraparvals);
                logL(i) = feval(flike, data, model, parnames, parvals);
            end

            % now scale the parameters, so that uniform parameters range from 0->1,
            % and Gaussian parameters have a mean of zero and unit standard deviation
            for i=1:Nlive
                livepoints(i,:) = scale_parameters(prior, livepoints(i,:));
            end

            % initial tolerance
            tol = inf;

            % initial width of prior volume (from X_0=1 to X_1=exp(-1/N))
            logw = log(1 - exp(-1/Nlive));

            % initial log evidence (Z=0)
            logZ = -inf;

            % initial information
            H = 0;

            % initialize array of samples for posterior
            nest_samples = zeros(1,D+1);

            %%%%%%%%%%%%%%%%
            % some initial values if MultiNest sampling is used
            h = 1.1; % h values from bottom of p. 1605 of Feroz and Hobson
            FS = h; % start FS at h, so ellipsoidal partitioning is done first time
            K = 1; % start with one cluster of live points

            % get maximum likelihood
            logLmax = max(logL);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % initialize iteration counter
            j = 1;

            %figure;

            % MAIN LOOP
            while tol > tolerance || j <= Nlive

                % expected value of true remaining prior volume X
                VS = exp(-j/Nlive);

                % find minimum of likelihoods
                [logLmin, idx] = min(logL);

                % set the sample to the minimum value
                nest_samples(j,:) = [livepoints(idx, :) logLmin];

                % get the log weight (Wt = L*w)
                logWt = logLmin + logw;

                % save old evidence and information
                logZold = logZ;
                Hold = H;

                % update evidence, information, and width
                logZ = logplus(logZ, logWt);
                H = exp(logWt - logZ)*logLmin + ...
                    exp(logZold - logZ)*(Hold + logZold) - logZ;
                %logw = logw - logt(Nlive);
                logw = logw - 1/Nlive;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if Nmcmc > 0

                    % do MCMC nested sampling

                    % get the Cholesky decomposed covariance of the live points
                    % (do every 100th iteration - CAN CHANGE THIS IF REQUIRED)
                    if mod(j-1, 100) == 0
                        % NOTE that for numbers of parameters >~10 covariances are often
                        % not positive definite and cholcov will have "problems".
                        %cholmat = cholcov(propscale*cov(livepoints));

                        % use modified Cholesky decomposition, which works even for
                        % matrices that are not quite positive definite
                        % from http://infohost.nmt.edu/~borchers/ldlt.html
                        % (via http://stats.stackexchange.com/questions/6364
                        % /making-square-root-of-covariance-matrix-positive-definite-matlab
                        cv = cov(livepoints);
                        [l, d] = mchol(propscale*cv);
                        cholmat = l.'*sqrt(d);

                        %plot3(livepoints(:,1), livepoints(:,2), livepoints(:,3), 'r.');
                        %drawnow();
                    end

                    % draw a new sample using mcmc algorithm
                    [livepoints(idx, :), logL(idx)] = draw_mcmc(livepoints, cholmat, ...
                        logLmin, prior, data, flike, model, Nmcmc, Nsloppy, ...
                        covfrac, diffevfrac, walkfrac, stretchfrac, parnames, ...
                        extraparvals);

                else

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % do MultiNest nested sampling

                    % separate out ellipsoids
                    if FS >= h
                        % NOTE: THIS CODE IS GUARANTEED TO RUN THE 1ST TIME THROUGH
                        % calculate optimal ellipsoids
                        [Bs, mus, VEs, ns] = optimal_ellipsoids(livepoints, VS);
                        K = length(VEs); % number of ellipsoids (subclusters)

                    else
                        % simply rescale the bounding ellipsoids
                        for k=1:K
                            scalefac = max([1 (exp(-(j+1)/Nlive)*ns(k)/Nlive)/VEs(k)]);

                            % scale bounding matrix and volume
                            if scalefac ~= 1
                                Bs((k-1)*D+1:k*D,:) = Bs((k-1)*D+1:k*D,:)*scalefac^(2/D);
                                VEs(k) = scalefac*VEs(k);
                            end
                        end

                    end

                    if DEBUG && D==2
                        % plot 2-dimensionsal live points and bounding ellipses
                        plot_2d_livepoints_with_ellipses(livepoints, Bs, mus);
                    end

                    % calculate ratio of volumes (FS>=1) and cumulative fractional volume
                    Vtot = sum(VEs);
                    FS = Vtot/VS;
                    fracvol = cumsum(VEs)/Vtot;

                    % draw a new sample using multinest algorithm
                    [livepoints(idx, :), logL(idx)] = draw_multinest(fracvol, ...
                        Bs, mus, logLmin, prior, data, flike, model, ...
                        parnames, extraparvals);

                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % update maximum likelihood if appropriate
                if logL(idx) > logLmax
                    logLmax = logL(idx);
                end

                % work out tolerance for stopping criterion
                tol = logplus(logZ, logLmax - (j/Nlive)) - logZ;

                % display progress (optional)
                if verbose
                    fprintf(1, 'log(Z): %.5e, tol = %.5e, K = %d, iteration = %d\n', ...
                        logZ, tol, K, j);
                end

                % update counter
                j = j+1;

            end

            % sort the remaining points (in order of likelihood) and add them on to
            % the evidence
            [logL_sorted, isort] = sort(logL);
            livepoints_sorted = livepoints(isort, :);

            for i=1:Nlive
                logZ = logplus(logZ, logL_sorted(i) + logw);
            end

            % append the additional livepoints to the nested samples
            nest_samples = [nest_samples; livepoints_sorted logL_sorted];

            % rescale the samples back to their true ranges
            for i=1:length(nest_samples)
                nest_samples(i,1:end-1) = ...
                    rescale_parameters(prior, nest_samples(i,1:end-1));
            end

            % convert nested samples into posterior samples - nest2pos assumes that the
            % final column in the sample chain is the log likelihood
            post_samples = nest2pos(nest_samples, Nlive);
        end
        function [par_mean,par_std] = posteriors(post_samples, wp, parnames, opts)
            %% function posteriors(post_samples, wp)
            %
            % This function displays joint and/or marginalized posterior
            % distributions, calculates summary statistics, etc.
            %
            % post_samples is an Nx(npars+2) array where N is the number of
            % posterior samples, npars is the number of parameters, and the last two
            % columns in this array contain the values of logL and logPosterior
            % (=prior weighted likelihood/evidence).
            %
            % wp is a vector containing the parameter(s) for which you want to create
            % the posterior e.g., wp=1 will provide a posterior only on the 1st
            % parameter; wp=[1 3] will provide a 2D posterior on the 1st and 3rd
            % parameters. For the moment wp can only contain a maximum of two params.
            %
            % parnames is a cell array of string names corresponding to the
            % parameters specified by wp.  Eg., parnames = {'b', 'slope'}
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            arguments
                post_samples double
                wp double
                parnames {mustBeText}
                opts.show_plot logical = true
            end

            npars = size(post_samples,2)-2;
            lwp = length(wp);
            lparnames = length(parnames);
            nbins = 100; % number of bins for historgram plots

            % check the length of the vector containing which parameters
            if lwp > npars || max(wp) > npars || lwp < 1 || lwp~=lparnames % lwp > 2 
                error('Error... parameters not properly specified!');
            end

            if lwp == 1
                % calculate mean and std deviation
                par_mean = mean(post_samples(:,wp));
                par_std = std(post_samples(:,wp));
                if opts.show_plot
                    fprintf('Parameter %s: Mean = %f, stddev = %f\n', ...
                        parnames{1}, par_mean, par_std);
                end
            end

            if opts.show_plot
                if lwp == 1
                    % make histogram plot
                    figure()
                    histogram(post_samples(:,wp), nbins);
                    xlabel(parnames{1},'fontsize',16,'FontWeight','bold');
                elseif lwp == 2
                    % make 2-d histogram plot
                    p1 = wp(1);
                    p2 = wp(2);
                    edges1 = linspace(min(post_samples(:,p1)), max(post_samples(:,p1)), nbins);
                    edges2 = linspace(min(post_samples(:,p2)), max(post_samples(:,p2)), nbins);
                    histmat = hist2(post_samples(:,p1), post_samples(:,p2), edges1, edges2);
    
                    figure()
                    imagesc(edges1, edges2, transpose(histmat));
                    colormap(viridis)
                    set(gca,'YDir','normal')
                    xlabel(parnames{1},'fontsize',16,'FontWeight','bold');
                    ylabel(parnames{2},'fontsize',16,'FontWeight','bold');
                elseif lwp > 2
                    % make tiled 2-d histogram plots
                    figure()
                    tiledlayout(lwp, lwp, TileSpacing="tight");
                    for p2 = wp
                        for p1 = flip(wp)
                            ax = nexttile;

                            edges1 = linspace(min(post_samples(:,p1)), max(post_samples(:,p1)), nbins);
                            edges2 = linspace(min(post_samples(:,p2)), max(post_samples(:,p2)), nbins);
                            histmat = hist2(post_samples(:,p1), post_samples(:,p2), edges1, edges2);

                            imagesc(edges1, edges2, transpose(histmat));
                            colormap(viridis)
                            set(ax,'YDir','normal')
                            ax.Visible = "off";
                            if p2 == wp(end)
                                ax.XAxis.Visible = "on";
                                xlabel(parnames{p1},'fontsize',12,'FontWeight','normal');
                            end
                            if p1 == wp(end)
                                ax.YAxis.Visible = "on";
                                ylabel(parnames{p2},'fontsize',12,'FontWeight','normal');
                            end
                        end
                    end
                end
            end
        end
    end
    
    %% PROTECTED
    
    properties (Access = protected)
        logZ_
        nest_samples_
        post_samples_
        product_
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
