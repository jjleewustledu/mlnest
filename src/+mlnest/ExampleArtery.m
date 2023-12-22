classdef ExampleArtery < handle & mlsystem.IHandle
    %% line1
    %  line2
    %  
    %  Created 20-Dec-2023 23:05:16 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlnest/src/+mlnest.
    %  Developed on Matlab 23.2.0.2459199 (R2023b) Update 5 for MACA64.  Copyright 2023 John J. Lee.
    
    properties
        artery
        Data        
        measurement
        measurement_sigma
        times_sampled
    end

    methods
        function ic = simulate(this, product, opts)
            arguments
                this mlnest.ExampleArtery
                product struct
                opts.filepath {mustBeFolder} = this.artery.filepath
                opts.fileprefix {mustBeTextScalar} = ""
            end
            if isemptytext(opts.fileprefix)
                opts.fileprefix = mlpipeline.Bids.adjust_fileprefix( ...
                    this.artery.fileprefix, ...
                    post_proc=stackstr(use_dashes=true));
            end
            
            pr = this.prior;
            parnames = [pr(:,1); fields(this.Data)];
            parvals = [num2cell(ascol(product.ks)); struct2cell(this.Data)];
            img = this.signalmodel(this.times_sampled, parnames, parvals);
            ic = copy(this.artery);
            ic = ic.selectImagingTool(img=img);
            ic.filepath = opts.filepath;
            ic.fileprefix = opts.fileprefix;
        end
    end

    methods (Static)
        function this = create(opts)
            arguments
                opts.artery mlfourd.ImagingContext2
                opts.tracer {mustBeTextScalar} = "Unknown"
            end

            this = mlnest.ExampleArtery();

            this.artery = mlpipeline.ImagingMediator.ensureFiniteImagingContext(opts.artery);
            ifc = this.artery.imagingFormat;
            this.times_sampled = asrow(ifc.json_metadata.timesMid);
            this.measurement = asrow(double(ifc.img));
            M0 = max(this.measurement);
            this.measurement_sigma = 0.05*M0*ones(size(this.measurement));

            this.Data.M0 = M0;
            this.Data.N = floor(this.times_sampled(end) - this.times_sampled(1)) + 1; % missing ts(1) -> bugs elsewhere?
            this.Data.tracer = opts.tracer;
        end
        function p = prior()
            p = {'a', 'uniform', 0.5, 2, 'fixed'; ...
                 'b', 'uniform', 0.05, 0.15, 'fixed'; ...
                 'p', 'uniform', 0.25, 3, 'fixed'; ...
                 't0', 'uniform', 0, 30, 'fixed'; ...
                 'ss', 'uniform', 0, 0.2, 'fixed'; ...
                 'g', 'uniform', 0.01, 0.2, 'fixed'};
        end
        function qs = signalmodel(t, parnames, parvals)

            import mlaif.ArteryLee2021Model.slide

            a = cell2mat(parvals(strcmp(parnames, 'a')));
            b = cell2mat(parvals(strcmp(parnames, 'b')));
            p = cell2mat(parvals(strcmp(parnames, 'p')));
            t0 = cell2mat(parvals(strcmp(parnames, 't0')));
            ss_frac = cell2mat(parvals(strcmp(parnames, 'ss')));
            g = cell2mat(parvals(strcmp(parnames, 'g')));
            M0 = cell2mat(parvals(strcmp(parnames, 'M0')));

            if (t(1) >= t0) 
                t_ = t - t0;
                k_ = t_.^a .* exp(-(b*t_).^p);
                ss_ = 1 - exp(-g*t_);
                qs = (1 - ss_frac)*k_ + ss_frac*ss_;
            else % k is complex for t - t0 < 0
                t_ = t - t(1);
                k_ = t_.^a .* exp(-(b*t_).^p);
                ss_ = 1 - exp(-g*t_);
                qs = (1 - ss_frac)*k_ + ss_frac*ss_;
                qs = slide(qs, t, t0 - t(1));
            end
            qs = qs/max(qs); % \in [0 1]
            assert(all(imag(qs) == 0))

            tindices = floor(t - t(1)) + 1;
            tindices(tindices > length(qs)) = [];
            qs = qs(tindices);
            try
                qs = M0*qs;
            catch ME
                handexcept(ME)
            end
        end
    end

    %% PROTECTED

    methods (Static, Access=protected)
        function this = ExampleArtery()
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
