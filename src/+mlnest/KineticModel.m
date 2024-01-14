classdef (Abstract) KineticModel < handle & mlsystem.IHandle
    %% line1
    %  line2
    %  
    %  Created 13-Jan-2024 14:47:12 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlnest/src/+mlnest.
    %  Developed on Matlab 23.2.0.2459199 (R2023b) Update 5 for MACA64.  Copyright 2024 John J. Lee.
    
    methods (Abstract, Static)
        create()
        signalmodel()
    end

    properties 
        input_func % ImagingContext2
        Data        
        scanner % ImagingContext2

        measurement % double, expected by MultiNest
        measurement_sigma % "
        times_sampled % "
    end

    properties (Dependent)
        artery % alas needed by MultiNest
    end

    methods %% GET
        function g = get.artery(this)
            g = this.input_func;
        end
    end

    methods
        function [signal,ideal] = simulate(this, product, opts)
            %% returns ImagingContext2

            arguments
                this mlnest.KineticModel
                product struct
                opts.filepath {mustBeFolder} = this.scanner.filepath
                opts.fileprefix {mustBeTextScalar} = ""
            end
            if isemptytext(opts.fileprefix)
                opts.fileprefix = mlpipeline.Bids.adjust_fileprefix( ...
                    this.scanner.fileprefix, ...
                    post_proc=stackstr(3, use_dashes=true));
            end
            fp_signal = mlpipeline.Bids.adjust_fileprefix( ...
                opts.fileprefix, ...
                new_mode="signal");
            
            pr = this.prior(this.Data);
            parnames = [pr(:,1); fields(this.Data)];
            parvals = [num2cell(ascol(product.mean_ks)); struct2cell(this.Data)];
            signal_numer = this.signalmodel(ascol(this.times_sampled), parnames, parvals);

            ic = copy(this.scanner);

            ic.json_metadata.mean_ks = product.mean_ks;
            ic.json_metadata.std_ks = product.std_ks;
            ic.json_metadata.logZ = product.logZ;
            ic.json_metadata.loss = product.loss;
            ic.json_metadata.M0 = this.Data.M0;
            ic.json_metadata.model_kind = this.Data.model_kind;
            ic.json_metadata.preferred_names = this.preferred_names;
            ic.json_metadata.Prior = this.prior(this.Data);
            ic.json_metadata.tracer = this.Data.tracer;

            signal = copy(ic);
            signal = signal.selectImagingTool(img=asrow(signal_numer));
            signal.filepath = opts.filepath;
            signal.fileprefix = fp_signal;

            ideal = [];
        end
    end

    methods (Static)
        function ds = datestr()
            ds = string(datetime("now", Format="yyyyMMddHHmmss"));
        end 

        function t = Data2t(Data)
            %% return timesStart as col

            t = (Data.timesMid(1)-Data.taus(1)/2):(Data.timesMid(end)+Data.taus(end)/2)';
        end
        function [vec,T] = ensureRow(vec)
            if (~isrow(vec))
                vec = vec';
                T = true;
                return
            end
            T = false; 
        end
    end

    %% PROTECTED

    methods (Access = protected)
        function this = KineticModel()
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
