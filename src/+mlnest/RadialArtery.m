classdef RadialArtery < handle & mlnest.Artery
    %% line1
    %  line2
    %  
    %  Created 22-Dec-2023 23:17:27 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlnest/src/+mlnest.
    %  Developed on Matlab 23.2.0.2459199 (R2023b) Update 5 for MACA64.  Copyright 2023 John J. Lee.
    
    properties
    end

    methods
    end

    methods (Static)
        function this = create(opts)
            arguments
                opts.artery mlfourd.ImagingContext2
                opts.kernel {mustBeNumeric} = []
                opts.measurement_sigma {mustBeNumeric} = 1410 % Twilite II
                opts.model_kind {mustBeTextScalar} = "4bolus"
                opts.tracer {mustBeTextScalar} = "Unknown"
            end

            switch opts.model_kind
                case "4bolus"
                    this = mlnest.RadialArtery4bolus();
                case "3bolus"
                    this = mlnest.RadialArtery3bolus();
                case "2bolus"
                    this = mlnest.RadialArtery2bolus();
                otherwise
                    error("mlnest:ValueError", stackstr())
            end

            this.artery = opts.artery;
            ifc = this.artery.imagingFormat;
            this.times_sampled = asrow(ifc.json_metadata.timesMid);
            this.measurement = asrow(double(ifc.img));
            M0 = max(this.measurement);
            this.measurement_sigma = opts.measurement_sigma*ones(size(this.measurement));

            j = this.artery.json_metadata;
            if isempty(opts.kernel)
                cath = mlswisstrace.Catheter_DT20190930(tracer=opts.tracer, model_kind="nomodel");
                opts.kernel = ascol(cath.kernel(Nt=length(this.times_sampled)+1));
            end
            this.Data.kernel = opts.kernel;
            this.Data.M0 = M0;
            this.Data.model_kind = opts.model_kind;
            this.Data.taus = j.taus;
            this.Data.timesMid = j.timesMid;
            this.Data.t = mlnest.Artery.Data2t(this.Data);
            this.Data.tracer = opts.tracer;
        end

        function vec = apply_dispersion(vec_, Data)
            %% Args:            
            %      vec double, has 1 Hz sampling
            %      Data struct
            %  Returns:
            %      vec double, sampled at opts.Data.timesMid

            arguments
                vec_ double 
                Data struct = []
            end      
            if isempty(Data)
                vec = vec_;
                return
            end
            vec = conv(vec_, Data.kernel);
            vec = vec(1:length(vec_)-1);
        end
    end

    %% PROTECTED

    methods (Access = protected)
        function this = RadialArtery()
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
