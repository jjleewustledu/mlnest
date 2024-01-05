classdef (Abstract) Boxcar < handle & mlnest.Artery
    %% line1
    %  line2
    %  
    %  Created 22-Dec-2023 23:17:43 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlnest/src/+mlnest.
    %  Developed on Matlab 23.2.0.2459199 (R2023b) Update 5 for MACA64.  Copyright 2023 John J. Lee.
    
    properties
    end

    methods
    end

    methods (Static)
        function this = create(opts)
            arguments
                opts.artery mlfourd.ImagingContext2
                opts.measurement_sigma {mustBeNumeric} = 7867/sqrt(11) % CCIRRadMeas 20210421, for 0.05*218 voxels in centerline_on_pet
                opts.model_kind {mustBeTextScalar} = "4bolus"
                opts.tracer {mustBeTextScalar} = "Unknown"
            end

            switch opts.model_kind
                case "4bolus"
                    this = mlnest.Boxcar4bolus();
                case "3bolus"
                    this = mlnest.Boxcar3bolus();
                case "2bolus"
                    this = mlnest.Boxcar2bolus();
                otherwise
                    error("mlnest:ValueError", stackstr())
            end

            this.artery = mlpipeline.ImagingMediator.ensureFiniteImagingContext(opts.artery);
            %this.artery.save();
            ifc = this.artery.imagingFormat;
            this.times_sampled = asrow(ifc.json_metadata.timesMid);
            this.measurement = asrow(double(ifc.img));
            M0 = max(this.measurement);
            this.measurement_sigma = opts.measurement_sigma*ones(size(this.measurement));

            j = this.artery.json_metadata;
            this.Data.M0 = M0;
            this.Data.model_kind = opts.model_kind;
            this.Data.taus = j.taus;
            this.Data.timesMid = j.timesMid;
            this.Data.t = mlnest.Artery.Data2t(this.Data);
            this.Data.tracer = opts.tracer;
        end

        function vec = apply_boxcar(vec_, Data)
            %% Args:            
            %      vec double_, has 1 Hz sampling
            %      Data struct
            %  Returns:
            %      vec double, sampled at opts.Data.timesMid

            arguments
                vec_ double 
                Data struct = []
            end      
            if isempty(Data)
                return
            end
            timesMid = Data.timesMid;
            taus = Data.taus;
            times0 = timesMid - taus/2;
            timesF = timesMid + taus/2;

            vec_sampled = NaN(size(timesMid));
            for vi = 1:length(timesMid)
                s = times0(vi) + 1;
                s1 = timesF(vi);
                %s1 = min(timesF(vi), length(vec_));
                vec_sampled(vi) = mean(vec_(s:s1));
            end
            vec = vec_sampled;
        end
    end

    %% PROTECTED
    
    methods (Access = protected)
        function this = Boxcar()
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
