classdef (Abstract) Ichise < handle & mlnest.KineticModel
    %% line1
    %  line2
    %  
    %  Created 13-Jan-2024 14:46:40 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlnest/src/+mlnest.
    %  Developed on Matlab 23.2.0.2459199 (R2023b) Update 5 for MACA64.  Copyright 2024 John J. Lee.
    
    properties
    end

    methods
    end

    methods (Static)
        function this = create(opts)
            arguments
                opts.input_func mlfourd.ImagingContext2 % interpolated, sec, -> Bq/mL
                opts.scanner mlfourd.ImagingContext2 % sampled, sec -> Bq/mL
                opts.scanner_sigma {mustBeNumeric} = 20000 % Bq/mL
                opts.model_kind {mustBeTextScalar} = "vasc"
                opts.tracer {mustBeTextScalar} = "Unknown"
            end

            switch opts.model_kind
                case "vasc"
                    this = mlnest.IchiseVasc();
                case "novasc"
                    this = mlnest.IchiseNoVasc();
                otherwise
                    error("mlnest:ValueError", stackstr())
            end

            this.input_func = opts.input_func;
            this.scanner = opts.scanner;
            this.times_sampled = opts.scanner.json_metadata.timesMid;
            M0 = max(double(this.scanner));

            j = opts.scanner.json_metadata;
            this.Data.M0 = M0;
            this.Data.model_kind = opts.model_kind;
            this.Data.taus = ascol(j.taus); % col
            this.Data.timesMid = ascol(j.timesMid); % col
            %this.Data.t = mlnest.KineticModel.Data2t(this.Data); % col, timesStart
            this.Data.tracer = opts.tracer;
            this.Data.scanner_sampled = ascol(double(opts.scanner));
            this.Data.plasma_interpolated = ascol(double(opts.input_func)); 

            % expected by MultiNest; copying for performance
            this.measurement = double(this.scanner);
            this.measurement_sigma = opts.scanner_sigma*ones(size(this.measurement));
        end

        function vec = apply_ichise(vec_, Data)
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
                vec = vec_;
                return
            end
            % timesMid = Data.timesMid;
            % taus = Data.taus;
            % times0 = timesMid - taus/2;
            % timesF = timesMid + taus/2;

            vec = vec_;
        end
    end

    %% PROTECTED
    
    methods (Access = protected)
        function this = Ichise()
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
