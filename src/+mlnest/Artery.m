classdef (Abstract) Artery < handle & mlsystem.IHandle
    %% line1
    %  line2
    %  
    %  Created 22-Dec-2023 23:26:53 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlnest/src/+mlnest.
    %  Developed on Matlab 23.2.0.2459199 (R2023b) Update 5 for MACA64.  Copyright 2023 John J. Lee.
        
    methods (Abstract, Static)
        create()
        signalmodel()
    end

    properties 
        artery
        Data        
        measurement
        measurement_sigma
        times_sampled
    end

    methods
        function [signal,ideal] = simulate(this, product, opts)
            %% returns ImagingContext2

            arguments
                this mlnest.Artery
                product struct
                opts.filepath {mustBeFolder} = this.artery.filepath
                opts.fileprefix {mustBeTextScalar} = ""
            end
            if isemptytext(opts.fileprefix)
                opts.fileprefix = mlpipeline.Bids.adjust_fileprefix( ...
                    this.artery.fileprefix, ...
                    post_proc=stackstr(3, use_dashes=true));
            end
            fp_signal = mlpipeline.Bids.adjust_fileprefix( ...
                opts.fileprefix, ...
                new_mode="signal");
            fp_ideal = mlpipeline.Bids.adjust_fileprefix( ...
                opts.fileprefix, ...
                new_mode="ideal");
            
            pr = this.prior(this.Data);
            parnames = [pr(:,1); fields(this.Data)];
            parvals = [num2cell(ascol(product.mean_ks)); struct2cell(this.Data)];
            [signal_numer,ideal_numer] = this.signalmodel(ascol(this.times_sampled), parnames, parvals);

            ic = copy(this.artery);

            ic.json_metadata.mean_ks = product.mean_ks;
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

            ideal = copy(ic);
            ideal = ideal.selectImagingTool(img=asrow(ideal_numer));
            ideal.json_metadata.taus = ones(size(this.Data2t(this.Data)));
            ideal.json_metadata.timesMid = this.Data2t(this.Data);
            ideal.json_metadata.times = this.Data2t(this.Data);
            ideal.filepath = opts.filepath;
            ideal.fileprefix = fp_ideal;
        end
    end

    methods (Static)
        function ds = datestr()
            ds = string(datetime("now", Format="yyyyMMddHHmmss"));
        end 

        function t = Data2t(Data)
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
        function p = prior_from_mlaif(~)
            % no peaks to distributions for k2, k9-12

            pm = mlaif.ArteryLee2021Model.preferredMap();
            keys = ascol(pm.keys);
            lkeys = length(keys);
            mins = cell(size(keys));
            maxs = cell(size(keys));
            for kidx = 1:lkeys
                mins{kidx} = pm(keys{kidx}).min;
                maxs{kidx} = pm(keys{kidx}).max;
            end
            types = repmat({'uniform'}, [lkeys, 1]);
            behaviors = repmat({'fixed'}, [lkeys, 1]);

            p = [keys, types, mins, maxs, behaviors];
            p = natsortrows(p, [], logical([1,0,0,0,0]));
        end   
        function conc = slide(conc, t, Dt)
            %% SLIDE slides discretized function conc(t) to conc(t - Dt);
            %  Dt > 0 will slide conc(t) towards later times t.
            %  Dt < 0 will slide conc(t) towards earlier times t.
            %  It works for inhomogeneous t according to the ability of interp1 to interpolate.
            %  It may not preserve information according to the Nyquist-Shannon theorem.  
            
            import mlnest.Artery.ensureRow;
            
            [conc,trans] = ensureRow(conc);
            t            = ensureRow(t);
            
            tspan = t(end) - t(1);
            tinc  = t(2) - t(1);
            t_    = [(t - tspan - tinc) t];   % prepend times
            conc_ = [zeros(size(conc)) conc]; % prepend zeros
            conc_(isnan(conc_)) = 0;
            conc  = interp1(t_, conc_, t - Dt); % interpolate onto t shifted by Dt; Dt > 0 shifts to right
            
            if (trans)
                conc = conc';
            end
        end
    end

    %% PROTECTED

    methods (Access = protected)
        function this = Artery()
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
