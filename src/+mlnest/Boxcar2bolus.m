classdef Boxcar2bolus < handle & mlnest.Boxcar
    %% line1
    %  line2
    %  
    %  Created 27-Dec-2023 23:08:34 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlnest/src/+mlnest.
    %  Developed on Matlab 23.2.0.2459199 (R2023b) Update 5 for MACA64.  Copyright 2023 John J. Lee.
    
    properties (Dependent)
        preferred_names
    end

    methods %% GET
        function g = get.preferred_names(~)
            g = ["\alpha-1"; ...
                "1/\beta"; ...
                "p"; ...
                "t_0"; ...
                "frac. ss"; ...
                "amplitude"; ...
                "1/\gamma"];
        end
    end

    methods
        function this = Boxcar2bolus(varargin)
            this = this@mlnest.Boxcar(varargin{:});
        end
    end

    methods (Static)
        function p = prior(~)
            %% no peaks to distributions for k2, k9-12

            p(1,:) = {'k1', 'uniform', eps, 10, 'fixed'}; % \alpha-1 for OO
            p(2,:) = {'k2', 'jeffreys', 0.1, 10, 'fixed'}; % 1/\beta for OO
            p(3,:) = {'k3', 'uniform', 0.5, 3, 'fixed'}; % p
            p(4,:) = {'k4', 'jeffreys', 5, 30, 'fixed'}; % t_0
            p(5,:) = {'k5', 'uniform', 0, 0.2, 'fixed'}; % frac. ss
            p(6,:) = {'k6', 'uniform', 1, 2, 'fixed'}; % amplitude
            p(7,:) = {'k7', 'jeffreys', 5, 100, 'fixed'}; % 1/\gamma
        end
        function [signal,ideal] = signalmodel(times_sampled, parnames, parvals)
            %% returns numerical signal on times_sampled; ideal on Data.t

            ks = zeros(1, 7);
            ks(1) = cell2mat(parvals(strcmp(parnames, 'k1')));
            ks(2) = cell2mat(parvals(strcmp(parnames, 'k2')));
            ks(3) = cell2mat(parvals(strcmp(parnames, 'k3'))); % stretching exp p
            ks(4) = cell2mat(parvals(strcmp(parnames, 'k4')));
            ks(5) = cell2mat(parvals(strcmp(parnames, 'k5')));
            ks(6) = cell2mat(parvals(strcmp(parnames, 'k6')));
            ks(7) = cell2mat(parvals(strcmp(parnames, 'k7')));

            % from Data ~ extra params 
            Data.M0 = cell2mat(parvals(strcmp(parnames, 'M0')));
            Data.model_kind = cell2str(parvals(strcmp(parnames, 'model_kind'))); % needed by solution()
            Data.taus = cell2mat(parvals(strcmp(parnames, 'taus'))); % needed by apply_boxcar()
            Data.timesMid = cell2mat(parvals(strcmp(parnames, 'timesMid'))); % needed by apply_boxcar()
            Data.t = ascol(cell2mat(parvals(strcmp(parnames, 't')))); % needed by solution*() 

            % similar to mlaif.BoxcarModel.sampled(), followed by 
            % mlpet.TracerSimulAnneal.plot()
            A = ks(6);
            qs_ = A*mlnest.Boxcar2bolus.solution(ks, Data); % \in [0,1], 1 Hz sampling
            qs = mlnest.Boxcar.apply_boxcar(qs_, Data); % sampled to Data.timesMid
            A_qs = 1/max(qs);
            signal = Data.M0*A_qs*qs; assert(length(signal) == length(times_sampled))
            ideal = Data.M0*A_qs*ascol(qs_); assert(length(ideal) == length(Data.t))
        end
        function qs = solution(ks, Data)
            %% returns col vec required by MultiNest.

            %% stretched gamma distribution + rising steadystate

            import mlnest.Artery.slide

            t = Data.t;
            t0 = ks(4);
            a = ks(1);
            b = 1/ks(2);
            g = 1/ks(7);
            ss_frac = ks(5);
            
            if (t(1) >= t0) 
                t_ = t - t0;
                k_ = t_.^a .* exp(-(b*t_).^p);
                k_ = k_/max(k_); % else overwhelmed by ss_
                ss_ = 1 - exp(-g*t_);
                qs = (1 - ss_frac)*k_ + ss_frac*ss_;
            else % k is complex for t - t0 < 0
                t_ = t - t(1);
                k_ = t_.^a .* exp(-(b*t_).^p);
                k_ = k_/max(k_); % else overwhelmed by ss_
                ss_ = 1 - exp(-g*t_);
                qs = (1 - ss_frac)*k_ + ss_frac*ss_;
                qs = slide(qs, t, t0 - t(1));
            end
            
            qs = real(qs);
            qs(qs < 0) = 0;
            qs = qs/max(qs);
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
