classdef IchiseVasc < handle & mlnest.Ichise
    %% line1
    %  line2
    %  
    %  Created 13-Jan-2024 15:07:54 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlnest/src/+mlnest.
    %  Developed on Matlab 23.2.0.2459199 (R2023b) Update 5 for MACA64.  Copyright 2024 John J. Lee.
    
    properties (Dependent)
        preferred_names
    end

    methods %% GET
        function g = get.preferred_names(~)
            g = ["1/K_1"; ...
                "1/k_2"; ...
                "1/k_3"; ...
                "1/k_4"; ...
                "V_P"; ...
                "V_N + V_S"; ...
                "amplitude"];
        end
    end

    methods
        function this = IchiseVasc(varargin)
            this = this@mlnest.Ichise(varargin{:});
        end
    end
    
    methods (Static)
        function p = prior(~)
            % no blocking
            % p(1,:) = {'k1', 'uniform', 0.05, 20, 'fixed'}; % 1/K_1
            % p(2,:) = {'k2', 'uniform', 20, 2e4, 'fixed'}; % 1/k_2
            % p(3,:) = {'k3', 'uniform', 20, 2e4, 'fixed'}; % 1/k_3
            % p(4,:) = {'k4', 'uniform', 20, 2e4, 'fixed'}; % 1/k_4
            % p(5,:) = {'k5', 'uniform', 1, 10, 'fixed'}; % V_P
            % p(6,:) = {'k6', 'uniform', 10, 100, 'fixed'}; % V_N + V_S
            % p(7,:) = {'k7', 'uniform', 0.8, 1.2, 'fixed'}; % amplitude

            % blocking
            p(1,:) = {'k1', 'uniform', 0.05, 20, 'fixed'}; % 1/K_1
            p(2,:) = {'k2', 'uniform', 5, 1e3, 'fixed'}; % 1/k_2
            p(3,:) = {'k3', 'uniform', 5, 1e3, 'fixed'}; % 1/k_3
            p(4,:) = {'k4', 'uniform', 5, 1e3, 'fixed'}; % 1/k_4
            p(5,:) = {'k5', 'uniform', 1, 10, 'fixed'}; % V_P
            p(6,:) = {'k6', 'uniform', 10, 100, 'fixed'}; % V_N + V_S
            p(7,:) = {'k7', 'uniform', 0.8, 1.2, 'fixed'}; % amplitude
        end
        function [signal,ideal] = signalmodel(times_sampled, parnames, parvals)
            %% returns numerical signal on times_sampled; ideal on Data.t

            ks = zeros(1, size(mlnest.IchiseVasc.prior(), 1));
            ks(1) = cell2mat(parvals(strcmp(parnames, 'k1')));
            ks(2) = cell2mat(parvals(strcmp(parnames, 'k2')));
            ks(3) = cell2mat(parvals(strcmp(parnames, 'k3'))); 
            ks(4) = cell2mat(parvals(strcmp(parnames, 'k4')));
            ks(5) = cell2mat(parvals(strcmp(parnames, 'k5')));
            ks(6) = cell2mat(parvals(strcmp(parnames, 'k6')));
            ks(7) = cell2mat(parvals(strcmp(parnames, 'k7')));

            % from Data ~ extra params 
            Data.M0 = cell2mat(parvals(strcmp(parnames, 'M0')));
            Data.model_kind = cell2str(parvals(strcmp(parnames, 'model_kind'))); % needed by solution()
            Data.taus = cell2mat(parvals(strcmp(parnames, 'taus'))); % needed by apply_boxcar()
            Data.timesMid = cell2mat(parvals(strcmp(parnames, 'timesMid'))); % needed by apply_boxcar()
            %Data.t = ascol(cell2mat(parvals(strcmp(parnames, 't')))); % needed by solution*() 
            Data.scanner_sampled = ascol(cell2mat(parvals(strcmp(parnames, 'scanner_sampled'))));
            Data.plasma_interpolated = ascol(cell2mat(parvals(strcmp(parnames, 'plasma_interpolated'))));

            % similar to mlaif.BoxcarModel.sampled(), followed by 
            % mlpet.TracerSimulAnneal.plot()
            CT = mlnest.IchiseVasc.solution(ks, Data); % \in [0,1], 1 Hz sampling
            signal = Data.M0*CT; assert(length(signal) == length(times_sampled))
            ideal = [];
        end
        function C_T = solution(ks, Data)
            %% k1 ~ K1, k5 ~ VP, regional plasma volume-of-distrib., k6 ~ VN + VS
            %  internally rescale to sec -> muCi/cc

            meas_samp = Data.scanner_sampled/37000;
            times_samp = Data.timesMid;
            plasma_interpolated = Data.plasma_interpolated/37000;
            Np = length(plasma_interpolated);
            dt = 1; % sec

            K1 = 1/ks(1);
            VP = ks(5);
            Vstar = sum(ks(5:6));
            g1 = Vstar/(ks(2)*ks(4));
            g2 = -1/(ks(2)*ks(4));
            g3 = -(1/ks(2) + 1/ks(3) + 1/ks(4));
            g4star = K1;
            g5 = VP;

            C_T = zeros(size(times_samp));
            for tidx = 2:length(times_samp)
                m_ = meas_samp(1:tidx);
                t_ = times_samp(1:tidx);
                times__ = (0:dt:min(times_samp(tidx), Np - 1))'; % integration interval
                plasma__ = plasma_interpolated(times__+1); % integration interval
                
                int3 = trapz(t_, m_); 
                int4 = trapz(times__, plasma__);
                int2 = 0.5*trapz(t_, cumtrapz(t_, m_));
                int1 = 0.5*trapz(times__, cumtrapz(times__, plasma__));
                
                C_T(tidx) = g1*int1 + g2*int2 + g3*int3 + g4star*int4 + ...
                    g5*plasma_interpolated(times__(end)+1);
            end
            C_T(C_T < 0) = 0;
            C_T = ks(7)*C_T/max(C_T); % rescale with amplitude
        end
    end

    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
