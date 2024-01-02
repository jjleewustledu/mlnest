classdef RadialArtery4bolus < handle & mlnest.RadialArtery
    %% line1
    %  line2
    %  
    %  Created 28-Dec-2023 13:28:56 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlnest/src/+mlnest.
    %  Developed on Matlab 23.2.0.2459199 (R2023b) Update 5 for MACA64.  Copyright 2023 John J. Lee.
    
    properties (Dependent)
        preferred_names
    end

    methods %% GET
        function g = get.preferred_names(~)
            g = ["\alpha-1"; ...
                "1/\beta"; ...
                "p"; ...
                "\delta p"; ...
                "t_0"; ...
                "frac. ss"; ...
                "frac. bolus_2"; ...
                "delay bolus_2"; ...
                "frac. recirc."; ...
                "delay recirc."; ...
                "amplitude"; ...
                "1/\gamma"];
        end
    end

    methods
        function this = RadialArtery4bolus(varargin)
            this = this@mlnest.RadialArtery(varargin{:});
        end
    end

    methods (Static)
        function p = prior(~)
            %% no peaks to distributions for k2, k9-12

            p(1,:) = {'k1', 'uniform', eps, 10, 'fixed'}; % \alpha-1 for OO
            p(2,:) = {'k2', 'jeffreys', 0.1, 10, 'fixed'}; % 1/\beta for OO
            p(3,:) = {'k3', 'uniform', 0.5, 3, 'fixed'}; % p
            p(4,:) = {'k4', 'uniform', -1, 0, 'fixed'}; % \delta p
            p(5,:) = {'k5', 'jeffreys', 5, 30, 'fixed'}; % t_0
            p(6,:) = {'k6', 'uniform', 0, 0.2, 'fixed'}; % frac. ss
            p(7,:) = {'k7', 'uniform', 0.2, 0.8, 'fixed'}; % frac. bolus_2
            p(8,:) = {'k8', 'jeffreys', 5, 20, 'fixed'}; % delay bolus_2
            p(9,:) = {'k9', 'uniform', 0.2, 0.4, 'fixed'}; % frac. recirc.
            p(10,:) = {'k10', 'jeffreys', 5, 20, 'fixed'}; % delay recirc.
            p(11,:) = {'k11', 'uniform', 0.8, 1.2, 'fixed'}; % amplitude
            p(12,:) = {'k12', 'jeffreys', 5, 100, 'fixed'}; % 1/\gamma
        end
        function [signal,ideal] = signalmodel(times_sampled, parnames, parvals)
            %% requires removal of baseline

            ks = zeros(1, 12);
            ks(1) = cell2mat(parvals(strcmp(parnames, 'k1')));
            ks(2) = cell2mat(parvals(strcmp(parnames, 'k2')));
            ks(3) = cell2mat(parvals(strcmp(parnames, 'k3'))); % stretching exp p
            ks(4) = cell2mat(parvals(strcmp(parnames, 'k4')));
            ks(5) = cell2mat(parvals(strcmp(parnames, 'k5')));
            ks(6) = cell2mat(parvals(strcmp(parnames, 'k6')));
            ks(7) = cell2mat(parvals(strcmp(parnames, 'k7')));
            ks(8) = cell2mat(parvals(strcmp(parnames, 'k8')));
            ks(9) = cell2mat(parvals(strcmp(parnames, 'k9')));
            ks(10) = cell2mat(parvals(strcmp(parnames, 'k10')));
            ks(11) = cell2mat(parvals(strcmp(parnames, 'k11')));
            ks(12) = cell2mat(parvals(strcmp(parnames, 'k12')));

            % from Data ~ extra params 
            Data.M0 = cell2mat(parvals(strcmp(parnames, 'M0')));
            Data.kernel = cell2mat(parvals(strcmp(parnames, 'kernel'))); % needed by apply_dispersion()
            Data.model_kind = cell2str(parvals(strcmp(parnames, 'model_kind'))); % needed by solution()
            Data.t = ascol(cell2mat(parvals(strcmp(parnames, 't')))); % needed by solution*() 

            % similar to mlaif.RadialArteryLee2024Model.sampled(), followed by 
            % mlpet.TracerSimulAnneal.plot()
            A = ks(11);
            qs_ = A*mlnest.RadialArtery4bolus.solution(ks, Data); % \in [0,1], 1 Hz sampling
            qs = mlnest.RadialArtery.apply_dispersion(qs_, Data);
            A_qs = 1/max(qs);            
            signal = Data.M0*A_qs*qs; assert(length(signal) == length(times_sampled))
            ideal = Data.M0*A_qs*ascol(qs_); assert(length(ideal) == length(Data.t))
        end
        function qs = solution(ks, Data)
            %% returns col vec required by MultiNest.

            % 2x stretched gamma distribution + recirc stretched gamma distribution + rising steadystate; 
            %  forcing p2 = p - dp2 < p, to be more dispersive

            import mlnest.RadialArtery4bolus.solution_3bolus;
            import mlnest.Artery.slide;
            t = Data.t;
            bolus2_frac = ks(7);
            bolus2_delay = ks(8);
            
            qs1 = solution_3bolus(ks, Data);
            qs1 = qs1/max(qs1);
            qs2 = solution_3bolus(ks, Data);
            qs2 = qs2/max(qs2);   
            qs2 = slide(qs2, t, bolus2_delay);  

            qs = (1 - bolus2_frac)*qs1 + bolus2_frac*qs2;
            qs = qs/max(qs); % \in [0 1] 
        end
        function qs = solution_3bolus(ks, Data)
            %% stretched gamma distribution + recirc stretched gamma distribution + rising steadystate; 
            %  forcing p2 = p - dp2 < p, to be more dispersive

            import mlnest.RadialArtery4bolus.prior
            import mlnest.RadialArtery4bolus.solution_1bolus
            import mlnest.RadialArtery4bolus.solution_2bolus
            import mlnest.Artery.slide

            t = Data.t;
            p1 = ks(3);
            prior_ = prior(Data);
            p2 = max(prior_{3,3}, ks(3) + ks(4));
            recirc_frac = ks(9);
            recirc_delay = ks(10);
            
            qs1 = solution_2bolus(ks, Data, p1);
            qs_recirc = solution_1bolus(ks, Data, p2);
            qs_recirc = slide(qs_recirc, t, recirc_delay);
            qs = (1 - recirc_frac)*qs1 + recirc_frac*qs_recirc;
            qs = qs/max(qs); % \in [0 1] 
        end
        function qs = solution_2bolus(ks, Data, p)
            %% stretched gamma distribution + rising steadystate

            import mlnest.Artery.slide
            t = Data.t;
            t0 = ks(5);
            a = ks(1);
            b = 1/ks(2);
            g = 1/ks(12);
            ss_frac = ks(6);
            
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
        function qs = solution_1bolus(ks, Data, p)
            %% stretched gamma distribution

            import mlnest.Artery.slide
            t = Data.t;
            t0 = ks(5);
            a = ks(1);
            b = 1/ks(2);
            
            if (t(1) >= t0) 
                t_ = t - t0;
                qs = t_.^a .* exp(-(b*t_).^p);
            else % k is complex for t - t0 < 0
                t_ = t - t(1);
                qs = t_.^a .* exp(-(b*t_).^p);
                qs = slide(qs, t, t0 - t(1));
            end
            
            qs = real(qs);
            qs(qs < 0) = 0;
            qs = qs/max(qs);
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
