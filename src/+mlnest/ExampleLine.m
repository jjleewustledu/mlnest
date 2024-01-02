classdef ExampleLine < handle & mlsystem.IHandle
    %% line1
    %  line2
    %  
    %  Created 19-Dec-2023 14:35:40 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlnest/src/+mlnest.
    %  Developed on Matlab 23.2.0.2459199 (R2023b) Update 5 for MACA64.  Copyright 2023 John J. Lee.
    
    properties
        Data
        measurement
        measurement_sigma
        times_sampled
    end

    properties (Dependent)
        preferred_names
    end

    methods %% GET
        function g = get.preferred_names(~)
            g = ["m"; "b"];
        end
    end

    methods
        function this = ExampleLine()
            this.Data = struct();

            data = readdata_line();
            this.times_sampled = data{1}(5:end);
            this.measurement = data{2}(5:end);
            this.measurement_sigma = data{3}(5:end);
        end
        function [s,s1] = simulate(this)
            arguments
                this mlnest.ExampleLine
            end
            
            s = [];
            s1 = s;
        end
    end

    methods (Static)
        function p = prior(~)
            p = {'m', 'uniform', 0, 10, 'fixed'; ...
                 'b', 'uniform', 0, 100, 'fixed'};
        end
        function s = signalmodel(x, parnames, parvals)
            s = line_model(x, parnames, parvals);
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
