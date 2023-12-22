classdef Test_MultiNest < matlab.unittest.TestCase
    %% line1
    %  line2
    %  
    %  Created 19-Dec-2023 16:16:53 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlnest/test/+mlnest_unittest.
    %  Developed on Matlab 23.2.0.2459199 (R2023b) Update 5 for MACA64.  Copyright 2023 John J. Lee.
    
    properties
        deriv_ho_pet_pth
        testObj
    end
    
    methods (Test)
        function test_afun(this)
            import mlnest.*
            this.assumeEqual(1,1);
            this.verifyEqual(1,1);
            this.assertEqual(1,1);
        end
        function test_idif(this)
            ic = mlfourd.ImagingContext2( ...
                fullfile(this.deriv_ho_pet_pth, ...
                "sub-108293_ses-20210421152358_trc-ho_proc-MipIdif-finite-deconv_idif.nii.gz"));
            %figure;
            plot(ic)
        end
        function test_ExampleArtery(this)
            tic
            ic = mlfourd.ImagingContext2( ...
                fullfile(this.deriv_ho_pet_pth, ...
                "sub-108293_ses-20210421152358_trc-ho_proc-MipIdif-finite-deconv_idif.nii.gz"));
            obj = mlnest.MultiNest(context=mlnest.ExampleArtery.create(artery=ic));
            obj = obj.solve(signal_model=@mlnest.ExampleArtery.signalmodel, Nlive=100); % 10 -> 0.8 sec; 100 -> 12 sec; 500 -> 892 sec
            disp(obj)
            disp(obj.ks())
            disp(obj.loss(signal_model=@mlnest.ExampleArtery.signalmodel))
            obj.plot_posteriors();
            ic = obj.simulate();
            plot(ic);
            toc
        end
        function test_ExampleLine(this)
            tic
            obj = this.testObj.solve(signal_model=@line_model, Nlive=100); % 100 -> 0.4 sec; 1000 -> 185 sec
            disp(obj)
            disp(obj.ks())
            disp(obj.loss(signal_model=@line_model))
            obj.plot_posteriors();
            toc
        end
    end
    
    methods (TestClassSetup)
        function setupMultiNest(this)
            import mlnest.*
            this.deriv_ho_pet_pth = fullfile( ...
                getenv("SINGULARITY_HOME"), "CCIR_01211", "derivatives", "sub-108293", "ses-20210421152358", "pet");
            this.testObj_ = MultiNest(context=mlnest.ExampleLine());
        end
    end
    
    methods (TestMethodSetup)
        function setupMultiNestTest(this)
            this.testObj = this.testObj_;
            this.addTeardown(@this.cleanTestMethod)
        end
    end
    
    properties (Access = private)
        testObj_
    end
    
    methods (Access = private)
        function cleanTestMethod(this)
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
