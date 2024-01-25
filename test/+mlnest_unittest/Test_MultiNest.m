classdef Test_MultiNest < matlab.unittest.TestCase
    %% line1
    %  line2
    %  
    %  Created 19-Dec-2023 16:16:53 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlnest/test/+mlnest_unittest.
    %  Developed on Matlab 23.2.0.2459199 (R2023b) Update 5 for MACA64.  Copyright 2023 John J. Lee.
    
    properties
        deriv_ho_pet_pth
        deriv_oo_pet_pth
        testObj
    end
    
    methods (Test)
        function test_afun(this)
            import mlnest.*
            this.assumeEqual(1,1);
            this.verifyEqual(1,1);
            this.assertEqual(1,1);
        end
        function test_ctor(this)
            disp(this.testObj)
        end
        function test_idif(this)

            return

            ic = mlfourd.ImagingContext2( ...
                fullfile(this.deriv_oo_pet_pth, ...
                "sub-108293_ses-20210421150523_trc-oo_proc-MipIdif_idif.nii.gz"));
            %figure;
            disp(ic)
            disp(ic.json_metadata)
            plot(ic)
        end
        function test_deconv_idif(this)

            return

            ic = mlfourd.ImagingContext2( ...
                fullfile(this.deriv_ho_pet_pth, ...
                "sub-108293_ses-20210421152358_trc-ho_proc-MipIdif-finite-deconv_idif.nii.gz"));
            %figure;
            plot(ic)
        end
        function test_Boxcar_ho(this)
            ic = mlfourd.ImagingContext2( ...
                fullfile(this.deriv_ho_pet_pth, ...
                "sub-108293_ses-20210421152358_trc-ho_proc-MipIdif_idif.nii.gz"));
            boxcar = mlnest.Boxcar.create(artery=ic, model_kind="3bolus");
            cd(ic.filepath)

            tic % <--------------

            obj = mlnest.MultiNest(context=boxcar);
            obj.filepath = ic.filepath;
            obj = obj.solve( ...
                signal_model=@boxcar.signalmodel, ...
                verbose=false, ...
                Nlive=200, ...
                Nmcmc=0); 

            % 3bolus logZ ~ -1295
            % 4bolus logZ ~ -1638

            disp(obj)
            fprintf("Multinest.product: \n"); 
            disp(obj.product)
            fprintf("ks(): %g\n", obj.ks())
            fprintf("loss: %s\n", obj.loss())
            obj.plot_posteriors(singles=true);
            obj.save();

            toc % <-------------- Nlive = 55 -> 245 s for slide_slow; 115-158 s for slide_fast; 926 s for slide_int
                %                 Nlive = 200 >> 12 h

            [sig,idl] = boxcar.simulate(obj.product);
            figure;
            plot_over_figure(boxcar.artery, 'o', MarkerSize=12); hold on;
            plot_over_figure(sig, '--', LineWidth=2);
            plot_over_figure(idl, '-', LineWidth=1.5); hold off;
            ylabel("activity (Bq/mL)")
            xlabel("time (s)")
            saveFigures(stackstr+"_Nlive=200_Nmcmc=0", closeFigure=false);
            close('all');            
        end
        function test_Boxcar_oo(this)
            tic            
            ic = mlfourd.ImagingContext2( ...
                fullfile(this.deriv_oo_pet_pth, ...
                "sub-108293_ses-20210421150523_trc-oo_proc-MipIdif_idif.nii.gz"));

            cd(ic.filepath)

            boxcar = mlnest.Boxcar.create(artery=ic, model_kind="4bolus");
            obj = mlnest.MultiNest(context=boxcar);
            obj.filepath = ic.filepath;
            obj = obj.solve( ...
                signal_model=@boxcar.signalmodel, ...
                verbose=false, ...
                Nlive=55, ...
                Nmcmc=0);                   

            % 4bolus logZ ~ -912
            % 3bolus logZ ~ -985

            disp(obj)
            fprintf("Multinest.product: \n"); 
            disp(obj.product)
            fprintf("ks(): %g\n", obj.ks())
            fprintf("loss: %s\n", obj.loss())
            obj.plot_posteriors(singles=true);
            obj.save();

            [sig,idl] = boxcar.simulate(obj.product);
            figure;
            plot_over_figure(boxcar.artery, 'o', MarkerSize=12); hold on;
            plot_over_figure(sig, '--', LineWidth=2);
            plot_over_figure(idl, '-', LineWidth=1.5); hold off;
            ylabel("activity (Bq/mL)")
            xlabel("time (s)")
            saveFigures(stackstr+"_Nlive=55_Nmcmc=0", closeFigure=false);
            close('all');
            toc
        end
        function test_Boxcar_oo_mcmc(this)
            tic
            ic = mlfourd.ImagingContext2( ...
                fullfile(this.deriv_oo_pet_pth, ...
                "sub-108293_ses-20210421150523_trc-oo_proc-MipIdif_idif.nii.gz"));

            filepath = fullfile(ic.filepath, stackstr());
            ensuredir(filepath);
            cd(filepath);

            boxcar = mlnest.Boxcar.create(artery=ic, model_kind="4bolus");
            obj = mlnest.MultiNest(context=boxcar);
            obj.filepath = filepath;
            obj = obj.solve( ...
                signal_model=@boxcar.signalmodel, ...
                verbose=false, ...
                Nlive=800, ...
                Nmcmc=3200); 
            
            % Nlive:  25 -> 2.8 sec; 50 -> 45 sec; 55 -> 131 sec; 60 -> __
            % sec; 70 -> __ sec; 100 ~ 3000 sec; 200 -> Inf sec;
            % after ~40h, tol~4, tolerance=0.1, j=3180, Nlive=20, forced
            % stop.  nested_sampler() main loops ~ __ min.

            % (Nlive,Nmcmc) ~ (25,25) -> 2.0 sec;
            %               ~ (25,50) -> 1.4 sec;
            %               ~ (25,100) -> 1.5 sec;
            %               ~ (25,200) -> 1.5 sec;
            %               ~ (25,400) -> 2.8 sec;
            %               ~ (25,800) -> 4.8 sec;
            %               ~ (25,1600) -> 7.2 sec;
            %               ~ (25,3200) -> 13 sec;
            %               ~ (25,6400) -> 24 sec;
            %               ~ (50,25) -> 1.6 sec;
            %               ~ (50,50) -> 1.8 sec;
            %               ~ (50,100) -> 2.3 sec;
            %               ~ (50,200) -> 3.0 sec;
            %               ~ (50,400) -> 3.0 sec;
            %               ~ (50,800) -> 7.7 sec;
            %               ~ (100,50) -> 2.9 sec;
            %               ~ (100,100) -> 3.3 sec;
            %               ~ (100,200) -> 5.0 sec;
            %               ~ (100,400) -> 7.4 sec;
            %               ~ (100,800) -> 13.8 sec;
            %               ~ (200,100) -> 6.0 sec;
            %               ~ (200,200) -> 8.8 sec;
            %               ~ (200,400) -> 14.7 sec;
            %               ~ (200,800) -> 25.4 sec;
            %               ~ (200,1600) -> 48 sec;
            %               ~ (400,200) -> 17.0 sec; too few chains
            %               ~ (400,400) -> 27.2 sec; *
            %               ~ (400,800) -> 49 sec; 
            %               ~ (400,1600) -> 94 sec; overly correlated samples?
            %               ~ (400,3200) -> 191 sec; overly correlated samples?
            %               ~ (800,400) -> 59 sec; overly correlated samples?
            %               ~ (800,800) -> 99 sec; overly correlated samples?
            %               ~ (800,1600) -> 192 sec;
            %               ~ (800,3200) -> 370 sec;
                        
            disp(obj)
            fprintf("Multinest.product: \n"); 
            disp(obj.product)
            fprintf("ks(): %g\n", obj.ks())
            fprintf("loss: %s\n", obj.loss())
            obj.plot_posteriors(singles=true);
            obj.save();

            [sig,idl] = boxcar.simulate(obj.product);
            figure;
            plot_over_figure(boxcar.artery, 'o', MarkerSize=12); hold on;
            plot_over_figure(sig, '--', LineWidth=2);
            plot_over_figure(idl, '-', LineWidth=1.5); hold off;
            ylabel("activity (Bq/mL)")
            xlabel("time (s)")
            saveFigures(stackstr+"_Nlive=800_Nmcmc=3200")
            toc
        end
        function test_RadialArtery_oo_mcmc(this)
            tic
            ifc = mlfourd.ImagingFormatContext2( ...
                fullfile(getenv("HOME"), ...
                "MATLAB-Drive", "mlkinetics", "data", "sourcedata", "sub-108293", 'ses-20210421', "pet", ...
                "sub-108293_ses-20210421150523_trc-oo_proc-TwiliteKit-do-make-input-func-nomodel_inputfunc.nii.gz"));  

            baseline = mean(ifc.json_metadata.baselineActivityDensity);
            ifc.img = ifc.img - baseline;
            ifc.img(ifc.img < 0) = 0;
            ic = mlfourd.ImagingContext2(ifc);

            filepath = fullfile(ic.filepath, stackstr());
            ensuredir(filepath);
            cd(filepath);

            radart = mlnest.RadialArtery.create(artery=ic, model_kind="4bolus");
            obj = mlnest.MultiNest(context=radart);
            obj.filepath = filepath;
            obj = obj.solve( ...
                signal_model=@radart.signalmodel, ...
                verbose=false, ...
                Nlive=800, ...
                Nmcmc=3200); 

            boxcar = mlnest.RadialArtery.create(artery=ic, model_kind="4bolus");
            obj = mlnest.MultiNest(context=boxcar);
            obj.filepath = filepath;
            obj = obj.solve( ...
                signal_model=@boxcar.signalmodel, ...
                verbose=false, ...
                Nlive=800, ...
                Nmcmc=3200); 
            
            % Nlive:  25 -> 2.8 sec; 50 -> 45 sec; 55 -> 131 sec; 60 -> __
            % sec; 70 -> __ sec; 100 ~ 3000 sec; 200 -> Inf sec;
            % after ~40h, tol~4, tolerance=0.1, j=3180, Nlive=20, forced
            % stop.  nested_sampler() main loops ~ __ min.

            % (Nlive,Nmcmc) ~ (25,25) -> 2.0 sec;
            %               ~ (25,50) -> 1.4 sec;
            %               ~ (25,100) -> 1.5 sec;
            %               ~ (25,200) -> 1.5 sec;
            %               ~ (25,400) -> 2.8 sec;
            %               ~ (25,800) -> 4.8 sec;
            %               ~ (25,1600) -> 7.2 sec;
            %               ~ (25,3200) -> 13 sec;
            %               ~ (25,6400) -> 24 sec;
            %               ~ (50,25) -> 1.6 sec;
            %               ~ (50,50) -> 1.8 sec;
            %               ~ (50,100) -> 2.3 sec;
            %               ~ (50,200) -> 3.0 sec;
            %               ~ (50,400) -> 3.0 sec;
            %               ~ (50,800) -> 7.7 sec;
            %               ~ (100,50) -> 2.9 sec;
            %               ~ (100,100) -> 3.3 sec;
            %               ~ (100,200) -> 5.0 sec;
            %               ~ (100,400) -> 7.4 sec;
            %               ~ (100,800) -> 13.8 sec;
            %               ~ (200,100) -> 6.0 sec;
            %               ~ (200,200) -> 8.8 sec;
            %               ~ (200,400) -> 14.7 sec;
            %               ~ (200,800) -> 25.4 sec;
            %               ~ (200,1600) -> 48 sec;
            %               ~ (400,200) -> 17.0 sec; too few chains
            %               ~ (400,400) -> 27.2 sec; *
            %               ~ (400,800) -> 49 sec; 
            %               ~ (400,1600) -> 94 sec; overly correlated samples?
            %               ~ (400,3200) -> 191 sec; overly correlated samples?
            %               ~ (800,400) -> 59 sec; overly correlated samples?
            %               ~ (800,800) -> 99 sec; overly correlated samples?
            %               ~ (800,1600) -> 192 sec;
            %               ~ (800,3200) -> 370 sec;
                        
            disp(obj)
            fprintf("Multinest.product: \n"); 
            disp(obj.product)
            fprintf("ks(): %g\n", obj.ks())
            fprintf("loss: %s\n", obj.loss())
            obj.plot_posteriors(singles=true);
            obj.save();

            [sig,idl] = boxcar.simulate(obj.product);
            figure;
            plot_over_figure(boxcar.artery, 'o', MarkerSize=12); hold on;
            plot_over_figure(sig, '--', LineWidth=2);
            plot_over_figure(idl, '-', LineWidth=1.5); hold off;
            ylabel("activity (Bq/mL)")
            xlabel("time (s)")
            saveFigures(stackstr+"_Nlive=800_Nmcmc=3200")
            toc
        end
        function test_RadialArtery_ho(this)
            tic
            ifc = mlfourd.ImagingFormatContext2( ...
                fullfile(getenv("HOME"), ...
                "MATLAB-Drive", "mlkinetics", "data", "sourcedata", "sub-108293", 'ses-20210421', "pet", ...
                "sub-108293_ses-20210421152358_trc-ho_proc-TwiliteKit-do-make-input-func-nomodel_inputfunc.nii.gz")); 

            cd(ifc.filepath)

            baseline = mean(ifc.json_metadata.baselineActivityDensity);
            ifc.img = ifc.img - baseline;
            ifc.img(ifc.img < 0) = 0;
            ic = mlfourd.ImagingContext2(ifc);

            radart = mlnest.RadialArtery.create(artery=ic, model_kind="3bolus");
            obj = mlnest.MultiNest(context=radart);
            obj.filepath = ifc.filepath;
            obj = obj.solve( ...
                signal_model=@radart.signalmodel, ...
                verbose=false, ...
                Nlive=55, ...
                Nmcmc=0); 

            % 4bolus logZ ~ -3364
            % 3bolus logZ ~ -3319

            disp(obj)
            fprintf("Multinest.product: \n"); 
            disp(obj.product)
            fprintf("ks(): %g\n", obj.ks())
            fprintf("loss: %s\n", obj.loss())
            obj.plot_posteriors(singles=true);
            obj.save();

            [sig,idl] = radart.simulate(obj.product);
            figure;
            plot_over_figure(radart.artery, 'o', MarkerSize=12); hold on;
            plot_over_figure(sig, '--', LineWidth=2);
            plot_over_figure(idl, '-', LineWidth=1.5); hold off;
            ylabel("activity (Bq/mL)")
            xlabel("time (s)")
            saveFigures(stackstr+"_Nlive=55_Nmcmc=0", closeFigure=false);
            close('all');
            toc
        end
        function test_RadialArtery_oo(this)
            tic
            ifc = mlfourd.ImagingFormatContext2( ...
                fullfile(getenv("HOME"), ...
                "MATLAB-Drive", "mlkinetics", "data", "sourcedata", "sub-108293", 'ses-20210421', "pet", ...
                "sub-108293_ses-20210421150523_trc-oo_proc-TwiliteKit-do-make-input-func-nomodel_inputfunc.nii.gz"));  

            cd(ifc.filepath)

            baseline = mean(ifc.json_metadata.baselineActivityDensity);
            ifc.img = ifc.img - baseline;
            ifc.img(ifc.img < 0) = 0;
            ic = mlfourd.ImagingContext2(ifc);

            radart = mlnest.RadialArtery.create(artery=ic, model_kind="4bolus");
            obj = mlnest.MultiNest(context=radart);
            obj.filepath = ifc.filepath;
            obj = obj.solve( ...
                signal_model=@radart.signalmodel, ...
                verbose=false, ...
                Nlive=55, ...
                Nmcmc=0); 

            % 4bolus logZ ~ -6651; Nive=55, Nmcmc=0, elapsed time ~ 1012 sec
            % 3bolus logZ ~ -35302

            disp(obj)
            fprintf("Multinest.product: \n"); 
            disp(obj.product)
            fprintf("ks(): %g\n", obj.ks())
            fprintf("loss: %s\n", obj.loss())
            obj.plot_posteriors(singles=true);
            obj.save();

            [sig,idl] = radart.simulate(obj.product);
            figure;
            plot_over_figure(radart.artery, 'o', MarkerSize=12); hold on;
            plot_over_figure(sig, '--', LineWidth=2);
            plot_over_figure(idl, '-', LineWidth=1.5); hold off;
            ylabel("activity (Bq/mL)")
            xlabel("time (s)")
            saveFigures(stackstr+"_Nlive=55_Nmcmc=0", closeFigure=false);
            close('all');
            toc
        end
        function test_ExampleArtery(this)

            return

            tic
            ic = mlfourd.ImagingContext2( ...
                fullfile(this.deriv_ho_pet_pth, ...
                "sub-108293_ses-20210421152358_trc-ho_proc-MipIdif-finite-deconv_idif.nii.gz"));
            obj = mlnest.MultiNest(context=mlnest.ExampleArtery.create(artery=ic));
            obj = obj.solve(signal_model=@mlnest.ExampleArtery.signalmodel, Nlive=50); % 10 -> 0.8 sec; 100 -> 12 sec; 500 -> 892 sec
            disp(obj)
            disp(obj.ks())
            disp(obj.loss())
            obj.plot_posteriors(singles=false);
            ic = obj.simulate();
            plot(ic);
            toc
        end
        function test_ExampleLine(this)

            return

            tic
            obj = this.testObj.solve(signal_model=@line_model, Nlive=100); % 100 -> 0.4 sec; 1000 -> 185 sec
            disp(obj)
            disp(obj.ks())
            disp(obj.loss())
            obj.plot_posteriors();
            toc
        end
        function test_RadialArtery(this)
        end
    end
    
    methods (TestClassSetup)
        function setupMultiNest(this)
            import mlnest.*
            this.deriv_ho_pet_pth = fullfile( ...
                getenv("SINGULARITY_HOME"), "CCIR_01211", "derivatives", "sub-108293", "ses-20210421152358", "pet");
            this.deriv_oo_pet_pth = fullfile( ...
                getenv("SINGULARITY_HOME"), "CCIR_01211", "derivatives", "sub-108293", "ses-20210421150523", "pet");
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
