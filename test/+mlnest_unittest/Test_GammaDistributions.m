classdef Test_GammaDistributions < matlab.unittest.TestCase
	%% TEST_GAMMADISTRIBUTIONS 

	%  Usage:  >> results = run(mlnest_unittest.Test_GammaDistributions)
 	%          >> result  = run(mlnest_unittest.Test_GammaDistributions, 'test_dt')
 	%  See also:  file:///Applications/Developer/MATLAB_R2014b.app/help/matlab/matlab-unit-test-framework.html

	%  $Revision$
 	%  was created 29-Jan-2020 17:57:40 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlnest/test/+mlnest_unittest.
 	%% It was developed on Matlab 9.7.0.1261785 (R2019b) Update 3 for MACI64.  Copyright 2020 John Joowon Lee.
 	
	properties
        map
 		registry
 		testObj
 	end

	methods (Test)
		function test_afun(this)
 			import mlnest.*;
 			this.assumeEqual(1,1);
 			this.verifyEqual(1,1);
 			this.assertEqual(1,1);
        end
        function test_measurement(this)
            figure
            plot(this.testObj_.Measurement)
        end
        function test_run(this)
            this.testObj_.run(this.testObj_)
        end
	end

 	methods (TestClassSetup)
		function setupGammaDistributions(this)
 			import mlnest.*;
            this.map = containers.Map;
            this.map('a')        = struct('min',   0.1,  'max',  10);
            this.map('b')        = struct('min',   0.01, 'max',  20);
            this.map('p')        = struct('min',   0.2,  'max',   2);
            this.map('t0')       = struct('min', -20,    'max',  20);
            this.map('w')        = struct('min',  -1,    'max',   1);
 			this.testObj_ = GammaDistributions( ...
                'paramMap', this.map, ...
                'timeInterpolants', 0:1:99, ...
                'sigma0', 0.001, ...
                'modelName', 'GeneralizedGammaDistributionP');
            Obj.a = 0.2;
            Obj.b = 6;
            Obj.p = 1/3;
            Obj.t0 = 10;
            Obj.w = 0.05;
            this.testObj_.Measurement = this.testObj_.estimatorGenGammaP(Obj);
 		end
	end

 	methods (TestMethodSetup)
		function setupGammaDistributionsTest(this)
 			this.testObj = this.testObj_;
 			this.addTeardown(@this.cleanTestMethod);
            rng('default')
 		end
	end

	properties (Access = private)
 		testObj_
 	end

	methods (Access = private)
		function cleanTestMethod(this)
 		end
	end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

