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
            
            % # iterates = 3000
            % Evidence:  ln(Z) = -70.164236 +/- 0.523950
            % Information:  H = 27.452321 nats = 39.605328 bits
            % 	a = 0.218827 +/- 0.000067
            % 	b = 5.645965 +/- 0.000131
            % 	p = 0.329690 +/- 0.000002
            % 	t0 = 10.000007 +/- 0.001019
            % 	w = -0.000020 +/- 0.000000
            % logL(k = 3000) = -42.709194
            % logWt(k = 3000) = -77.309360            
        end
        
        function test_Obj2uniform2native(this)
            o.logL = -10;
            o.logWt = -1;
            o.a = 10;
            o.b = 10;
            o.p = 2;
            o.t0 = 0;
            o.w = -1e-3;

            u = this.testObj.Obj2uniform(o);
            this.verifyEqual(u.a, 1)
            this.verifyEqual(u.b, 1)
            this.verifyEqual(u.p, 1)
            this.verifyEqual(u.t0, 0)
            this.verifyEqual(u.w, 0)

            n = this.testObj.Obj2native(u);
            this.verifyEqual(n.a, o.a)
            this.verifyEqual(n.b, o.b)
            this.verifyEqual(n.p, o.p)
            this.verifyEqual(n.t0, o.t0)
            this.verifyEqual(n.w, o.w)
        end
        function test_Obj2vec(this)
            o.logL = -10;
            o.logWt = -1;
            o.a = 10;
            o.b = 20;
            o.p = 2;
            o.t0 = -20;
            o.w = -1;
            
            disp(fields(o))
            
            v = this.testObj.Obj2vec(o);
            this.verifyEqual(v, [-10 -1 10 20 2 -20 -1])
        end
	end

 	methods (TestClassSetup)
		function setupGammaDistributions(this)
 			import mlnest.*;
            this.map = containers.Map;
            this.map('a')        = struct('min',   0.1,  'max',  10,    'init', 0.2);
            this.map('b')        = struct('min',   1,    'max',  10,    'init', 5);
            this.map('p')        = struct('min',   0.2,  'max',   2,    'init', 1/3);
            this.map('t0')       = struct('min',   0,    'max',  20,    'init', 10);
            this.map('w')        = struct('min',  -1e-3, 'max',   1e-3, 'init', 0);
 			this.testObj_ = GammaDistributions( ...
                'paramMap', this.map, ...
                'timeInterpolants', 0:1:99, ...
                'sigma0', 0.001, ...
                'modelName', 'GeneralizedGammaDistributionP');
            Obj.a = 0.2;
            Obj.b = 5;
            Obj.p = 1/3;
            Obj.t0 = 10;
            Obj.w = 0;
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

