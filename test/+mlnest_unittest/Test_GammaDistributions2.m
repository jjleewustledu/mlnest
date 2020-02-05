classdef Test_GammaDistributions2 < matlab.unittest.TestCase
	%% TEST_GAMMADISTRIBUTIONS2

	%  Usage:  >> results = run(mlnest_unittest.Test_GammaDistributions2)
 	%          >> result  = run(mlnest_unittest.Test_GammaDistributions2, 'test_dt')
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
            % Evidence:  ln(Z) = -56.593567 +/- 0.536939
            % Information:  H = 28.830401 nats = 41.593477 bits
            % 	a = 0.234955 +/- 0.000001
            % 	b = 6.957048 +/- 0.000005
            % 	p = 0.323471 +/- 0.000000
            % 	t0 = 9.988227 +/- 0.000008
            % 	w = 0.111852 +/- 0.000002
            % logL(k = 3000) = -27.763164
            % logWt(k = 3000) = -62.363330  
            % 54.1621 seconds testing time.
        end
        function test_runGenGammaDist(this)            
 			import mlnest.*;
            this.map = containers.Map;
            this.map('a')        = struct('min',   0.1,  'max',  10,    'init', 0.2);
            this.map('b')        = struct('min',   1,    'max',  10,    'init', 5);
            this.map('p')        = struct('min',   0.2,  'max',   2,    'init', 1/3);
            this.map('t0')       = struct('min',   0,    'max',  20,    'init', 10);
            %this.map('w')        = struct('min',  -1e-3, 'max',   1e-3, 'init', 0);
 			this.testObj = GammaDistributions2( ...
                'paramMap', this.map, ...
                'timeInterpolants', 0:1:99, ...
                'sigma0', 0.001, ...
                'modelName', 'GeneralizedGammaDistribution');
            Obj.a = 0.2;
            Obj.b = 5;
            Obj.p = 1/3;
            Obj.t0 = 10;
            %Obj.w = 0;
            this.testObj.Measurement = this.testObj_.estimatorGenGamma(Obj);
            this.testObj.run(this.testObj)
            
            % # iterates = 3000
            % Evidence:  ln(Z) = -52.165618 +/- 0.527786
            % Information:  H = 27.855811 nats = 40.187440 bits
            % 	a = 0.224177 +/- 0.000141
            % 	b = 6.068307 +/- 0.000471
            % 	p = 0.327203 +/- 0.000015
            % 	t0 = 9.994158 +/- 0.000320
            % logL(k = 3000) = -24.274314
            % logWt(k = 3000) = -58.874480  
            % 45.4871 seconds testing time.
        end
        function test_runGammaDistP(this)            
 			import mlnest.*;
            this.map = containers.Map;
            this.map('a')        = struct('min',   0.1,  'max',  5,    'init', 0.2*(1 + randn()));
            this.map('b')        = struct('min',   0.01, 'max',  5,    'init', 0.2*(1 + randn()));
            %this.map('p')        = struct('min',   0.2,  'max',   2,    'init', 1/3);
            this.map('t0')       = struct('min',   0,    'max',  20,    'init', 10*(1 + randn()));
            this.map('w')        = struct('min',  -1,    'max',   1,    'init', 0.1*(1 + randn()));
 			this.testObj = GammaDistributions2( ...
                'paramMap', this.map, ...
                'timeInterpolants', 0:1:99, ...
                'sigma0', 0.001, ...
                'modelName', 'GammaDistributionP');
            Obj.a = 0.2;
            Obj.b = 0.2;
            %Obj.p = 1/3;
            Obj.t0 = 10;
            Obj.w = 0.1;
            this.testObj.Measurement = this.testObj_.estimatorGammaP(Obj);
                        
            this.testObj.run(this.testObj)
                      
        end
        function test_runGammaDist(this)            
 			import mlnest.*;
            this.map = containers.Map;
            this.map('a')        = struct('min',   0.1,  'max',   5,    'init', 0.2*(1 + randn()));
            this.map('b')        = struct('min',   0.01, 'max',   5,    'init', 0.2*(1 + randn()));
            %this.map('p')        = struct('min',   0.2,  'max',   2,    'init', 1/3);
            this.map('t0')       = struct('min',   0,    'max',  20,    'init', 10*(1 + randn()));
            %this.map('w')        = struct('min',  -1,    'max',   1,    'init', 0.1);
 			this.testObj = GammaDistributions2( ...
                'paramMap', this.map, ...
                'timeInterpolants', 0:1:99, ...
                'sigma0', 0.001, ...
                'modelName', 'GammaDistribution');
            Obj.a = 0.2;
            Obj.b = 0.2;
            %Obj.p = 1/3;
            Obj.t0 = 10;
            %Obj.w = 0.1;
            this.testObj.Measurement = this.testObj_.estimatorGamma(Obj);
            
            this.testObj.run(this.testObj)
                      
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
            this.map('a')        = struct('min',   0.1,  'max',  10,    'init', 0.   *(1 + randn()));
            this.map('b')        = struct('min',   1,    'max',  10,    'init', 5    *(1 + randn()));
            this.map('p')        = struct('min',   0.2,  'max',   2,    'init', (1/3)*(1 + randn()));
            this.map('t0')       = struct('min',   0,    'max',  20,    'init', 10   *(1 + randn()));
            this.map('w')        = struct('min',  -1,    'max',   1,    'init', 0.1  *(1 + randn()));
 			this.testObj_ = GammaDistributions2( ...
                'paramMap', this.map, ...
                'timeInterpolants', 0:1:99, ...
                'sigma0', 1e-3, ...
                'modelName', 'GeneralizedGammaDistributionP');
            Obj.a = 0.2;
            Obj.b = 5;
            Obj.p = 1/3;
            Obj.t0 = 10;
            Obj.w = 0.1;
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

