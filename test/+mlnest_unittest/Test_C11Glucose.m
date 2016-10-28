classdef Test_C11Glucose < matlab.unittest.TestCase
	%% TEST_C11GLUCOSE 

	%  Usage:  >> results = run(mlnest_unittest.Test_C11Glucose)
 	%          >> result  = run(mlnest_unittest.Test_C11Glucose, 'test_dt')
 	%  See also:  file:///Applications/Developer/MATLAB_R2014b.app/help/matlab/matlab-unit-test-framework.html

	%  $Revision$
 	%  was created 20-Oct-2016 21:18:44
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlnest/test/+mlnest_unittest.
 	%% It was developed on Matlab 9.0.0.341360 (R2016a) for MACI64.
 	

	properties
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
	end

 	methods (TestClassSetup)
		function setupC11Glucose(this)
 			import mlnest.*;
 			this.testObj_ = C11Glucose;
 		end
	end

 	methods (TestMethodSetup)
		function setupC11GlucoseTest(this)
 			this.testObj = this.testObj_;
 			this.addTeardown(@this.cleanFiles);
 		end
	end

	properties (Access = private)
 		testObj_
 	end

	methods (Access = private)
		function cleanFiles(this)
 		end
	end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

