classdef Test_Lighthouse < matlab.unittest.TestCase 
	%% TEST_LIGHTHOUSE  

	%  Usage:  >> results = run(mlnest_unittest.Test_Lighthouse)
 	%          >> result  = run(mlnest_unittest.Test_Lighthouse, 'test_dt')
 	%  See also:  file:///Applications/Developer/MATLAB_R2014b.app/help/matlab/matlab-unit-test-framework.html

	%  $Revision$ 
 	%  was created $Date$ 
 	%  by $Author$,  
 	%  last modified $LastChangedDate$ 
 	%  and checked into repository $URL$,  
 	%  developed on Matlab 8.5.0.197613 (R2015a) 
 	%  $Id$ 

	properties 
 		testObj 
 	end 

	methods (Test) 
 		function test_run(this) 			
            import mlnest.*;
 			this.testObj = Lighthouse.run;
 		end 
 	end 

 	methods (TestClassSetup) 
 		function setupLighthouse(this) 
 		end 
 	end 

 	methods (TestClassTeardown) 
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
 end 

