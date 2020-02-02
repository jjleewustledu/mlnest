classdef Test_Linear < matlab.unittest.TestCase 
	%% TEST_LINEAR  

	%  Usage:  >> results = run(mlnest_unittest.Test_Linear)
 	%          >> result  = run(mlnest_unittest.Test_Linear, 'test_dt')
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
        iterations = 3
 	end 

	methods (Test) 
 		function test_run(this) 
            for iter = 1:this.iterations
                this.testObj.run(this.testObj)
            end
            
            % # iterates ~ MAX = 500
            % # sampling particles = 25.000000
            % MCMC_Counter = 50.000000
            % STEP_Initial = 0.200000
            % Stopping criteria = 1.269903
            % Evidence:  ln(Z) = -25.306110 +/- 0.793706
            % Information:  H = 15.749232 nats = 22.721339 bits
            % Model:
            % 	intercept = -20.012062 +/- 0.000000
            % 	slope = 2.000717 +/- 0.000000
            % 	sigma0 = 0.001000
            % 	sampled logL(k = 500) = -8.998522
            % 	sampled logWt(k = 500) = -32.197331
            % 11.1031 seconds testing time.
    
 		end 
 		function test_Linear2_run(this) 
 			obj = mlnest.Linear2;
            obj.MAX = 500;
            obj.n = 25;
            obj.MCMC_Counter = 50;
            obj.STEP_Initial = 0.2;            
            obj.sigma0 = 1e-3;
            for iter = 1:this.iterations
                obj.run(obj)
            end
 		end 
 	end 

 	methods (TestClassSetup) 
 		function setupLinear(this)   
            rng('default')         
 			this.testObj = mlnest.Linear;
            this.testObj.MAX = 500;
            this.testObj.n = 25;
            this.testObj.MCMC_Counter = 50;
            this.testObj.STEP_Initial = 0.2;            
            this.testObj.sigma0 = 1e-3;
 		end 
 	end 

 	methods (TestClassTeardown) 
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
 end 

