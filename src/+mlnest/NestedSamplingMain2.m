classdef NestedSamplingMain2 < handle & mlnest.NestedSamplingMain
	%% NESTEDSAMPLINGMAIN2 implements main.c from "Data Analysis:  A Bayesian Tutorial, Second Edition"
    %  by D.S. Sivia and J. Skilling, section 9.2.4.  
    %  (GNU General Public License software, (C) Sivia and Skilling 2006)

	%  $Revision$ 
 	%  was created $Date$ 
 	%  by $Author$,  
 	%  last modified $LastChangedDate$ 
 	%  and checked into repository $URL$,  
 	%  developed on Matlab 8.4.0.150421 (R2014b) 
 	%  $Id$ 
    
    methods        
  		function this = NestedSamplingMain2(app) 
 			%% NESTEDSAMPLINGMAIN2
            %  @param app is an application object implementing mlnest.IApply
            %  @return this
 			%  Usage:  this = NestedSamplingMain2(mlnest.IApply object)
            %  Internally:
            %     Object Obj[n];        // Collection of n objects
            %     Object Samples[MAX];  // Objects stored for posterior results
            %     double logwidth;      // ln(width in prior mass)
            %     double logLstar;      // ln(Likelihood constraint)
            %     double H    = 0.0;    // Information, initially 0
            %     double logZ =-DBL_MAX;// ln(Evidence Z, initially 0)
            %     double logZnew;       // Updated logZ
            %     int    i;             // Object counter
            %     int    copy;          // Duplicated object
            %     int    worst;         // Worst object
            %     int    nest;          // Nested sampling iteration count
            
            this = this@mlnest.NestedSamplingMain(app);
        end 
    end 
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

