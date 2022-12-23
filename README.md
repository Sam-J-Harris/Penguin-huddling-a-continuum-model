# Penguin huddling: a continuum model
Supplementary material accompanying the paper "Penguin huddling: a continuum model". The attached MATLAB code outputs the propagation of a penguin huddle boundary and other related free boundary problems.


## Map:
It's often hard to follow which code is which, so I am providing a "map" of where to look.

	Start at penguin_RUN
	(This is the only code you need to run.)

	penguin_RUN -> penguin_ode_solve ; penguin_error; penguin_plots ; penguin_plots_steady_shape ; penguin_plots_heat_flux
	
		penguin_ode_solve -> penguin_initial_shape ; odeqns ; sigsolve ; AAA_LS_solve
		
			odeqns -> SOCcalc
			
				SOCcalc -> sigsolve ; AAA_LS_solve ; adzeta
					
					AAA_LS_solve -> VAevald
	
		penguin_plots_steady_shape ; penguin_plots_heat_flux -> centrepoly


	To recreate the figures from the paper, run the code: penguin_huddling_a_continuum_model

All other functions used are either from MATLAB or the chebfun package.

## Function Glossary:
A brief description of what each function does (in some kind of "running" order).

## References:
For full reference list, see the paper Penguin huddling: a continuum model. References used in the attached code are listed below.

Costa & Trefethen 2021

Ladd et al. 2020

Rycroft and Bazant 2016

T. A. Driscoll, N. Hale, L. N. Trefethen, Chebfun Guide. Pafnuty Press, Oxford, 2014;
see also www.chebfun.org

Dallaston and McCue 2016

