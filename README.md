## Penguin-huddling-a-continuum-model
Supplementary material accompanying the paper "Penguin huddling: a continuum model". The attached MATLAB code outputs the propagation of a penguin huddle boundary and other related free boundary problems.


# Map:
It's often hard to follow which code is which, so I am providing a "map" of where to look.

	Start at penguin_RUN. 
	(This is the only code you need to run.)

	penguin_RUN -> penguin_ode_solve ; penguin_plots ; penguin_steady_shape_plot.
	
		penguin_ode_solve -> penguin_initial_shape ; odeqns.
		
			odeqns -> SOCcalc
			
				SOCcalc -> sigsolve ; AAA_LS_solve ; adzeta
	
		penguin_steady_shape_plot -> centrepoly
