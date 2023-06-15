# Penguin huddling: a continuum model
Supplementary material accompanying [1]. The attached MATLAB code outputs the propagation of a penguin huddle boundary and other related free boundary problems. Requires MATLAB and the chebfun package -- see [6] below.

## Map:
It's often hard to follow which code is which, so I am providing a "map" of where to look.

	Start at A1_penguin_RUN
	(This is the only code you need to run.)

	A1_penguin_RUN -> penguin_ode_solve ; penguin_error; penguin_plots ; penguin_plots_steady_shape ; penguin_plots_heat_flux
	
		penguin_ode_solve -> penguin_initial_shape ; odeqns ; sigsolve ; AAA_LS_solve
		
			odeqns -> SOCcalc
			
				SOCcalc -> sigsolve ; AAA_LS_solve ; adzeta
					
					AAA_LS_solve -> VAevald
	
		penguin_plots_steady_shape ; penguin_plots_heat_flux; penguin_plots_aspect_ratio -> centrepoly

All other functions used are either from MATLAB or the chebfun package.

## Function Glossary:
A brief description of what each function does (in some kind of "running" order).

A1_penguin_RUN: run this code to obtain the evolution of the free boundary.

penguin_huddling_a_continuum_model: run this code to recreate the figures from [1] (excluding figure 1). 

penguin_ode_solve: produces the free boundary evolution and the canonical and physical heat fluxes of the final free boundary shape.

penguin_error: returns the root mean-squared error and area error for the free boundary shapes. 

penguin_plots: plots the free boundary evolution.

penguin_plots_steady_shape: plots the steady shape (and compares this steady shapes from other experiments).

penguin_plots_heat_flux: three plots of heat flux data across the steady shape.

penguin_initial_shape: returns Laurent coefficients and their derivatives of desired initial boundary shape -- see [2]. 

odeqns: the system of ODEs to be solved by decic/ode15i.

SOCcalc: calculates the exterior, interior and area conserving effects of the Polubarinova-Galin type equation -- see [1].

sigsolve: solves the integral equation governing the exterior effect sigma(theta) -- see [3].

AAA_LS_solve: uses AAA-least squares algorithm to return F'(z) where u=Re(F) is harmonic in the interior -- see [1], [4], [5].

adzeta: find |dz/dzeta|, for use with the "integrate" function.

VAevald: the function VAeval from chebfun returning only the basis vector for derivatives -- see [6].

penguin_plots_aspect_ratio: plots the aspect ratio (AR) of the steady shape vs a parameter, either Pe or beta.

centrepoly: centres a polygon around the origin.

## References:
For full reference list, see [1]. References used in the attached code are listed below.

[1]	Harris, S.J., McDonald, N.R. (2023) "Penguin Huddling: A Continuum Model". Acta Appl. Math. 185, 7. https://doi.org/10.1007/s10440-023-00578-2.

[2]     Rycroft, C. H., & Bazant, M. Z. (2016). "Asymmetric collapse by dissolution or melting in a uniform flow". Proc. R. Soc. A, 472(2185), 20150531.

[3]     Ladd, A. J., Yu, L., & Szymczak, P. (2020). "Dissolution of a cylindrical disk in Hele-Shaw flow: a conformal-mapping approach". J. Fluid Mech., 903.

[4]     Costa, S., & Trefethen, L. N. (2021). "AAA-least squares rational approximation and solution of Laplace problems". arXiv preprint arXiv:2107.01574.

[5]     Nakatsukasa Y., Sete 0., Trefethen L.N. (2018) "The AAA algorithm for rational approximation". SIAM J. Sci. Comp. 40, A1494-A1522.

[6]	T. A. Driscoll, N. Hale, L. N. Trefethen, Chebfun Guide. Pafnuty Press, Oxford, 2014; see also www.chebfun.org.


