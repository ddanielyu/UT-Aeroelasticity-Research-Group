Run main_monte_carlo.m

This code performs Markov Chain Monte Carlo simulations of a ball falling from 
100 m in the presence of wind and gravity. 

The drag area is 0.14 m2 and initial velocity at 10 m height is 5 mph (2.2m/s).

The hourly wind velocity data is obtained from Austin airport data at 10 m.

A statistical profile at 10 Hz is obtained using Markov Chain method and
extrapolated to other altitudes using log-law atmospheric wind profile.

The equations of drag & gravity are simultaneously solved to obtain trajectory.

Ref [1]: Papaefthymiou, G. and Klockl, B., 2008. MCMC for wind power simulation
Ref [2]: weather data from https://www.wunderground.com/history/daily/us/tx/austin/KAUS/date/2021-9-1 through 2021-9-30

Images of sample results is attached.