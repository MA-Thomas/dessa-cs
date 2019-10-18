The simulator is run from the file "run_spatial_simulation.m".

0. Initialize data structures
1. Call the function/file: "preCompute_propensityCDF_curves.m" before the simulation begins to evaluate the integrated propensity curves
2. Initialize the queue with events for all particles. Function/file: "initializeRxQueue______v2_alt.m"
3. Main event loop

In each iteration of the main event loop, the min wait time event is selected from the queue.  
The script "script__unimolecular_update______v2_alt.m" is called in the case of a unimolecular event.
The script "script__bimolecular_update______v2_alt.m" is called in the case of a bimolecular event.
The scripts "Part1_script__positionOnly_update______v2_alt.m" and "Part2_script__positionOnly_update______v2_alt.m" are called in the case of a single (or multiple successive) position only update 
event(s). 
	If multiple events of this type, all are run together.

After running "run_spatial_simulation", run "plot_MichaelisMenten.m" to visualize the resulting populations as a function of time.