# pde-solver

Solve the elliptic Poisson partial differential equation using finite difference methods.

<pre>
.
└── main.py                   				# Main program which contains all the code run for the project. These might take a while
└── main_plots.py              				# Generates all the plots used in the report from stored data
└── modules                   				# Contains main classes
│   ├── ceramic_cooling.py                 	# Class which sets up the ceramic-microprocessor structure
│   ├── heat_structure.py           		# Class which sets up the microprocessor with heat sink structure
│   ├── simulated_annealing.py         		# Simulated annealing method to perform discrete optimisation
│	├── automated_runs						# Folder containing all the automated runs such as sweeps
│	│	├── forced_convection_runs.py 		# Sweep performed for forced convection to try to find parameters to cool below 80 degrees
│	│	├── natural_convection_runs.py 		# Sweep performed for natural convection to find best params to cool microprocessor
│	│	└── scale_runs.py 					# Sweep performed to see how microprocessor convergence temperature varied with step size and scale
│	├── plots 								# Folder containing the modules which generate all the plots for the report
│	│	├── forced_convection_sweep.py 		# Generates plots from the forced convection sweeps
│	│	├── natural_convection_sweep.py 	# Generates plots from the natural convection sweeps
│	│	├── microprocessor_ceramic_plots.py # Generates plots for microprocessor-ceramic natural convection cooling
│	│	└── simulated_annealing_analysis.py # Generates plots from the natural convection simulated annealing run
│  	└── utils
│		└── solvers.py 						# Contains the functions for different iterative solvers - Jacobi, Gauss-Seidel, SOR, Red-Black SOR
├── tests
│	└── test_solvers.py 					# Unit tests which validates the solver works as expected
└── data 									# Contains data files used to generate plots to cut down wait time
</pre>

### Set up
1. Ensure that Python version: 3.7.0 is installed
2. Install project requirements:
<pre>
$ pip install -r requirements.txt
</pre>
3. Run unit tests
<pre>
$ cd into project folder
$ pytest
</pre>
4. Generate plots
<pre>
$ python main_plots.py
</pre>
5. Run code (will take about 10 hours to run everything, including sweeps and simulated annealing optimiser)
<pre>
$ python main.py
</pre>

### Usage of classes

1. HeatStructure - class which sets up a heat sink structure atop the microprocessor
	Accepts 9 parameters: 
	-	scale (scale of 1 means each mesh point is 1mm)
	-	b (fin separation)
	-	c (fin thickness)
	-	f_h (fin height)
	-	n_fins (number of fins)
	-	conv_ratio = 1E-6
	-	convection_type = "forced"
	-	solver = jacobi_solver
	-	wind_speed = 20

Example usage:
<pre>
hs = HeatStructure(2, 5, 1, 30, 5, conv_ratio=1E-6, convection_type="natural", solver=solv.red_black_SOR)
T, n = hs.solve_mesh()
</pre>

2. Microprocessor - class which sets up a lone microprocessor with a ceramic case
	Accepts 3 parameters:
	- 	scale
	- 	conv_ratio
	-	solver

Example usage:
<pre>
mp = Microprocessor(5, conv_ratio=1E-8, solver=solv.jacobi_solver)
T, n = mp.solve_mesh()
</pre>

### Simulated Annealing

1. Specify start parameters, e.g. x0 = [2, 2, 30, 30]
2. Specify start annealing temperature, e.g. T = 1.0
3. Specify minimum annealing temperature, e.g. T_min = 0.00001
4. Specify the factor alpha which determines how fast T falls

<pre>
simulate_annealing(x0, T, T_min, alpha)	
</pre>
