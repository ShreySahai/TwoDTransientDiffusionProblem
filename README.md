# 2DHeatTransfer

A python library that simulates two-dimensional transient heat conduction in a square plate using ADI finite difference method


## Getting started 

The repository can be cloned and the library `TwoDHeat.py` can be simply imported for use

### Prerequisites 

The library uses the following two python libraries

* numpy
* matplotlib

### Problem statement

A square plate of length 2L and thermal diffusivity α with initial temperature Tin is suddenly subjected to temperature T0

The problem is a diffusion type problem and the governing equation is given by 

![Equation 1](/Images/1.png)

with IC, T = Tin at t = 0, ∀ x,y
and BCs, T = T0 at x = L and y = L for ∀ t > 0

After appropriate non-dimensionalisation the equation can be given as

![Equation 2](/Images/2.png)

Alternating Directional Implicit finite difference method is used to calculate the temperature 

Implicit FDM formulation in x direction for the first half time step results in

![Equation 3&4](/Images/34.png)

where θ*  denotes θ at time step t= p+½ 

These equations don’t incorporates boundary cases and image point technique is used for such cases

Further FDM over y direction results in

![Equation 5&6](/Images/56.png)

## Usage

The object to the class has to be passed with a dictionary of parameters for the simulation
Objects's `solve()` function is to be called for the animation to start 

```python
from 	TwoDHeatTransfer	import	TwoDHeatTransfer as problem
parameters ={
			"Length"		: 1.,			#Length of the plate
			"T0"			: 300.,			#Surrounding Temperature
			"Tin"			: 100.,			#Initial Temperature
			"alpha"			: 0.5e-6,		#Thermal diffusivity
			"Nelements"		: 20,			#Number of elements along the half plate
			"timestep"		: 0.01,			#Value of dt(dimensionless time step) 
			"Niteration"	: 150,			#Number of iterations over 
} 
p1 = problem(parameters)					#initialise object with these parameters
p1.solve()				
```
This is one snapshot of the animation 

![Snapshot taken during Animation](/Images/Example.png)






