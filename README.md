Vehicular Messaging Simulation in Transportation Networks
===================================

V2V technologies bridge two infrastructures: communication and transportation. These infrastructures are interconnected and interdependent. To capture this inter-dependence, which may vary in time and space, we propose a new methodology for modeling information propagation between V2V-enabled vehicles. The model is based on a continuous-time Markov chain which is shown to converge, under appropriate conditions, to a set of clustered epidemiological differential equations. The fraction of vehicles which have received a message, as a function of space and time may be obtained as a solution of these differential equations, which can be solved efficiently, independently of the number of vehicles. 

The V2X_simulations.R file contains both 1) the code for information propagation simulation based on a continuous-time Markov chain and 2) the code for numerically solving the corresponding set of clustered epidemiological differential equations.

Once you open the file in RStudio, please fold the code (go to the menu bar and select: Edit > Folding > Collapse All) for easy navigation and then unfold the section of interest by clicking the caret.

### Recommended citation
Please cite our paper below if you use this code.

J. Kim, S. Sarkar, S. S. Venkatesh, M. S. Ryerson, and D. Starobinski, 2019. An Epidemiological Diffusion Framework for Vehicular Messaging in General Transportation Networks. (In press) Transportation Research Part B: Methodological

And respective BibTeX entry:

@article{KSVRS18,  
  title={An Epidemiological Diffusion Framework for Vehicular Messaging in General Transportation Networks},  
  author={Kim, Jungyeol and Sarkar, Saswati and Venkatesh, Santosh S and Ryerson, Megan Smirti and Starobinski, David},  
  journal={(In press) Transportation Research Part B: Methodological},  
  year={2019},  
  publisher={Elsevier}  
}

### Table of Contents 
The code is written in the following structure:
- R Packages we use
- Transportation network
  - Setting a size of grid topology
  - Creating directed mobility network and corresponding adjacency metrix
  - Temporal variation of traffic density and routing: setting routing probability pertaining to vehicle movement
- Information propagation simulation based on a continuous-time Markov chain
  - Creating function: simulations for information propagation in V2V-enabled transportation network
  - Running multiple simulations to calculate ensemble average
- Solution of the set of clustered epidemiological differential equations
  - Automatic generation of the set of differential equations
  - Creating function: set of ordinary differential equations
  - Getting numerical solution of the set of the differential equations
- Geographical representation of traffic density and information propagation for both simulation results and model solutions
