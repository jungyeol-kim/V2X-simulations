Vehicular Messaging Simulation in Transportation Networks
===================================

V2V technologies bridge two infrastructures: communication and transportation. These infrastructures are interconnected and interdependent. To capture this inter-dependence, which may vary in time and space, we propose a new methodology for modeling information propagation between V2V-enabled vehicles. The model is based on a continuous-time Markov chain which is shown to converge, under appropriate conditions, to a set of clustered epidemiological differential equations. The fraction of vehicles which have received a message, as a function of space and time may be obtained as a solution of these differential equations, which can be solved efficiently, independently of the number of vehicles.

**Objective**: Our goal is to model the spread of V2V messages and obatin the fraction of vehicles which have received a message in arbitrary transportation networks, as a function of space and time, using a set of clustered epidemiological differential equations.

### Recommended citation
Please cite our paper below if you use this code.

J. Kim, S. Sarkar, S. S. Venkatesh, M. S. Ryerson, and D. Starobinski, 2019. An Epidemiological Diffusion Framework for Vehicular Messaging in General Transportation Networks. Transportation Research Part B: Methodological, 131, pp.160-190.


### Table of Contents 
The manual is written in the following structure:
- R Packages we use
- Introduction
  - Recommended citation
  - Requirements and installation
    - Input files (required)
      - File 1: mobility-network.csv
      - File 2: communication-network.csv
      - File 3: initial-condition.csv
  - Step-by-step instructions
  - Case study: Grid road topology
- Mobility and Communicatoin Network
  - Importing edge list of mobility network
  - Importing edge list of communication network
  - Defining a neighborhood of a cluster
- Generation of Differential Equations
- Solving Differential Equations
  - Initial condition
  - Solution of differential equations
- Generating Figures
  - Fraction of overall informed vehicles over time
  - Fraction of informed and non-informed vehicles per cluster over time
