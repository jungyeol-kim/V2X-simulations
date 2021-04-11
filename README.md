Vehicular Messaging Simulation in Transportation Networks
===================================

V2V technologies bridge two infrastructures: communication and transportation. These infrastructures are interconnected and interdependent. To capture this inter-dependence, which may vary in time and space, we propose a new methodology for modeling information propagation between V2V-enabled vehicles. The model is based on a continuous-time Markov chain which is shown to converge, under appropriate conditions, to a set of clustered epidemiological differential equations. The fraction of vehicles which have received a message, as a function of space and time may be obtained as a solution of these differential equations, which can be solved efficiently, independently of the number of vehicles.

**Objective**: Our goal is to model the spread of V2V messages and obatin the fraction of vehicles which have received a message in arbitrary transportation networks, as a function of space and time, using a set of clustered epidemiological differential equations.

### Recommended citation
Please cite our paper below if you use this code.

*J. Kim, S. Sarkar, S. S. Venkatesh, M. S. Ryerson, and D. Starobinski, 2019. An Epidemiological Diffusion Framework for Vehicular Messaging in General Transportation Networks. Transportation Research Part B: Methodological, 131, pp.160-190.*



### Step-by-step instructions
#### Please find a "readme.pdf" file illustrating the usage of the code
For first time users, here is the quick version of the instructions. Follow these steps:

1. Open RStudio.
2. Set working directory.
    - Go to Tools > Change Working Dir... menu (Session > Set Working Directory on a mac), or use setwd() function. eg., setwd("/Users/jy/V2Vproject").
    - Check that the “home” directory for your project is the working directory of our current R process by typing getwd() in the Console.
3. Keep all the files (required input csv files, R scripts) associated with a project located together in the working directory.
4. Open the script file V2Vsimulation.R in the editor.
    - Go to File > Open > V2Vsimulation.R or type file.edit('V2Vsimulation.R') in the Console.
5. Setting parameters values. (set values in lines 10,11,12 of the script file)
    - total.clusters: the total number of clusters
    - sim.time: end time of the time sequence
    - step.size: step, increment of time
    - e.g., if step.size=1 and sim.time=100, the result will be generated every 1 time unit, from 0 to 100 units.
6. Check output
    - Type source("V2Vsimulation.R", echo = TRUE) in the Console.
7. Check output of this software
    - The output of this software consist of csv file representing the fraction of vehicles which have received a message, as a function of space and time.
    - The output files are saved in the the current working directory.
    - fraction_of_informed_vehicles_per_cluster.csv: the fraction of informed vehicles per clusters over time.
    - fraction_of_non_formed_vehicles_per_cluster.csv: the fraction of non-informed vehicles per clusters over time.
    - Most users will want to create graphs using the output files. Refer to the last section of this manual.


### Table of Contents 
The software manual is written in the following structure:
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
