# Induced_Seismicity

The structure and function of each file is as follows:

•	Pfs.py – (Probability of Fault Slip) Package written by David Healy, contains the python code required to perform many of the calculations.

•	Sga.py – (Structured geology Algorithms) Package written by David Healy, contains the python code for geological functions (calculate shear stress, convert cartesian to spherical coordinates etc.)

•	K-Test.ipynb – Python notebook to show k-test validity of using a Monte-Carlo Von-Mises sampler

•	WorkedExample1/2 – Python notebook showing Healy’s results using the original implementation

•	WorkedExample3 – Python notebook containing the updated python implementation with the changes discussed in this paper. 

•	Points.txt – txt file outputted by C++ code of all 5000 angular datapoints using Von-Mises Sampler

•	Pfs.cpp & sga.cpp – Healy’s Python code refactored to C++, contains Von-Mises Monte-Carlo Sampler. Pfs.cpp contains the main function to output a points.txt file.
