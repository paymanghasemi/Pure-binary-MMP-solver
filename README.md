# Pure-binary-MMP-algorithm
An Algorithm for solving Pure Binary Linear Maximum Multiplicative Programs 

This is an algorithm to solve a class of Pure Binary Linear Maximum Multiplicative Programs of format:

![Images](Images/problem.jpg)

where n represents the number of  binary decision variables. Also, D is a p * n matrix, d is a p-vector, A is an m * n matrix, and b is an m-vector.

This project is a Netbeans IDE 8.2 C++ project which was written in Linux (Ubuntu).

To compile the code, CPLEX (default version 12.7) must be installed on your computer. The default directory for CPLEX used is /opt/CPLEX/12.7/. Changing the directory of CPLEX to your preferred directory can be done either in the Makefile or through Netbeans. If you would like to do it in the Makefile you should go to nbproject/Makefile-Debug.mk and nbproject/Makefile-Release.mk and change all instances of /opt/CPLEX/12.7/ to your preferred directory. If you would like to do it through Netbeans, you can open the project in Netbeans and right click on the name of the project and choose Properties. You can then change the directory in the Include Directories box which is located in the C++ Compiler sub-menu. Moreover, you should also change the directory in the Libraries box which is located in the Linker sub-menu.

# Data Files
The instances used in the computational experiments of this algorithm are available in https://goo.gl/otpJwX.

For further instances, each data file should be written as a CPLEX LP file as follows,

![Images](Images/instance.jpg)


where the first p x's are the objective functions respresenting the y1, y2,..., and yp of the problem, and the rest of x are the decision variables.

Our algorithm is only capable of solving instances with only binary variables. So, in order to implement the algorithm on pure integer instances, convert all the integer variables to binaries.


# Compiling and Running

Compiling the project can be done using the Terminal by going to the directory in which the project is saved and typing ”make clean” followed by ”make” (you can also compile through Netbeans).

An instance can be solved by typing 

./Algorithm <*address*>/<instance*> <*number of objective*>

where instance* is the original .lp file of the instance.

For better understanging, we have provided a folder named Instances in the algorithm folder. In instance folder, 1.lp is the problems original file, which has 4 objective functions. In order to solve this instance, one can use the code

make clean

make

./Algorithm instance/1.lp 4

# Supporting and Citing

The software in this ecosystem was developed as part of academic research. If you would like to help support it, please star the repository as such metrics may help us secure funding in the future. We would be grateful if you could cite:

[Ghasemi Saghand, P., Charkhgard, H., A Criterion Space Search Algorithm for Binary Linear Maximum Multiplicative Programs]
