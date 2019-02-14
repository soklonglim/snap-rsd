========================================================================
    Rumor Source Detection RSD
========================================================================
The application implements:
   -an infected directed graph using Independent Cascade and record infected 
   node at each time step from input graph
   -reverse dissemination method using Wavefront model to calculate suspected
   node (source of infection) 

For more details and motivation what this code is trying to achive see
“Maximizing the Spread of Influence through a Social Network” by D. Kempe, 
J. Kleinberg, E. Tardos, 2003.
“Rumor Source Identification in Social Networks with Time-varying Topology” by
J. Jiang, S. Wen, S. Yu, Y. Xiang, W. Zhou, 2015.
“Approximating the k-Minimum Distance Rumor Source Detection in Online Social Networks” by
S. Lim, J. Hao, Z. Lu, X. Zhang, Z. Zhang, 2018.


The code works under Windows with Visual Studio or Cygwin with GCC,
Mac OS X, Linux and other Unix variants with GCC. Make sure that a
C++ compiler is installed on the system. Visual Studio project files
and makefiles are provided. For makefiles, compile the code with
"make all".

/////////////////////////////////////////////////////////////////////////////
Input files:
   -CollegeMsg.txt
   -email-Eu-core.txt

/////////////////////////////////////////////////////////////////////////////
Parameters:
   -i:Input graph (tab separated list of source-node-id destination-node-id weight-value) 
   -o:Output file name (default:'demo')
   -t:Infection threshold (infected probability), can be fixed or random 

/////////////////////////////////////////////////////////////////////////////
Usage:

./rsd -i:../email-Eu-core.txt -o:output.txt -t:0.1
