# Welcome to OM-MADE

OM-MADE is an open-source progam written in Python v3. It is designed to simulate one-dimensional solute transport in multiple exchanging conduits and storage zones.

If you use OM-MADE, please cite the following reference : 
"Tinet A.J., Collon P., Philippe C., Dewaide L., Hallet V. (Subm.). OM-MADE: an open-source program to simulate one-dimensional solute transport in multiple exchanging conduits and storage zones. "

This paper details all technical aspects of OM-MADE, the hypothesis, the equations, the numerical schemes and the validations made against analytical and numerical solutions.

#Installation

This package needs a Python version >=3.3 . Furthermore it uses numpy and matplotlib.

We recommand to use an integrated environment of development like Spyder that you can easily download here : https://pythonhosted.org/spyder/installation.html or with the Python Scientific Distribution Anaconda : https://www.anaconda.com/download/
This solution will automatically install all packages (numpy and matplotlib) that are required. 

#Running a simulation

The user can use the "OMMADE_generic_Main" file to launch the program: Open it in Spyder, and click on the 
green arrow to execute the file. In this case, input parameter files should have been created 
in a specific folder and the corresponding path indicated in the _INPUTFILE.txt file. 
It would generate output files in the same folder => look at the folder in an Explorer to see the output files.

The user can also create his own main program, to directly adapt it to his problem
and generate the output files he desires. matplotlib could be used in this case
to directly vizualise the results in the Python development environment. Several examples 
of adapted main programs are thus provided in the folders "validations" and "Furfooz_TracerTest3".

#Authors

The OM-MADE main developpers are: 
  - Dr. Anne-Julie Tinet (anne-julie.tinet@univ-lorraine.fr) - Associate Professor - Université de Lorraine, CNRS, GeoRessouces laboratory
  - Dr. Pauline Collon ( pauline.collon@univ-lorraine.fr) - Associate Professor - Université de Lorraine, CNRS, GeoRessouces laboratory
  
  #License
  
  The current license of the software is BSD 3-Clause "New" or "Revised" License
