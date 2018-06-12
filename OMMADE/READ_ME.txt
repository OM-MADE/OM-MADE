# ==============================================
#
#      GENERAL INFORMATION FOR USING OM-MADE
#
# ==============================================

This program has been written in Python 3.6. It requires the python libraries: numpy and matplotlib.
To easily install all these features for non-used user, we advice to directly download a free integrated
Python development environment like Spyder, created by Pierre Raybaut, maintained by the 
Spyder Project Contributors (https://github.com/spyder-ide/spyder/blob/master/AUTHORS).


In the folder "Codes_OMMADE", the files classDataPoint, classParameters, readData and timeLoops
should normally NOT be MODIFIED by a user. They consitutes the core of the code. All the main programs 
that have been written, call these files by relative path. Thus, the relative tree view SHOULD NOT be modified.

The user can use the "OMMADE_generic_Main" file to launch the program: Open it in Spyder, and click on the 
green arrow to execute the file. In this case, input parameter files should have been created 
in a specific folder and the corresponding path indicated in the _INPUTFILE.txt file. 
It would generate output files in the same folder => look at the folder in an Explorer to see the output files.

The user can also create his own main program, to directly adapt it to his problem
and generate the output files he desires. matplotlib could be used in this case
to directly vizualise the results in the Python development environment. Several examples 
of adapted main programs are thus provided in the folders "validations" and "Furfooz_TracerTest3".

For further information, one can refer to the associated paper : 
Tinet, Collon, Philippe, Dewaide, Hallet, (XXXXX) OM-MADE: an open-source program to simulate one-dimensional
 solute transport in multiple exchanging conduits and storage zones. Computers & Geosciences, Vol XX, pp. XXX-XXX.