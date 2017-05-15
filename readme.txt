This folder contains the main programs and routines of the project
valid on 15/12/2016
updated until 04/2017

Updated nuclear reaction rates are available in the folder 
[Accretion_disk/nuclear_rates_updated]
which contains a program to plot the burning rates 
in function of temperature

-----------------------------------------------------------------------------------
The user only needs to change parameters (initial composition of the disk,
 	accretion rate, black hole mass...) in the input file 
------------------------------------------------------------------------------------
List of "include" files

- includeall.h is the only one the user needs to care about. 
	It declares the functions of "functions.cpp" and includes "include_rates.h" 
	as well as mathematical functions and C++ print functions.

- include_rates.h contains the declaration of the nuclear burning rates.

------------------------------------------------------------------------------------
List of subroutines:

- functions_v3.cpp contains all the functions used in the main programs, 
	mainly regarding the accretion disk structure under Shakura-Sunavey regime
	and a model of ADAF to compare to Shakura-Sunyaev
	-> contains matrix A to be inverted for complete alpha network
	works on 14/11/2016 - updated until 04/2017

- functions_v4.cpp update of functions_v3.cpp for adding pp chain network.
	did not finish.
	do not know if it works.

- nuclear_rates.cpp contains the nuclear burning rates to use. 
	They are numbered "rateX()" X = 0, 1, ... with X increasing with the nuclei
	sizes or X = C12C12, O16O15, C12O16 etc. contains the two nuclei interacting

------------------------------------------------------------------------------------
List of the programs:

- abundances_v3.cpp latest version, works on 15/12/2016, updated 04/2017
	integrates the radial composition of the disk and the total yields 
	of NS during accretion onto the black hole.
	to be fed with input file (e.g. input.txt)

- rhoTplane.cpp, works on 15/12/2016
	computes:
	path of a gas blob through the disk in the rho-T plane
	path of the center of a star in the rho-T plane through its evolution
	path of a blob of gas in the rho-T plane during the Big Bang
	conditions necessary for nuclear burning in the rho-T plane
	data exported in 4 separate output files.

- rhoTplane_v3.cpp, works on 15/12/2016
	same as rhoTplane.cpp, updated with possiblity of BH spin
	with parameter 0 < a < 1.

- mmdotplane.cpp, works on 15/12/2016
	computes for various ranges of BH mass and accretion rates
	the greatest disk temperatures and density (inner most)
	the accretion time scale (viscous)
	the burning time scale for triple alpha and C12+alpha

- mmdotplane_v3.cpp, works on 02/2017
	same as mmdotplane.cpp, updated with possiblity of BH spin
	with parameter 0 < a < 1.


------------------------------------------------------------------------------------
Download LAPACK++ PACKAGE
sudo apt-get install build-essential liblapack-dev libblas-dev checkinstall
Download the version lapackpp-2.5.4.tar.gz from 
http://sourceforge.net/projects/lapackpp/files/ 
and unpack it
In the directory of the lapackpp-2.5.4, type
./configure
make
make install
sudo checkinstall

------------------------------------------------------------------------------------
Compilation in the directory of the code:
1. compile the files and link to objects:
g++ -c abundances_v3.cpp functions_v3.cpp nuclear_rate_v3.cpp -I /usr/local/include/lapackpp
2. compile objects to make program
g++ abundance_v3.o  functions_v3.o nuclear_rates.o -o program_name -llapackpp
3. give path to shared library
export LD_LIBRARY_PATH=/usr/local/lib

------------------------------------------------------------------------------------
Run:
./program_name input.txt

------------------------------------------------------------------------------------
description of input file: input.txt
--------------------------
M_BH	M_WD
Mdot	alpha	spin
n_iso	h

He
C
O
Ne
(... other alpha elements...)
Ni

Output_name
-------------------------
where: 
M_BH = black hole mass in solar mass
M_WD = white dwarf mass in solar mass
Mdot = accretion rate in Eddington units
alpha = viscosity parameter
spin = normalised BH spin (between 0 and 1)
n_iso = number of isotopes in the network (can only be 13 in current version!!!)
h = max time-step, about 0.01 s (visc timescale inner disk, increase with M_BH)
He = mass fraction of helium
Output_name = head name of output file (don't write ".txt", code does it)

===================================================================================
LIMITATIONS (explained in thesis text)
===================================================================================
- WD mass must be < 1.4 Msun
- BH mass must be < 300 Msun
- Accretion rate limited (Fig 5.6 thesis)
- Sum of mass fractions for initial disk composition must be 1


