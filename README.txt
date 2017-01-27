  __  __  __  __  __  __  __  __  __  __  __  __  __  __  __  __  __  __   
 /_/ /_/ /_/ /_/ /_/ /_/ /_/ /_/ /_/ /_/ /_/ /_/ /_/ /_/ /_/ /_/ /_/ /_/

            RNANR - RNA Non-Redundant - program to compute 
          and sample non-redundantly RNA socondary structures
        

           __  __  __  __  __  __  __  __  __  __  __  __  __  
          /_/ /_/ /_/ /_/ /_/ /_/ /_/ /_/ /_/ /_/ /_/ /_/ /_/
             ______    __    _    ____    __    _   ______
            |  __  \  |  \  | |  / __ \  |  \  | | |  __  \
            | |__| |  |   \ | | / |__| \ |   \ | | | |__| |
            |  _   /  | |\ \| | |  __  | | |\ \| | |  _   /
            | | \  \  | | \   | | |  | | | | \   | | | \  \
            |_|  \__\ |_|  \__| |_|  |_| |_|  \__| |_|  \__\
           __  __  __  __  __  __  __  __  __  __  __  __  __  
          /_/ /_/ /_/ /_/ /_/ /_/ /_/ /_/ /_/ /_/ /_/ /_/ /_/
	  
	  		

WHAT IS IT ?
------------

RNANR is a C program that reads an RNA sequence and computes all
possible locally optimal secondary structures for this sequence (or
optionnally just counting them). The folding model is based on the
Turner model. It can be refined through a series of parameters. See
the section "OPTIONAL PARAMETERS - TOPOLOGY OF THE SECONDARY STRUCTURE".

One of features of RNANR is also the ability to sample locally optimal
secondary structures without redundancy. The core of the functionality
is based on the computation of Boltzmann probability of each substructure
and generating the complete structure based on them. To know more about it,
please refer to the section "HOW DO I PERFORM SAMPLING?"


HOW TO INSTALL IT ?
-------------------

Go to the "src" directory, and type 'make'. This will create an executable
file "locopt".


HOW TO RUN IT ?
---------------

You can try the program with the file "test_file.fa" provided in the
samplel directory.

> RNANR -i ../sample/test_file.fa 

The input file is an RNA sequence in FASTA format. DNA alphabet
(A,C,G,T) is also accepted. 

By default, results are displayed on the standard output. With the
"test_file.fa" and default paramater values, you should get 1106
locally optimal secondary structures, that are represented in
bracket-dot format. The number at the end of each line is the total
number of base pairs in the structure.

HOW DO I PERFORM SAMPLING ?
---------------------------

There are three ways to perform sampling using optional parameters:

-s : Performs a non-redundant sampling with number of samples equal to
     the number entered following this parameteres. This value has to
     be given. Note that if higher number of samples than existing 
     structures is asked for, the program will return an error message 
     specifying so.

-z : Performs a non-redundant sampling, but unlike the previous option,
	 the number of samples corresponds to attained coverage (the sum of
	 Boltzmann probabilities of returned samples). Has to have value 
	 between 0 and 1.

-s <int> -r : Similar to option -s, but with redundancy, and consequently
			  without upper limit on number of structures returned.


HOW TO GET SUPPLEMENTARY INFORMATION ABOUT ALL PARAMETERS ?
-----------------------------------------------------------

The parameters are described more in detail below this paragraph. 
However, you can always print the manual included in the program itself
by simply typing:

> RNANR -h 

The help will be also printed if incorrect parameters or their values 
were entered.  


OPTIONAL PARAMETER - COUNTING
-----------------------------

-c : Count the number of locally optimal secondary structures, instead
     of building them.

OPTIONAL PARAMETERS - TOPOLOGY OF THE SECONDARY STRUCTURE
---------------------------------------------------------

The algorithm of locopt includes a series of combinatorial parameters,
that allow to drastically reduce search space and to generate more
realistic structures.

-a : Minimum helix length
     The minimum number of consecutive base pairs in a helix. 
     Default is 3. 

-k : Minimum hairpin loop length. 
     The minimum number of single stranded bases between any two bases
     that are paired together. Default is 3.

-l : Maximum hairpin loop length
     The maximum number of single stranded bases between any two bases 
     that are paired togethe. Default is 12.

-m : Maximum loop length
     The maximum number  of single stranded bases between any two 
     paired bases that do not form a terminal loop. Default is 7.

-n : Maximum span
     The maximum number of positions for a base pair to span. Default 
     is the length of the input sequence. 

-p : Minimum percentage of paired positions. Default is 0. Not
     compatible with option -c. 
     
-q : Maximum number of branches within multiloop. Default option is the
	 length of entered sequence, which basically means there is no limit.

For example, the secondary structure below fullfils the following parameters. 

       ((((.(((..((...))))).....(((....)))))))

Minimum helix length: 2
Maximum stemloop length: 4 
Minimum stemloop length: 3
Maximum loop length: 5 
Maximum span: 37
Minimum percentage of paired positions: 64
Maximum degree of multiloop: 2


OPTIONAL PARAMETERS - SET OF BASE PAIRS
---------------------------------------

By default, the set of all possible base pairs is composed of all A-U,
C-G and G-U base pairs present in the RNA sequence. It is also
possible to specify your own set of base pairs, and/or to combine both
solutions.

In this case, you should provide the program with a text file that
contains all the base pairs in dot-bracket format. This file can
contain as many lines as necessary. 

-d : This option allows the user to specify his/her own set of base
     pairs. The set of locally optimal secondary structures will be
     built solely on this set of base pairs.

-e : This option enables to specify additional base pairs, in
     complement to A-U, C-G and G-U base pairs present in the RNA
     sequence.

OPTIONAL PARAMETERS - TEMPERATURE SCALING
-----------------------------------------

The program has two ways of influencing the temperature at which energy 
computation and sampling is done:

- t : Changes the global temperature of the system to the value entered
      in Celsuis. Default value is 37°C. 

- b : Scaling of temperature that does not affect the free energy of
      computed structures and only intervenes during the computation
      of Boltzmann probabilities. Has positive values. Higher values can
      be used to randomize sampling. Default is 1.0

OPTIONAL PARAMETERS - OUTPUT FILE
---------------------------------

By default, the result is displayed on the standart output. 

-o : Name of the output file that will contain all locally optimal 
     secondary structures. Not compatible with option -c or sampling 
     options
     
OPTIONAL PARAMETERS - COMPUTING DIFFERENT VALUES
------------------------------------------------

Besides the count of all locally optimal secondary structures, there
are two options that can be called and return a global information about
the structures of entered sequence:

- f : Computes partition function of all locally optimal secondary 
      structures and exits.
      
- u : Computes the number of all flat structures and exits. 

KNOWN BUGS : The sampling starts to send invalid samples once their 
----------   Boltzmann probability starts getting low. This is caused
             by numerical imprecision. Fix is on their way.


AUTHOR :  Helene Touzet (helene.touzet@univ-lille1.fr) 
------    Juraj Michalik (juraj.michalik@inria.fr)
		  Yann Ponty	

Please send comments and bug reports to the address above.
		  
  __  __  __  __  __  __  __  __  __  __  __  __  __  __  __  __  __  __   
 /_/ /_/ /_/ /_/ /_/ /_/ /_/ /_/ /_/ /_/ /_/ /_/ /_/ /_/ /_/ /_/ /_/ /_/

