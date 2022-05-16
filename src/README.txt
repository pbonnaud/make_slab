 1.78    * Modify the group for better taking into account the dynamics of 
           different groups 
         * Add a set_lib keyword in order to let the user define its own location 
           for library files


 1.77    * Modify the code for generating tethered particles for wall atoms in      X 
           moving wall simulaitons

 1.76    * Add a keyword that modifies properties of the insert keyword             X
         * Modified the insertion routine to correctly assign corss interactions    X
           for pair potentials
         * Modify the code for properly preparing the lammps input file             X

 1.75    * Modify the code for inserting molecules from a library through the usage X
           of keywords in the input file

 1.74    * Modify the code for generating molecular structure (duplication) through X
           the usage of keywords in the input file

 1.73    * Modify the program to generate several polymer molecules from pre-       X
           existing monomers in the xyz-format
         * Move the code part looking for the location of R1 and R1 in monomers in  X
           a separated routine

 1.72    * Modify the program to insert xyz files at defined locations of centers   X 
           of mass or defined location of geometric centers

 1.71    * MODIFY THE ROUTINE GENERATING POLYMER MOLECULES IN ORDER TO PUT A LIMIT  X 
           IN THE WHILE LOOP AVOIDING INFINITE LOOPS 

 1.70    * MODIFY THE PROGRAM TO GENERATE POLYMER MOLECULES FROM PRE-EXISTING       X
           MONOMERS IN THE XYZ FORMAT 
         * MODIFY THE PROGRAM TO GENERATE MOLECULAR CONFIGURATIONS WITH THE AMBER   X
           FORCE FIELD

 1.69    * MODIFY THE PROGRAM TO READ MOLECULE WITH PARAMETERS TAKEN FROM THE       X 
           PADUA'S FORCE FIELD
         * MERGE SAME MOLECULES WITH DIFFERENT ATOM IDs                             X
         /!\ Check the building of dihedral angles some parameters are not written

 1.68    *

 1.67    * MODIFY THE ROUTINE THAT REMOVES MOLECULES - DEFINED MORE PRECISE 
           CRITERIA IN ORDER TO BE ABLE TO DISCRIMINATE TWO MOLECULES WITH THE 
           SAME MOLECULAR WEIGHT

 1.66    * MODIFY THE PROGRAM TO INSERT MOLECULES WITH THE LAMMPS TEMPLATE IN A     X 
           GIVEN REGION OF THE SIMULATION BOX
         * CORRECTED A BUG RELATED TO THE INSERTION OF MOLECULES IN A NEW           X 
           SIMULATION BOX WITHOUT DEFINING SPECIFIC REGIONS

 1.65    * CORRECT A BUG IN THE INSERTION OF DUMMY PARTICLES AND EXTEND THE         X
           ALGORITHM TO INSERT DUMMY PARTICLES IN AN EXISTING, DENSE MOLECULAR 
           CONFIGURATION
         * CORRECTION OF A BUG WHEN READING THE PARAMETERS OF A TEMPLATED MOLECULE  X 
           WHILE REMOVING HOST MOLECULES IN THE MOLECULAR STRUCTURE

 1.64    * MODIFY THE PROGRAM TO INSERT DUMMY PARTICLES IN THE SIMULATION BOX       X

 1.63    * MODIFY THE PROGRAM IN ORDER TO BE ABLE TO REMOVE MOLECULES FROM THE      X 
           INPUT MOLECULAR STRUCTURE

 1.62    * MODIFY THE PROGRAM TO EXIT THE LOOP WHEN TWO ATOMS ARE TOO CLOSE IN THE  X 
           	ROUTINE OF RANDOM INSERTION OF MOLECULES 

 1.61    * MODIFY THE PROGRAM TO CREATE A TEST WHEN CHOSING THE LOCATION OF THE     X
           CENTER OF MASS IN ORDER TO SAVE TIME IN THE INSERTION OF MOLECULES

 1.60    * MODIFY THE PROGRAM TO CREATE A NEW CONFIGURATION FROM SIMULATION DATA    X

 1.59    * UPDATE THE SCRIPT TO GENERATE AUTOMATICALLY THE MAKE FILE                X

 1.58    * MODIFY THE PROGRAM TO ALLOCATE AUTOMATICALLY THE SIZE OF ARRAYS          X

 1.57    * MODIFY THE PROGRAM TO CREATE A MODULE CONTAINING THE SIZE OF ALL         X 
           THE ARRAYS

 1.56    * MODIFY THE PROGRAM TO CORRECT A BUG                                      X
         * ADD MODULES TO BETTER MANAGE ARRAYS                                      X

 1.55    * MODIFY THE PROGRAM TO CREATE AN INITIAL CONFIGURATION WITH FLEXIBLE      X
           KEROGEN MOLECULES

 1.54    * MODIFY THE PROGRAM TO READ MOLECULES IN THE LIBRARY THAT ARE IN THE      X
           LAMMPS FILE FORMAT

 1.53    * MODIFY THE PROGRAM TO BUILD A MOLECULAR CONFIGURATION FROM A LIBRARY OF  X
           MOLECULES

 1.52    *                                                                          X

 1.51    * MODIFY THE PROGRAM TO TAKE A CONFIGURATION COMING TO A PROGRAM           X
           DIFFERENT FROM PIBITO AND APPLY MODIFICATIONS ON IT 

 1.50    *                                                                          X

 1.49    * MODIFICATION OF INITIAL CONFIGURATION READING                            X
         * MODIFY SET_SUPERCELL TO CREATE A PLATELET PARTICLE (CYLINDER)            X
         * ADD DUPLICATION OF SUBSTRATE IN APPLY_REPLICA_CREATION                   X

 1.48    * READ CONFIGURATION COMMING FROM THE AMERICAN MINERALOGIST CRYSTAL        X
         * STRUCTURE DATABASE

 1.47    * SKIP READING CENTER OF MASSES WHEN THEY ARE NOT WRITTEN                  X 

 1.46    * MODIFY PROGRAM TO APPLY A  SINGLE TRANSLATION TO TWO EXISTING PARTICLES  X
         * MODIFY WRITE_FIELD.F90 TO WRITE THE UPDATED CENTER OF MASS
           !!! NEED TO CHECK SET_SUPERCELL FOR THE NEX UDE  !!!
           !!! THERE SHOULD BE A PROBLEM                    !!!

 1.45    * TAKE THE DIRECTION GIVEN BY THE CENTER OF MASS OF THE TWO PARTICULES     X
           AS THE DIRECTION OF TRANSLATION
