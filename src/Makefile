# Master makefile for POSTLAMMPS                                        
# Author: P. Bonnaud SEPTEMBER 2017                                     
#                                                                       
#=======================================================================
# Define default settings                                               
#=======================================================================
                                                                        
.SUFFIXES: .F90                                                         
.SUFFIXES: .f90                                                         
                                                                        
#=======================================================================
# Define object files                                                   
                                                                        
OBJS = ALLOCATE_CONFIG_ARRAYS.o ALLOCATE_CONFIG_ARRAY_OLD.o ALLOCATE_CONFIG_GLOBAL_ARRAYS.o \
ALLOCATE_CONFIG_LIBRARY_ARRAYS.o ALLOCATE_DUPLICATE_ARRAYS.o ALLOCATE_INSERT_LIBRARY_ARRAYS.o \
ALLOCATE_LIBRARY_ARRAYS.o ALLOCATE_PARAMETER_ARRAY_OLD.o ALLOCATE_PARAMETER_GLOBAL_ARRAYS.o \
APPLY_BUILD_CELLBOX.o APPLY_ION_INSERTION.o APPLY_MODIFICATION.o \
APPLY_PBC.o APPLY_PLATELET.o APPLY_REBUILD_SILICA_CHAINS.o \
APPLY_REBUILD_WATER.o APPLY_REBUILD_WATER_Hw.o APPLY_REPLICA_CREATION.o \
APPLY_TRANSLATION.o BUILD_INPUT_KEYWORD_LIST.o BUILD_POLYMER_CHAINS.o \
BUILD_PRIMITIVE_CELL.o BUILD_SUPER_CELL.o CHECK_DISTANCES.o \
CHECK_ELECTRONEUTRALITY.o CHECK_FILE_DUPLICATE.o CHECK_FILE_LIBRARY.o \
CHECK_INSERT_FILE_LIBRARY.o CHECK_INTRAMOLECULAR_BONDS.o CHECK_OVERLAPPING_ATOMS.o \
CHECK_SILICACHAINS.o CLOSING_PROGRAM.o COMPUTE_CONFIG_TOPOLOGY.o \
CONVERT_INT_TO_CHAR.o CONVERT_REAL_TO_CHAR.o DEALLOCATE_CONFIG_ARRAY_OLD.o \
DEALLOCATE_CONFIG_GLOBAL_ARRAYS.o DEALLOCATE_INPUT_KEYWORD_LIST.o DEALLOCATE_INSERT_COMMAND.o \
DEALLOCATE_INSERT_LIBRARY_ARRAYS.o DEALLOCATE_LIBRARY_ARRAYS.o DEALLOCATE_LMP_INPUT_COMMAND.o \
DEALLOCATE_PARAMETER_ARRAY_OLD.o DEALLOCATE_PARAMETER_GLOBAL_ARRAYS.o DEALLOCATE_REGIONS.o \
DUPLICATE_ASSIGN_FORCE_FIELD.o DUPLICATE_BUILD_SUPER_CELL.o DUPLICATE_READ_XYZ_TEMPLATE.o \
FIND_ATOM_ELEMENTID.o FIND_DIRECTORY.o FIND_MONOMER_RADICALS.o \
GENERATE_CONFIG_PREFIX.o INITIALIZATION_RANDOM_GENERATOR.o INSERT_LIBRARY_MOLECULE.o \
INSERT_SET_INTERATOMIC_POTENTIALS.o LIBRARY_INITIALIZATION_ARRAYS.o LMP_INPUT_UPDATE_GROUP_UNION.o \
MAKE_PROGRAM_TITLE.o MAKE_ROUTINE_TITLE.o MAKE_SUPCELLCONF.o \
MAKE_SWITCH_ATOMIDS.o MAKE_TETHERED_ATOMS.o MAKE_TRANSLATION_MOLECULEIDS.o \
MANA_ERROR.o RANDOM_DELETION_MOLECULE.o RANDOM_INSERTION_DUMMY_PARTICLE.o \
RANDOM_INSERTION_MOLECULE.o READ_CONFIG.o READ_DATA.o \
READ_FIELD.o READ_FILE_INFO_TEMPLATE.o READ_FILE_LIBRARY.o \
READ_INPUT.o READ_INTERATOMIC_POTENTIALS_TEMPLATE.o READ_LAMMPS_CONFIG.o \
READ_LAMMPS_CONFIG_LIBRARY.o READ_LAMMPS_TEMPLATE.o READ_XYZ_TEMPLATE.o \
SEEK_FOR_CREATE_BOX_COMMAND.o SEEK_FOR_INSERT_COMMAND.o SEEK_FOR_INSERT_MODIFY_COMMAND.o \
SEEK_FOR_LMP_INPUT_COMMAND.o SEEK_FOR_REGIONS.o SEEK_FOR_SET_LIB_COMMAND.o \
SET_ATOM_PROPERTIES.o SET_BOND_PROPERTIES.o SET_CENTROSYM.o \
SET_CONFIG.o SET_PASSAGE.o SET_PBC.o \
SET_ROTATION_MATRIX.o SET_SLAB_KEYWORDS.o SET_TRANSLATION.o \
SLAB.o TEMPLATE_FILE_SET_CONFIG_NAT.o TOOLS.o \
UPDATE_CONFIG_ARRAY.o UPDATE_CONFIG_ARRAY_OLD.o UPDATE_PARAMETER_ARRAY.o \
UPDATE_PARAMETER_ARRAY_OLD.o UPDATE_SIZE_OF_ARRAYS.o WRITE_CONFIG.o \
WRITE_DATA.o WRITE_ERROR_MESSAGE.o WRITE_FIELD.o \
WRITE_FINAL_CONFIG.o WRITE_LAMMPS_CONFIG.o WRITE_LAMMPS_INPUT.o \
WRITE_XYZ_CONFIG.o module_config.o module_data_in.o \
module_duplicate.o module_inserts.o module_library.o \
module_lmp_input.o module_osef.o module_physical_constants.o \
module_regions.o module_simulation_box.o module_size_arrays.o \
module_slab_keywords.o
                                                                      
#=====================================================================
                                                                      
MODULE = module_config.mod module_data_in.mod module_duplicate.mod\
module_inserts.mod module_library.mod module_lmp_input.mod\
module_osef.mod module_physical_constants.mod module_regions.mod\
module_simulation_box.mod module_size_arrays.mod module_slab_keywords.mod
                                                                      
#=====================================================================
# Define targets                                                      
                                                                      
all: SLAB_AMD.exe                                                     
                                                                      
SLAB_AMD.exe: $(MODULE) $(OBJS)                                       
	ifort -O3 -o SLAB_AMD.exe $(OBJS)                                   
                                                                      
module_config.mod: $(MODULE) module_config.F90
	ifort -c module_config.F90
                                                                      
module_data_in.mod: $(MODULE) module_data_in.F90
	ifort -c module_data_in.F90
                                                                      
module_duplicate.mod: $(MODULE) module_duplicate.F90
	ifort -c module_duplicate.F90
                                                                      
module_inserts.mod: $(MODULE) module_inserts.F90
	ifort -c module_inserts.F90
                                                                      
module_library.mod: $(MODULE) module_library.F90
	ifort -c module_library.F90
                                                                      
module_lmp_input.mod: $(MODULE) module_lmp_input.F90
	ifort -c module_lmp_input.F90
                                                                      
module_osef.mod: $(MODULE) module_osef.F90
	ifort -c module_osef.F90
                                                                      
module_physical_constants.mod: $(MODULE) module_physical_constants.F90
	ifort -c module_physical_constants.F90
                                                                      
module_regions.mod: $(MODULE) module_regions.F90
	ifort -c module_regions.F90
                                                                      
module_simulation_box.mod: $(MODULE) module_simulation_box.F90
	ifort -c module_simulation_box.F90
                                                                      
module_size_arrays.mod: $(MODULE) module_size_arrays.F90
	ifort -c module_size_arrays.F90
                                                                      
module_slab_keywords.mod: $(MODULE) module_slab_keywords.F90
	ifort -c module_slab_keywords.F90
                                                                      
                                                                      
.f90.o:                                                               
	ifort -c $<                                                         
                                                                      
.F90.o:                                                               
	ifort -c $<                                                         
                                                                      
clean:                                                                
	rm -f *.o *.mod                                                     
