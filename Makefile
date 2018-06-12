FC = gfortran
#FC = ifort
#FFLAGS = -O3 -ipo
FFLAGS = -g

OBJS = 	\
	exact.o \
	flux.o \
	boundaries.o \
	output.o \
	parameters.o \
	conversion.o

TVDOBJS = \
	evolve.o \
	compute_flux.o
#	rhs.o \
#	reconstruct.o

MODS = 	real_type_mod.o \
	system_parameters_mod.o \
	grid_parameters_mod.o \
	method_parameters_mod.o \
	data_mod.o \
	grid_mod.o \
        reconstruct_mod.o \
        flux_splitting_mod.o \
	rhs_mod.o \
	sources_mod.o \
	box_mod.o \
	generic_list_mod.o \
	box_list_mod.o \
	grid_list_mod.o \
	level_mod.o
#	gnode_mod.o \
#	glist_mod.o \
#	level_mod.o \
#	box_mod.o \
#	boxnode_mod.o \
#	boxlist_mod.o
#	parameters_mod.o

#godunov: godunov.o $(OBJS) $(MODS)
#	$(FC) $(FFLAGS) -o godunov godunov.o $(OBJS) $(MODS)

#tvd: tvd.o $(TVDOBJS) $(OBJS) $(MODS)
#	$(FC) $(FFLAGS) -o tvd tvd.o $(TVDOBJS) $(OBJS) $(MODS)

amr: amr.o $(TVDOBJS) $(OBJS) $(MODS)
	$(FC) $(FFLAGS) -o amr amr.o $(TVDOBJS) $(OBJS) $(MODS)

clean:
	rm -f *.o *.mod

cleandata:
	rm -f *.dat *.errors

parameters_mod.o: real_type_mod.o
system_parameters_mod.o: real_type_mod.o
grid_parameters_mod.o: real_type_mod.o
method_parameters_mod.o: real_type_mod.o
data_mod.o: real_type_mod.o method_parameters_mod.o
grid_mod.o: real_type_mod.o system_parameters_mod.o data_mod.o box_mod.o box_list_mod.o
rhs_mod.o: grid_mod.o sources_mod.o
sources_mod.o: real_type_mod.o system_parameters_mod.o
grid_list_mod.o: grid_mod.o box_list_mod.o generic_list_mod.o
level_mod.o: grid_list_mod.o box_list_mod.o
box_list_mod.o: box_mod.o generic_list_mod.o
flux_splitting_mod.o: real_type_mod.o data_mod.o method_parameters_mod.o

$(OBJS): %.o: $(MODS)
$(TVDOBJS): %.o: $(MODS)

amr.o: amr.f90 $(MODS)
tvd.o: tvd.f90 $(MODS)
godunov.o: godunov.f90 $(MODS)

%.o: %.f90
	$(FC) -c $(FFLAGS) $< -o $@
