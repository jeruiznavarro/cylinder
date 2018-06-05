program Molecular_Dynamics_Simulator
	implicit none!all variables must be declared
	character(11),parameter::version='0.7.4.00001'!version of the program
	integer,parameter::double_kind=kind(0.d0)!to increase portability between compilers, the double precision kind is set in a parameter variable here
	integer(double_kind),parameter::maxbinrad=10000!maximum number of radial bins
	integer(double_kind),parameter::maxbinhist=100!maximum number of bins for the velocity histogram
	integer(double_kind),parameter::maxcyl=100!maximum number of cylinders
	integer(double_kind),parameter::maxpartic=50000!maximum number of particles in the system
	integer(double_kind),parameter::maxpartcell=1000!maximum number of particles in a single cell
	integer(double_kind),parameter::maxcells=20!maximum number of cells per dimension
	!real(double_kind),parameter::kboltz=1.d0!Boltzmann's constant
	real(double_kind),parameter::pi=dacos(-1.d0)!number pi
	logical::andersen=.false.!if true the andersen thermostat will be active
	logical::berendsen=.false.!if true the berendsen thermostat will be active
	logical::cylinders=.false.!if true, independent measurements will be taken in the different cylindrical intervals of the system
	logical::extforce=.false.!if true, an external force will be applied to all free particles in the system
	logical::firststep=.true.!after the first call to this subroutine is completed, this will be false
	logical::fixed(maxpartic)=.false.!this is true for fixed particles and false for free ones that can move
	logical::gaussdistribution=.false.!if false, a uniform distribution for the velocities will be used instead of a gaussian one
	logical::gnuplot=.false.!if true, the program will output data to be visualized in gnuplot
	logical::halfdensity=.false.!if true, the fluid in the porous system will have half of its particles erased to achieve a half of the previous density
	logical::interaction(maxpartic,maxpartic)=.false.!this matrix stores which interactions have been calculated in each timestep so they don't have to be computed twice
	logical::measuring=.false.!if true, measurements will be taken
	logical::measurediffusion=.false.!if true, measurements of the mean square displacement of the particles will be taken
	logical::measureprofile=.false.!if true, measurements of the velocity profile in a cylindrical area will be taken
	logical::measureraddist=.false.!if true, the radial distribution function will be measured
	logical::measurevhist=.false.!if true, a histogram of the components of the velocities of the particles will be taken each time a measurement is performed
	logical::noinitialinterac=.false.!if false, the initial interactions between all the particles won't be computed
	logical::normaloutput=.false.!if true, measurements of the different energies, temperature and pressure will be taken
	logical::nosehoover=.false.!if true, the Nosé-Hoover thermostat will be used to sample the NVT ensenmble
	logical::periodicspeed=.false.!if true, no matter how fast the particles move, the periodic boundary conditions will apply, but this might make some bugs much harder to detect (and uses a tiny bit more resources)
	logical::press3d=.false.!if true, a measurement of the pressure will be taken in each cell of the system, temp3d is needed to be true as well
	logical::progress=.true.!if true, the progress of the simulation will be printed out in the screen if the gnuplot output is not active
	logical::removedrift=.false.!if true, the average momentum of the system will be eliminated so there is no global drift
	logical::samemass=.false.!if true, all the particles have the same mass
	logical::saving=.false.!if true, the state of the system will be saved
	logical::separateoutput=.false.!if true, separate files with the current date and time will be used for the output of the measurements
	logical::temp3d=.false.!if true, a measurement of the temperature will be taken in each cell of the system
	logical::timetomeasure=.false.!if true, a measurement of the properties of the system will be performed
	logical::title=.false.!if at least one particle is fixed, a special title will be displayed in each gnuplot frame
	logical::xyz=.false.!if true, visualization data for VMD will be kept in a file
	character(8)::datechar=''!this character variable stores the current date for the names of the output files
	character(2)::element=''!this is used for the VMD output, to tell the program which element the atoms belong to
	character(3)::struc=''!this variable will tell the program which structure should be generated, fcc for cubic face centered, fcy is the same but has a cylinder of free particles, scc for simple cubic and por for a porous matrix
	character(10)::timechar=''!this character variable stores the current time for the names of the output files
	integer(double_kind)::amntfix=0!amount of fixed particles in the simulation
	integer(double_kind)::amntfree=0!amount of free particles in the simulation
	integer(double_kind)::amnttot=0!amount of total particles in the simulation
	integer(double_kind)::binhist=0!number of bins for the velocity histograms
	integer(double_kind)::binrad=0!number of bins for the radial distribution
	integer(double_kind)::distbins(maxbinrad)=0!these bins will count how many particles are in the divisions of the (0,L/2] interval in which the radial distribution will be measured
	integer(double_kind)::cells(maxpartcell,maxcells,maxcells,maxcells)=0!lists containing the tags of the particles in each cell, the first index are the particles, the second one's the x axis, the third one's y and the last one's z
	integer(double_kind)::cvcount=0!this counter will keep track of how many measurements are used to calculate the heat capacity so the time average can be computed
	integer(double_kind)::cvcountcy(maxcyl)=0!this counter will keep track of how many measurements are used to calculate the heat capacity for each cylinder so the time average can be computed
	integer(double_kind)::difftag(maxcyl)=0!this variable is the tag of the particle that will be followed to measure diffusion properties in each cylinder
	integer(double_kind)::iterations=0!number of iterations carried out by the program
	integer(double_kind)::itermeas=0!after this amount of iterations, measurements will start
	integer(double_kind)::maxiter=0!maxmimum number of iterations
	integer(double_kind)::meascounter=0!number of measurements taken
	integer(double_kind)::meascylin=0!this is the number of cylindrical intervals in which there will be independent measurements of variables
	integer(double_kind)::measiter=0!number of iterations between measurements
	integer(double_kind)::numcells(3)=0!number of cells for each dimension (x=1, y=2, z=3)
	integer(double_kind)::occup(maxcells,maxcells,maxcells)=0!number of particles in each cell, occupancy lists
	integer(double_kind)::presscount=0!this counter keeps track of how many timesteps are used to obtain the pressure
	integer(double_kind)::pores=0!amount of pores in the system
	integer(double_kind)::radcount=0!this counts how many timesteps are used in the measurement of the radial distribution
	integer(double_kind)::seed(8)=0!seed for the random number generator
	real(double_kind)::avec(3,maxpartic)=0.d0!accelerations for each particle (ax->(1,:), ay->(2,:), az->(3,:)) are here
	real(double_kind)::cellsize(3)=0.d0!size of the cells for each dimension (x=1, y=2, z=3)
	real(double_kind)::cellvolume=0.d0!volume of a single cell
	real(double_kind)::convfactorenergy=0.d0!conversion factor from energy simulation units to SI units
	real(double_kind)::convfactorlength=0.d0!conversion factor from length simulation units to SI units
	real(double_kind)::convfactormass=0.d0!conversion factor from mass simulation units to SI units
	real(double_kind)::convfactorpressure=0.d0!conversion factnormaloutput=.false.!if true, measurements of the different energies, temperature, pressure and mean square displacement will be takenor from pressure simulation units to SI units
	real(double_kind)::convfactortemperature=0.d0!conversion factor from temperature simulation units to SI units
	real(double_kind)::convfactortime=0.d0!conversion factor from time simulation units to SI units
	real(double_kind)::cv=0.d0!heat capacity of the system
	real(double_kind)::cvcy(maxcyl)=0.d0!heat capacity of each cylinder
	real(double_kind)::dens=0.d0!number density of the system
	real(double_kind)::denscells(maxcells,maxcells,maxcells)=0.d0!number density of each cell
	real(double_kind)::displ(3)=0.d0!displacement vector to calculate the appropiate distance between particles with the periodic boundary conditions
	real(double_kind)::endsim=0.d0!this variable is used to calculate how long the simulation takes to finish
	real(double_kind)::ekine=0.d0!kinetic energy of the system
	real(double_kind)::ekinecells(maxcells,maxcells,maxcells)=0.d0!kinetic energy of each cell
	real(double_kind)::ekinecy(maxcyl)=0.d0!kinetic energy of each cylinder in the system
	real(double_kind)::epot=0.d0!potential energy of the system
	real(double_kind)::epotcy(maxcyl)=0.d0!potential energy of each cylinder
	real(double_kind)::epsi=0.d0!depth parameter of the Lennard-Jones interaction
	real(double_kind)::fvec(3,maxpartic)=0.d0!forces for each particle (fx->(1,:), fy->(2,:), fz->(3,:)) are here
	real(double_kind)::intrad=0.d0!interaction radius, cut-off for the potential, beyond this distance there is no interaction
	real(double_kind)::l2=0.d0!this is half of the smallest side of the simulation box and is used as the limit to look for particles when calculating the radial distribution
	real(double_kind)::latcon=0.d0!lattice constant
	real(double_kind)::mass=0.d0!generic mass of the particles
	real(double_kind)::mvec(maxpartic)=0.d0!if the masses of the particles are not the same they can be stored here
	real(double_kind)::measintv=0.d0!interval of time between measurements
	real(double_kind)::nosehooverconstant=0.d0!coupling constant of the Nosé-Hoover thermostat
	real(double_kind)::outforce(3)=0.d0!external force applied to all the particles in the system
	real(double_kind)::porevolume=0.d0!the volume occupied by the pored if they were non-overlapping
	real(double_kind)::porosity=0.d0!relative fraction of non-solid space in the system
	real(double_kind)::potentialshift=0.d0!shift of the interaction potential
	real(double_kind)::press=0.d0!pressure of the system
	real(double_kind)::presscells(maxcells,maxcells,maxcells)=0.d0!pressure in each cell
	real(double_kind)::profile(maxcyl)=0.d0!the velocity profile as a function of radius is kept here
	real(double_kind)::radius=0.d0!radius of the non solid zone of the system where particles can flow
	real(double_kind)::rvec(3,maxpartic)=0.d0!positions for each particle (x->(1,:), y->(2,:), z->(3,:)) are here
	real(double_kind)::rvecini(3,maxpartic)=0.d0!initial positions for each particle (x->(1,:), y->(2,:), z->(3,:)) are here
	real(double_kind)::rvecnonper(3,maxpartic)=0.d0!positions for each particle (x->(1,:), y->(2,:), z->(3,:)) without periodic boundary conditions (to calculate mean square displacement) are here
	real(double_kind)::sigm=0.d0!size parameter of the Lennard-Jones interaction
	real(double_kind)::siz(3)=0.d0!dimensions (x=1, y=2, z=3) of the simulation box
	real(double_kind)::squaredispl=0.d0!mean square displacement used to calculate the diffusion constant
	real(double_kind)::sumekine=0.d0!the sum of the kinetic energy over time is kept here
	real(double_kind)::sumekinesqur=0.d0!the sum of the kinetic energy squared over time is kept here
	real(double_kind)::sumekinecy(maxcyl)=0.d0!the sum of the kinetic energy over time for each cylinder is kept here
	real(double_kind)::sumekinecysqur(maxcyl)=0.d0!the sum of the kinetic energy squared over time for each cylinder is kept here
	real(double_kind)::startsim=0.d0!this variable is used to keep track of when the simulation started
	real(double_kind)::temp=0.d0!temperature of the system
	real(double_kind)::tempbath=0.d0!temperature of the heat bath
	real(double_kind)::tempcells(maxcells,maxcells,maxcells)=0.d0!temperature in each cell
	real(double_kind)::tempmax=0.d0!maximum temperature during equilibrium after termalization
	real(double_kind)::tempmean=0.d0!mean temperature
	real(double_kind)::tempmin=1.d100!minimum temperature during equilibrium after termalization
	real(double_kind)::tempcy(maxcyl)=0.d0!temperature of each cylinder in the system
	real(double_kind)::threshold=0.d0!this variable is used to test if a particle will be affected by the andersen thermostat
	real(double_kind)::time=0.d0!time for which the simulation has been running
	real(double_kind)::timetosave=0.d0!the state of the system will be saved after this point in time
	real(double_kind)::tmax=0.d0!time when the simulation should stop
	real(double_kind)::tmeas=0.d0!when time is greater or equal than this, measurements will start
	real(double_kind)::trelax=0.d0!relaxation time for the andersen and berendsen thermostats
	real(double_kind)::tstep=0.d0!timestep of the simulation
	real(double_kind)::volume=0.d0!volume of the system
	real(double_kind)::vvec(3,maxpartic)=0.d0!velocities for each particle (vx->(1,:), vy->(2,:), vz->(3,:)) are here


	call input()!reading initial parameters
	if(progress.and.(.not.gnuplot))then
		call cpu_time(startsim)!the simulation starts
		call date_and_time(date=datechar,time=timechar,values=seed)!getting the date and the time
		write(*,*) 'Simulation started on ',datechar//' at '//timechar!showing work in progress
	end if
	call startup()!setting up the system
	do!simulation loop
		if(iterations>maxiter)exit!the program ends when it reaches the maximum time
		if(((mod(iterations,measiter)==0).and.(iterations>=itermeas)).and.(measuring.and.(.not.gnuplot)))then
			timetomeasure=.true.!checking if a measurment should be taken
		else
			timetomeasure=.false.!checking if a measurment should be taken
		end if
		if(nosehoover)then
			call nose_hoover()!using the Nosé-Hoover thermostat to carry out all the important parts of the timestep
		else
			call half_leap(.false.)!updating the velocities of the particles without the thermostat, because it should only be applied to the speeds once per timestep
			call move_system()!updating the positions of the particles
			call update_particles_in_cells()!after the particles have moved, they may have entered a different cell and this needs to be checked and taken care of
			call calculate_forces()!the Verlet integrator requires to update the speeds again after moving the particle so the forces need to be calculated again
			call half_leap(.true.)!updating the velocities of the particles with the thermostat (if it's active)
		end if
		if(saving.and.(dabs(time-timetosave)<=tstep*1.d-1))call save_state()!storing the generated structure in the apropiate file
		if(xyz)call xyz_output()!storing the output for VMD visualization
		if(cylinders.and.timetomeasure)then
			call measuring_with_cylinders()!with a certain frequency measurements are taken in each radial interval of the cylindrical fluid cavity of the system
		else if(timetomeasure)then
			call normal_measuring()!with a certain frequency measurements are performed
		end if
		time=time+tstep!time advances
		iterations=iterations+1!one more iteration
		firststep=.false.!this is no longer true after the first iteration
		if((progress.and.(.not.gnuplot)).and.(mod(iterations,maxiter/100)==0))then
			call date_and_time(date=datechar,time=timechar,values=seed)!getting the date and the time
			write(*,*) nint(time*1.d2/tmax,double_kind),'% completed on ',datechar//' at '//timechar!showing work in progress
		end if
	end do
	if(measuring)call output()!finishing up the writing to file of the results
	if(progress.and.(.not.gnuplot))then
		call cpu_time(endsim)!the simulation ends
		call date_and_time(date=datechar,time=timechar,values=seed)!getting the date and the time
		write(*,*) 'Total execution time: ',endsim-startsim,' finished on ',datechar//' at '//timechar!showing work in progress
	end if
	stop!simulation is over



contains



	subroutine input()!reading initial parameters
		open(10,file='input.dat',action='read')!all the initial parameters are in this file
		read(10,*) andersen!if true the andersen thermostat will be active
		read(10,*) berendsen!if true the berendsen thermostat will be active
		read(10,*) binhist!number of bins for the velocity histograms
		read(10,*) binrad!number of bins for the radial distribution
		read(10,*) convfactorenergy!conversion factor from energy simulation units to SI units
		read(10,*) convfactorlength!conversion factor from length simulation units to SI units
		read(10,*) convfactormass!conversion factor from mass simulation units to SI units
		read(10,*) convfactorpressure!conversion factor from pressure simulation units to SI units
		read(10,*) convfactortemperature!conversion factor from temeprature simulation units to SI units
		read(10,*) convfactortime!conversion factor from time simulation units to SI units
		read(10,*) cylinders!if true, independent measurements will be taken in the different cylindrical intervals of the system
		read(10,*) element!this is used for the VMD output, to tell the program which element the atoms belong to
		read(10,*) epsi!depth parameter of the Lennard-Jones interaction
		read(10,*) extforce!if true, an external force will be applied to all free particles in the system
		read(10,*) gaussdistribution!if false, a uniform distribution for the velocities will be used instead of a gaussian one
		read(10,*) gnuplot!if true, the program will output data to be visualized in gnuplot
		read(10,*) halfdensity!if true, the fluid in the porous system will have half of its particles erased to achieve a half of the previous density
		read(10,*) intrad!interaction radius
		read(10,*) latcon!lattice constant
		read(10,*) mass!generic mass of the particles
		read(10,*) meascylin!this is the number of cylindrical intervals in which there will be independent measurements of variables
		read(10,*) measintv!interval of time between measurements
		read(10,*) measuring!if true, measurements will be taken
		read(10,*) measurediffusion!if true, measurements of the mean square displacement of the particles will be taken
		read(10,*) measureprofile!if true, measurements of the velocity profile in a cylindrical area will be taken
		read(10,*) measureraddist!if true, the radial distribution function will be measured
		read(10,*) measurevhist!if true, a histogram of the components of the velocities of the particles will be taken each time a measurement is performed
		read(10,*) noinitialinterac!if false, the initial interactions between all the particles won't be computed
		read(10,*) normaloutput!if true, measurements of the different energies, temperature, pressure and mean square displacement will be taken
		read(10,*) nosehoover!if true, the Nosé-Hoover thermostat will be used to sample the NVT ensenmble
		read(10,*) nosehooverconstant!coupling constant of the Nosé-Hoover thermostat
		read(10,*) outforce(1)!x component of the external force applied to all the particles in the system
		read(10,*) outforce(2)!y component of the external force applied to all the particles in the system
		read(10,*) outforce(3)!z component of the external force applied to all the particles in the system
		read(10,*) periodicspeed!if true, no matter how fast the particles move, the periodic boundary conditions will apply, but this might make some bugs much harder to detect (and uses a tiny bit more resources)
 		read(10,*) pores!amount of pores in the system
 		read(10,*) press3d!if true, a measurement of the pressure will be taken in each cell of the system, temp3d is needed to be true as well
 		read(10,*) progress!if true, the progress of the simulation will be printed out in the screen if the gnuplot output is not active
		read(10,*) radius!radius of the non solid zone of the system where particles can flow
		read(10,*) removedrift!if true, the average momentum of the system will be eliminated so there is no global drift
		read(10,*) samemass!if true, all the particles have the same mass
		read(10,*) saving!if true, the state of the system will be saved
		read(10,*) separateoutput!if true, separate files with the current date and time will be used for the output of the measurements
		read(10,*) sigm!size parameter of the Lennard-Jones interaction
		read(10,*) siz(1)!size in the x dimension of the simulation box
		read(10,*) siz(2)!size in the y dimension of the simulation box
		read(10,*) siz(3)!size in the z dimension of the simulation box
		read(10,*) struc!this variable will tell the program which structure should be generated, fcc for cubic face centered, fcy is the same but has a cylinder of free particles, scc for simple cubic, por for a porous matrix and two for a system with just particles near the centre (useful for debugging)
		read(10,*) temp!temperature of the system
		read(10,*) tempbath!temperature of the heat bath
		read(10,*) temp3d!if true, a measurement of the temperature will be taken in each cell of the system
		read(10,*) timetosave!the state of the system will be saved after this point in time
		read(10,*) tmax!time when the simulation should stop
		read(10,*) tmeas!when time is greater or equal than this, measurements will start
		read(10,*) trelax!relaxation time for the andersen and berendsen thermostats
		read(10,*) tstep!timestep of the simulation
		read(10,*) xyz!if true, visualization data for VMD will be kept in a file
		close(10)!disconnecting the "intput.dat" file
		if((cylinders.or.measurediffusion.or.measureprofile.or.measureraddist.or.measurevhist.or.normaloutput.or.press3d.or.separateoutput.or.temp3d).and.(.not.measuring))measuring=.true.!making sure that there are no conflicts with these logical variables
		if((berendsen.and.andersen).or.(nosehoover.and.andersen).or.(berendsen.and.nosehoover).or.((berendsen.and.andersen).and.nosehoover))stop "More than one thermostat active at once, aborting."!this conflict is not allowed, so the program is terminated
		if((meascylin>maxcyl).and.(cylinders))stop 'More cylinders than allowed, the variable "maxcyl" should be at least as big as "meascylin", aborting.'!this conflict is not allowed, so the program is terminated
		if(progress.and.gnuplot)stop 'Two conflicting variables are true at the same time: "gnuplot" and "progress", aborting.'
		return
	end subroutine input



	subroutine startup()!the system to be simulated is set up here
		call date_and_time(date=datechar,time=timechar,values=seed)!creating a seed using the clock of the computer
		call dran_ini(1000*seed(8)+3*seed(7)*seed(6)/10)!initializing the random number generator
		if(separateoutput)then
			if(xyz)open(20,file='./data/position_output_'//datechar//'_'//timechar//'.xyz',action='write',status='new')!if the positions are to be visualized later they must be saved to this file
			if(normaloutput)open(30,file='./data/energy_temperature_pressure_diffusion_output_'//datechar//'_'//timechar//'.dat',action='write',status='new')!this is where the time evolution of these magnitudes will be kept
			if(measureraddist)open(40,file='./data/radial_distribution_'//datechar//'_'//timechar//'.dat',action='write',status='new')!and the radial distribution function values are here
			if(measurevhist)open(60,file='./data/velocity_histogram_'//datechar//'_'//timechar//'.gnu',action='write',status='new')!the velocity histogram is kept here in a format such that it can be visualized as an animation in gnuplot
			if(temp3d)open(70,file='./data/temperature_3d_output_'//datechar//'_'//timechar//'.dat',action='write',status='new')!the 3d data of the temperature
			if(press3d.and.temp3d)open(80,file='./data/pressure_3d_output_'//datechar//'_'//timechar//'.dat',action='write',status='new')!the 3d data of the pressure
			if(measureprofile)open(90,file='./data/velocity_profile_'//datechar//'_'//timechar//'.dat',action='write',status='new')!the velocity profile data is saved here
		else
			if(xyz)open(20,file='./data/position_output.xyz',action='write')!if the positions are to be visualized later they must be saved to this file
			if(normaloutput)open(30,file='./data/energy_temperature_pressure_diffusion_output.dat',action='write')!this is where the time evolution of these magnitudes will be kept
			if(measureraddist)open(40,file='./data/radial_distribution.dat',action='write')!and the radial distribution function values are here
			if(measurevhist)open(60,file='./data/velocity_histogram.gnu',action='write')!the velocity histogram is kept here in a format such that its evolution can be visualiazed as an animation in gnuplot
			if(temp3d)open(70,file='./data/temperature_3d_output.dat',action='write')!the 3d data of the temperature
			if(press3d.and.temp3d)open(80,file='./data/pressure_3d_output.dat',action='write')!the 3d data of the pressure
			if(measureprofile)open(90,file='./data/velocity_profile.dat',action='write')!the velocity profile data is saved here
		end if
		if(xyz)then
			write(20,*) amnttot!the xyz format needs to have the amount of particles in the first line
			write(20,*) 'This is an optional comment, but it must be here for the file to be readable by the visualizing program.'!and this comment is mandatory
		end if
		if(normaloutput)write(30,*) ' #			time					total energy			potential energy		kinetic energy			temperature						pressure				mean square displacement/6'!magnitudes in this file
		if(measurevhist)then
			write(60,*) "plot '-' lc rgb "//'"#FFFFFF"'!this will allow to resize the window in gnuplot to get a plot with better resolution, first it's necessary to plot something random and disposable in the graph, so plotting a single white point is the best idea
			write(60,*) '0.5 0.5'!the coordinates of said point
			write(60,*) 'e'!end of input for the point
			write(60,*) 'set term'!now gnuplot will ask the user to select a terminal, pressing enter once, or maybe twice, will allow the user to resize the plotting window to the desired size
			write(60,*) 'set key top right'!placing the title in a place where it won't be much of a nuissance
			write(60,*) 'set grid'!the grid makes it easier to see the data
			write(60,*) 'set boxwidth 1 relative'!this is the width of the histogram columns
			write(60,*) 'set style fill transparent solid 0.5 noborder'!and this is the style of the columns, so they are half opaque and the grid can be visible behind them
		end if
		if(temp3d)then
			write(70,*) 'set grid'!the grid makes it easier to see the data
			write(70,*) 'set ticslevel 0'!setting the bottom of the graph at z=0
		end if
		if(press3d.and.temp3d)then
			write(80,*) 'set grid'!the grid makes it easier to see the data
			write(80,*) 'set ticslevel 0'!setting the bottom of the graph at z=0
		end if
		if(xyz)call flush(20)!forcing the program to write stuff
		if(normaloutput)call flush(30)!forcing the program to write stuff
		if(measurevhist)call flush(60)!forcing the program to write stuff
		if(temp3d)call flush(70)!forcing the program to write stuff
		if(press3d.and.temp3d)call flush(80)!forcing the program to write stuff
		if(gnuplot)then
			write(*,*) 'unset key'!this just gets in the way
			write(*,*) 'set border 4095'!gnuplot initial configuration, displaying the borders of the simulation box
			write(*,*) 'set ticslevel 0'!gnuplot initial configuration, setting the bottom of the simulation box at z=0
			write(*,*) 'set ztics mirror'!gnuplot initial configuration, setting the grid in the vertical axes
			write(*,*) 'set xrange [0:',siz(1),']'!gnuplot initial configuration, x size of the simulation box
			write(*,*) 'set yrange [0:',siz(2),']'!gnuplot initial configuration, y size of the simulation box
			write(*,*) 'set zrange [0:',siz(3),']'!gnuplot initial configuration, z size of the simulation box
			write(*,*) "splot '-' w p pt 7 lc rgb variable"!new gnuplot frame
		end if
		call calculate_some_constants()!getting some values for some variables
		if(struc=='fcc')then
			call create_fcc_structure()!face centered cubic structure
		else if(struc=='fcy')then
			call create_fcy_structure()!face centered cubic structure with a cylinder of free particles
		else if(struc=='por')then
			call create_porous_structure()!porous structure
		else if(struc=='scc')then
			call create_sc_structure()!simple cubic structure
		else if(struc=='two')then
			call create_two_particles()!just_two_particles_to_debug
		else
			stop 'No defined structure, aborting.'!the program needs to know which system should be generated
		end if
		if(gnuplot)write(*,*) 'e'!end of gnuplot frame
		call create_cell_lists()!cell lists are needed to calculate the interactions and to speed up the simulation
		if(measurevhist)call velocity_histogram()!getting a first snapshot of the velocity distribution
		if(removedrift)call remove_average_momentum()!the system should't have a global net movement so it is removed here if there's any (previously loaded systems are usually driftless,that's why this function can be turned off)
		if(noinitialinterac)call calculate_forces()!before the first integration iteration can be carried out, the forces need to be calculated
		return
	end subroutine startup



	subroutine calculate_some_constants()!just making the startup subroutine cleaner
		l2=dmin1(siz(1),siz(2),siz(3))*5.d-1!the value of L/2 is chosen to be half of the smallest side of the simulation box
		itermeas=nint(tmeas/tstep,double_kind)!iteration threshold after which there will be measurements
		maxiter=nint(tmax/tstep,double_kind)!maximum number of iterations in the simulation
		measiter=nint(measintv/tstep,double_kind)!number of iterations between measurements
		nosehooverconstant=1.d0/nosehooverconstant!this constant is always dividing, so it's better to calculate this and just multiply all the time
		if(samemass.and.extforce)then
			outforce(1)=outforce(1)/mass!if the masses are all the same, the force has to be converted to an acceleration for all three components
			outforce(2)=outforce(2)/mass!if the masses are all the same, the force has to be converted to an acceleration for all three components
			outforce(3)=outforce(3)/mass!if the masses are all the same, the force has to be converted to an acceleration for all three components
		end if
		potentialshift=4.d0*(intrad**(-6)-1.d0)*intrad**(-6)!potential shift so that the potential energy is a continious function after the cut-off (if constants are needed, multiply by epsi)
		threshold=tstep/trelax!this is the value that the threshold for the andersen thermostat should have
		volume=siz(1)*siz(2)*siz(3)!volume of the system
		return
	end subroutine calculate_some_constants



	subroutine create_fcc_structure()!a face centered cubic structure will be generated here
		integer(double_kind)::i1=0!x dimension counter
		integer(double_kind)::i2=0!y dimension counter
		integer(double_kind)::i3=0!z dimension counter
		integer(double_kind)::j=1!particle counter
		integer(double_kind)::k=0!dimension counter
		integer(double_kind)::l=0!auxiliary particle counter
		real(double_kind)::dran_g!gaussian random number
		real(double_kind)::dran_u!uniform random number
		do i1=0,nint(siz(1)/latcon,double_kind)-1
			do i2=0,nint(siz(2)/latcon,double_kind)-1
				do i3=0,nint(siz(3)/latcon,double_kind)-1
					rvec(1,j)=latcon*dble(i1)!x position for the first particle of the crystalline cell
					rvec(2,j)=latcon*dble(i2)!y position for the first particle of the crystalline cell
					rvec(3,j)=latcon*dble(i3)!z position for the first particle of the crystalline cell
					rvec(1,j+1)=latcon*(5.d-1+dble(i1))!x position for the second particle of the crystalline cell
					rvec(2,j+1)=latcon*(5.d-1+dble(i2))!y position for the second particle of the crystalline cell
					rvec(3,j+1)=latcon*dble(i3)!z position for the second particle of the crystalline cell
					rvec(1,j+2)=latcon*dble(i1)!x position for the third particle of the crystalline cell
					rvec(2,j+2)=latcon*(5.d-1+dble(i2))!y position for the third particle of the crystalline cell
					rvec(3,j+2)=latcon*(5.d-1+dble(i3))!z position for the third particle of the crystalline cell
					rvec(1,j+3)=latcon*(5.d-1+dble(i1))!x position for the fourth particle of the crystalline cell
					rvec(2,j+3)=latcon*dble(i2)!y position for the fourth particle of the crystalline cell
					rvec(3,j+3)=latcon*(5.d-1+dble(i3))!z position for the fourth particle of the crystalline cell
					if(gnuplot)then
						do k=0,3
							write(*,*) rvec(1,j+k),rvec(2,j+k),rvec(3,j+k),'0x000000'!data feed to gnuplot
						end do
					end if
					if(samemass)then
						if(gaussdistribution)then
							do l=0,3
								do k=1,3
									vvec(k,j+l)=dsqrt(temp/mass)*dran_g()!random gaussian components for the velocities (if constants needed use kboltz*temp/mass instead of temp/mass)
									rvecini(k,j+l)=rvec(k,j+l)!storing the initial position of the particles
								end do
							end do
						else
							do l=0,3
								do k=1,3
									vvec(k,j+l)=3.d0*dsqrt(temp/mass)*(dran_u()-5.d-1)!random uniform components for the velocities (if constants needed use kboltz*temp/mass instead of temp/mass)
									rvecini(k,j+l)=rvec(k,j+l)!storing the initial position of the particles
								end do
							end do
						end if
					else
						mvec(j)=mass*dsin(dble(j))!setting the mass of the first particle
						mvec(j+1)=mass*dsin(dble(j))!setting the mass of the second particle
						mvec(j+2)=mass*dsin(dble(j))!setting the mass of the third particle
						mvec(j+3)=mass*dsin(dble(j))!setting the mass of the fourth particle
						if(gaussdistribution)then
							do l=0,3
								do k=1,3
									vvec(k,j+l)=dsqrt(temp/mvec(j+l))*dran_g()!random gaussian components for the velocities (if constants needed use kboltz*temp/mass instead of temp/mvec(j+l))
									rvecini(k,j+l)=rvec(k,j+l)!storing the initial position of the particles
								end do
							end do
						else
							do l=0,3
								do k=1,3
									vvec(k,j+l)=3.d0*dsqrt(temp/mvec(j+l))*(dran_u()-5.d-1)!random uniform components for the velocities (if constants needed use kboltz*temp/mass instead of temp/mvec(j+l))
									rvecini(k,j+l)=rvec(k,j+l)!storing the initial position of the particles
								end do
							end do
						end if
					end if
					j=j+4!next particles
				end do
			end do
		end do
		amnttot=j-1!total amount of particles
		if(amnttot>maxpartic)then
			write(*,*) 'Not enough space for so many particles, the variable "maxpartic" should be at least as big as "amnttot", ',amnttot,' in this case. Aborting.'!avoiding segmentation faults
			stop!terminating execution
		end if
		amntfree=amnttot!just in case, to avoid problems since there are no fixed particles in this situation
		dens=amnttot/volume!this is the number density of all the particles
		return
	end subroutine create_fcc_structure



	subroutine create_fcy_structure()!a face centered cubic structure with a cylinder of free particles will be generated here
		integer(double_kind)::i=0!particle counter
		integer(double_kind)::j=0!dimension counter
		call load_state()!the pores are created in an already thermalized system that is in a fluid, preferably liquid, state, so the particles are loaded from a previously prepared file
		amntfree=amnttot!right now all the particles are free
		if(gnuplot)write(*,*) "set title 'Fixed particles are red and free ones are black.'"!information about the representation of free and fixed particles
		title=.true.!this same title has to be displayed on each frame so setting this variable to true will make sure that happens
		do i=1,amnttot
			do j=1,3
				rvecini(j,i)=rvec(j,i)!storing the initial position of the particles
			end do
			if(dsqrt((rvec(1,i)-siz(1)*5.d-1)**2+(rvec(2,i)-siz(2)*5.d-1)**2)>=radius)then
				fixed(i)=.true.!fixing particles inside the pores
				if(gnuplot)write(*,*) rvec(1,i),rvec(2,i),rvec(3,i),'0xFF0000'!data feed to gnuplot for fixed particles
				do j=1,3
					vvec(j,i)=0.d0!now the particle is fixed so the velocity is not needed anymore
					avec(j,i)=0.d0!same with the acceleration (if constants needed replace avec with fvec)
				end do
				amntfix=amntfix+1!one more fixed particle
				amntfree=amntfree-1!one less free particle
			else
				if(gnuplot)write(*,*) rvec(1,i),rvec(2,i),rvec(3,i),'0x000000'!data feed to gnuplot for free particles
			end if
		end do
		if(halfdensity)call halving_density()!reducing the amount of particles by half
		dens=amnttot/volume!this is the number density of all the particles
		return
	end subroutine create_fcy_structure



	subroutine create_porous_structure()!a structure composed of a fcc crystal with a matrix of empty cylinders is generated here
		integer(double_kind)::i=0!particle counter
		integer(double_kind)::j=0!pore counter
		integer(double_kind)::k=0!dimension counter
		real(double_kind)::dran_u!uniform random number
		real(double_kind)::porematrix(3,maxcyl)!the coordinates (x->(1,:), y->(2,:), z->(3,:)) of the centres of the pores will be kept here
		real(double_kind)::poreradius(maxcyl)!the radius of each pore
		do j=1,pores
			do k=1,3
				porematrix(k,j)=dran_u()*siz(k)!random positions for the pores
			end do
			poreradius(j)=dran_u()+5.874d0!random radii for the pores
			porevolume=porevolume+4.d0*pi*poreradius(j)/3.d0!adding another non-overlapping volume to the total
		end do
		call load_state()!the pores are created in an already thermalized system that is in a fluid, preferably liquid, state, so the particles are loaded from a previously prepared file
		amntfree=amnttot!right now all the particles are free
		if(gnuplot)write(*,*) "set title 'Fixed particles are red and free ones are black.'"!information about the representation of free and fixed particles
		title=.true.!this same title has to be displayed on each frame so setting this variable to true will make sure that happens
		do i=1,amnttot
			do j=1,pores
				if(dsqrt((rvec(1,i)-porematrix(1,j))**2+(rvec(2,i)-porematrix(2,j))**2+(rvec(3,i)-porematrix(3,j))**2)<=poreradius(j))then
					fixed(i)=.true.!fixing particles inside the pores
					amntfix=amntfix+1!one more fixed particle
					amntfree=amntfree-1!one less free particle
					exit!this particle is already fixed, moving onto the next one
				end if
			end do
			if(gnuplot)then
				if(fixed(i))then
					write(*,*) rvec(1,i),rvec(2,i),rvec(3,i),'0xFF0000'!data feed to gnuplot for fixed particles
				else
					write(*,*) rvec(1,i),rvec(2,i),rvec(3,i),'0x000000'!data feed to gnuplot for free particles
				end if
			end if
		end do
		if(halfdensity)call halving_density()!reducing the amount of particles by half
		porosity=amntfree/amnttot!approximately finding out how much free space is there in the system, exploiting the fact that the particles are homogeneously distributed in a fluid state
		dens=amnttot/volume!this is the number density of all the particles
		return
	end subroutine create_porous_structure



	subroutine create_sc_structure()!a simple cubic structure will be generated here
		integer(double_kind)::i1=0!x dimension counter
		integer(double_kind)::i2=0!y dimension counter
		integer(double_kind)::i3=0!z dimension counter
		integer(double_kind)::j=1!particle counter
		integer(double_kind)::k=0!dimension counter
		integer(double_kind)::vecindx(3)=0!here the indeces relevant to each dimensions are kept, for example the indeces of the cells in the crystalline lattice
		real(double_kind)::dran_g!gaussian random number
		real(double_kind)::dran_u!uniform random number
		if(gnuplot)write(*,*) "set title 'Fixed particles are red and free ones are black.'"!information about the representation of free and fixed particles
		title=.true.!this same title has to be displayed on each frame so setting this variable to true will make sure that happens
		do i1=0,nint(siz(1)/latcon,double_kind)-1
			vecindx(1)=i1!index of the current cell in the first dimension
			do i2=0,nint(siz(2)/latcon,double_kind)-1
				vecindx(2)=i2!index of the current cell in the second dimension
				do i3=0,nint(siz(3)/latcon,double_kind)-1
					vecindx(3)=i3!index of the current cell in the third dimension
					do k=1,3
						rvec(k,j)=latcon*dble(vecindx(k))!position in the solid lattice of the material, free particles will leave these as they move
					end do
					if(dsqrt((rvec(1,j)-5.d-1*siz(1))**2+(rvec(2,j)-5.d-1*siz(1))**2)>=radius)then
						fixed(j)=.true.!the particles outside of this cylinder are fixed
						amntfix=amntfix+1!counting the number of fixed particles
						if(gnuplot)write(*,*) rvec(1,j),rvec(2,j),rvec(3,j),'0xFF0000'!data feed to gnuplot
					else
						if(samemass)then
							if(gaussdistribution)then
								do k=1,3
									vvec(3,j)=dsqrt(temp)*dran_g()!random gaussian components for the velocities of the free particles (if constants needed use kboltz*temp/mass instead of temp)
									rvecini(k,j)=rvec(k,j)!storing the initial position of the particles
								end do
							else
								do k=1,3
									vvec(3,j)=dsqrt(temp)*(dran_u()-5.d-1)!random gaussian components for the velocities of the free particles (if constants needed use kboltz*temp/mass instead of temp)
									rvecini(k,j)=rvec(k,j)!storing the initial position of the particles
								end do
							end if
							amntfree=amntfree+1!counting the number of free particles
						else
							mvec(j)=mass*dsin(dble(j))!setting the mass of the particles
							if(gaussdistribution)then
								do k=1,3
									vvec(3,j)=dsqrt(temp/mvec(j))*dran_g()!random gaussian components for the velocities of the free particles (if constants needed use kboltz*temp/mass instead of temp/mvec(j+l))
									rvecini(k,j)=rvec(k,j)!storing the initial position of the particles
								end do
							else
								do k=1,3
									vvec(3,j)=dsqrt(temp/mvec(j))*(dran_u()-5.d-1)!random gaussian components for the velocities of the free particles (if constants needed use kboltz*temp/mass instead of temp/mvec(j+l))
									rvecini(k,j)=rvec(k,j)!storing the initial position of the particles
								end do
							end if
							amntfree=amntfree+1!counting the number of free particles
						end if
					end if
					j=j+1!next particle
					if(gnuplot)write(*,*) rvec(1,j),rvec(2,j),rvec(3,j),'0x000000'!data feed to gnuplot
				end do
			end do
		end do
		amnttot=j-1!total amount of particles
		if(amnttot>maxpartic)then
			write(*,*) 'Not enough space for so many particles, the variable "maxpartic" should be at least as big as "amnttot", ',amnttot,' in this case. Aborting.'!avoiding segmentation faults
			stop!terminating execution
		end if
		dens=amntfree/(siz(3)*pi*radius**2)!this is the number density of the free particles
		return
	end subroutine create_sc_structure



	subroutine create_two_particles()!simple system for debugging interactions, cell updating and so on
		integer(double_kind)::i=0!dimension counter
		real(double_kind)::dran_g!gaussian random number
		do i=1,3
			rvec(i,1)=siz(i)*5.d-1-0.5d0!coordinates for the first particle
			rvec(i,2)=siz(i)*5.d-1+0.5d0!coordinates for the second particle
			vvec(i,1)=dsqrt(temp/mass)*dran_g()!random gaussian components for the velocities (if constants needed use kboltz*temp/mass instead of temp/mvec(j+l))
			vvec(i,2)=dsqrt(temp/mass)*dran_g()!random gaussian components for the velocities (if constants needed use kboltz*temp/mass instead of temp/mvec(j+l))
		end do
		amnttot=2!total amount of particles
		amntfree=amnttot!just in case, to avoid problems since there are no fixed particles in this situation
		dens=amnttot/volume!this is the number density of all the particles
		return
	end subroutine create_two_particles



	subroutine halving_density()!if a less dense fluid is needed this subroutine is used
		integer(double_kind)::i=0!particle counter
		integer(double_kind)::k=0!dimension counter
		real(double_kind)::dran_u!uniform random number
		do
			i=1!starting particle counter
			if((.not.fixed(i)).and.(dran_u()<5.d-1))then
				if(samemass)then
					do k=1,3
						rvec(k,i)=rvec(k,amnttot)!destroying half of the particles to create a half-density fluid
						rvec(k,amntfree)=0.d0!erasing old unnecessary information
						vvec(k,i)=vvec(k,amnttot)!destroying half of the particles to create a half-density fluid
						vvec(k,amntfree)=0.d0!erasing old unnecessary information
						avec(k,i)=avec(k,amnttot)!destroying half of the particles to create a half-density fluid
						avec(k,amntfree)=0.d0!erasing old unnecessary information
					end do
				else
					do k=1,3
						rvec(k,i)=rvec(k,amnttot)!destroying half of the particles to create a half-density fluid
						rvec(k,amntfree)=0.d0!erasing old unnecessary information
						vvec(k,i)=vvec(k,amnttot)!destroying half of the particles to create a half-density fluid
						vvec(k,amntfree)=0.d0!erasing old unnecessary information
						fvec(k,i)=fvec(k,amnttot)!destroying half of the particles to create a half-density fluid
						fvec(k,amntfree)=0.d0!erasing old unnecessary information
					end do
				end if
				fixed(i)=fixed(amnttot)!destroying half of the particles to create a half-density fluid
				fixed(amnttot)=.false.!erasing old unnecessary information
				if(.not.samemass)then
					mvec(i)=mvec(amnttot)
					mvec(amnttot)=0.d0!erasing old unnecessary information
				end if
				amntfree=amntfree-1!one less free particle
				amnttot=amnttot-1!which also needs to be substracted from the total
			else
				rvecini(k,i)=rvec(k,i)!storing the initial position of the particles
			end if
			i=i+1!next particle
			if(i>=amnttot)exit!once all the particles have been checked for deletion, the loop ends
		end do
		return
	end subroutine halving_density



	subroutine create_cell_lists()!cell lists are used to speed up the simulation, see the article in wikipedia for more general details
		integer(double_kind)::i=0!particle counter
		integer(double_kind)::inds(3)=0!indeces of the cell in which a particle is
		integer(double_kind)::j=0!dimension counter
		do j=1,3
			numcells(j)=floor(siz(j)/intrad,double_kind)!computing how many cells each dimension has
			if(numcells(j)>maxcells)then
				write(*,*) 'Not enough space for so many cells, the variable "maxcells" should be at least as big as "numcells", ',numcells(j),' in this case. Aborting.'!avoiding segmentation faults
				stop!terminating execution
			end if
			cellsize(j)=siz(j)/dble(numcells(j))!calculating how big the cells are
		end do
		cellvolume=cellsize(1)*cellsize(2)*cellsize(3)!computing the volume of a single cell
		!$omp parallel do
		do i=1,amnttot
			do j=1,3
				inds(j)=floor(rvec(j,i)/cellsize(j),double_kind)+1!particle j is the cell which has these indeces
			end do
			occup(inds(1),inds(2),inds(3))=occup(inds(1),inds(2),inds(3))+1!updating the occupancy list with the new particle
			cells(occup(inds(1),inds(2),inds(3)),inds(1),inds(2),inds(3))=i!adding the new particle to the apropiate cell
		end do
		!$omp end parallel do
		return
	end subroutine create_cell_lists



	subroutine remove_average_momentum()!to avoid the system having a drift, the average momentum must be substracted to all the particles in the system
		integer(double_kind)::i=0!particle counter
		integer(double_kind)::j=0!dimension counter
		real(double_kind)::mmntsum(3)=0.d0!the global momentum will be stored here
		if(samemass)then
			do i=1,amnttot
				if(fixed(i))cycle!fixed particles don't contribute to the average momentum because they have no velocity
				do j=1,3
					mmntsum(j)=mmntsum(j)+vvec(j,i)!calculating the global momentum if the particles have all the same mass (if constants needed multiply vvec by mass)
				end do
			end do
			mmntsum=mmntsum/dble(amnttot)!dividing by the number of particles so that the final result is the average momentum of the system
			do i=1,amnttot
				if(fixed(i))cycle!the velocity of fixed particles shouldn't change
				do j=1,3
					vvec(j,i)=vvec(j,i)-mmntsum(j)!substracting the average momentum to the momentum of each particle and converting it back to velocity (if constants needed divide mmntsum by mass)
				end do
			end do
		else
			do i=1,amnttot
				if(fixed(i))cycle!fixed particles don't contribute to the average momentum because they have no velocity
				do j=1,3
					mmntsum(j)=mmntsum(j)+vvec(j,i)*mvec(i)!calculating the global momentum if the mass of each particle is different
				end do
			end do
			mmntsum=mmntsum/dble(amnttot)!dividing by the number of particles so that the final result is the average momentum of the system
			do i=1,amnttot
				if(fixed(i))cycle!the velocity of fixed particles shouldn't change
				do j=1,3
					vvec(j,i)=vvec(j,i)-mmntsum(j)/mvec(i)!substracting the average momentum to the momentum of each particle and converting it back to velocity
				end do
			end do
		end if
		return
	end subroutine remove_average_momentum



	subroutine half_leap(thermostat)!the verlet integrator requires to update the velocities in two steps, the subroutine distinguishes when the particles have the same mass and when they don't
		logical,intent(in)::thermostat!this variable controls if the thermostat should be applied or not, because the function is called two times per time step
		integer(double_kind)::i=0!x dimension cell counter
		integer(double_kind)::j=0!y dimension cell counter
		integer(double_kind)::k=0!z dimension cell counter
		integer(double_kind)::m=0!dimension counter
		integer(double_kind)::n=0!particle counter in the cell lists
		integer(double_kind)::tag=0!particle tag
		real(double_kind)::dran_g!gaussian random number
		real(double_kind)::dran_u!uniform random number
		real(double_kind)::andersenfactor=0.d0!to avoid computing the factor of the andersen thermostat many times it's calculated only once at the beginning and stored here
		real(double_kind)::berendsenfactor=0.d0!to avoid computing the factor of the berendsen thermostat many times it's calculated only once at the beginning and stored here
		real(double_kind)::ekineberendsen=0.d0!instantaneous kinetic energy for the berendsen temperature
		real(double_kind)::tempberendsen=0.d0!this temperature is updated every timestep to allow the berendsen thermostat to work properly
		if(firststep)tempberendsen=temp!setting the initial value for the instantaneous temperature
		if(andersen.and.firststep)andersenfactor=dsqrt(tempbath)!to avoid computing this constant value many times it's obtained here once in the first step of the simulation
		if(berendsen)berendsenfactor=dsqrt(1.d0+(tempbath/tempberendsen-1.d0)*threshold)!this is the amount that will multiply the components of the velocities in the berendsen thermostat
		if(berendsen.and.(.not.thermostat))ekineberendsen=0.d0!resetting the instantaneous kinetic energy in the first call of the time step to the subroutine
		!$omp parallel do
		do k=1,numcells(3)
			do j=1,numcells(2)
				do i=1,numcells(1)
					do n=1,occup(i,j,k)
						tag=cells(n,i,j,k)!this value is going to be called from memory many times so it will be kept here for convenience
						if(fixed(tag))cycle!fixed particles don't move
						if(samemass)then
							if((.not.berendsen).and.(.not.andersen))then
								do m=1,3
									vvec(m,tag)=vvec(m,tag)+avec(m,tag)*tstep*5.d-1!updating speeds without a thermostat if all the particles have the same mass (if constants needed divide by mass and replace avec with fvec)
								end do
							else if(berendsen.and.thermostat)then
								do m=1,3
									vvec(m,tag)=berendsenfactor*(vvec(m,tag)+avec(m,tag)*tstep*5.d-1)!updating speeds with the berendsen thermostat (if constants needed divide by mass and replace avec with fvec)
								end do
							else if(berendsen.and.(.not.thermostat))then
								do m=1,3
									ekineberendsen=ekineberendsen+vvec(m,tag)**2!computing the instantaneous kinetic energy
								end do
							else if((andersen.and.thermostat).and.(dran_u()<threshold))then
								do m=1,3
									vvec(m,tag)=andersenfactor*dran_g()!updating speeds with the andersen thermostat
								end do
							end if
						else
							if((.not.berendsen).and.(.not.andersen))then
								do m=1,3
									vvec(m,tag)=vvec(m,tag)+fvec(m,tag)*tstep*5.d-1/mvec(tag)!updating speeds without a thermostat if the particles have different masses
								end do
							else if(berendsen.and.thermostat)then
								do m=1,3
									vvec(m,tag)=berendsenfactor*(vvec(m,tag)+fvec(m,tag)*tstep*5.d-1/mvec(tag))!updating speeds with the berendsen thermostat
								end do
							else if(berendsen.and.(.not.thermostat))then
								ekineberendsen=ekineberendsen+mvec(tag)*(vvec(1,tag)**2+vvec(2,tag)**2+vvec(3,tag)**2)!computing the instantaneous kinetic energy
							else if((andersen.and.thermostat).and.(dran_u()<threshold))then
								do m=1,3
									vvec(m,tag)=andersenfactor*dran_g()/dsqrt(mvec(tag))!updating speeds with the andersen thermostat
								end do
							end if
						end if
					end do
				end do
			end do
		end do
		!$omp end parallel do
		if(berendsen.and.(.not.thermostat))tempberendsen=(ekineberendsen)/(3.d0*dble(amntfree))!instantaneous temperature (if constants needed add kboltz in the denominator and mass in the numerator when the particles have the same mass)
		return
	end subroutine half_leap



	subroutine nose_hoover()!this is the Nosé-Hoover integrator, it contains all the parts of a timestep since it calls the functions to move the particle, calculate the forces and so on
		integer(double_kind)::i=0!x dimension cell counter
		integer(double_kind)::j=0!y dimension cell counter
		integer(double_kind)::k=0!z dimension cell counter
		integer(double_kind)::m=0!dimension counter
		integer(double_kind)::n=0!particle counter in the cell lists
		integer(double_kind)::tag=0!particle tag
		real(double_kind)::nosehooverfactor=0.d0!this value will control the intensity of the friction force applied to the particles
		real(double_kind)::ekinenosehoover=0.d0!instantaneous kinetic energy for the thermostat
		real(double_kind)::ekinenosehooverhalfstep=0.d0!instantaneous kinetic energy for the thermostat at half of the step
		real(double_kind)::vvechalfstep(3,maxpartic)=0.d0!velocities for each particle (vx->(1,:), vy->(2,:), vz->(3,:)) calculated at half of the step are here
		ekinenosehoover=0.d0!resetting the instantaneous kinetic energy
		ekinenosehooverhalfstep=0.d0!resetting the instantaneous kinetic energy at half of the step
		if(samemass)then
			!$omp parallel do
			do k=1,numcells(3)
				do j=1,numcells(2)
					do i=1,numcells(1)
						do n=1,occup(i,j,k)
							tag=cells(n,i,j,k)!this value is going to be called from memory many times so it will be kept here for convenience
							if(fixed(tag))cycle!fixed particles don't move
							do m=1,3
								ekinenosehoover=ekinenosehoover+vvec(m,tag)**2!calculating the kinetic energy for the timestep
								vvechalfstep(m,tag)=vvec(m,tag)+(avec(m,tag)-nosehooverfactor*vvec(m,tag))*tstep*5.d-1!updating speeds (if constants needed divide by mass and replace avec with fvec)
								ekinenosehooverhalfstep=ekinenosehooverhalfstep+vvechalfstep(m,tag)**2!and calculating the kinetic energy at half of the step
							end do
						end do
					end do
				end do
			end do
			!$omp end parallel do
			ekinenosehoover=ekinenosehoover*mass!multiplying by the necessary constants, it's not necessary to divide by 2 because it can be factored out when calculating the nosehooverfactor variable
			ekinenosehooverhalfstep=ekinenosehooverhalfstep*mass!multiplying by the necessary constants, it's not necessary to divide by 2 because it can be factored out when calculating the nosehooverfactor variable
			call move_system()!updating the positions of the particles
			call update_particles_in_cells()!after the particles have moved, they may have entered a different cell and this needs to be checked and taken care of
			call calculate_forces()!the forces need to be calculated again so that the speed may be updated
			nosehooverfactor=nosehooverfactor+(ekinenosehoover-dble(3*amntfree+1)*tempbath)*tstep*nosehooverconstant!computing the scaling factor for the friction force (if constants needed multiply tempbath by kboltz)
			nosehooverfactor=nosehooverfactor+(ekinenosehooverhalfstep-dble(3*amntfree+1)*tempbath)*tstep*nosehooverconstant!same as before, this is at the next timestep and the previous one was at half of the step (if constants needed multiply tempbath by kboltz)
			!$omp parallel do
			do k=1,numcells(3)
				do j=1,numcells(2)
					do i=1,numcells(1)
						do n=1,occup(i,j,k)
							tag=cells(n,i,j,k)!this value is going to be called from memory many times so it will be kept here for convenience
							if(fixed(tag))cycle!fixed particles don't move
							do m=1,3
								vvec(m,tag)=(vvechalfstep(m,tag)+avec(m,tag)*tstep*5.d-1)/(1.d0+nosehooverfactor*tstep*5.d-1)!updating speeds (if constants needed divide by mass and replace avec with fvec)
							end do
						end do
					end do
				end do
			end do
			!$omp end parallel do
		else
			!$omp parallel do
			do k=1,numcells(3)
				do j=1,numcells(2)
					do i=1,numcells(1)
						do n=1,occup(i,j,k)
							tag=cells(n,i,j,k)!this value is going to be called from memory many times so it will be kept here for convenience
							if(fixed(tag))cycle!fixed particles don't move
							ekinenosehoover=ekinenosehoover+mvec(tag)*(vvec(1,tag)**2+vvec(2,tag)**2+vvec(3,tag)**2)!calculating the kinetic energy for the timestep
							do m=1,3
								vvechalfstep(m,tag)=vvec(m,tag)+(fvec(m,tag)/mvec(tag)-nosehooverfactor*vvec(m,tag))*tstep*5.d-1!updating speeds
							end do
							ekinenosehooverhalfstep=ekinenosehooverhalfstep+mvec(tag)*(vvechalfstep(1,tag)**2+vvechalfstep(2,tag)**2+vvechalfstep(3,tag)**2)!and calculating the kinetic energy at half of the step
						end do
					end do
				end do
			end do
			!$omp end parallel do
			call move_system()!updating the positions of the particles
			call update_particles_in_cells()!after the particles have moved, they may have entered a different cell and this needs to be checked and taken care of
			call calculate_forces()!the forces need to be calculated again so that the speed may be updated
			nosehooverfactor=nosehooverfactor+(ekinenosehoover-dble(3*amntfree+1)*tempbath)*tstep*nosehooverconstant!computing the scaling factor for the friction force (if constants needed multiply tempbath by kboltz)
			nosehooverfactor=nosehooverfactor+(ekinenosehooverhalfstep-dble(3*amntfree+1)*tempbath)*tstep*nosehooverconstant!same as before, this is at the next timestep and the previous one was at half of the step (if constants needed multiply tempbath by kboltz)
			!$omp parallel do
			do k=1,numcells(3)
				do j=1,numcells(2)
					do i=1,numcells(1)
						do n=1,occup(i,j,k)
							tag=cells(n,i,j,k)!this value is going to be called from memory many times so it will be kept here for convenience
							if(fixed(tag))cycle!fixed particles don't move
							do m=1,3
								vvec(m,tag)=(vvechalfstep(m,tag)+fvec(m,tag)*tstep*5.d-1/mvec(tag))/(1.d0+nosehooverfactor*tstep*5.d-1)!updating speeds
							end do
						end do
					end do
				end do
			end do
			!$omp end parallel do
		end if
		return
	end subroutine nose_hoover



	subroutine move_system()!this is where the actual movement of the particles takes place
		integer(double_kind)::i=0!x dimension cell counter
		integer(double_kind)::j=0!y dimension cell counter
		integer(double_kind)::k=0!z dimension cell counter
		integer(double_kind)::m=0!dimension counter
		integer(double_kind)::n=0!particle counter in the cell lists
		integer(double_kind)::tag=0!particle tag
		if(gnuplot.and.title)write(*,*) "set title 'Fixed particles are red and free ones are black.'"
		if(gnuplot)write(*,*) "splot '-' w p pt 7 lc rgb variable"!new gnuplot frame
		!$omp parallel do
		do k=1,numcells(3)
			do j=1,numcells(2)
				do i=1,numcells(1)
					do n=1,occup(i,j,k)
						tag=cells(n,i,j,k)!this value is going to be called from memory many times so it will be kept here for convenience
						if(fixed(tag))then
							if(gnuplot)write(*,*) rvec(1,tag),rvec(2,tag),rvec(3,tag),'0xFF0000'!data feed to gnuplot
							cycle!fixed particles don't move
						end if
						do m=1,3
							rvec(m,tag)=rvec(m,tag)+vvec(m,tag)*tstep!moving the particle
							if(.not.periodicspeed)then
								if(rvec(m,tag)<0.d0)then
									rvec(m,tag)=rvec(m,tag)+siz(m)!if the particle has reached the lower end of one dimension of the simulation box it should be at the other side of it
									rvecnonper(m,tag)=rvecnonper(m,tag)-siz(m)!tracking the movement along the periodic boundary conditions
								else if(rvec(m,tag)>=siz(m))then
									rvec(m,tag)=rvec(m,tag)-siz(m)!if the particle has reached the higher end of one dimension of the simulation box it should be at the other side of it
									rvecnonper(m,tag)=rvecnonper(m,tag)+siz(m)!tracking the movement along the periodic boundary conditions
								end if
							else
								if(rvec(m,tag)<0.d0)then
									rvec(m,tag)=rvec(m,tag)-dble(floor(rvec(m,tag)/siz(m),double_kind))*siz(m)!same as with .not.periodicspeed, but it will work no matter how far away the particle moves
									rvecnonper(m,tag)=rvecnonper(m,tag)+dble(floor(rvec(m,tag)/siz(m),double_kind))*siz(m)!same as with .not.periodicspeed, but it will work no matter how far away the particle moves
								else if(rvec(m,tag)>=siz(m))then
									rvec(m,tag)=rvec(m,tag)-dble(floor(rvec(m,tag)/siz(m),double_kind))*siz(m)!same as with .not.periodicspeed, but it will work no matter how far away the particle moves
									rvecnonper(m,tag)=rvecnonper(m,tag)+dble(floor(rvec(m,tag)/siz(m),double_kind))*siz(m)!same as with .not.periodicspeed, but it will work no matter how far away the particle moves
								end if
							end if
						end do
						if(gnuplot)write(*,*) rvec(1,tag),rvec(2,tag),rvec(3,tag),'0x000000'!data feed to gnuplot
					end do
				end do
			end do
		end do
		!$omp end parallel do
		if(gnuplot)write(*,*) 'e'!end of gnuplot frame
		return
	end subroutine move_system



	subroutine update_particles_in_cells()!the particles may be in a different cell after moving, so the lists need to be updated
		integer(double_kind)::i=0!x dimension cell counter
		integer(double_kind)::inew=0!x dimension, new-cell counter
		integer(double_kind)::j=0!y dimension cell counter
		integer(double_kind)::jnew=0!y dimension, new-cell counter
		integer(double_kind)::k=0!z dimension cell counter
		integer(double_kind)::knew=0!z dimension, new-cell counter
		integer(double_kind)::n=0!particle counter in the cell lists
		integer(double_kind)::tag=0!particle tag
		!$omp parallel do
		do k=1,numcells(3)
			do j=1,numcells(2)
				do i=1,numcells(1)
					n=1!starting particle counter at the beginning of the current list
					tag=cells(n,i,j,k)!this value is going to be called from memory many times so it will be kept here for convenience
					if(tag==0)cycle!this cell is empty, moving onto the next one
					do
						tag=cells(n,i,j,k)!this value is going to be called from memory many times so it will be kept here for convenience
						if(fixed(tag))then
							n=n+1!next particle
							if(n>occup(i,j,k))exit!no more particles in this cell
							cycle!fixed particles don't need updating
						end if
						inew=floor(rvec(1,tag)/cellsize(1),double_kind)+1!computing the x index of the cell in which the particle is located after moving
						jnew=floor(rvec(2,tag)/cellsize(2),double_kind)+1!computing the y index of the cell in which the particle is located after moving
						knew=floor(rvec(3,tag)/cellsize(3),double_kind)+1!computing the z index of the cell in which the particle is located after moving
						if((i/=inew).or.(j/=jnew).or.(k/=knew))then!after the particle has moved, the cell in which it is contained may have changed so it needs to be updated
							occup(inew,jnew,knew)=occup(inew,jnew,knew)+1!updating the new occupancy list with the particle
							cells(occup(inew,jnew,knew),inew,jnew,knew)=tag!adding the particle to the corresponding new cell
							cells(n,i,j,k)=cells(occup(i,j,k),i,j,k)!taking away the updated particle from the old cell, its place is occupied by the last one in that cell
							cells(occup(i,j,k),i,j,k)=0!erasing old information
							occup(i,j,k)=occup(i,j,k)-1!updating the old occupancy list
						end if
						n=n+1!next particle
						if(n>occup(i,j,k))exit!no more particles in this cell
					end do
				end do
			end do
		end do
		!$omp end parallel do
		return
	end subroutine update_particles_in_cells



	subroutine calculate_forces()!obtaining the values for the forces of the particles
		integer(double_kind)::i1=0!x dimension cell counter for the first particle
		integer(double_kind)::i2=0!x dimension cell counter for the second particle
		integer(double_kind)::i3=0!x dimension auxiliary cell counter
		integer(double_kind)::j1=0!y dimension cell counter for the first particle
		integer(double_kind)::j2=0!y dimension cell counter for the second particle
		integer(double_kind)::j3=0!y dimension auxiliary cell counter
		integer(double_kind)::k1=0!z dimension cell counter for the first particle
		integer(double_kind)::k2=0!z dimension cell counter for the second particle
		integer(double_kind)::k3=0!z dimension auxiliary cell counter
		integer(double_kind)::n1=0!particle counter in the cells, for the first one
		integer(double_kind)::n2=0!particle counter in the cells, for the second one
		integer(double_kind)::tag1=0!tag for the first particle
		integer(double_kind)::tag2=0!tag for the second particle
		!$omp parallel do
		if(samemass.and.(.not.extforce))then
			do n1=1,amnttot
				do n2=1,3
					avec(n2,n1)=0.d0!setting accelerations to zero (if constants needed replace avec with fvec)
				end do
			end do
		else if(samemass.and.extforce)then
			do n1=1,amnttot
				do n2=1,3
					avec(n2,n1)=outforce(n2)!setting accelerations to the external acceleration, it has already been divided by the mass (if constants needed replace avec with fvec)
				end do
			end do
		else if((.not.samemass).and.extforce)then
			do n1=1,amnttot
				do n2=1,3
					fvec(n2,n1)=outforce(n2)!setting forces to to the external force
				end do
			end do
		else if((.not.samemass).and.(.not.extforce))then
			do n1=1,amnttot
				do n2=1,3
					fvec(n2,n1)=0.d0!setting forces to zero
				end do
			end do
		end if
		do n1=1,amnttot
			do n2=1,amnttot
				interaction(n2,n1)=.false.!the interactions must be reset so they can be taken into account again in the new iteration
			end do
		end do
		!$omp end parallel do
		!$omp parallel do
		do k1=1,numcells(3)
			do j1=1,numcells(2)
				do i1=1,numcells(1)
					do n1=1,occup(i1,j1,k1)
					tag1=cells(n1,i1,j1,k1)!this value is going to be called from memory many times so it will be kept here for convenience
						do k3=k1-1,k1+1
							call periodic_interaction(3_8,k3,k2)!checking for interactions beyond the periodic boundary conditions in the z dimension (the _8 is to indicate an 8 byte integer constant)
							do j3=j1-1,j1+1
								call periodic_interaction(2_8,j3,j2)!checking for interactions beyond the periodic boundary conditions in the y dimension (the _8 is to indicate an 8 byte integer constant)
								do i3=i1-1,i1+1
									call periodic_interaction(1_8,i3,i2)!checking for interactions beyond the periodic boundary conditions in the x dimension (the _8 is to indicate an 8 byte integer constant)
									do n2=1,occup(i2,j2,k2)
										tag2=cells(n2,i2,j2,k2)!this value is going to be called from memory many times so it will be kept here for convenience
										!$omp flush (interaction)
										if((tag1==tag2).or.(interaction(tag1,tag2)).or.(interaction(tag2,tag1)))cycle!a particle shouldn't interact with itself or more than needed with any other one
										call interacting(i1,j1,k1,tag1,tag2)!if none of the previous conflicts happen, the interaction is computed
									end do
								end do
							end do
						end do
					end do
				end do
			end do
		end do
		!$omp end parallel do
		presscount=presscount+1!this counter keeps track of how many timesteps are used to sample the pressure
		return
	end subroutine calculate_forces



	subroutine periodic_interaction(dimen,inindex,outindex)!in the case that some cell is beyond the limits, this will apply the periodic boundary conditions to the interaction between the particles
		integer(double_kind),intent(in)::dimen!dimension in which the subroutine is carried out, 1 is x-axis, 2 is y-axis and 3 is z-axis
		integer(double_kind),intent(in)::inindex!this number runs through the range of possible values for neighbouring cells
		integer(double_kind),intent(out)::outindex!and this is the resulting index of that cell
		if(inindex==0)then
			displ(dimen)=-siz(dimen)!if the neighbouring cell is at the end of the simulation box in the x direction the displacement vector has this value
			outindex=numcells(dimen)!and this is the actual index of that cell
		!else if(inindex==numcells(dimen)+1)then
		else if(inindex>numcells(dimen))then
			displ(dimen)=siz(dimen)!if the neighbouring cell is at the beginning of the simulation box in the x direction the displacement vector has this value
			outindex=1!and this is the actual index of that cell
		else
			displ(dimen)=0.d0!if the case is not one the former, the displacement vector is null
			outindex=inindex!and the actual index of the cell doesn't change
		end if
		return
	end subroutine periodic_interaction



	subroutine interacting(i,j,k,tag1,tag2)!this is where the actual interaction happens
		integer(double_kind),intent(in)::i!x dimension cell counter
		integer(double_kind),intent(in)::j!y dimension cell counter
		integer(double_kind),intent(in)::k!z dimension cell counter
		integer(double_kind)::l=0!dimension counter
		integer(double_kind),intent(in)::tag1!tag for the first particle
		integer(double_kind),intent(in)::tag2!tag for the second particle
		real(double_kind)::acceleration=0.d0!modulus of the acceleration vector
		real(double_kind)::dist=0.d0!distance between particles
		real(double_kind)::dist_2=0.d0!distance between particles squared
		real(double_kind)::dist_6=0.d0!distance between particles raised to the minus sixth power
		real(double_kind)::distvec(3)=0.d0!vectorial distance between particles
		!$omp critical
		interaction(tag1,tag2)=.true.!accounting for this interaction
		!$omp end critical
		do l=1,3
			distvec(l)=rvec(l,tag1)-rvec(l,tag2)-displ(l)!these are the components of the vectorial distance between the particles
		end do
		dist=dsqrt(distvec(1)**2+distvec(2)**2+distvec(3)**2)!and this is the distance between the two particles (if constants needed, add sigm/ before the square root)
		if(dist>intrad)return!if the distance between the particles is bigger than the interacion radius there is no need to continue because this interaction is very small and can be disregarded
		dist_2=dist**2!keeping this value in the memory to avoid calculating it many times
		dist_6=dist**(-6)!keeping this value in the memory to avoid calculating it many times
		if(fixed(tag2).or.fixed(tag1))then!the interaction can be slightly different between fixed and free particles
			acceleration=2.4d1*dist_6*(2.d0*dist_6-1.d0)/dist_2!this is the module of the acceleration between a free and a fixed particle
		else
			acceleration=2.4d1*dist_6*(2.d0*dist_6-1.d0)/dist_2!and this is the module of the acceleration between two free particles
		end if
		if(samemass)then
			if((.not.fixed(tag1)).and.(.not.fixed(tag2)))then
				do l=1,3
					avec(l,tag1)=avec(l,tag1)+distvec(l)*acceleration!and these are the components of the acceleration vector for the first particle (if constants needed replace avec with fvec)
					avec(l,tag2)=avec(l,tag2)-distvec(l)*acceleration!the same for the second particle (if constants needed replace avec with fvec)
				end do
			else if(fixed(tag2).and.(.not.fixed(tag1)))then
				do l=1,3
					avec(l,tag1)=avec(l,tag1)+distvec(l)*acceleration!and these are the components of the acceleration vector for the first particle (if constants needed replace avec with fvec)
				end do
			else if(fixed(tag1).and.(.not.fixed(tag2)))then
				do l=1,3
					avec(l,tag2)=avec(l,tag2)-distvec(l)*acceleration!and these are the components of the acceleration vector for the first particle (if constants needed replace avec with fvec)
				end do
			end if
		else
			if((.not.fixed(tag1)).and.(.not.fixed(tag2)))then
				do l=1,3
					fvec(l,tag1)=fvec(l,tag1)+distvec(l)*acceleration*mvec(tag1)!and these are the components of the force vector for the first particle
					fvec(l,tag2)=fvec(l,tag2)-distvec(l)*acceleration*mvec(tag2)!the same for the second particle
				end do
			else if(fixed(tag2).and.(.not.fixed(tag1)))then
				do l=1,3
					fvec(l,tag1)=fvec(l,tag1)+distvec(l)*acceleration*mvec(tag1)!and these are the components of the force vector for the first particle
				end do
			else if(fixed(tag1).and.(.not.fixed(tag2)))then
				do l=1,3
					fvec(l,tag2)=fvec(l,tag2)-distvec(l)*acceleration*mvec(tag2)!and these are the components of the force vector for the first particle
				end do
			end if
		end if
		if((timetomeasure.or.(iterations==0)).and.(cylinders.and.normaloutput))then
			!#####HERE MAYBE THERE SHOULD BE A CONTRIBUTION TO THE POTENTIAL ENERGY FOR PARTICLES LOCATED EXCLUSIVELY IN EACH CYLINDER
			epotcy(1)=epot+4.d0*dist_6*(2.d0*dist_6-1.d0)-potentialshift
		else if((timetomeasure.or.(iterations==0)).and.normaloutput)then
			epot=epot+4.d0*dist_6*(2.d0*dist_6-1.d0)-potentialshift!this is the potential energy
		end if
		if(samemass)then
			if(normaloutput)press=press+dist_2*acceleration*mass!to calculate the pressure, the sum over all interactions must be carried out
			if(temp3d.and.press3d)presscells(i,j,k)=presscells(i,j,k)+dist_2*acceleration*mass!same as above for the current cell
		else
			if(normaloutput)press=press+dist_2*acceleration*mvec(tag1)!to calculate the pressure, the sum over all interactions must be carried out
			if(temp3d.and.press3d)presscells(i,j,k)=presscells(i,j,k)+dist_2*acceleration*mvec(tag1)!same as above for the current cell
		end if
		return
	end subroutine interacting



	subroutine normal_measuring()!here the different measurements are performed
		integer(double_kind)::freecounter=0!the amount of free particles per cell is kept here
		integer(double_kind)::i=0!x dimension cell counter
		integer(double_kind)::j=0!y dimension cell counter
		integer(double_kind)::k=0!z dimension cell counter
		integer(double_kind)::n=0!particle counter
		integer(double_kind)::tag=0!particle tag
		if(normaloutput)then
			ekine=0.d0!setting the kinetic energy to zero to calculate the new values
			if(samemass)then
				do k=1,numcells(3)
					do j=1,numcells(2)
						do i=1,numcells(1)
							ekinecells(i,j,k)=0.d0!resetting the values of the kinetic energy of each cell
							freecounter=0!setting the counter to zero for the next cell
							do n=1,occup(i,j,k)
								tag=cells(n,i,j,k)!this value is going to be called from memory many times so it will be kept here for convenience
								if(fixed(tag))cycle!fixed particles don't contribute to the kinetic energy
								ekinecells(i,j,k)=ekinecells(i,j,k)+vvec(1,tag)**2+vvec(2,tag)**2+vvec(3,tag)**2!measuring the kinetic energy for the current cell
								freecounter=freecounter+1!another free particle
							end do
							ekine=ekine+ekinecells(i,j,k)!measuring the global kinetic energy
							if(temp3d)then
								ekinecells(i,j,k)=ekinecells(i,j,k)*5.d-1!don't forget the 1/2 (if constants are needed multiply by mass)
								tempcells(i,j,k)=(2.d0*ekinecells(i,j,k))/(3.d0*dble(freecounter))!temperature in the current cell (if constants needed add kboltz in the denominator)
								if(press3d)then
									denscells(i,j,k)=dble(freecounter)/cellvolume!and same with the number density
									presscells(i,j,k)=denscells(i,j,k)*tempcells(i,j,k)+presscells(i,j,k)/(3.d0*cellvolume*dble(presscount))!final result for the pressure in the current cell (if constants needed multiply temp by kboltz)
								end if
							end if
						end do
					end do
				end do
				ekine=ekine*mass*5.d-1!don't forget the 1/2
			else
				do k=1,numcells(3)
					do j=1,numcells(2)
						do i=1,numcells(1)
							ekinecells(i,j,k)=0.d0!resetting the values of the kinetic energy of each cell
							freecounter=0!setting the counter to zero for the next cell
							do n=1,occup(i,j,k)
								tag=cells(n,i,j,k)!this value is going to be called from memory many times so it will be kept here for convenience
								if(fixed(tag))cycle!fixed particles don't contribute to the kinetic energy
								ekinecells(i,j,k)=ekinecells(i,j,k)+mvec(tag)*(vvec(1,tag)**2+vvec(2,tag)**2+vvec(3,tag)**2)!measuring the kinetic energy for the current cell if the particles have different masses
								freecounter=freecounter+1!another free particle
							end do
							ekine=ekine+ekinecells(i,j,k)!measuring the global kinetic energy
							if(temp3d)then
								ekinecells(i,j,k)=ekinecells(i,j,k)*5.d-1!don't forget and the 1/2
								tempcells(i,j,k)=(2.d0*ekinecells(i,j,k))/(3.d0*dble(freecounter))!temperature in the current cell (if constants needed add kboltz in the denominator)
								if(press3d)then
									denscells(i,j,k)=dble(freecounter)/cellvolume!and same with the number density
									presscells(i,j,k)=denscells(i,j,k)*tempcells(i,j,k)+presscells(i,j,k)/(3.d0*cellvolume*dble(presscount))!final result for the pressure in the current cell (if constants needed multiply temp by kboltz)
								end if
							end if
						end do
					end do
				end do
				ekine=ekine*5.d-1!don't forget the 1/2
			end if
			sumekine=sumekine+ekine!summing the kinetic energy
			sumekinesqur=sumekinesqur+ekine**2!and the squared kinetic energy
			cvcount=cvcount+1!another measurement for the heat capacity
			temp=(2.d0*ekine)/(3.d0*dble(amntfree))!and with the kinetic energy the temperature is measured too (if constants needed add kboltz in the denominator)
			if((iterations>=itermeas).or.(iterations>0))then
				tempmax=dmax1(temp,tempmax)!recursively finding the maximum temperature
				tempmean=tempmean+temp!keeping track of the mean temperature
				tempmin=dmin1(temp,tempmin)!recursively finding the minimum temperature
			end if
			press=dens*temp+press/(3.d0*volume*dble(presscount))!final result for the pressure (if constants needed multiply temp by kboltz)
			if(iterations==0)epot=epot*5.d-1!in the initial iteration the potential energy is counted twice since the function needs to be called when the system is set up in startup()
			write(30,*) time_to_SI(time),energy_to_SI(epot+ekine),energy_to_SI(epot),energy_to_SI(ekine),temperature_to_SI(temp),pressure_to_SI(press),area_to_SI(squaredispl)/(6.d0*amnttot)!writing data
			epot=0.d0!setting the potential energy to zero so it can be calculated again
			press=0.d0!discarding the old value of the pressure so the new one can be computed in the next steps
			presscount=0!setting the counter to zero for the next measurement of the presure
			call flush(30)!forcing the program to write stuff
		end if
		if(measurevhist.and.(.not.(iterations==0)))call velocity_histogram()!obtaining another snapshot of the velocity distribution
		if(measureraddist)call radial_distribution()!getting a new sample of the radial distribution
		if(measureprofile)call velocity_profile()!
		if(measurediffusion)call diffusion()!measuring the diffusion of the particles
		if(temp3d)then
			write(70,*) 'set title "The time is: ',time!the title of the plot is used to keep track of the time
			write(70,*) "splot '-' u 1:2:3:4 lc palette title 'Temperature'"!frame start
			do k=1,numcells(3)
				do j=1,numcells(2)
					do i=1,numcells(1)
						write(70,*) (dble(i)-5.d-1)*cellsize(1),(dble(j)-5.d-1)*cellsize(2),(dble(k)-5.d-1)*cellsize(3),tempcells(i,j,k)!writing data for spatially distributed magnitudes
					end do
				end do
			end do
			write(70,*) 'e'!end of frame
			write(70,*) 'set term'!pause between frames, press enter twice for the next frame
			call flush(70)!forcing the program to write stuff
		end if
		if(press3d.and.temp3d)then
			write(80,*) 'set title "The time is: ',time!the title of the plot is used to keep track of the time
			write(80,*) "splot '-' u 1:2:3:4 lc palette title 'Pressure'"!frame start
			do k=1,numcells(3)
				do j=1,numcells(2)
					do i=1,numcells(1)
						write(80,*) (dble(i)-5.d-1)*cellsize(1),(dble(j)-5.d-1)*cellsize(2),(dble(k)-5.d-1)*cellsize(3),presscells(i,j,k)!writing data for spatially distributed magnitudes
						presscells(i,j,k)=0.d0!setting the pressure in the current cell to zero so it can be calculated again
					end do
				end do
			end do
			write(80,*) 'e'!end of frame
			write(80,*) 'set term'!pause between frames, press enter twice for the next frame
			call flush(80)!forcing the program to write stuff
		end if
		meascounter=meascounter+1!another measurement
		return
	end subroutine normal_measuring



	subroutine diffusion()!here measurements of the square displacement of the particles are taken
		integer(double_kind)::i=0!particle counter
		integer(double_kind)::j=0!dimension counter
		real(double_kind)::tmp=0.d0!temporary variable to calculate the distance travelled by the particles
		squaredispl=0.d0!same for the mean square displacement
		do i=1,amnttot
			if(fixed(i))cycle!fixed particles are irrelevant for diffusion processes
			do j=1,3
				tmp=tmp+(rvec(j,i)-rvecini(j,i)+rvecnonper(j,i))**2!calculating the distance travelled by the particles
			end do
			squaredispl=squaredispl+dsqrt(tmp)!computing the mean square displacement needed for the diffusion
			tmp=0.d0!resetting the temporary variable
		end do
		return
	end subroutine diffusion



	subroutine velocity_profile()!the velocity profile in a cylindrical section can be measured with this subroutine
		integer(double_kind)::i=0!particle counter
		integer(double_kind)::tag=0!tag for the interval in which a particle is located inside the cylinder
		real(double_kind)::distance=0.d0!distance to the centre of the cylinder
		do i=1,amnttot
			if(fixed(i))cycle!fixed particles are irrelevant here
			distance=dsqrt((rvec(1,i)-siz(1)*5.d-1)**2+(rvec(2,i)-siz(2)*5.d-1)**2)!calculating the distance to the centre
			tag=nint(dble(meascylin)*distance/radius,double_kind)!this is where the particle is
			if(tag<1)tag=1!just in case the nint function fails or the particles are in a weird place
			if(tag>meascylin)tag=meascylin!just in case the nint function fails or the particles are in a weird place
			profile(tag)=profile(tag)+vvec(3,i)!adding the velocities of the particles to the current
		end do
		return
	end subroutine velocity_profile



	subroutine measuring_with_cylinders()!here the different measurements are performed in each cylinder, whether all the particles have the same mass or not
		integer(double_kind)::i=0!particle counter
		integer(double_kind)::j=0!cylinder counter
		do j=1,meascylin
			ekinecy(j)=0.d0!setting the kinetic energy of the cylinders to zero to calculate the new values
		end do
		if(samemass)then
			do i=1,amnttot
				if(.not.fixed(i))then
					j=floor(dble(meascylin)*dsqrt(rvec(1,i)**2+rvec(2,i)**2)/radius,double_kind)!finding in which cylinder the particle is located
					if(j>meascylin)j=meascylin!just in case the particle is slightly beyond the radius of the cylinder
					ekinecy(j)=ekinecy(j)+vvec(1,i)**2+vvec(2,i)**2+vvec(3,i)**2!measuring the kinetic energy per cylinder (if constants are needed multiply by mass)
				end if
			end do
			do j=1,meascylin
				ekinecy(j)=ekinecy(j)*mass*5.d-1!don't forget the 1/2
				sumekinecy(j)=sumekinecy(j)+ekinecy(j)!summing the kinetic energy per cylinder
				sumekinecysqur(j)=sumekinecysqur(j)+ekinecy(j)**2!and the squared kinetic energy per cylinder
				cvcountcy(j)=cvcountcy(j)+1!another measurement for the heat capacity per cylinder
				tempcy(j)=(2.d0*ekinecy(j))/(3.d0*amntfree)!and with the kinetic energy the temperature is measured too per cylinder (if constants needed add kboltz in the denominator)
			end do
		else
			do i=1,amnttot
				if(.not.fixed(i))then
					j=floor(dble(meascylin)*dsqrt(rvec(1,i)**2+rvec(2,i)**2)/radius,double_kind)!finding in which cylinder the particle is located
					if(j>meascylin)j=meascylin!just in case the particle is slightly beyond the radius of the cylinder
					ekinecy(j)=ekinecy(j)+mvec(i)*(vvec(1,i)**2+vvec(2,i)**2+vvec(3,i)**2)!measuring the kinetic energy per cylinder
				end if
			end do
			do j=1,meascylin
				ekinecy(j)=ekinecy(j)*5.d-1!don't forget the 1/2
				sumekinecy(j)=sumekinecy(j)+ekinecy(j)!summing the kinetic energy per cylinder
				sumekinecysqur(j)=sumekinecysqur(j)+ekinecy(j)**2!and the squared kinetic energy per cylinder
				cvcountcy(j)=cvcountcy(j)+1!another measurement for the heat capacity per cylinder
				tempcy(j)=(2.d0*ekinecy(j))/(3.d0*amntfree)!and with the kinetic energy the temperature is measured too per cylinder (if constants needed add kboltz in the denominator)
			end do
		end if
		if(measurevhist)call velocity_histogram()!getting a sample of the velocities
		if(measurediffusion)call diffusion_with_cylinders()!obtaining another diffusion measurement
		return
	end subroutine measuring_with_cylinders



	subroutine diffusion_with_cylinders()!here measurements of the square displacement of the particles are taken for each cylinder of the system
		real(double_kind)::dran_u!uniform random number
		do
			difftag(1)=nint(dble(amntfree)*dran_u(),double_kind)!choosing randomly the particle that will be used to measure diffusion properties
			if(difftag(1)>0)exit!just in case the tag is 0 or negative
		end do
		return
	end subroutine diffusion_with_cylinders



	subroutine radial_distribution()!computing the radial distribution
		integer(double_kind)::bintag=0!tag of the current bin
		integer(double_kind)::i=0!particle counter
		integer(double_kind)::tag=0!tag for the first particle
		real(double_kind)::dist=1.d10!distance between particles
		real(double_kind)::tmp=0.d0!temporary storage for the recursive search of the minimum distance to the centre of the simulation box
		if(firststep)then
			do i=1,amnttot
				tmp=dsqrt((rvec(1,i)-siz(1)*5.d-1)**2+(rvec(2,i)-siz(2)*5.d-1)**2+(rvec(3,i)-siz(3)*5.d-1)**2)!distance between the particle i and the centre of the simulation box
				if(tmp<dist)then
					dist=tmp!here the minimum is found recursively
					tag=i!and the tag is kept for later use
				end if
			end do
		end if
		do i=1,amnttot
			dist=dsqrt((rvec(1,tag)-rvec(1,i))**2+(rvec(2,tag)-rvec(2,i))**2+(rvec(3,tag)-rvec(3,i))**2)!distance between particles
			if(dist<=l2)then
				bintag=nint(dist*dble(binrad)/l2,double_kind)!this is the bin in which the particle is
				if(bintag>binrad)bintag=binrad!just in case the nint function fails
				if(bintag<1)bintag=1!just in case the particles are too close
				distbins(bintag)=distbins(bintag)+1!counting another particle in the appropiate bin
			end if
		end do
		radcount=radcount+1!one more sample for the measurement
		return
	end subroutine radial_distribution



	subroutine velocity_histogram()!computing the velocity histogram
		character(26)::timetitle=''!this is used to convert double precision reals to characters
		integer(double_kind)::histogram(3,maxbinhist)=0!in this vector the values of the velocity histogram will be kept for each component (x->(1,:), y->(2,:), z->(3,:))
		integer(double_kind)::i=0!particle counter
		integer(double_kind)::j=0!dimension counter
		integer(double_kind)::temph(3)=0!temporary placeholder (x->(1), y->(2), z->(3))
		real(double_kind)::maxv(3)=0.d0!maximum value found for the components of the velocity vector (x->(1), y->(2), z->(3))
		real(double_kind)::minv(3)=0.d0!minimum value found for the components of the velocity vector (x->(1), y->(2), z->(3))
		real(double_kind)::sizeofbins(3)=0.d0!how big the bins are for each dimension (x->(1), y->(2), z->(3))
		do j=1,3
			maxv(j)=-1.d100!resetting maximum
			minv(j)=1.d100!resetting minimum
		end do
		do i=1,binhist
			do j=1,3
				histogram(j,i)=0!resetting the histogram
			end do
		end do
		do i=1,amnttot
			do j=1,3
				maxv(j)=dmax1(maxv(j),vvec(j,i))!recursively finding the maximum value of the modules
				minv(j)=dmin1(minv(j),vvec(j,i))!recursively finding the minimum value of the modules
			end do
		end do
		do j=1,3
			sizeofbins(j)=(maxv(j)-minv(j))/dble(binhist)!finding the size of the bins
		end do
		do i=1,amnttot
			do j=1,3
				temph(j)=ceiling((vvec(j,i)-minv(j))/sizeofbins(j),double_kind)!finding to which bin this value belongs
				if(temph(j)>binhist)temph(j)=binhist!just in case the ceiling function fails
				if(temph(j)<1)temph(j)=1!just in case the ceiling function fails
				histogram(j,temph(j))=histogram(j,temph(j))+1!updating the histogram
			end do
		end do
		write(timetitle,*) time!converting the time to a character variable
		write(60,*) 'set title "The time is:'//timetitle//'"'!the title of the plot is used to keep track of time
		write(60,*) "plot '-' u 1:2:3 w boxes lc rgb variable title 'Histogram of velocity components: Red -> x, Green -> y, Blue -> z'"!frame start
		write(60,*) dmin1(minv(1),minv(2),minv(3))-dmin1(sizeofbins(1),sizeofbins(2),sizeofbins(3))*5.d-1,0,'0xFFFFFF'!writing histogram data to file
		do i=1,binhist
			write(60,*) (dble(i)-5.d-1)*sizeofbins(1)+minv(1),histogram(1,i),'0xFF0000'!writing histogram data to file
			write(60,*) (dble(i)-5.d-1)*sizeofbins(2)+minv(2),histogram(2,i),'0x00FF00'!writing histogram data to file
			write(60,*) (dble(i)-5.d-1)*sizeofbins(3)+minv(3),histogram(3,i),'0x0000FF'!writing histogram data to file
		end do
		write(60,*) dmax1(maxv(1),maxv(2),maxv(3))+dmax1(sizeofbins(1),sizeofbins(2),sizeofbins(3))*5.d-1,0,'0xFFFFFF'!writing histogram data to file
		write(60,*) 'e'!end of frame
		write(60,*) 'set term'!pause between frames, press enter twice for the next frame
		call flush(60)!forcing the program to write stuff
		return
	end subroutine velocity_histogram



	subroutine output()!after the simulation is done, some more stuff needs to be computed and the files need to be closed
		integer(double_kind)::i=0!cylinder and bin counter
		if(cylinders.and.normaloutput)then
			do i=1,meascylin
				cvcy(i)=(3.d0)/((4.d0*amntfree*(sumekinecy(i)**2/dble(cvcountcy(i))-sumekinecysqur(i))/(dble(cvcountcy(i))*3.d0*tempcy(i)**2))-2.d0)!computing the final value of the heat capacity (if constants needed multiply by kboltz both 3.d0's)
				write(30,*) '# The heat capacity in cylinder',i,'is:',cvcy(i)!writing data
			end do
			write(30,*) '# The mean temperature is:',temperature_to_SI(tempmean)/dble(meascounter),', the maximum temperature is:',temperature_to_SI(tempmax),', and the minimum temperature is:',temperature_to_SI(tempmin)!writing data
			write(30,*) '# The timestep is: ',tstep,' and the x-y-z dimensions of the system are: ',siz(1),siz(2),siz(3)!keeping the time and the size of the system in the file in case they are needed later
			call flush(30)!forcing the program to write stuff
		else
			if(normaloutput)then
				cv=(3.d0)/((4.d0*amnttot*(sumekine**2/dble(cvcount)-sumekinesqur)/(dble(cvcount)*3.d0*(tempmean/dble(meascounter))**2))-2.d0)!computing the final value of the heat capacity (if constants needed multiply by kboltz both 3.d0's)
				write(30,*) '# The heat capacity is:',cv,', the mean temperature is:',temperature_to_SI(tempmean)/dble(meascounter),', the maximum temperature is:',temperature_to_SI(tempmax),', and the minimum temperature is:',temperature_to_SI(tempmin)!writing data
				if(struc=='por')write(30,*) '# The porosity is:',porosity,', the non-overlapping porosity is:',porevolume,', and the ratio between them is:',porosity/porevolume!writing porosity data
				write(30,*) '# The timestep is: ',tstep,' and the x-y-z dimensions of the system are: ',siz(1),siz(2),siz(3)!keeping the time and the size of the system in the file in case they are needed later
				call flush(30)!forcing the program to write stuff
			end if
			if(measureraddist)then
				do i=1,binrad
					if(i==1)distbins(i)=0!there seems to be a hard to find bug that creates huge values in the first bin of the radial distribution when they should zero or almost zero, so the program is forced to show that value
					write(40,*) (dble(i)+5.d-1)*l2/dble(binrad),dble(distbins(i))/(dble(radcount)*dens*(4.d0*pi*(l2/dble(binrad))**3)*(dble(i)**2+dble(i)+1.d0/3.d0))!writing the data of the radial distribution
				end do
				call flush(40)!forcing the program to write stuff
			end if
			if(measureprofile)then
				do i=1,meascylin
					write(90,*) (dble(i)-5.d-1)*radius/dble(meascylin),profile(i)/dble(meascounter)!writing the data of the velocity profile
				end do
				call flush(90)!forcing the program to write stuff
			end if
		end if
		if(xyz)close(20)!disconnecting the "position_output.xyz" file
		if(normaloutput)close(30)!disconnecting the "energy_temperature_pressure_diffusion_output.dat" file
		if(measureraddist)close(40)!disconnecting the "radial_distribution.dat" file
		if(measurevhist)close(60)!disconnecting the "velocity_histogram.gnu" file
		if(temp3d)close(70)!disconnecting the "temperature_3d_output.dat" file
		if(temp3d.and.press3d)close(80)!disconnecting the "pressure_3d_output.dat" file
		if(measureprofile)close(90)!disconnecting the "velocity_profile.dat" file
		return
	end subroutine output



	subroutine save_state()!if the current state is to be preserved it can be saved to a file with this subroutine
		integer(double_kind)::i=0!particle counter
		if(separateoutput)then
			open(50,file='./data/saved_state_'//datechar//'_'//timechar//'.dat',action='write',status='new',form='unformatted',access='stream')!the state of the system is saved here
		else
			open(50,file='./data/saved_state.dat',action='write',form='unformatted',access='stream')!the state of the system is saved here
		end if
		write(50) amnttot!writing out the total amount of particles in the simulation
		if(samemass)then
			do i=1,amnttot
				write(50) rvec(1,i),rvec(2,i),rvec(3,i),vvec(1,i),vvec(2,i),vvec(3,i),avec(1,i),avec(2,i),avec(3,i),fixed(i)!writing positions, velocities and forces when the masses are all equal (if constants needed replace avec with fvec)
			end do
		else
			do i=1,amnttot
				write(50) rvec(1,i),rvec(2,i),rvec(3,i),vvec(1,i),vvec(2,i),vvec(3,i),fvec(1,i),fvec(2,i),fvec(3,i),mvec(i),fixed(i)!writing positions, velocities and forces when the masses are different
			end do
		end if
		call flush(50)!forcing the program to write stuff
		close(50)!disconnecting the "saved_state.dat" file
		return
	end subroutine save_state



	subroutine load_state()!and to retrieve a saved state this is the subroutine needed
		logical::fileexists=.false.!this variable will be used to see if the file to load the system from exists
		integer(double_kind)::i=0!particle counter
		inquire(file='./data/saved_state.dat',exist=fileexists)!testing if the file is there
		if(fileexists)then
			open(50,file='./data/saved_state.dat',action='read',form='unformatted',access='stream')!opening the saved_state file
		else
			stop 'No file named "saved_state.dat" exists, the system must be loaded from that specific file, aborting.'!no file to load the system from, check if there is a file named saved_state_(date)_(time).dat
		end if
		read(50) amnttot!reading in the total amount of particles in the simulation
		if(amnttot>maxpartic)stop 'Not enough space for so many particles, the variable "maxpartic" should be at least as big as "amnttot", aborting.'!avoiding segmentation faults
		if(samemass)then
			do i=1,amnttot
				read(50) rvec(1,i),rvec(2,i),rvec(3,i),vvec(1,i),vvec(2,i),vvec(3,i),avec(1,i),avec(2,i),avec(3,i),fixed(i)!reading positions, velocities and forces when the masses are all equal (if constants needed replace avec with fvec)
			end do
		else
			do i=1,amnttot
				read(50) rvec(1,i),rvec(2,i),rvec(3,i),vvec(1,i),vvec(2,i),vvec(3,i),fvec(1,i),fvec(2,i),fvec(3,i),mvec(i),fixed(i)!reading positions, velocities and forces when the masses are different
			end do
		end if
		close(50)!disconnecting the "saved_state.dat" file
		return
	end subroutine load_state



	subroutine xyz_output()!output for VMD
		integer(double_kind)::i=0!particle counter
		do i=1,amnttot
			write(20,*) element,rvec(1,i),rvec(2,i),rvec(3,i)!writing positions
		end do
		call flush(20)!forcing the program to write stuff
		return
	end subroutine xyz_output



	real(double_kind) function length_to_SI(x)!this functions converts distances in units of the simulation to metres
		real(double_kind),intent(in)::x!the value to convert
		length_to_SI=x*convfactorlength!the conversion
	end function length_to_SI



	real(double_kind) function length_from_SI(x)!this functions converts distances in metres to units of the simulation
		real(double_kind),intent(in)::x!the value to convert
		length_from_SI=x/convfactorlength!the conversion
	end function length_from_SI



	real(double_kind) function area_to_SI(x)!this functions converts areas in units of the simulation to squared metres
		real(double_kind),intent(in)::x!the value to convert
		area_to_SI=x*convfactorlength**2!the conversion
	end function area_to_SI



	real(double_kind) function area_from_SI(x)!this functions converts areas in squared metres to units of the simulation
		real(double_kind),intent(in)::x!the value to convert
		area_from_SI=x/convfactorlength**2!the conversion
	end function area_from_SI



	real(double_kind) function volume_to_SI(x)!this functions converts volumes in units of the simulation to cubic metres
		real(double_kind),intent(in)::x!the value to convert
		volume_to_SI=x*convfactorlength**3!the conversion
	end function volume_to_SI



	real(double_kind) function volume_from_SI(x)!this functions converts volumes in cubic metres to units of the simulation
		real(double_kind),intent(in)::x!the value to convert
		volume_from_SI=x/convfactorlength**3!the conversion
	end function volume_from_SI



	real(double_kind) function time_to_SI(x)!this functions converts times in units of the simulation to seconds
		real(double_kind),intent(in)::x!the value to convert
		time_to_SI=x*convfactortime!the conversion
	end function time_to_SI



	real(double_kind) function time_from_SI(x)!this functions converts times in seconds to units of the simulation
		real(double_kind),intent(in)::x!the value to convert
		time_from_SI=x/convfactortime!the conversion
	end function time_from_SI



	real(double_kind) function mass_to_SI(x)!this functions converts masses in units of the simulation to kilograms
		real(double_kind),intent(in)::x!the value to convert
		mass_to_SI=x*convfactormass!the conversion
	end function mass_to_SI



	real(double_kind) function mass_from_SI(x)!this functions converts masses in kilograms to units of the simulation
		real(double_kind),intent(in)::x!the value to convert
		mass_from_SI=x/convfactormass!the conversion
	end function mass_from_SI



	real(double_kind) function energy_to_SI(x)!this functions converts energies in units of the simulation to joules
		real(double_kind),intent(in)::x!the value to convert
		energy_to_SI=x*convfactorenergy!the conversion
	end function energy_to_SI



	real(double_kind) function energy_from_SI(x)!this functions converts enegies in joules to units of the simulation
		real(double_kind),intent(in)::x!the value to convert
		energy_from_SI=x/convfactorenergy!the conversion
	end function energy_from_SI



	real(double_kind) function temperature_to_SI(x)!this functions converts temperatures in units of the simulation to joules
		real(double_kind),intent(in)::x!the value to convert
		temperature_to_SI=x*convfactortemperature!the conversion
	end function temperature_to_SI



	real(double_kind) function temperature_from_SI(x)!this functions converts temperatures in joules to units of the simulation
		real(double_kind),intent(in)::x!the value to convert
		temperature_from_SI=x/convfactortemperature!the conversion
	end function temperature_from_SI



	real(double_kind) function pressure_to_SI(x)!this functions converts pressures in units of the simulation to pascals
		real(double_kind),intent(in)::x!the value to convert
		pressure_to_SI=x*convfactorpressure!the conversion
	end function pressure_to_SI



	real(double_kind) function pressure_from_SI(x)!this functions converts pressure in pascals to units of the simulation
		real(double_kind),intent(in)::x!the value to convert
		pressure_from_SI=x/convfactorpressure!the conversion
	end function pressure_from_SI



end program Molecular_Dynamics_Simulator