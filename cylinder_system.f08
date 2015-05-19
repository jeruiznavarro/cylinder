program Molecular_Dynamics_in_3D!main program
	implicit none!all variables must be declared
	integer(8),parameter::maxbinrad=10000!maximum number of radial bins
	integer(8),parameter::maxbinhist=100!maximum number of bins for the velocity histogram
	integer(8),parameter::maxcyl=100!maximum number of cylinders
	integer(8),parameter::maxpartic=50000!maximum number of particles in the system
	integer(8),parameter::maxpartcell=1000!maximum number of particles in a single cell
	integer(8),parameter::maxcells=20!maximum number of cells per dimension
	real(8),parameter::kboltz=1.38d-8!Boltzmann's constant
	real(8),parameter::pi=dacos(-1.d0)!number pi
	logical::andersen=.false.!if true the andersen thermostat will be active
	logical::berendsen=.false.!if true the berendsen thermostat will be active
	logical::cylinders=.false.!if true, independent measurements will be taken in the different cylindrical intervals of the system
	logical::firststep=.true.!after the first call to this subroutine is completed, this will be false
	logical::fixed(maxpartic)=.false.!this is true for fixed particles and false for free ones that can move
	logical::gaussdistribution=.true.!if false, a uniform distribution for the velocities will be used instead of a gaussian one
	logical::gnuplot=.false.!if true, the program will output data to be visualized in gnuplot
	logical::halfdensity=.false.!if true, the fluid in the porous system will have half of its particles erased to achieve a half of the previous density
	logical::interaction(maxpartic,maxpartic)=.false.!this matrix stores which interactions have been calculated in each timestep so they don't have to be computed twice
	logical::measuring=.true.!if true, measurements will be taken
	logical::measureraddist=.false.!if true, the radial distribution function will be measured
	logical::measurevhist=.false.!if true, a histogram of the components of the velocities of the particles will be taken each time a measurement is performed
	logical::noinitialinterac=.false.!if false, the initial interactions between all the particles won't be computed
	logical::periodicspeed=.false.!if true, no matter how fast the particles move, the periodic boundary conditions will apply, but this might make some bugs much harder to detect (and uses a tiny bit more resources)
	logical::removedrift=.false.!if true, the average momentum of the system will be eliminated so there is no global drift
	logical::samemass=.false.!if true, all the particles have the same mass
	logical::saving=.false.!if true, the state of the system will be saved
	logical::separateoutput=.false.!if true, separate files with the current date and time will be used for the output of the measurements
	logical::timetomeasure=.false.!if true, a measurement of the properties of the system will be performed
	character(8)::datechar='--------'!this character variable stores the current date for the names of the output files
	character(3)::struc='---'!this variable will tell the program which structure should be generated, fcc for cubic face centered, fcy is the same but has a cylinder of free particles and scc for simple cubic and por for a porous matrix
	character(10)::timechar='----------'!this character variable stores the current time for the names of the output files
	integer(8)::amntfix=0!amount of fixed particles in the simulation
	integer(8)::amntfree=0!amount of free particles in the simulation
	integer(8)::amnttot=0!amount of total particles in the simulation
	integer(8)::binhist=0!number of bins for the velocity histograms
	integer(8)::binrad=0!number of bins for the radial distribution
	integer(8)::distbins(maxbinrad)=0!these bins will count how many particles are in the divisions of the (0,L/2] interval in which the radial distribution will be measured
	integer(8)::cells(maxpartcell,maxcells,maxcells,maxcells)=0!lists containing the tags of the particles in each cell, the first index are the particles, the second one's the x axis, the third one's y and the last one's z
	integer(8)::cvcount=0!this counter will keep track of how many measurements are used to calculate the heat capacity so the time average can be computed
	integer(8)::cvcountcy(maxcyl)=0!this counter will keep track of how many measurements are used to calculate the heat capacity for each cylinder so the time average can be computed
	integer(8)::difftag(maxcyl)=0!this variable is the tag of the particle that will be followed to measure diffusion properties in each cylinder
	integer(8)::iterations=0!number of iterations carried out by the program
	integer(8)::maxiter=0!maxmimum number of iterations
	integer(8)::meascounter=0!number of measurements taken
	integer(8)::meascylin=0!this is the number of cylindrical intervals in which there will be independent measurements of variables
	integer(8)::measiter=0!number of iterations between measurements
	integer(8)::numcells(3)=0!number of cells for each dimension (x=1, y=2, z=3)
	integer(8)::occup(maxcells,maxcells,maxcells)=0!number of particles in each cell, occupancy lists
	integer(8)::presscount=0!this counter keeps track of how many timesteps are used to obtain the pressure
	integer(8)::pores=0!amount of pores in the system
	integer(8)::radcount=0!this counts how many timesteps are used in the measurement of the radial distribution
	real(8)::cellsize(3)=0.d0!size of the cells for each dimension (x=1, y=2, z=3)
	real(8)::cellvolume=0.d0!volume of a single cell
	real(8)::cv=0.d0!heat capacity of the system
	real(8)::cvcy(maxcyl)=0.d0!heat capacity of each cylinder
	real(8)::dens=0.d0!number density of the system
	real(8)::denscells(maxcells,maxcells,maxcells)=0.d0!number density of each cell
	real(8)::displ(3)=0.d0!displacement vector to calculate the appropiate distance between particles with the periodic boundary conditions
	real(8)::ekine=0.d0!kinetic energy of the system
	real(8)::ekinecells(maxcells,maxcells,maxcells)=0.d0!kinetic energy of each cell
	real(8)::ekinecy(maxcyl)=0.d0!kinetic energy of each cylinder in the system
	real(8)::epot=0.d0!potential energy of the system
	real(8)::epsi=0.d0!depth parameter of the Lennard-Jones interaction
	real(8)::fvec(3,maxpartic)=0.d0!forces for each particle (fx->(1,:), fy->(2,:), fz->(3,:)) are here
	real(8)::intrad=0.d0!interaction radius, cut-off for the potential, beyond this distance there is no interaction
	real(8)::l2=0.d0!this is half of the smallest side of the simulation box and is used as the limit to look for particles when calculating the radial distribution
	real(8)::latcon=0.d0!lattice constant
	real(8)::mass=0.d0!generic mass of the particles
	real(8)::mvec(maxpartic)=0.d0!if the masses of the particles are not the same they can be stored here
	real(8)::measintv=0.d0!interval of time between measurements
	real(8)::porevolume=0.d0!the volume occupied by the pored if they were non-overlapping
	real(8)::porosity=0.d0!relative fraction of non-solid space in the system
	real(8)::potentialshift=0.d0!shift of the interaction potential
	real(8)::press=0.d0!pressure of the system
	real(8)::presscells(maxcells,maxcells,maxcells)=0.d0!pressure in each cell
	real(8)::radius=0.d0!radius of the non solid zone of the system where particles can flow
	real(8)::rvec(3,maxpartic)=0.d0!positions for each particle (x->(1,:), y->(2,:), z->(3,:)) are here
	real(8)::rvecini(3,maxpartic)=0.d0!initial positions for each particle (x->(1,:), y->(2,:), z->(3,:)) are here
	real(8)::rvecnonper(3,maxpartic)=0.d0!positions for each particle (x->(1,:), y->(2,:), z->(3,:)) without periodic boundary conditions (to calculate mean square displacement) are here
	real(8)::sigm=0.d0!size parameter of the Lennard-Jones interaction
	real(8)::siz(3)=0.d0!dimensions (x=1, y=2, z=3) of the simulation box
	real(8)::squaredispl=0.d0!mean square displacement used to calculate the diffusion constant
	real(8)::sumekine=0.d0!the sum of the kinetic energy over time is kept here
	real(8)::sumekinesqur=0.d0!the sum of the kinetic energy squared over time is kept here
	real(8)::sumekinecy(maxcyl)=0.d0!the sum of the kinetic energy over time for each cylinder is kept here
	real(8)::sumekinecysqur(maxcyl)=0.d0!the sum of the kinetic energy squared over time for each cylinder is kept here
	real(8)::temp=0.d0!temperature of the system
	real(8)::tempbath=0.d0!temperature of the heat bath
	real(8)::tempcells(maxcells,maxcells,maxcells)=0.d0!temperature in each cell
	real(8)::tempmax=0.d0!maximum temperature during equilibrium after termalization
	real(8)::tempmean=0.d0!mean temperature
	real(8)::tempmin=1.d100!minimum temperature during equilibrium after termalization
	real(8)::tempcy(maxcyl)=0.d0!temperature of each cylinder in the system
	real(8)::threshold=0.d0!this variable is used to test if a particle will be affected by the andersen thermostat
	real(8)::time=0.d0!time for which the simulation has been running
	real(8)::timetosave=0.d0!the state of the system will be saved after this point in time
	real(8)::tmax=0.d0!time when the simulation should stop
	real(8)::tmeas=0.d0!when time is greater or equal than this, measurements will start
	real(8)::trelax=0.d0!relaxation time for the andersen and berendsen thermostats
	real(8)::tstep=0.d0!timestep of the simulation
	real(8)::volume=0.d0!volume of the system
	real(8)::vvec(3,maxpartic)=0.d0!velocities for each particle (vx->(1,:), vy->(2,:), vz->(3,:)) are here



	call input()!reading initial parameters
	call startup()!setting up the system
	do!simulation loop
		if(iterations>maxiter)exit!the program ends when it reaches the maximum time
		if(firststep.and.(.not.gnuplot))write(*,*) 'Simulation started at: ',ctime(time8())!showing work in progress
		call half_leap(.false.)!updating the velocities of the particles without the thermostat, because it should only be applied to the speeds once per timestep
		call move_system()!updating the positions of the particles
		call update_particles_in_cells()!after the particles have moved, they may have entered a different cell and this needs to be checked and taken care of
		call calculate_forces()!the Verlet integrator requires to update the speeds again after moving the particle so the forces need to be calculated again
		call half_leap(.true.)!updating the velocities of the particles with the thermostat (if it's active)
		if(saving.and.(time>=timetosave))call save_state()!storing the generated structure in the apropiate file
		if((mod(iterations,measiter)==0).and.(measuring.and.(.not.gnuplot)))then
			timetomeasure=.true.!checking if a measurment should be taken
		else
			timetomeasure=.false.!checking if a measurment should be taken
		end if
		if(timetomeasure)then
			call normal_measuring()!with a certain frequency measurements are performed
		else if(cylinders.and.timetomeasure)then
			call measuring_with_cylinders()!with a certain frequency measurements are taken in each radial interval of the cylindrical fluid cavity of the system
		end if
		time=time+tstep!time advances
		iterations=iterations+1!one more iteration
		firststep=.false.!this is no longer true after the first iteration
		if((.not.gnuplot).and.(mod(iterations,maxiter/100)==0))write(*,*) nint(time*1.d2/tmax,8),'% completed at: ',ctime(time8())!showing work in progress
	end do
	call output()!finishing up writing the results



contains



	subroutine input()!reading initial parameters
		open(10,file='input.dat',action='read')!all the parameters are in this file
		read(10,*) andersen!if true the andersen thermostat will be active
		read(10,*) berendsen!if true the berendsen thermostat will be active
		read(10,*) binhist!number of bins for the velocity histograms
		read(10,*) binrad!number of bins for the radial distribution
		read(10,*) cylinders!if true, independent measurements will be taken in the different cylindrical intervals of the system
		read(10,*) epsi!depth parameter of the Lennard-Jones interaction
		read(10,*) gaussdistribution!if false, a uniform distribution for the velocities will be used instead of a gaussian one
		read(10,*) gnuplot!if true, the program will output data to be visualized in gnuplot
		read(10,*) halfdensity!if true, the fluid in the porous system will have half of its particles erased to achieve a half of the previous density
		read(10,*) intrad!interaction radius
		read(10,*) latcon!lattice constant
		read(10,*) mass!generic mass of the particles
		read(10,*) meascylin!this is the number of cylindrical intervals in which there will be independent measurements of variables
		read(10,*) measintv!interval of time between measurements
		read(10,*) measuring!if true, measurements will be taken
		read(10,*) measureraddist!if true, the radial distribution function will be measured
		read(10,*) measurevhist!if true, a histogram of the components of the velocities of the particles will be taken each time a measurement is performed
		read(10,*) noinitialinterac!if false, the initial interactions between all the particles won't be computed
		read(10,*) periodicspeed!if true, no matter how fast the particles move, the periodic boundary conditions will apply, but this might make some bugs much harder to detect (and uses a tiny bit more resources)
 		read(10,*) pores!amount of pores in the system
		read(10,*) radius!radius of the non solid zone of the system where particles can flow
		read(10,*) removedrift!if true, the average momentum of the system will be eliminated so there is no global drift
		read(10,*) samemass!if true, all the particles have the same mass
		read(10,*) saving!if true, the state of the system will be saved
		read(10,*) separateoutput!if true, separate files with the current date and time will be used for the output of the measurements
		read(10,*) sigm!size parameter of the Lennard-Jones interaction
		read(10,*) siz(1)!size in the x dimension of the simulation box
		read(10,*) siz(2)!size in the y dimension of the simulation box
		read(10,*) siz(3)!size in the z dimension of the simulation box
		read(10,*) struc!this variable will tell the program which structure should be generated, fcc for cubic face centered, fcy is the same but has a cylinder of free particles and scc for simple cubic and por for a porous matrix
		read(10,*) temp!temperature of the system
		read(10,*) tempbath!temperature of the heat bath
		read(10,*) timetosave!the state of the system will be saved after this point in time
		read(10,*) tmax!time when the simulation should stop
		read(10,*) tmeas!when time is greater or equal than this, measurements will start
		read(10,*) trelax!relaxation time for the andersen and berendsen thermostats
		read(10,*) tstep!timestep of the simulation
		close(10)!disconnecting the "intput.dat" file
		if(berendsen.and.andersen)stop "Both thermostats can't be active at the the same time, aborting."!this conflict is not allowed, so the program is terminated
		if((meascylin>maxcyl).and.(cylinders))stop 'More cylinders than allowed, the variable "maxcyl" should be at least as big as "meascylin", aborting.'!this conflict is not allowed, so the program is terminated
		return
	end subroutine input



	subroutine startup()!the system to be simulated is set up here
		integer(8)::seed(8)=0!seed for the random number generator
		call date_and_time(date=datechar,time=timechar,values=seed)!creating a seed using the clock of the computer
		call dran_ini(1000*seed(8)+3*seed(7)*seed(6)/10)!initializing the random number generator
		if(separateoutput)then
			open(20,file='./data/position_output_'//datechar//'_'//timechar//'.dat',action='write',status='new')!if the positions are to be visualized later they must be saved to this file
			open(30,file='./data/energy_temperature_pressure_diffusion_output_'//datechar//'_'//timechar//'.dat',action='write',status='new')!this is where the time evolution of these magnitudes will be kept
			open(40,file='./data/radial_distribution_'//datechar//'_'//timechar//'.dat',action='write',status='new')!and the radial distribution function values are here
			open(60,file='./data/velocity_histogram_'//datechar//'_'//timechar//'.dat',action='write',status='new')!the velocity histogram is kept here in a format such that its evolution can be visualiazed as an animation in gnuplot
			open(70,file='./data/temperature_pressure_3d_output_'//datechar//'_'//timechar//'.dat',action='write',status='new')!the 3d data of the
		else
			open(20,file='./data/position_output.dat',action='write')!if the positions are to be visualized later they must be saved to this file
			open(30,file='./data/energy_temperature_pressure_diffusion_output.dat',action='write')!this is where the time evolution of these magnitudes will be kept
			open(40,file='./data/radial_distribution.dat',action='write')!and the radial distribution function values are here
			open(60,file='./data/velocity_histogram.dat',action='write')!the velocity histogram is kept here in a format such that its evolution can be visualiazed as an animation in gnuplot
			open(70,file='./data/temperature_pressure_3d_output.dat',action='write')!the 3d data of the
		end if
		write(20,*) amnttot!the xyz format needs to have the amount of particles in the first line
		write(20,*) 'This is an optional comment, but it must be here for the file to be readable by the visualizing program.'!and this comment is mandatory
		write(30,*) '#	time	total energy	potential energy	kinetic energy	temperature	pressure	mean square displacement/6'!magnitudes in this file
		write(60,*) 'set key top right'!placing the title in a place where it won't be much of a nuissance
		write(60,*) 'set grid'!the grid makes it easier to read
		write(60,*) 'set boxwidth 1 relative'!this is the width of the histogram columns
		write(60,*) 'set style fill transparent solid 0.5 noborder'!and this is the style of the columns, so they are half opaque and the grid  can be visible behind them
		call flush(20)!forcing the program to write stuff
		call flush(30)!forcing the program to write stuff
		call flush(60)!forcing the program to write stuff
		if(gnuplot)then
			write(*,*) 'unset key'!this just gets in the way
			write(*,*) 'set border 4095'!gnuplot initial configuration, displaying the borders of the simulation box
			write(*,*) 'set ticslevel 0'!gnuplot initial configuration, setting the bottom of the simulation box at z=0
			write(*,*) 'set ztics mirror'!gnuplot initial configuration, setting the grid in the vertical axes
			write(*,*) 'set xrange [0:',siz(1),']'!gnuplot initial configuration, x size of the simulation box
			write(*,*) 'set yrange [0:',siz(2),']'!gnuplot initial configuration, y size of the simulation box
			write(*,*) 'set zrange [0:',siz(3),']'!gnuplot initial configuration, z size of the simulation box
			write(*,*) "splot '-' w p pt 7 lc 0"!new gnuplot frame
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
		else
			stop 'No defined structure, aborting.'!the program needs to know which system should be generated
		end if
		if(gnuplot)write(*,*) 'e'!end of gnuplot frame
		call create_cell_lists()!cell lists are needed to calculate the interactions and to speed up the simulation
		if(measurevhist)call velocity_histogram()!getting a first snapshot of the velocity distribution
		if(removedrift)call remove_average_momentum()!the system should't have a global net movement so it is removed here if there's any (this system is already driftless in a porous system because they are loaded from a previous one)
		if(noinitialinterac)call calculate_forces()!before the first integration iteration can be carried out, the forces need to be calculated
		return
	end subroutine startup



	subroutine calculate_some_constants()!just making the startup subroutine cleaner
		l2=dmin1(siz(1),siz(2),siz(3))*5.d-1!the value of L/2 is chosen to be half of the smallest side of the simulation box
		measiter=nint(measintv/tstep,8)!number of iterations between measurements
		maxiter=nint(tmax/tstep,8)!maximum number of iterations in the simulation
		potentialshift=4.d0*epsi*(intrad**(-6)-1.d0)*intrad**(-6)!potential shift so that the potential energy is a continious function after the cut-off
		threshold=tstep/trelax!this is the value that the threshold for the andersen thermostat should have
		volume=siz(1)*siz(2)*siz(3)!volume of the system
		return
	end subroutine calculate_some_constants



	subroutine create_fcc_structure()!a face centered cubic structure will be generated here
		integer(8)::i1=0!x dimension counter
		integer(8)::i2=0!y dimension counter
		integer(8)::i3=0!z dimension counter
		integer(8)::j=1!particle counter
		integer(8)::k=0!dimension counter
		integer(8)::l=0!auxiliary particle counter
		real(8)::dran_g!gaussian random number
		real(8)::dran_u!uniform random number
		do i1=0,nint(siz(1)/latcon,8)-1
			do i2=0,nint(siz(2)/latcon,8)-1
				do i3=0,nint(siz(3)/latcon,8)-1
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
					if(gaussdistribution)then
						do l=0,3
							do k=1,3
								vvec(k,j+l)=dsqrt(kboltz*temp/mass)*dran_g()!random gaussian components for the velocities
								rvecini(k,j+l)=rvec(k,j+l)!storing the initial position of the particles
							end do
						end do
					else
						do l=0,3
							do k=1,3
								vvec(k,j+l)=3.d0*dsqrt(kboltz*temp/mass)*(dran_u()-5.d-1)!random uniform components for the velocities
								rvecini(k,j+l)=rvec(k,j+l)!storing the initial position of the particles
							end do
						end do
					end if
					if(.not.samemass)then
						mvec(j)=mass*dsin(dble(j))!setting the mass of the first particle
						mvec(j+1)=mass*dsin(dble(j))!setting the mass of the second particle
						mvec(j+2)=mass*dsin(dble(j))!setting the mass of the third particle
						mvec(j+3)=mass*dsin(dble(j))!setting the mass of the fourth particle
					end if
					j=j+4!next particles
					if(gnuplot)then
						do k=0,3
							write(*,*) rvec(1,j+k),rvec(2,j+k),rvec(3,j+k)!data feed to gnuplot
						end do
					end if
				end do
			end do
		end do
		amnttot=j-1!total amount of particles
		if(amnttot>maxpartic)then
			write(*,*) 'Not enough space for so many particles, the variable "maxpartic" should be at least as big as "amnttot", ',amnttot,' in this case.'
			stop 'Aborting.'!avoiding segmentation faults
		end if
		amntfree=amnttot!just in case, to avoid problems since there are no fixed particles in this situation
		dens=amnttot/volume!this is the number density of all the particles
		return
	end subroutine create_fcc_structure



	subroutine create_fcy_structure()!a face centered cubic structure with a cylinder of free particles will be generated here
		integer(8)::i=0!particle counter
		integer(8)::j=0!dimension counter
		call load_state()!the pores are created in an already thermalized system that is in a fluid, preferably liquid, state, so the particles are loaded from a previously prepared file
		amntfree=amnttot!right now all the particles are free
		do i=1,amnttot
			do j=1,3
				rvecini(j,i)=rvec(j,i)!storing the initial position of the particles
			end do
			if(dsqrt((rvec(1,i)-siz(1)*5.d-1)**2+(rvec(2,i)-siz(2)*5.d-1)**2)>=radius)then
				fixed(i)=.true.!fixing particles inside the pores
				do j=1,3
					vvec(j,i)=0.d0!now the particle is fixed so the velocity is not needed anymore
					fvec(j,i)=0.d0!same with the force
				end do
				amntfix=amntfix+1!one more fixed particle
				amntfree=amntfree-1!one less free particle
			end if
		end do
		dens=amnttot/volume!this is the number density of all the particles
		return
	end subroutine create_fcy_structure



	subroutine create_porous_structure()!a structure composed of a fcc crystal with a matrix of empty cylinders is generated here
		integer(8)::i=0!particle counter
		integer(8)::j=0!pore counter
		integer(8)::k=0!dimension counter
		real(8)::dran_u!uniform random number
		real(8)::porematrix(3,maxcyl)!the coordinates (x->(1,:), y->(2,:), z->(3,:)) of the centres of the pores will be kept here
		real(8)::poreradius(maxcyl)!the radius of each pore
		do j=1,pores
			do k=1,3
				porematrix(k,j)=dran_u()*siz(k)!random positions for the pores
			end do
			poreradius(j)=dran_u()+2.d0!random radii for the pores
			porevolume=porevolume+4.d0*pi*poreradius(j)/3.d0!adding another non-overlapping volume to the total
		end do
		call load_state()!the pores are created in an already thermalized system that is in a fluid, preferably liquid, state, so the particles are loaded from a previously prepared file
		amntfree=amnttot!right now all the particles are free
		do i=1,amnttot
			do j=1,pores
				if(dsqrt((rvec(1,i)-porematrix(1,j))**2+(rvec(2,i)-porematrix(2,j))**2+(rvec(3,i)-porematrix(3,j))**2)<=poreradius(j))then
					fixed(i)=.true.!fixing particles inside the pores
					amntfix=amntfix+1!one more fixed particle
					amntfree=amntfree-1!one less free particle
				else if(fixed(i))then
					cycle!this particle is already fixed, moving onto the next one
				end if
			end do
		end do
		if(halfdensity)then
			do
				i=1!starting particle counter
				if((.not.fixed(i)).and.(dran_u()<5.d-1))then
					do k=1,3
						rvec(k,i)=rvec(k,amnttot)!destroying half of the particles to create a half-density fluid
						rvec(k,amntfree)=0.d0!erasing old unnecessary information
						vvec(k,i)=vvec(k,amnttot)!destroying half of the particles to create a half-density fluid
						vvec(k,amntfree)=0.d0!erasing old unnecessary information
						fvec(k,i)=fvec(k,amnttot)!destroying half of the particles to create a half-density fluid
						fvec(k,amntfree)=0.d0!erasing old unnecessary information
					end do
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
		end if
		porosity=amntfree/amnttot!approximately finding out how much free space is there in the system, exploiting the fact that the particles are homogeneously distributed in a fluid state
		dens=amnttot/volume!this is the number density of all the particles
		return
	end subroutine create_porous_structure



	subroutine create_sc_structure()!a simple cubic structure will be generated here
		integer(8)::i1=0!x dimension counter
		integer(8)::i2=0!y dimension counter
		integer(8)::i3=0!z dimension counter
		integer(8)::j=1!particle counter
		integer(8)::k=0!dimension counter
		integer(8)::vecindx(3)=0!here the indeces relevant to each dimensions are kept, for example the indeces of the cells in the crystalline lattice
		real(8)::dran_g!gaussian random number
		real(8)::dran_u!uniform random number
		do i1=0,nint(siz(1)/latcon,8)-1
			vecindx(1)=i1!index of the current cell in the first dimension
			do i2=0,nint(siz(2)/latcon,8)-1
				vecindx(2)=i2!index of the current cell in the second dimension
				do i3=0,nint(siz(3)/latcon,8)-1
					vecindx(3)=i3!index of the current cell in the third dimension
					do k=1,3
						rvec(k,j)=latcon*dble(vecindx(k))!position in the solid lattice of the material, free particles will leave these as they move
					end do
					if(dsqrt((rvec(1,j)-5.d-1*siz(1))**2+(rvec(2,j)-5.d-1*siz(1))**2)>=radius)then
						fixed(j)=.true.!the particles outside of this cylinder are fixed
						amntfix=amntfix+1!counting the number of fixed particles
					else
						if(gaussdistribution)then
							do k=1,3
								vvec(3,j)=dsqrt(kboltz*temp/mass)*dran_g()!random gaussian components for the velocities of the free particles
								rvecini(k,j)=rvec(k,j)!storing the initial position of the particles
							end do
						else
							do k=1,3
								vvec(3,j)=dsqrt(kboltz*temp/mass)*(dran_u()-5.d-1)!random gaussian components for the velocities of the free particles
								rvecini(k,j)=rvec(k,j)!storing the initial position of the particles
							end do
						end if
						amntfree=amntfree+1!counting the number of free particles
					end if
					if(.not.samemass)then
						mvec(j)=mass*dsin(dble(j))!setting the mass of the particles
					end if
					j=j+1!next particle
					if(gnuplot)write(*,*) rvec(1,j),rvec(2,j),rvec(3,j)!data feed to gnuplot
				end do
			end do
		end do
		amnttot=j-1!total amount of particles
		if(amnttot>maxpartic)then
			write(*,*) 'Not enough space for so many particles, the variable "maxpartic" should be at least as big as "amnttot", ',amnttot,' in this case.'
			stop 'Aborting.'!avoiding segmentation faults
		end if
		dens=amntfree/(siz(3)*pi*radius**2)!this is the number density of the free particles
		return
	end subroutine create_sc_structure



	subroutine create_cell_lists()!cell lists are used to speed up the simulation, see the article in wikipedia for more general details
		integer(8)::i=0!particle counter
		integer(8)::inds(3)=0!indeces of the cell in which a particle is
		integer(8)::j=0!dimension counter
		do j=1,3
			numcells(j)=floor(siz(j)/intrad,8)!computing how many cells each dimension has
			if(numcells(j)>maxcells)then
				write(*,*) 'Not enough space for so many cells, the variable "maxcells" should be at least as big as "numcells", ',numcells(j),' in this case.'
				stop 'Aborting.'!avoiding segmentation faults
			end if
			cellsize(j)=siz(j)/dble(numcells(j))!calculating how big the cells are
		end do
		cellvolume=cellsize(1)*cellsize(2)*cellsize(3)!computing the volume of a single cell
		do i=1,amnttot
			do j=1,3
				inds(j)=floor(rvec(j,i)/cellsize(j),8)+1!particle j is the cell which has these indeces
			end do
			occup(inds(1),inds(2),inds(3))=occup(inds(1),inds(2),inds(3))+1!updating the occupancy list with the new particle
			cells(occup(inds(1),inds(2),inds(3)),inds(1),inds(2),inds(3))=i!adding the new particle to the apropiate cell
		end do
		return
	end subroutine create_cell_lists



	subroutine remove_average_momentum()!to avoid the system having a drift, the average momentum must be substracted to all the particles in the system
		integer(8)::i=0!particle counter
		integer(8)::j=0!dimension counter
		real(8)::mmntsum(3)=0.d0!the global momentum will be stored here
		if(samemass)then
			do i=1,amnttot
				if(fixed(i))cycle!fixed particles don't contribute to the average momentum because they have no velocity
				do j=1,3
					mmntsum(j)=mmntsum(j)+vvec(j,i)*mass!calculating the global momentum if the particles have all the same mass
				end do
			end do
			mmntsum=mmntsum/dble(amnttot)!dividing by the number of particles so that the final result is the average momentum of the system
			do i=1,amnttot
				if(fixed(i))cycle!the velocity of fixed particles shouldn't change
				do j=1,3
					vvec(j,i)=vvec(j,i)-mmntsum(j)/mass!substracting the average momentum to the momentum of each particle and converting it back to velocity
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
		integer(8)::i=0!x dimension cell counter
		integer(8)::j=0!y dimension cell counter
		integer(8)::k=0!z dimension cell counter
		integer(8)::m=0!dimension counter
		integer(8)::n=0!particle counter in the cell lists
		integer(8)::tag=0!particle tag
		real(8)::dran_g!gaussian random number
		real(8)::dran_u!uniform random number
		real(8)::andersenfactor=0.d0!to avoid computing the factor of the andersen thermostat many times it's calculated only once at the beginning and stored here
		real(8)::berendsenfactor=0.d0!to avoid computing the factor of the berendsen thermostat many times it's calculated only once at the beginning and stored here
		real(8)::ekineberendsen=0.d0!instantaneous kinetic energy for the berendsen temperature
		real(8)::tempberendsen=0.d0!this temperature is updated every timestep to allow the berendsen thermostat to work properly
		if(firststep)tempberendsen=temp!setting the initial value for the instantaneous temperature
		if(andersen)andersenfactor=dsqrt(kboltz*tempbath/mass)!this is the amount that will multiply the components of the velocities in the andersen thermostat if the masses are the same
		if(berendsen)berendsenfactor=dsqrt(1.d0+(tempbath/tempberendsen-1.d0)*threshold)!this is the amount that will multiply the components of the velocities in the berendsen thermostat
		if(berendsen.and.(.not.thermostat))ekineberendsen=0.d0!resetting the instantaneous kinetic energy in the first call of the time step to the subroutine
		do k=1,numcells(3)
			do j=1,numcells(2)
				do i=1,numcells(1)
					do n=1,occup(i,j,k)
						tag=cells(n,i,j,k)!this value is going to be called from memory many times so it will be kept here for convenience
						if(fixed(tag))cycle!fixed particles don't move
						if(samemass)then
							if((.not.berendsen).and.(.not.andersen))then
								do m=1,3
									vvec(m,tag)=vvec(m,tag)+fvec(m,tag)*tstep*5.d-1/mass!updating speeds without a thermostat if all the particles have the same mass
								end do
							else if(berendsen.and.thermostat)then
								do m=1,3
									vvec(m,tag)=berendsenfactor*(vvec(m,tag)+fvec(m,tag)*tstep*5.d-1/mass)!updating speeds with the berendsen thermostat
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
									vvec(m,tag)=andersenfactor*dsqrt(mass/mvec(tag))*dran_g()!updating speeds with the andersen thermostat
								end do
							end if
						end if
					end do
				end do
			end do
		end do
		if(berendsen.and.(.not.thermostat))then
			if(samemass)then
				tempberendsen=(ekineberendsen*mass)/(3.d0*dble(amntfree)*kboltz)!instantaneous temperature
			else
				tempberendsen=(ekineberendsen)/(3.d0*dble(amntfree)*kboltz)!instantaneous temperature
			end if
		end if
		return
	end subroutine half_leap



	subroutine move_system()!this is where the actual movement of the particles takes place
		integer(8)::i=0!x dimension cell counter
		integer(8)::j=0!y dimension cell counter
		integer(8)::k=0!z dimension cell counter
		integer(8)::m=0!dimension counter
		integer(8)::n=0!particle counter in the cell lists
		integer(8)::tag=0!particle tag
		if(gnuplot)write(*,*) "splot '-' w p pt 7 lc 0"!new gnuplot frame
		do k=1,numcells(3)
			do j=1,numcells(2)
				do i=1,numcells(1)
					do n=1,occup(i,j,k)
						tag=cells(n,i,j,k)!this value is going to be called from memory many times so it will be kept here for convenience
						if(fixed(tag))cycle!fixed particles don't move
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
									rvec(m,tag)=rvec(m,tag)-dble(ceiling(rvec(m,tag)/siz(m),8))*siz(m)!same as with .not.periodicspeed, but it will work no matter how far away the particle moves
									rvecnonper(m,tag)=rvecnonper(m,tag)+dble(ceiling(rvec(m,tag)/siz(m),8))*siz(m)!same as with .not.periodicspeed, but it will work no matter how far away the particle moves
								else if(rvec(m,tag)>=siz(m))then
									rvec(m,tag)=rvec(m,tag)-dble(floor(rvec(m,tag)/siz(m),8))*siz(m)!same as with .not.periodicspeed, but it will work no matter how far away the particle moves
									rvecnonper(m,tag)=rvecnonper(m,tag)+dble(floor(rvec(m,tag)/siz(m),8))*siz(m)!same as with .not.periodicspeed, but it will work no matter how far away the particle moves
								end if
							end if
						end do
						if(gnuplot)write(*,*) rvec(1,tag),rvec(2,tag),rvec(3,tag)!data feed to gnuplot
					end do
				end do
			end do
		end do
		if(gnuplot)write(*,*) 'e'!end of gnuplot frame
		return
	end subroutine move_system



	subroutine update_particles_in_cells()!the particles may be in a different cell after moving, so the lists need to be updated
		integer(8)::i=0!x dimension cell counter
		integer(8)::inew=0!x dimension, new-cell counter
		integer(8)::j=0!y dimension cell counter
		integer(8)::jnew=0!y dimension, new-cell counter
		integer(8)::k=0!z dimension cell counter
		integer(8)::knew=0!z dimension, new-cell counter
		integer(8)::n=0!particle counter in the cell lists
		integer(8)::tag=0!particle tag
		do k=1,numcells(3)
			do j=1,numcells(2)
				do i=1,numcells(1)
					n=1!starting particle counter at the beginning of the current list
					do
						tag=cells(n,i,j,k)!this value is going to be called from memory many times so it will be kept here for convenience
						if(fixed(tag))cycle!fixed particles don't need updating
						inew=floor(rvec(1,tag)/cellsize(1),8)+1!computing the x index of the cell in which the particle is located after moving
						jnew=floor(rvec(2,tag)/cellsize(2),8)+1!computing the y index of the cell in which the particle is located after moving
						knew=floor(rvec(3,tag)/cellsize(3),8)+1!computing the z index of the cell in which the particle is located after moving
						if((i/=inew).or.(j/=jnew).or.(k/=knew))then!after the particle has moved, the cell in which it is contained may have changed so it needs to be updated
							occup(inew,jnew,knew)=occup(inew,jnew,knew)+1!updating the new occupancy list with the particle
							cells(occup(inew,jnew,knew),inew,jnew,knew)=tag!adding the particle to the corresponding new cell
							cells(n,i,j,k)=cells(occup(i,j,k),i,j,k)!taking away the updated particle from the old cell, its place is occupied by the last one in that cell
							cells(occup(i,j,k),i,j,k)=0!erasing old information
							occup(i,j,k)=occup(i,j,k)-1!updating the old occupancy list
						end if
						n=n+1!next particle
						if(n>occup(i,j,k))exit
					end do
				end do
			end do
		end do
		return
	end subroutine update_particles_in_cells



	subroutine calculate_forces()!obtaining the values for the forces of the particles
		integer(8)::i1=0!x dimension cell counter for the first particle
		integer(8)::i2=0!x dimension cell counter for the second particle
		integer(8)::i3=0!x dimension auxiliary cell counter
		integer(8)::j1=0!y dimension cell counter for the first particle
		integer(8)::j2=0!y dimension cell counter for the second particle
		integer(8)::j3=0!y dimension auxiliary cell counter
		integer(8)::k1=0!z dimension cell counter for the first particle
		integer(8)::k2=0!z dimension cell counter for the second particle
		integer(8)::k3=0!z dimension auxiliary cell counter
		integer(8)::n1=0!particle counter in the cells, for the first one
		integer(8)::n2=0!particle counter in the cells, for the second one
		integer(8)::tag1=0!tag for the first particle
		integer(8)::tag2=0!tag for the second particle
		do n1=1,amnttot
			do n2=1,3
				fvec(n2,n1)=0.d0!setting forces to zero
			end do
			do n2=1,amnttot
				interaction(n2,n1)=.false.!the interactions must be reset so they can be taken into account again in the new iteration
			end do
		end do
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
		presscount=presscount+1!this counter keeps track of how many timesteps are used to sample the pressure
		return
	end subroutine calculate_forces



	subroutine periodic_interaction(dimen,inindex,outindex)!in the case that some cell is beyond the limits, this will apply the periodic boundary conditions to the interaction between the particles
		integer(8),intent(in)::dimen!dimension in which the subroutine is carried out, 1 is x-axis, 2 is y-axis and 3 is z-axis
		integer(8),intent(in)::inindex!this number runs through the range of possible values for neighbouring cells
		integer(8),intent(out)::outindex!and this is the resulting index of that cell
		if(inindex==0)then
			displ(dimen)=-siz(dimen)!if the neighbouring cell is at the end of the simulation box in the x direction the displacement vector has this value
			outindex=numcells(dimen)!and this is the actual index of that cell
		else if(inindex==numcells(dimen)+1)then
			displ(dimen)=siz(dimen)!if the neighbouring cell is at the beginning of the simulation box in the x direction the displacement vector has this value
			outindex=1!and this is the actual index of that cell
		else
			displ(dimen)=0.d0!if the case is not one the former, the displacement vector is null
			outindex=inindex!and the actual index of the cell doesn't change
		end if
		return
	end subroutine periodic_interaction



	subroutine interacting(i,j,k,tag1,tag2)!this is where the actual interaction happens
		integer(8),intent(in)::i!x dimension cell counter
		integer(8),intent(in)::j!y dimension cell counter
		integer(8),intent(in)::k!z dimension cell counter
		integer(8)::l=0!dimension counter
		integer(8),intent(in)::tag1!tag for the first particle
		integer(8),intent(in)::tag2!tag for the second particle
		real(8)::dist=0.d0!distance between particles
		real(8)::distvec(3)=0.d0!vectorial distance between particles
		real(8)::force=0.d0!modulus of the force vector
		interaction(tag1,tag2)=.true.!accounting for this interaction
		do l=1,3
			distvec(l)=rvec(l,tag1)-rvec(l,tag2)-displ(l)!these are the components of the vectorial distance between the particles
		end do
		dist=dsqrt(distvec(1)**2+distvec(2)**2+distvec(3)**2)!and this is the distance between the two particles
		if(fixed(tag2).or.fixed(tag1))then!the interaction can be slightly different between fixed and free particles
			force=2.4d1*epsi*(sigm/dist)**6*(2.d0*(sigm/dist)**6-1.d0)/dist**2!this is the module of the force between a free and a fixed particle
		else
			force=2.4d1*epsi*(sigm/dist)**6*(2.d0*(sigm/dist)**6-1.d0)/dist**2!and this is the module of the force between two free particles
		end if
		if((.not.fixed(tag1)).and.(.not.fixed(tag2)))then
			do l=1,3
				fvec(l,tag1)=fvec(l,tag1)+distvec(l)*force!and these are the components of the force vector for the first particle
				fvec(l,tag2)=fvec(l,tag2)-distvec(l)*force!the same for the second particle
			end do
		else if(fixed(tag2).and.(.not.fixed(tag1)))then
			do l=1,3
				fvec(l,tag1)=fvec(l,tag1)+distvec(l)*force!and these are the components of the force vector for the first particle
			end do
		else if(fixed(tag1).and.(.not.fixed(tag2)))then
			do l=1,3
				fvec(l,tag2)=fvec(l,tag2)-distvec(l)*force!and these are the components of the force vector for the first particle
			end do
		end if
		if(timetomeasure.and.cylinders)then
			!#####HERE MAYBE THERE SHOULD BE A CONTRIBUTION TO THE POTENTIAL ENERGY FOR PARTICLES LOCATED EXCLUSIVELY IN EACH CYLINDER
		else if(timetomeasure)then
			epot=epot+4.d0*epsi*(sigm/dist)**6*((sigm/dist)**6-1.d0)-potentialshift!this is the potential energy
		end if
		press=press+dist**2*force!to calculate the pressure, the sum over all interactions must be carried out
		presscells(i,j,k)=presscells(i,j,k)+dist**2*force!same as above for the current cell
		return
	end subroutine interacting



	subroutine normal_measuring()!here the different measurements are performed
		integer(8)::freecounter=0!the amount of free particles per cell is kept here
		integer(8)::i=0!x dimension cell counter
		integer(8)::j=0!y dimension cell counter
		integer(8)::k=0!z dimension cell counter
		integer(8)::n=0!particle counter
		integer(8)::tag=0!particle tag
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
						ekinecells(i,j,k)=ekinecells(i,j,k)*mass*5.d-1!don't forget the mass and the 1/2
						tempcells(i,j,k)=(2.d0*ekinecells(i,j,k))/(3.d0*dble(freecounter)*kboltz)!temperature in the current cell
						denscells(i,j,k)=dble(freecounter)/cellvolume!and same with the number density
						presscells(i,j,k)=denscells(i,j,k)*tempcells(i,j,k)*kboltz+presscells(i,j,k)/(3.d0*cellvolume*dble(presscount))!final result for the pressure
					end do
				end do
			end do
			ekine=ekine*mass*5.d-1!don't forget the mass and the 1/2
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
						ekinecells(i,j,k)=ekinecells(i,j,k)*5.d-1!don't forget the mass and the 1/2
						tempcells(i,j,k)=(2.d0*ekinecells(i,j,k))/(3.d0*dble(freecounter)*kboltz)!temperature in the current cell
						denscells(i,j,k)=dble(freecounter)/cellvolume!and same with the number density
						presscells(i,j,k)=denscells(i,j,k)*tempcells(i,j,k)*kboltz+presscells(i,j,k)/(3.d0*cellvolume*dble(presscount))!final result for the pressure
					end do
				end do
			end do
			ekine=ekine*5.d-1!don't forget the 1/2
		end if
		sumekine=sumekine+ekine!summing the kinetic energy
		sumekinesqur=sumekinesqur+ekine**2!and the squared kinetic energy
		cvcount=cvcount+1!another measurement for the heat capacity
		temp=(2.d0*ekine)/(3.d0*dble(amntfree)*kboltz)!and with the kinetic energy the temperature is measured too
		if(iterations>0)then
			tempmax=dmax1(temp,tempmax)!recursively finding the maximum temperature
			tempmean=tempmean+temp!keeping track of the mean temperature
			tempmin=dmin1(temp,tempmin)!recursively finding the minimum temperature
		end if
		press=dens*temp*kboltz+press/(3.d0*volume*dble(presscount))!final result for the pressure
		if(measurevhist)call velocity_histogram()!obtaining another snapshot of the velocity distribution
		if(measureraddist)call radial_distribution()!getting a new sample of the radial distribution
		call diffusion()!measuring the diffusion of the particles
		meascounter=meascounter+1!another measurement
		write(30,*) time,epot+ekine,epot,ekine,temp,press,squaredispl/(6.d0*amnttot)!writing data
		call flush(30)!forcing the program to write stuff
		do k=1,numcells(3)
			do j=1,numcells(2)
				do i=1,numcells(1)
					write(70,*) (dble(i)-5.d-1)*cellsize(1),(dble(j)-5.d-1)*cellsize(2),(dble(k)-5.d-1)*cellsize(3),tempcells(i,j,k),presscells(i,j,k),time!writing data for spatially distributed magnitudes
					call flush(70)!forcing the program to write stuff
					presscells(i,j,k)=0.d0!setting the pressure in the current cell to zero so it can be calculated again
				end do
			end do
		end do
		epot=0.d0!setting the potential energy to zero so it can be calculated again
		press=0.d0!discarding the old value of the pressure so the new one can be computed in the next steps
		presscount=0!setting the counter to zero for the next measurement of the presure
		return
	end subroutine normal_measuring



	subroutine diffusion()!here measurements of the square displacement of the particles are taken
		integer(8)::i=0!particle counter
		integer(8)::j=0!dimension counter
		real(8)::tmp=0.d0!temporary variable to calculate the distance travelled by the particles
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



	subroutine measuring_with_cylinders()!here the different measurements are performed in each cylinder, whether all the particles have the same mass or not
		integer(8)::i=0!particle counter
		integer(8)::j=0!cylinder counter
		do j=1,meascylin
			ekinecy(j)=0.d0!setting the kinetic energy of the cylinders to zero to calculate the new values
		end do
		if(samemass)then
			do i=1,amnttot
				if(.not.fixed(i))then
					j=floor(dble(meascylin)*dsqrt(rvec(1,i)**2+rvec(2,i)**2)/radius,8)!finding in which cylinder the particle is located
					if(j>meascylin)j=meascylin!just in case the particle is slightly beyond the radius of the cylinder
					ekinecy(j)=ekinecy(j)+mass*(rvec(1,i)**2+rvec(2,i)**2+rvec(3,i)**2)!measuring the kinetic energy per cylinder
				end if
			end do
			do j=1,meascylin
				ekinecy(j)=ekinecy(j)*5.d-1!don't forget the 1/2
				sumekinecy(j)=sumekinecy(j)+ekinecy(j)!summing the kinetic energy per cylinder
				sumekinecysqur(j)=sumekinecysqur(j)+ekinecy(j)**2!and the squared kinetic energy per cylinder
				cvcountcy(j)=cvcountcy(j)+1!another measurement for the heat capacity per cylinder
				tempcy(j)=(2.d0*ekinecy(j))/(3.d0*amntfree*kboltz)!and with the kinetic energy the temperature is measured too per cylinder
			end do
		else
			do i=1,amnttot
				if(.not.fixed(i))then
					j=floor(dble(meascylin)*dsqrt(rvec(1,i)**2+rvec(2,i)**2)/radius,8)!finding in which cylinder the particle is located
					if(j>meascylin)j=meascylin!just in case the particle is slightly beyond the radius of the cylinder
					ekinecy(j)=ekinecy(j)+mvec(i)*(rvec(1,i)**2+rvec(2,i)**2+rvec(3,i)**2)!measuring the kinetic energy per cylinder
				end if
			end do
			do j=1,meascylin
				ekinecy(j)=ekinecy(j)*5.d-1!don't forget the 1/2
				sumekinecy(j)=sumekinecy(j)+ekinecy(j)!summing the kinetic energy per cylinder
				sumekinecysqur(j)=sumekinecysqur(j)+ekinecy(j)**2!and the squared kinetic energy per cylinder
				cvcountcy(j)=cvcountcy(j)+1!another measurement for the heat capacity per cylinder
				tempcy(j)=(2.d0*ekinecy(j))/(3.d0*amntfree*kboltz)!and with the kinetic energy the temperature is measured too per cylinder
			end do
		end if
		if(measurevhist)call velocity_histogram()!getting a sample of the velocities
		call diffusion_with_cylinders()!obtaining another diffusion measurement
		return
	end subroutine measuring_with_cylinders



	subroutine diffusion_with_cylinders()!here measurements of the square displacement of the particles are taken for each cylinder of the system
		real(8)::dran_u!uniform random number
		do
			difftag(1)=nint(dble(amntfree)*dran_u(),8)!choosing randomly the particle that will be used to measure diffusion properties
			if(difftag(1)>0)exit!just in case the tag is 0 or negative
		end do
		return
	end subroutine diffusion_with_cylinders



	subroutine radial_distribution()!computing the radial distribution
		integer(8)::bintag=0!tag of the current bin
		integer(8)::i=0!particle counter
		integer(8)::tag=0!tag for the first particle
		real(8)::dist=1.d10!distance between particles
		real(8)::tmp=0.d0!temporary storage for the recursive search of the minimum distance to the centre of the simulation box
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
				bintag=nint(dist*dble(binrad)/l2,8)!this is the bin in which the particle is
				if(bintag>binrad)bintag=binrad!just in case the nint function fails
				if(bintag<1)bintag=1!just in case the particles are too close
				distbins(bintag)=distbins(bintag)+1!counting another particle in the appropiate bin
			end if
		end do
		radcount=radcount+1!one more sample for the measurement
		return
	end subroutine radial_distribution



	subroutine velocity_histogram()!computing the velocity histogram
		character(26)::timetitle='------------------------'!this is used to convert double precision reals to characters
		integer(8)::histogram(3,maxbinhist)=0!in this vector the values of the velocity histogram will be kept for each component (x->(1,:), y->(2,:), z->(3,:))
		integer(8)::i=0!particle counter
		integer(8)::j=0!dimension counter
		integer(8)::temph(3)=0!temporal placeholder (x->(1), y->(2), z->(3))
		real(8)::maxv(3)=0.d0!maximum value found for the components of the velocity vector (x->(1), y->(2), z->(3))
		real(8)::minv(3)=0.d0!minimum value found for the components of the velocity vector (x->(1), y->(2), z->(3))
		real(8)::sizeofbins(3)=0.d0!how big the bins are for each dimension (x->(1), y->(2), z->(3))
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
				temph(j)=ceiling((vvec(j,i)-minv(j))/sizeofbins(j),8)!finding to which bin this value belongs
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
		integer(8)::i=0!cylinder and bin counter
		if(cylinders)then
			do i=1,meascylin
				cvcy(i)=(3.d0*kboltz)/((4.d0*amntfree*(sumekinecy(i)**2/dble(cvcountcy(i))-sumekinecysqur(i))/(dble(cvcountcy(i))*3.d0*kboltz*tempcy(i)**2))-2.d0)!computing the final value of the heat capacity
				write(30,*) '# The heat capacity in cylinder',i,'is:',cvcy(i)!writing data
			end do
			write(30,*) '# And the mean temperature is:',tempmean!writing data
		else
			do i=1,binrad
				write(40,*) (dble(i)+5.d-1)*l2/dble(binrad),dble(distbins(i))/(dble(radcount)*dens*(4.d0*pi*(l2/dble(binrad))**3)*(dble(i)**2+dble(i)+1/3.d0))!writing the data of the radial distribution
			end do
			cv=(3.d0*kboltz)/((4.d0*amnttot*(sumekine**2/dble(cvcount)-sumekinesqur)/(dble(cvcount)*3.d0*kboltz*temp**2))-2.d0)!computing the final value of the heat capacity
			write(30,*) '# The heat capacity is:',cv,', the mean temperature is:',tempmean/dble(meascounter),', the maximum temperature is:',tempmax,', the minimum temperature is:',tempmin!writing data
			if(struc=='por')write(30,*) '# The porosity is:',porosity,', the non-overlapping porosity is:',porevolume,', and the ratio between them is:',porosity/porevolume!writing porosity data
		end if
		close(20)!disconnecting the "position_output.dat" file
		close(30)!disconnecting the "energy_temperature_pressure_diffusion_output.dat" file
		close(40)!disconnecting the "radial_distribution.dat" file
		close(60)!disconnecting the "velocity_histogram.dat" file
		close(70)!disconnecting the "temperature_pressure_3d_output.dat" file
		return
	end subroutine output



	subroutine save_state()!if the current state is to be preserved it can be saved to a file with this subroutine
		integer(8)::i=0!particle counter
		if(separateoutput)then
			open(50,file='./data/saved_state_'//datechar//'_'//timechar//'.dat',action='write',status='new',form='unformatted',access='stream')!the state of the system is saved here
		else
			open(50,file='./data/saved_state.dat',action='write',form='unformatted',access='stream')!the state of the system is saved here
		end if
		write(50) amnttot!writing out the total amount of particles in the simulation
		if(samemass)then
			do i=1,amnttot
				write(50) rvec(1,i),rvec(2,i),rvec(3,i),vvec(1,i),vvec(2,i),vvec(3,i),fvec(1,i),fvec(2,i),fvec(3,i)!writing positions, velocities and forces when the masses are all equal
			end do
			call flush(50)!forcing the program to write stuff
		else
			do i=1,amnttot
				write(50) rvec(1,i),rvec(2,i),rvec(3,i),vvec(1,i),vvec(2,i),vvec(3,i),fvec(1,i),fvec(2,i),fvec(3,i),mvec(i)!writing positions, velocities and forces when the masses are different
			end do
			call flush(50)!forcing the program to write stuff
		end if
		close(50)!disconnecting the "saved_state.dat" file
		return
	end subroutine save_state



	subroutine load_state()!and to retrieve a saved state this is the subroutine needed
		logical::fileexists=.false.!this variable will be used to see if the file to load the system from exists
		integer(8)::i=0!particle counter
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
				read(50) rvec(1,i),rvec(2,i),rvec(3,i),vvec(1,i),vvec(2,i),vvec(3,i),fvec(1,i),fvec(2,i),fvec(3,i)!reading positions, velocities and forces when the masses are all equal
			end do
		else
			do i=1,amnttot
				read(50) rvec(1,i),rvec(2,i),rvec(3,i),vvec(1,i),vvec(2,i),vvec(3,i),fvec(1,i),fvec(2,i),fvec(3,i),mvec(i)!reading positions, velocities and forces when the masses are different
			end do
		end if
		close(50)!disconnecting the "saved_state.dat" file
		return
	end subroutine load_state



	subroutine xyz_output()!output for VMD
		integer(8)::i=0!particle counter
		do i=1,amnttot
			write(20,*) 'Ar',rvec(1,i),rvec(2,i),rvec(3,i)!writing positions
		end do
		call flush(20)!forcing the program to write stuff
		return
	end subroutine xyz_output



	real(8) function length_to_SI(x)!this functions converts distances in units of the simulation to metres
		real(8),intent(in)::x!the value to convert
		length_to_SI=x*1.d-9!the conversion
	end function length_to_SI



	real(8) function length_from_SI(x)!this functions converts distances in metres to units of the simulation
		real(8),intent(in)::x!the value to convert
		length_from_SI=x*1.d9!the conversion
	end function length_from_SI



	real(8) function area_to_SI(x)!this functions converts areas in units of the simulation to squared metres
		real(8),intent(in)::x!the value to convert
		area_to_SI=x*1.d-18!the conversion
	end function area_to_SI



	real(8) function area_from_SI(x)!this functions converts areas in squared metres to units of the simulation
		real(8),intent(in)::x!the value to convert
		area_from_SI=x*1.d18!the conversion
	end function area_from_SI



	real(8) function volume_to_SI(x)!this functions converts volumes in units of the simulation to cubic metres
		real(8),intent(in)::x!the value to convert
		volume_to_SI=x*1.d-27!the conversion
	end function volume_to_SI



	real(8) function volume_from_SI(x)!this functions converts volumes in cubic metres to units of the simulation
		real(8),intent(in)::x!the value to convert
		volume_from_SI=x*1.d27!the conversion
	end function volume_from_SI



	real(8) function time_to_SI(x)!this functions converts times in units of the simulation to seconds
		real(8),intent(in)::x!the value to convert
		time_to_SI=x*1.d-15!the conversion
	end function time_to_SI



	real(8) function time_from_SI(x)!this functions converts times in seconds to units of the simulation
		real(8),intent(in)::x!the value to convert
		time_from_SI=x*1.d15!the conversion
	end function time_from_SI



	real(8) function mass_to_SI(x)!this functions converts masses in units of the simulation to kilograms
		real(8),intent(in)::x!the value to convert
		mass_to_SI=x*1.d-27!the conversion
	end function mass_to_SI



	real(8) function mass_from_SI(x)!this functions converts masses in kilograms to units of the simulation
		real(8),intent(in)::x!the value to convert
		mass_from_SI=x*1.d27!the conversion
	end function mass_from_SI



	real(8) function energy_to_SI(x)!this functions converts energies in units of the simulation to joules
		real(8),intent(in)::x!the value to convert
		energy_to_SI=x*1.d-15!the conversion
	end function energy_to_SI



	real(8) function energy_from_SI(x)!this functions converts enegies in joules to units of the simulation
		real(8),intent(in)::x!the value to convert
		energy_from_SI=x*1.d15!the conversion
	end function energy_from_SI



	real(8) function temperature_to_SI(x)!this functions converts temperatures in units of the simulation to joules
		real(8),intent(in)::x!the value to convert
		temperature_to_SI=x*1.d0!the conversion
	end function temperature_to_SI



	real(8) function temperature_from_SI(x)!this functions converts temperatures in joules to units of the simulation
		real(8),intent(in)::x!the value to convert
		temperature_from_SI=x*1.d0!the conversion
	end function temperature_from_SI



	real(8) function pressure_to_SI(x)!this functions converts pressures in units of the simulation to pascals
		real(8),intent(in)::x!the value to convert
		pressure_to_SI=x*1.d0!the conversion
	end function pressure_to_SI



	real(8) function pressure_from_SI(x)!this functions converts pressure in pascals to units of the simulation
		real(8),intent(in)::x!the value to convert
		pressure_from_SI=x*1.d0!the conversion
	end function pressure_from_SI



end program Molecular_Dynamics_in_3D
