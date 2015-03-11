module common_stuff
	implicit none
	save
	integer(8),parameter::maxbin=20!maximum number of radial bins
	integer(8),parameter::maxcyl=10!maximum number of cylinders
	integer(8),parameter::maxpartic=10000!maximum number of particles in the system
	integer(8),parameter::maxcells=20!maximum number of cells per dimension
	logical::andersen=.false.!if true the andersen thermostat will be active
	logical::berendsen=.false.!if true the berendsen thermostat will be active
	logical::cylinders=.false.!if true, independent measurements will be taken in the different cylindrical intervals of the system
	logical::fcc=.false.!if true a face centered cubic structure will be generated
	logical::fixed(maxpartic)=.false.!this is true for fixed particles and false for free ones that can move
	logical::gnuplot=.false.!if true, the program will output data to be visualized in gnuplot
	logical::samemass=.false.!if true, all the particles have the same mass
	logical::simpcub=.false.!if true, a simple cubic structure will be generated
	integer(8)::amntfix=0!amount of fixed particles in the simulation
	integer(8)::amntfree=0!amount of free particles in the simulation
	integer(8)::amnttot=0!amount of total particles in the simulation
	integer(8)::distbins(maxbin)=0!these bins will count how many particles are in the divisions of the (0,L/2] interval in which the radial distribution will be measured
	integer(8)::cells(maxpartic,maxcells,maxcells,maxcells)=0!lists containing the tags of the particles in each cell, the first index are the particles, the second one is the x axis, the third one is y and the last one is z
	integer(8)::cvcount=0!this counter will keep track of how many measurements are used to calculate the heat capacity so the time average can be computed
	integer(8)::cvcountcy(maxcyl)=0!this counter will keep track of how many measurements are used to calculate the heat capacity for each cylinder so the time average can be computed
	integer(8)::difftag=0!this variable is the tag of the particle that will be followed to measure diffusion properties
	integer(8)::meascylin=0!this is the number of cylindrical intervals in which there will be independent measurements of variables
	integer(8)::numcells(3)=0!number of cells for each dimension (x=1, y=2, z=3)
	integer(8)::occup(maxcells,maxcells,maxcells)=0!number of particles in each cell, occupancy lists
	integer(8)::presscount=0!this counter keeps track of how many timesteps are used to obtain the pressure
	integer(8)::radcount=0!this counts how many timesteps are used in the measurement of the radial distribution
	real(8)::cellsize(3)=0.d0!size of the cells for each dimension (x=1, y=2, z=3)
	real(8)::cv=0.d0!heat capacity of the system
	real(8)::cvcy(maxcyl)=0.d0!heat capacity of each cylinder
	real(8)::dens=0.d0!number density of the system
	real(8)::diffconst=0.d0!diffusion constant
	real(8)::ekine=0.d0!kinetic energy of the system
	real(8)::ekinecy(maxcyl)=0.d0!kinetic energy of each cylinder in the system
	real(8)::epot=0.d0!potential energy of the system
	real(8)::epsi=0.d0!depth parameter of the Lennard-Jones interaction
	real(8)::fvec(3,maxpartic)=0.d0!forces for each particle (fx->(1,:), fy->(2,:), fz->(3,:)) are here
	real(8)::intrad=0.d0!interaction radius
	real(8)::kboltz=1.38d-8!Boltzmann's constant
	real(8)::l2=0.d0!this is the L factor of the interval for the distance bins
	real(8)::latcon=0.d0!lattice constant
	real(8)::mass=0.d0!generic mass of the particles
	real(8)::meansqrdispl=0.d0!mean square displacement used to calculate the diffusion constant
	real(8)::mvec(maxpartic)=0.d0!if the masses of the particles are not the same they can be stored here
	real(8)::measintv=0.d0!interval of time between measurements
	real(8)::press=0.d0!pressure of the system
	real(8)::radius=0.d0!radius of the non solid zone of the system where particles can flow
	real(8)::rvec(3,maxpartic)=0.d0!positions for each particle (x->(1,:), y->(2,:), z->(3,:)) are here
	real(8)::rvecini(3,maxpartic)=0.d0!initial positions for each particle (x->(1,:), y->(2,:), z->(3,:)) are here
	real(8)::sigm=0.d0!size parameter of the Lennard-Jones interaction
	real(8)::siz(3)=0.d0!dimensions (x=1, y=2, z=3) of the simulation box
	real(8)::sumekine=0.d0!the sum of the kinetic energy over time is kept here
	real(8)::sumekinesqur=0.d0!the sum of the kinetic energy squared over time is kept here
	real(8)::sumekinecy(maxcyl)=0.d0!the sum of the kinetic energy over time for each cylinder is kept here
	real(8)::sumekinecysqur(maxcyl)=0.d0!the sum of the kinetic energy squared over time for each cylinder is kept here
	real(8)::time=0.d0!time for which the simulation has been running
	real(8)::temp=0.d0!temperature of the system
	real(8)::tempbath=0.d0!temperature of the heat bath
	real(8)::tempcy(maxcyl)=0.d0!temperature of each cylinder in the system
	real(8)::threshold=0.d0!this variable is used to test if a particle will be affected by the andersen thermostat
	real(8)::tmeas=0.d0!when time is greater or equal than this measurements will start
	real(8)::trelax=0.d0!relaxation time for the andersen and berendsen thermostats
	real(8)::tstep=0.d0!timestep of the simulation
	real(8)::tmax=0.d0!time when the simulation should stop
	real(8)::volum=0.d0!volumen of the system
	real(8)::vvec(3,maxpartic)=0.d0!velocities for each particle (vx->(1,:), vy->(2,:), vz->(3,:)) are here
end module common_stuff



program MD_3D_cylinders
	use common_stuff
	implicit none
	call input()!reading initial parameters
	call startup()!setting up the system
	do
		if(time>tmax)exit!the program ends when it reaches the maximum time
		call kinematics_and_dynamics()!moving the system
		if((int(dmod(time,measintv),8)==0).and.(time>=tmeas))then
			if(cylinders)then
				call measuring_with_cylinders()!with a certain frequency measurements are performed in each cylindrical section of the fluid cavity of the system
			else
				call measuring()!with a certain frequency measurements are performed
			end if
		end if
		time=time+tstep!time advances
	end do
	call output()!finishing up the writing of results
contains
	function length_to_SI(x)!this functions converts distances in units of the simulation to metres
		implicit none
		real(8)::length_to_SI!output value
		real(8),intent(in)::x!the value to convert
		length_to_SI=x*1.d-9!the conversion
	end function length_to_SI
	function length_from_SI(x)!this functions converts distances from metres to units of the simulation
		implicit none
		real(8)::length_from_SI!output value
		real(8),intent(in)::x!the value to convert
		length_from_SI=x*1.d9!the conversion
	end function length_from_SI
	function time_to_SI(x)!this functions converts times in units of the simulation to seconds
		implicit none
		real(8)::time_to_SI!output value
		real(8),intent(in)::x!the value to convert
		time_to_SI=x*1.d-15!the conversion
	end function time_to_SI
	function time_from_SI(x)!this functions converts times from seconds to units of the simulation
		implicit none
		real(8)::time_from_SI!output value
		real(8),intent(in)::x!the value to convert
		time_from_SI=x*1.d15!the conversion
	end function time_from_SI
	function mass_to_SI(x)!this functions converts masses in units of the simulation to kilograms
		implicit none
		real(8)::mass_to_SI!output value
		real(8),intent(in)::x!the value to convert
		mass_to_SI=x*1.d-27!the conversion
	end function mass_to_SI
	function mass_from_SI(x)!this functions converts masses from kilograms to units of the simulation
		implicit none
		real(8)::mass_from_SI!output value
		real(8),intent(in)::x!the value to convert
		mass_from_SI=x*1.d27!the conversion
	end function mass_from_SI
	function energy_to_SI(x)!this functions converts energies in units of the simulation to joules
		implicit none
		real(8)::energy_to_SI!output value
		real(8),intent(in)::x!the value to convert
		energy_to_SI=x*1.d-15!the conversion
	end function energy_to_SI
	function energy_from_SI(x)!this functions converts enegies from joules to units of the simulation
		implicit none
		real(8)::energy_from_SI!output value
		real(8),intent(in)::x!the value to convert
		energy_from_SI=x*1.d15!the conversion
	end function energy_from_SI
end program MD_3D_cylinders



subroutine input()
	use common_stuff
	implicit none
	open(10,file='input.dat',action='read')
	read(10,*) andersen,berendsen,cylinders,epsi,fcc,gnuplot,intrad,latcon,mass,meascylin,measintv,radius,samemass,sigm,simpcub,siz(1),siz(2),siz(3),temp,tempbath,tmeas,trelax,tstep,tmax!reading input
	close(10)
	tmax=tmax+tstep*1.d-1!to avoid problems with repeated real summations, tmax is slightly modified so the last step is not lost
	return
end subroutine input



subroutine startup()!the system to be simulated is set up here
	use common_stuff
	implicit none
	integer(8)::i1=0!x dimension counter
	integer(8)::i2=0!y dimension counter
	integer(8)::i3=0!z dimension counter
	integer(8)::j=1!particle counter
	integer(8)::k=0!dimension counter
	integer(8)::seed(8)=0!seed for the random number generator
	real(8)::dran_g!gaussian random number
	call date_and_time(values=seed)!creating a seed using the clock of the computer
	call dran_ini(0)!seed(8))!initializing the random number generator
	open(20,file='position_output.dat',action='write')
	open(30,file='energy_temperature_pressure_output.dat',action='write')
	open(40,file='diffusion.dat',action='write')
	radius=2.5d-1*dmin1(siz(1),siz(2))!radius based on the size of the system
	!radius=2.d1*2.d0**(1.d0/6.d0)*sigm!radius based on 20 equlibrium distances
	volum=siz(1)*siz(2)*siz(3)!volume of the system
	threshold=tstep/trelax!this is the value that the threshold for the andersen thermostat should have
	l2=dmin1(siz(1),siz(2),siz(3))*5.d-1!the value of l2 is chosen to be half of the smallest side of the simulation box
	if(gnuplot)then
		write(*,*) 'unset key'!this just gets in the way
		write(*,*) 'set border 4095'!gnuplot initial configuration, showing the borders of the simulation box
		write(*,*) 'set ticslevel 0'!gnuplot initial configuration, setting the bottom of the simulation box at z=0
		write(*,*) 'set ztics mirror'!gnuplot initial configuration, setting the grid in the vertical axes
		write(*,*) 'set xrange [0:',siz(1),']'!gnuplot initial configuration, x size of the simulation box
		write(*,*) 'set yrange [0:',siz(2),']'!gnuplot initial configuration, y size of the simulation box
		write(*,*) 'set zrange [0:',siz(3),']'!gnuplot initial configuration, z size of the simulation box
		write(*,*) "splot '-' w p pt 7 lc 0"!new gnuplot frame
	end if
	do i1=0,nint(siz(1)/latcon,8)-1
		do i2=0,nint(siz(2)/latcon,8)-1
			do i3=0,nint(siz(3)/latcon,8)-1
				if(simpcub)then!simple cubic structure
					rvec(1,j)=latcon*dble(i1)!x position in the solid lattice of the material, free particles will leave these as they move
					rvec(2,j)=latcon*dble(i2)!y position in the solid lattice of the material, free particles will leave these as they move
					rvec(3,j)=latcon*dble(i3)!z position in the solid lattice of the material, free particles will leave these as they move
					if(gnuplot)write(*,*) rvec(1,j),rvec(2,j),rvec(3,j)!data feed to gnuplot
					if(dsqrt((rvec(1,j)-5.d-1*siz(1))**2+(rvec(2,j)-5.d-1*siz(1))**2)>=radius)then
						fixed(j)=.true.!the particles outside of this cylinder are fixed
						amntfix=amntfix+1!counting the number of fixed particles
					else
						do k=1,3
							vvec(3,j)=dsqrt(kboltz*temp/mass)*dran_g()!free particles have a speed
						end do
						amntfree=amntfree+1!counting the number of free particles
					end if
					if(.not.samemass)then
						mvec(j)=mass!setting the mass of the particles
					end if
					j=j+1!next particle
				else if(fcc)then!face centered structure
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
					rvecini(1,j)=rvec(1,j)!storing the initial x position of the first particle
					rvecini(2,j)=rvec(2,j)!storing the initial y position of the first particle
					rvecini(3,j)=rvec(3,j)!storing the initial z position of the first particle
					rvecini(1,j+1)=rvec(1,j+1)!storing the initial x position of the second particle
					rvecini(2,j+1)=rvec(2,j+1)!storing the initial y position of the second particle
					rvecini(3,j+1)=rvec(3,j+1)!storing the initial z position of the second particle
					rvecini(1,j+2)=rvec(1,j+2)!storing the initial x position of the third particle
					rvecini(2,j+2)=rvec(2,j+2)!storing the initial y position of the third particle
					rvecini(3,j+2)=rvec(3,j+2)!storing the initial z position of the third particle
					rvecini(1,j+3)=rvec(1,j+3)!storing the initial x position of the fourth particle
					rvecini(2,j+3)=rvec(2,j+3)!storing the initial y position of the fourth particle
					rvecini(3,j+3)=rvec(3,j+3)!storing the initial z position of the fourth particle
					if(gnuplot)then
						write(*,*) rvec(1,j),rvec(2,j),rvec(3,j)!data feed to gnuplot
						write(*,*) rvec(1,j+1),rvec(2,j+1),rvec(3,j+1)!data feed to gnuplot
						write(*,*) rvec(1,j+2),rvec(2,j+2),rvec(3,j+2)!data feed to gnuplot
						write(*,*) rvec(1,j+3),rvec(2,j+3),rvec(3,j+3)!data feed to gnuplot
					end if
					do k=1,3
						vvec(k,j)=dsqrt(kboltz*temp/mass)*dran_g()!velocity of the first particle
						vvec(k,j+1)=dsqrt(kboltz*temp/mass)*dran_g()!velocity of the second particle
						vvec(k,j+2)=dsqrt(kboltz*temp/mass)*dran_g()!velocity of the third particle
						vvec(k,j+3)=dsqrt(kboltz*temp/mass)*dran_g()!velocity of the fourth particle
					end do
					if(.not.samemass)then
						mvec(j)=mass!setting the mass of the first particle
						mvec(j+1)=mass!setting the mass of the first particle
						mvec(j+2)=mass!setting the mass of the first particle
						mvec(j+3)=mass!setting the mass of the first particle
					end if
					j=j+4!next particles
				else if(simpcub.and.fcc)then
					write(*,*) 'Attempted to generate more than one structure, aborting.'
					stop!the program needs one of these two varaibles to be true
				else
					write(*,*) 'No crystalline structure defined, it must be either simple cubic or face centered cubic, stopping the program now, aborting.'
					stop!the program won't work without a defined structure
				end if
			end do
		end do
	end do
	if(gnuplot)write(*,*) 'e'!end of gnuplot frame
	if(simpcub)then
		amnttot=j-1!total amount of particles for the sc structure
	else if(fcc)then
		amnttot=j-1!total amount of particles for the fcc structure
		amntfree=amnttot!just in case to avoid problems since there are no fixed particles in this case
	end if
	call remove_average_momentum()
	write(20,*) amnttot!the xyz format needs to have the amount of particles in the first line
	write(20,*) 'This is an optional comment, but it must be here for the file to be readable by the visualizing program.'!and this comment is mandatory
	write(30,*) '#	time	total energy	potential energy	kinetic energy	temperature	pressure'!magnitudes in this file
	if(simpcub)then
		dens=amntfree/(siz(3)*dacos(-1.d0)*radius**2)!this is the number density of the free particles
	else if(fcc)then
		dens=amnttot/(siz(3)*dacos(-1.d0)*radius**2)!this is the number density of all the particles
	end if
	call create_cell_lists()!cells are needed to calculate forces
	call calculate_forces()!before the first integration iteration can be carried out, the forces need to be calculated
	return
end subroutine startup



subroutine remove_average_momentum()!to avoid the system having a drift, the average momentum must be substracted to all the particles in the system
	use common_stuff
	implicit none
	integer(8)::i=0!particle counter
	integer(8)::j=0!dimension counter
	real(8)::mmntsum(3)=0.d0!the global momentum will be stored here
	if(samemass)then
		do i=1,amnttot
			do j=1,3
				mmntsum(j)=mmntsum(j)+vvec(j,i)*mass!calculating the global momentum if the particles have all the same mass
			end do
		end do
		mmntsum=mmntsum/dble(amnttot)!dividing by the number of particles so that the final result is the average momentum of the system
		do i=1,amnttot
			do j=1,3
				vvec(j,i)=(vvec(j,i)*mass-mmntsum(j))/mass!substracting the average momentum to the momentum of each particle and converting it back to velocity
			end do
		end do
	else
		do i=1,amnttot
			do j=1,3
				mmntsum(j)=mmntsum(j)+vvec(j,i)*mvec(i)!calculating the global momentum if the mass of each particle is different
			end do
		end do
		mmntsum=mmntsum/dble(amnttot)!dividing by the number of particles so that the final result is the average momentum of the system
		do i=1,amnttot
			do j=1,3
				vvec(j,i)=(vvec(j,i)*mvec(i)-mmntsum(j))/mvec(i)!substracting the average momentum to the momentum of each particle and converting it back to velocity
			end do
		end do
	end if
	return
end subroutine remove_average_momentum



subroutine create_cell_lists()!cell lists are used to speed up the simulation, see the article in wikipedia for more general details
	use common_stuff
	implicit none
	integer(8)::i=0!dimension counter
	integer(8)::indeces(3)=0!indeces of the cell in which a particle is
	integer(8)::j=0!particle counter
	do i=1,3
		numcells(i)=floor(siz(i)/intrad,8)!computing how many cells each dimension has
		cellsize(i)=siz(i)/dble(numcells(i))!calculating how big the cells are
	end do
	do j=1,amnttot
		do i=1,3
			indeces(i)=floor(rvec(i,j)/cellsize(i),8)+1!particle j is the cell which has these indeces
		end do
		occup(indeces(1),indeces(2),indeces(3))=occup(indeces(1),indeces(2),indeces(3))+1!updating the occupancy list with the new particle
		cells(occup(indeces(1),indeces(2),indeces(3)),indeces(1),indeces(2),indeces(3))=j!adding the new particle to the apropiate cell
	end do
	return
end subroutine create_cell_lists!			#####BE CAREFUL, IF THE SIMULATION BOX IS NOT COMPLETELY ON THE NON NEGATIVE PARTS OF THE AXES SEGMENTATION FAULTS MAY APPEAR



subroutine kinematics_and_dynamics()!this is the integrator of the simulation, in this case the Verlet integrator
	use common_stuff
	implicit none
	call half_leap(.false.)!updating the velocities of the particles without the thermostat, because it should only be applied to the speeds one per timestep
	call move()!updating the positions of the particles
	call update_particle_in_cells()!after the particles have moved, they may have entered a different cell and this needs to be checked and taken care of
	call calculate_forces()!the Verlet integrator requires to update the speeds again after moving the particle so the forces need to be calculated again
	call half_leap(.true.)!updating the velocities of the particles with the thermostat
	return
end subroutine kinematics_and_dynamics



subroutine half_leap(thermostat)!the verlet integrator requires to update the velocities in two steps, the subroutine distinguishes when the particles have the same mass and when they don't
	use common_stuff
	implicit none
	logical,intent(in)::thermostat!this variable controls if the thermostat should be applied or not
	logical::andersenpart=.false.!if true, the particle will collide according to the thermostat rules
	integer(8)::i=0!x dimension cell counter
	integer(8)::j=0!y dimension cell counter
	integer(8)::k=0!z dimension cell counter
	integer(8)::m=0!dimension counter
	integer(8)::n=0!particle counter in the occupancy lists
	integer(8)::tag=0!particle tag
	real(8)::dran_g!gaussian random number
	real(8)::dran_u!uniform random number
	do k=1,numcells(3)
		do j=1,numcells(2)
			do i=1,numcells(1)
				do n=1,occup(i,j,k)
					tag=cells(n,i,j,k)!this value is going to be called from memory many times so it will be kept here for convenience
					if(fixed(tag))cycle!fixed particles don't move
					if((andersen.and.thermostat).and.(dran_u()<threshold))then
						andersenpart=.true.!this particle will suffer a collision in this case
					else
						andersenpart=.false.!and in this case there is no collision
					end if
					if(samemass)then
						if(berendsen.and.andersen)then
							stop "Both thermostats can't be active at the the same time, aborting."
						else if(berendsen.and.thermostat)then
							do m=1,3
								vvec(m,tag)=dsqrt(1.d0+(tempbath/temp-1.d0)*tstep/trelax)*(vvec(m,tag)+fvec(m,tag)*tstep*5.d-1/mass)!updating speeds with the berendsen thermostat
							end do
						else if((andersen.and.thermostat).and.andersenpart)then
							do m=1,3
								vvec(m,tag)=dsqrt(kboltz*tempbath/mass)*dran_g()!updating speeds with the andersen thermostat
							end do
						else
							do m=1,3
								vvec(m,tag)=vvec(m,tag)+fvec(m,tag)*tstep*5.d-1/mass!updating speeds without a thermostat
							end do
						end if
					else
						if(berendsen.and.andersen)then
							stop "Both thermostats can't be active at the the same time, aborting."
						else if(berendsen.and.thermostat)then
							do m=1,3
								vvec(m,tag)=dsqrt(1.d0+(tempbath/temp-1.d0)*tstep/trelax)*(vvec(m,tag)+fvec(m,tag)*tstep*5.d-1/mvec(tag))!updating speeds with the berendsen thermostat
							end do
						else if((andersen.and.thermostat).and.andersenpart)then
							do m=1,3
								vvec(m,tag)=dsqrt(kboltz*tempbath/mvec(tag))*dran_g()!updating speeds with the andersen thermostat
							end do
						else
							do m=1,3
								vvec(m,tag)=vvec(m,tag)+fvec(m,tag)*tstep*5.d-1/mvec(tag)!updating speeds without a thermostat
							end do
						end if
					end if
				end do
			end do
		end do
	end do
	return
end subroutine half_leap



subroutine move()!this is where the actual movement of the particles takes place
	use common_stuff
	implicit none
	integer(8)::i=0!x dimension cell counter
	integer(8)::j=0!y dimension cell counter
	integer(8)::k=0!z dimension cell counter
	integer(8)::m=0!dimension counter
	integer(8)::n=0!particle counter in the occupancy lists
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
					end do
					do m=1,3
						if(rvec(m,tag)<0.d0)then
							rvec(m,tag)=rvec(m,tag)+siz(m)!if the particle has reached the lower end of one dimension of the simulation box it should be at the other side of it (periodic boundary 							conditions)
						else if(rvec(m,tag)>=siz(m))then
							rvec(m,tag)=rvec(m,tag)-siz(m)!if the particle has reached the higher end of one dimension of the simulation box it should be at the other side of it (periodic boundary 								conditions)
						end if
					end do
					if(gnuplot)write(*,*) rvec(1,tag),rvec(2,tag),rvec(3,tag)!data feed to gnuplot
				end do
			end do
		end do
	end do
	if(gnuplot)write(*,*) 'e'!end of gnuplot frame
	return
end subroutine move



subroutine update_particle_in_cells()
	use common_stuff
	implicit none
	integer(8)::i1=0!x dimension, old-cell counter
	integer(8)::j1=0!y dimension, old-cell counter
	integer(8)::k1=0!z dimension, old-cell counter
	integer(8)::i2=0!x dimension, new-cell counter
	integer(8)::j2=0!y dimension, new-cell counter
	integer(8)::k2=0!z dimension, new-cell counter
	integer(8)::n=0!particle counter in the cell lists
	integer(8)::tag=0!particle tag
	do k1=1,numcells(3)
		do j1=1,numcells(2)
			do i1=1,numcells(1)
				n=1!starting particle counter at the beginning of the current list
				do
					tag=cells(n,i1,j1,k1)!this value is going to be called from memory many times so it will be kept here for convenience
					if(fixed(tag))cycle!fixed particles don't need updating
					i2=floor(rvec(1,tag)/cellsize(1),8)+1!computing the x index of the cell in which the particle is located after moving
					j2=floor(rvec(2,tag)/cellsize(2),8)+1!computing the y index of the cell in which the particle is located after moving
					k2=floor(rvec(3,tag)/cellsize(3),8)+1!computing the z index of the cell in which the particle is located after moving
					if((i1/=i2).or.(j1/=j2).or.(k1/=k2))then!after the particle has moved, the cell in which it is contained may have changed so it needs to be updated
						occup(i2,j2,k2)=occup(i2,j2,k2)+1!updating the new occupancy list with the particle
						cells(occup(i2,j2,k2),i2,j2,k2)=tag!adding the particle to the corresponding new cell
						cells(n,i1,j1,k1)=cells(occup(i1,j1,k1),i1,j1,k1)!taking away the updated particle from the old cell, its place is occupied by the last one in that cell
						cells(occup(i1,j1,k1),i1,j1,k1)=0!erasing old information
						occup(i1,j1,k1)=occup(i1,j1,k1)-1!updating the old occupancy list
					end if
					n=n+1!next particle
					if(n>occup(i1,j1,k1))exit
				end do
			end do
		end do
	end do
	return
end subroutine update_particle_in_cells



subroutine calculate_forces()
	use common_stuff
	implicit none
	logical::interaction(maxpartic,maxpartic)=.false.!this matrix stores which interactions have been calculated already
	logical::measure=.false.!if it is time to make a measurement this variable will be true
	integer(8)::i1=0!x dimension cell counter for the first particle
	integer(8)::i2=0!x dimension cell counter for the second particle
	integer(8)::i3=0!x dimension auxiliary cell counter
	integer(8)::j1=0!y dimension cell counter for the first particle
	integer(8)::j2=0!y dimension cell counter for the second particle
	integer(8)::j3=0!y dimension auxiliary cell counter
	integer(8)::k1=0!z dimension cell counter for the first particle
	integer(8)::k2=0!z dimension cell counter for the second particle
	integer(8)::k3=0!z dimension auxiliary cell counter
	integer(8)::l=0!dimension counter
	integer(8)::n1=0!particle counter in the cells, for the first one
	integer(8)::n2=0!particle counter in the cells, for the second one
	integer(8)::tag1=0!tag for the first particle
	integer(8)::tag2=0!tag for the second particle
	real(8)::displ(3)=0.d0!displacement vector to calculate the appropiate distance with the periodic boundary conditions
	real(8)::dist=0.d0!distance between particles
	real(8)::distvec(3)=0.d0!vectorial distance between particles
	real(8)::force=0.d0!modulus of the force vector
	if((nint(dmod(time,measintv),8)==0).and.(time>=tmeas))then
		measure=.true.!to avoid computing the above many times, this logical variable is set
		epot=0.d0!setting the potential energy to zero so it can be calculated
		presscount=0!setting the counter to zero for the next measurement of the presure
	end if
	do n1=1,amnttot
		do n2=1,3
			fvec(n2,n1)=0.d0!setting forces to zero
		end do
		do n2=1,amnttot
			interaction(n2,n1)=.false.!the interactions must be reset so they are taken into account
		end do
	end do
	do k1=1,numcells(3)
		do j1=1,numcells(2)
			do i1=1,numcells(1)
				do n1=1,occup(i1,j1,k1)
				tag1=cells(n1,i1,j1,k1)!this value is going to be called from memory many times so it will be kept here for convenience
					do k3=k1-1,k1+1
						if(k3==0)then
							displ(3)=-siz(3)!if the neighbouring cell is at the end of the simulation box in the x direction the displacement vector has this value
							k2=numcells(3)!and this is the actual index of that cell
						else if(k3==numcells(3)+1)then
							displ(3)=siz(3)!if the neighbouring cell is at the beginning of the simulation box in the x direction the displacement vector has this value
							k2=1!and this is the actual index of that cell
						else
							displ(3)=0.d0!if the case is not one the former, the displacement vector is null
							k2=k3!and the actual index of the cell doesn't change
						end if!there are 26 neighbouring cells, 3 per dimension, there must be a cycle to run through them and compute the appropiate displacement vector for each different case (the 							central cell is touching the edges of the simulation box or not)
						do j3=j1-1,j1+1
							if(j3==0)then
								displ(2)=-siz(2)!if the neighbouring cell is at the end of the simulation box in the y direction the displacement vector has this value
								j2=numcells(2)!and this is the actual index of that cell
							else if(j3==numcells(2)+1)then
								displ(2)=siz(2)!if the neighbouring cell is at the beginning of the simulation box in the y direction the displacement vector has this value
								j2=1!and this is the actual index of that cell
							else
								displ(2)=0.d0!if the case is not one the former, the displacement vector is null
								j2=j3!and the actual index of the cell doesn't change
							end if!there are 26 neighbouring cells, 3 per dimension, there must be a cycle to run through them and compute the appropiate displacement vector for each different case 								(the central cell is touching the edges of the simulation box or not)
							do i3=i1-1,i1+1
								if(i3==0)then
									displ(1)=-siz(1)!if the neighbouring cell is at the end of the simulation box in the z direction the displacement vector has this value
									i2=numcells(1)!and this is the actual index of that cell
								else if(i3==numcells(1)+1)then
									displ(1)=siz(1)!if the neighbouring cell is at the beginning of the simulation box in the z direction the displacement vector has this value
									i2=1!and this is the actual index of that cell
								else
									displ(1)=0.d0!if the case is not one the former, the displacement vector is null
									i2=i3!and the actual index of the cell doesn't change
								end if!there are 26 neighbouring cells, 3 per dimension, there must be a cycle to run through them and compute the appropiate displacement vector for each different 									case (the central cell is touching the edges of the simulation box or not)
								do n2=1,occup(i2,j2,k2)
									tag2=cells(n2,i2,j2,k2)!this value is going to be called from memory many times so it will be kept here for convenience
									if(tag1==tag2)cycle!a particle shouldn't interact with itself
									if((interaction(tag1,tag2)).or.(interaction(tag2,tag1)))cycle!if this interaction has already been accounted for there's no need to go on
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
									if(measure.and.cylinders)then
										! 			#####HERE MAYBE THERE SHOULD BE A CONTRIBUTION TO THE POTENTIAL ENERGY FOR PARTICLES LOCATED EXCLUSIVELY IN EACH CYLINDER
									else if(measure)then
										epot=epot+4.d0*epsi*(sigm/dist)**6*((sigm/dist)**6-1.d0)+9.06d-9!this is the potential energy
									end if
									press=press-force*dist!to calculate the pressure, the sum over all interactions must be carried out
									!press=press+(fvec(1,tag1)*distvec(1)+fvec(2,tag1)*distvec(2)+fvec(3,tag1)*distvec(3))!to calculate the pressure, the force per particle is needed 										#####THIS IS THE ACTUAL SCALAR PRODUCT, IF THE PREVIOUS STATEMENT CAUSES PROBLEMS OR DOESN'T WORK USE THIS ONE
								end do
							end do
						end do
					end do
				end do
			end do
		end do
	end do
	presscount=presscount+1!this counter keeps track of how many times the forces are sampled to calculate the pressure
	return
end subroutine calculate_forces



subroutine periodic_interaction(dimen,inindex,outindex,displ)!in the case that some cell is beyond the limits, this will apply the periodic boundary conditions to the interaction between the particles
	use common_stuff
	implicit none
	integer(8),intent(in)::dimen!dimension in which the subroutine is carried out, 1 is x-axis, 2 is y-axis and 3 is z-axis
	integer(8),intent(in)::inindex!this number runs through the range of possible values for neighbouring cells
	integer(8),intent(out)::outindex!and this is the resulting index of that cell
	real(8),intent(out)::displ!this is the value of the displacement vector in that dimension
	if(inindex==0)then
		displ=-siz(dimen)!if the neighbouring cell is at the end of the simulation box in the x direction the displacement vector has this value
		outindex=numcells(dimen)!and this is the actual index of that cell
	else if(inindex==numcells(dimen)+1)then
		displ=siz(dimen)!if the neighbouring cell is at the beginning of the simulation box in the x direction the displacement vector has this value
		outindex=1!and this is the actual index of that cell
	else
		displ=0.d0!if the case is not one the former, the displacement vector is null
		outindex=inindex!and the actual index of the cell doesn't change
	end if
	return
end subroutine periodic_interaction



subroutine interacting(measure,tag1,tag2,displ)!this is where the actual interaction happens
	use common_stuff
	implicit none
	logical,intent(in)::measure!if it's time to measure this will be true
	integer(8)::l=0!dimension counter
	integer(8),intent(in)::tag1!tag for the first particle
	integer(8),intent(in)::tag2!tag for the second particle
	real(8),intent(in)::displ(3)!displacement vector to calculate the appropiate distance with the periodic boundary conditions
	real(8)::dist=0.d0!distance between particles
	real(8)::distvec(3)=0.d0!vectorial distance between particles
	real(8)::force=0.d0!modulus of the force vector
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
	if(measure.and.cylinders)then
		! 			#####HERE MAYBE THERE SHOULD BE A CONTRIBUTION TO THE POTENTIAL ENERGY FOR PARTICLES LOCATED EXCLUSIVELY IN EACH CYLINDER
	else if(measure)then
		epot=epot+4.d0*epsi*(sigm/dist)**6*((sigm/dist)**6-1.d0)+9.06d-9!this is the potential energy
	end if
	press=press-force*dist!to calculate the pressure, the sum over all interactions must be carried out
	!press=press+(fvec(1,tag1)*distvec(1)+fvec(2,tag1)*distvec(2)+fvec(3,tag1)*distvec(3))!to calculate the pressure, the force per particle is needed 			#####THIS IS THE ACTUAL SCALAR PRODUCT, 	IF THE PREVIOUS STATEMENT CAUSES PROBLEMS OR DOESN'T WORK USE THIS ONE
	return
end subroutine interacting



subroutine measuring()!here the different measurements are performed
	use common_stuff
	implicit none
	integer(8)::i=0!particle counter
	ekine=0.d0!setting the kinetic energy to zero to calculate the new values
	meansqrdispl=0.d0!same for the mean square displacement
	if(samemass)then
		do i=1,amnttot
			if(.not.fixed(i))ekine=ekine+vvec(1,i)**2+vvec(2,i)**2+vvec(3,i)**2!measuring the kinetic energy if all the particles have the same mass
		end do
		ekine=ekine*mass*5.d-1!don't forget the mass and the 1/2
	else
		do i=1,amnttot
			if(.not.fixed(i))ekine=ekine+mvec(i)*(vvec(1,i)**2+vvec(2,i)**2+vvec(3,i)**2)!measuring the kinetic energy if the particles have different masses
		end do
		ekine=ekine*5.d-1!don't forget the 1/2
	end if
	sumekine=sumekine+ekine!summing the kinetic energy
	sumekinesqur=sumekinesqur+ekine**2!and the squared kinetic energy
	cvcount=cvcount+1!another measurement for the heat capacity
	temp=(2.d0*ekine)/(3.d0*amntfree*kboltz)!and with the kinetic energy the temperature is measured too
	press=dens*temp*kboltz+press/(3.d0*volum*dble(presscount))!final result for the pressure
	write(30,*) time,epot+ekine,epot,ekine,temp,press!writing data
	press=0.d0!discarding the old value of the pressure so the new one can be computed in the next steps
	do i=1,amnttot
		meansqrdispl=meansqrdispl+dsqrt((rvec(1,i)-revecini(1,i))**2+(rvec(2,i)-revecini(2,i))**2+(rvec(3,i)-revecini(3,i))**2)
	end do
	write(40,*) time,meansqrdispl/(6.d0*amnttot)!writing diffusion data
	return
end subroutine measuring



subroutine measuring_with_cylinders()!here the different measurements are performed in each cylinder, whether all the particles have the same mass or not
	use common_stuff
	implicit none
	integer(8)::i=0!particle counter
	integer(8)::j=0!cylinder counter
	do j=1,maxcyl
		ekinecy(j)=0.d0!setting the kinetic energy of the cylinders to zero to calculate the new values
	end do
	if(samemass)then
		do i=1,amnttot
			if(.not.fixed(i))then
				j=floor(1.d1*dsqrt(rvec(1,i)**2+rvec(2,i)**2)/radius,8)!finding in which cylinder the particle is located
				if(j>10)j=10!just in case the particle is slightly beyond the radius of the cylinder
				ekinecy(j)=ekinecy(j)+mass*(rvec(1,i)**2+rvec(2,i)**2+rvec(3,i)**2)!measuring the kinetic energy per cylinder
			end if
		end do
		do j=1,maxcyl
			ekinecy(j)=ekinecy(j)*5.d-1!don't forget the 1/2
			sumekinecy(j)=sumekinecy(j)+ekinecy(j)!summing the kinetic energy per cylinder
			sumekinecysqur(j)=sumekinecysqur(j)+ekinecy(j)**2!and the squared kinetic energy per cylinder
			cvcountcy(j)=cvcountcy(j)+1!another measurement for the heat capacity per cylinder
			tempcy(j)=(2.d0*ekinecy(j))/(3.d0*amntfree*kboltz)!and with the kinetic energy the temperature is measured too per cylinder
		end do
	else
		do i=1,amnttot
			if(.not.fixed(i))then
				j=floor(1.d1*dsqrt(rvec(1,i)**2+rvec(2,i)**2)/radius,8)!finding in which cylinder the particle is located
				if(j>10)j=10!just in case the particle is slightly beyond the radius of the cylinder
				ekinecy(j)=ekinecy(j)+mvec(i)*(rvec(1,i)**2+rvec(2,i)**2+rvec(3,i)**2)!measuring the kinetic energy per cylinder
			end if
		end do
		do j=1,maxcyl
			ekinecy(j)=ekinecy(j)*5.d-1!don't forget the 1/2
			sumekinecy(j)=sumekinecy(j)+ekinecy(j)!summing the kinetic energy per cylinder
			sumekinecysqur(j)=sumekinecysqur(j)+ekinecy(j)**2!and the squared kinetic energy per cylinder
			cvcountcy(j)=cvcountcy(j)+1!another measurement for the heat capacity per cylinder
			tempcy(j)=(2.d0*ekinecy(j))/(3.d0*amntfree*kboltz)!and with the kinetic energy the temperature is measured too per cylinder
		end do
	end if
	call diffusion_with_cylinders()
	return
end subroutine measuring_with_cylinders



subroutine diffusion_with_cylinders()
	use common_stuff
	implicit none
	real(8)::dran_u!uniform random number
	
	do
		difftag=nint(dble(amntfree)*dran_u(),8)!choosing randomly the particle that will be used to measure diffusion properties
		if(difftag>0)exit!just in case the tag is 0 or negative
	end do
	return
end subroutine diffusion_with_cylinders



subroutine radial_distribution()
	use common_stuff
	implicit none
	logical::firststep=.true.!after the first call to this subroutine is completed, this will be false
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
		firststep=.false.!the loop directly above is meant to be run only once
	end if
	do i=1,amnttot
		dist=dsqrt((rvec(1,tag)-rvec(1,i))**2+(rvec(2,tag)-rvec(2,i))**2+(rvec(3,tag)-rvec(3,i))**2)!distance between particles
		if(dist<=l2)then
			bintag=nint(dist*2.d1/l2,8)
			if(bintag<1)bintag=1!just in case the particles are too close
			distbins(bintag)=distbins(bintag)+1!counting another particle in the appropiate bin
		end if
	end do
	radcount=radcount+1!one more sample for the measurement
	return
end subroutine radial_distribution



subroutine output()!after the simulation, some more stuff needs to be computed and the files need to be closed
	use common_stuff
	implicit none
	integer(8)::i=0!cylinder and bin counter
	if(cylinders)then
		do i=1,maxcyl
			cvcy(i)=(3.d0*kboltz)/((4.d0*amntfree*(sumekinecy(i)**2/dble(cvcountcy(i))-sumekinecysqur(i))/(dble(cvcountcy(i))*3.d0*kboltz*tempcy(i)**2))-2.d0)!computing the final value of the heat 				capacity
			write(30,*) '# The heat capacity in cylinder',i,'is:',cvcy(i)!writing data
		end do
	else
		cv=(3.d0*kboltz)/((4.d0*amntfree*(sumekine**2/dble(cvcount)-sumekinesqur)/(dble(cvcount)*3.d0*kboltz*temp**2))-2.d0)!computing the final value of the heat capacity
		write(30,*) '# The heat capacity is:',cv!writing data
	end if
	close(20)
	close(30)
	close(40)
	return
end subroutine output



subroutine save_state()!if the current state is to be preserved it can be saved to a file with this subroutine
	use common_stuff
	implicit none
	integer(8)::i=0!particle counter
	open(50,file='saved_state.dat',action='write')
	write(50,*) amnttot
	if(samemass)then
		do i=1,amnttot
			write(50,*) rvec(1,i),rvec(2,i),rvec(3,i),vvec(1,i),vvec(2,i),vvec(3,i),fvec(1,i),fvec(2,i),fvec(3,i)!reading positions, velocities and forces when the masses are all equal
		end do
	else
		do i=1,amnttot
			write(50,*) rvec(1,i),rvec(2,i),rvec(3,i),vvec(1,i),vvec(2,i),vvec(3,i),fvec(1,i),fvec(2,i),fvec(3,i),mvec(i)!reading positions, velocities and forces when the masses are different
		end do
	end if
	close(50)
	return
end subroutine save_state



subroutine load_state()!and to retrieve a saved state this is the subroutine needed
	use common_stuff
	implicit none
	integer(8)::i=0!particle counter
	open(50,file='saved_state.dat',action='read')
	read(50,*) amnttot
	if(samemass)then
		do i=1,amnttot
			read(50,*) rvec(1,i),rvec(2,i),rvec(3,i),vvec(1,i),vvec(2,i),vvec(3,i),fvec(1,i),fvec(2,i),fvec(3,i)
		end do
	else
		do i=1,amnttot
			read(50,*) rvec(1,i),rvec(2,i),rvec(3,i),vvec(1,i),vvec(2,i),vvec(3,i),fvec(1,i),fvec(2,i),fvec(3,i),mvec(i)
		end do
	end if
	close(50)
	return
end subroutine load_state
