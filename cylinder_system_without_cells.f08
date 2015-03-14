module common_stuff
	implicit none
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
	real(8)::l2=0.d0!this is the L/2 factor of the interval for the distance bins
	real(8)::latcon=0.d0!lattice constant
	real(8)::mass=0.d0!generic mass of the particles
	real(8)::squaredispl=0.d0!mean square displacement used to calculate the diffusion constant
	real(8)::mvec(maxpartic)=0.d0!if the masses of the particles are not the same they can be stored here
	real(8)::measintv=0.d0!interval of time between measurements
	real(8)::potentialshift=0.d0!shift of the interaction potential
	real(8)::press=0.d0!pressure of the system
	real(8)::radius=0.d0!radius of the non solid zone of the system where particles can flow
	real(8)::rvec(3,maxpartic)=0.d0!positions for each particle (x->(1,:), y->(2,:), z->(3,:)) are here
	real(8)::rvecini(3,maxpartic)=0.d0!initial positions for each particle (x->(1,:), y->(2,:), z->(3,:)) are here
	real(8)::rvecnonper(3,maxpartic)=0.d0!positions for each particle (x->(1,:), y->(2,:), z->(3,:)) without periodic boundary conditions (to calculate mean square displacement) are here
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
	read(10,*) andersen,berendsen,cylinders,epsi,fcc,gnuplot,intrad,latcon,mass,meascylin,measintv,potentialshift,radius,samemass,sigm,simpcub,siz(1),siz(2),siz(3),temp,tempbath,tmeas,trelax,tstep,tmax
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
	integer(8)::vecindx(3)=0!here the indeces relevant to each dimensions are kept, for example the indeces of the cells in the crystalline lattice
	real(8)::dran_g!gaussian random number
	call date_and_time(values=seed)!creating a seed using the clock of the computer
	call dran_ini(0)!seed(8))!initializing the random number generator
	open(20,file='position_output.dat',action='write')!if the positions are to be visualized later they must be saved to this file
	open(30,file='energy_temperature_pressure_diffusion_output.dat',action='write')!this is where the time evolution of these magnitudes will be kept
	open(40,file='radial_distribution.dat',action='write')!and the radial distribution function values are here
	radius=2.5d-1*dmin1(siz(1),siz(2))!radius based on the size of the system
	!radius=2.d1*2.d0**(1.d0/6.d0)*sigm!radius based on 20 equlibrium distances
	volum=siz(1)*siz(2)*siz(3)!volume of the system
	threshold=tstep/trelax!this is the value that the threshold for the andersen thermostat should have
	l2=dmin1(siz(1),siz(2),siz(3))*5.d-1!the value of L/2 is chosen to be half of the smallest side of the simulation box
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
	do i1=0,nint(siz(1)/latcon,8)-1
		vecindx(1)=i1!index of the current cell in the first dimension
		do i2=0,nint(siz(2)/latcon,8)-1
			vecindx(2)=i2!index of the current cell in the second dimension
			do i3=0,nint(siz(3)/latcon,8)-1
				vecindx(3)=i3!index of the current cell in the third dimension
				if(simpcub.and.fcc)then
					write(*,*) 'Attempted to generate more than one structure, aborting.'
					stop!the program needs just one of these two varaibles to be true, it won't work properly with both at the same time
				else if(fcc)then!face centered cubic structure
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
					do k=1,3
						rvecini(k,j)=rvec(k,j)!storing the initial x position of the first particle
						rvecini(k,j+1)=rvec(k,j+1)!storing the initial x position of the second particle
						rvecini(k,j+2)=rvec(k,j+2)!storing the initial x position of the third particle
						rvecini(k,j+3)=rvec(k,j+3)!storing the initial x position of the fourth particle
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
					if(gnuplot)then
						do k=0,3
							write(*,*) rvec(1,j+k),rvec(2,j+k),rvec(3,j+k)!data feed to gnuplot
						end do
					end if
				else if(simpcub)then!simple cubic structure
					do k=1,3
						rvec(k,j)=latcon*dble(vecindx(k))!position in the solid lattice of the material, free particles will leave these as they move
					end do
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
					if(gnuplot)write(*,*) rvec(1,j),rvec(2,j),rvec(3,j)!data feed to gnuplot
				else
					write(*,*) 'No crystalline structure defined, it must be either simple cubic or face centered cubic, aborting.'
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
		amntfree=amnttot!just in case, to avoid problems since there are no fixed particles in this situation
	end if
	write(20,*) amnttot!the xyz format needs to have the amount of particles in the first line
	write(20,*) 'This is an optional comment, but it must be here for the file to be readable by the visualizing program.'!and this comment is mandatory
	write(30,*) '#	time	total energy	potential energy	kinetic energy	temperature	pressure	mean square displacement/6'!magnitudes in this file
	if(simpcub)then
		dens=amntfree/(siz(3)*dacos(-1.d0)*radius**2)!this is the number density of the free particles
	else if(fcc)then
		dens=amnttot/(siz(3)*dacos(-1.d0)*radius**2)!this is the number density of all the particles
	end if
	call remove_average_momentum()!the system should't have a global net movement so it is removed here if there's any
	call calculate_forces()!before the first integration iteration can be carried out, the forces need to be calculated
	return
contains
	subroutine remove_average_momentum()!to avoid the system having a drift, the average momentum must be substracted to all the particles in the system
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
end subroutine startup



subroutine kinematics_and_dynamics()!this is the integrator of the simulation, in this case the Verlet integrator
	use common_stuff
	implicit none
	integer(8)::i=0!first particle counter
	integer(8)::j=0!second particle counter
	integer(8)::k=0!dimension counter
	call half_leap(.false.)!updating the velocities of the particles without the thermostat, because it should only be applied to the speeds one per timestep
	call move()!updating the positions of the particles
	call calculate_forces()!the Verlet integrator requires to update the speeds again after moving the particle so the forces need to be calculated again
	call half_leap(.true.)!updating the velocities of the particles with the thermostat
	return
contains
	subroutine half_leap(thermostat)!the verlet integrator requires to update the velocities in two steps, the subroutine distinguishes when the particles have the same mass and when they don't
		implicit none
		logical,intent(in)::thermostat!this variable controls if the thermostat should be applied or not
		logical::andersenpart=.false.!if true, the particle will collide according to the thermostat rules
		real(8)::dran_g!gaussian random number
		real(8)::dran_u!uniform random number
		do i=1,amnttot
			if(fixed(i))cycle!fixed particles don't move
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
					vvec(m,i)=dsqrt(1.d0+(tempbath/temp-1.d0)*tstep/trelax)*(vvec(m,i)+fvec(m,i)*tstep*5.d-1/mass)!updating speeds with the berendsen thermostat
				end do
			else if((andersen.and.thermostat).and.andersenpart)then
				do m=1,3
					vvec(m,i)=dsqrt(kboltz*tempbath/mass)*dran_g()!updating speeds with the andersen thermostat
				end do
				else
					do m=1,3
						vvec(m,i)=vvec(m,i)+fvec(m,i)*tstep*5.d-1/mass!updating speeds without a thermostat
					end do
				end if
			else
				if(berendsen.and.andersen)then
					stop "Both thermostats can't be active at the the same time, aborting."
				else if(berendsen.and.thermostat)then
					do m=1,3
						vvec(m,i)=dsqrt(1.d0+(tempbath/temp-1.d0)*tstep/trelax)*(vvec(m,i)+fvec(m,i)*tstep*5.d-1/mvec(i))!updating speeds with the berendsen thermostat
					end do
				else if((andersen.and.thermostat).and.andersenpart)then
					do m=1,3
						vvec(m,i)=dsqrt(kboltz*tempbath/mvec(i))*dran_g()!updating speeds with the andersen thermostat
					end do
				else
					do m=1,3
						vvec(m,i)=vvec(m,i)+fvec(m,i)*tstep*5.d-1/mvec(i)!updating speeds without a thermostat
					end do
				end if
			end if
		end do
		return
	end subroutine half_leap
	subroutine move()!this is where the actual movement of the particles takes place
		implicit none
		if(gnuplot)write(*,*) "splot '-' w p pt 7 lc 0"!new gnuplot frame
		do i=1,numcells(1)
			if(fixed(i))cycle!fixed particles don't move
			do m=1,3
				rvec(m,i)=rvec(m,i)+vvec(m,i)*tstep!moving the particle
			end do
			do m=1,3
				if(rvec(m,i)<0.d0)then
					rvec(m,i)=rvec(m,i)+siz(m)!if the particle has reached the lower end of one dimension of the simulation box it should be at the other side of it
					rvecnonper(m,i)=rvecnonper(m,i)-siz(m)!tracking the movement along the periodic boundary conditions
				else if(rvec(m,i)>=siz(m))then
					rvec(m,i)=rvec(m,i)-siz(m)!if the particle has reached the higher end of one dimension of the simulation box it should be at the other side of it
					rvecnonper(m,i)=rvecnonper(m,i)+siz(m)!tracking the movement along the periodic boundary conditions
				end if
			end do
			if(gnuplot)write(*,*) rvec(1,i),rvec(2,i),rvec(3,i)!data feed to gnuplot
		end do
		if(gnuplot)write(*,*) 'e'!end of gnuplot frame
		return
	end subroutine move
end subroutine kinematics_and_dynamics




subroutine calculate_forces()
	use common_stuff
	implicit none
	logical::measure=.false.!if it is time to make a measurement this variable will be true
	integer(8)::i=0!first particle counter
	integer(8)::j=0!second particle counter
	integer(8)::k=0!dimension counter
	real(8)::dist=0.d0!distance between particles
	real(8)::distvec(3)=0.d0!vectorial distance between particles
	real(8)::force=0.d0!modulus of the force vector
	if((nint(dmod(time,measintv),8)==0).and.(time>=tmeas))then
		measure=.true.!to avoid computing the above many times, this logical variable is set
		epot=0.d0!setting the potential energy to zero so it can be calculated
		presscount=0!setting the counter to zero for the next measurement of the presure
	else
		measure=.false.!to make sure it stays false when it should
	end if
	do j=1,amnttot
		do k=1,3
			fvec(k,j)=0.d0!setting forces to zero
		end do
	end do
	do i=1,amnttot
		do j=i+1,amnttot
			if(i==j)cycle!a particle shouldn't interact with itself
			do k=1,3
				distvec(k)=rvec(k,i)-rvec(k,j)!these are the components of the vectorial distance between the particles
				if(disvec(k)>=cellsize(k)*5.d-1)then
					distvec(k)=rvec(k,i)-rvec(k,j)-cellsize(k)!if the particles are too far apart this is the actual component of the distance because of the periodic boundary conditions
				else if(disvec(k)>=-cellsize(k)*5.d-1)then
					distvec(k)=rvec(k,i)-rvec(k,j)+cellsize(k)!if the particles are too far apart this is the actual component of the distance because of the periodic boundary conditions
				end if
			end do
			dist=dsqrt(distvec(1)**2+distvec(2)**2+distvec(3)**2)!and this is the distance between the two particles
			if(fixed(j).or.fixed(i))then!the interaction can be slightly different between fixed and free particles
				force=2.4d1*epsi*(sigm/dist)**6*(2.d0*(sigm/dist)**6-1.d0)/dist**2!this is the module of the force between a free and a fixed particle
			else
				force=2.4d1*epsi*(sigm/dist)**6*(2.d0*(sigm/dist)**6-1.d0)/dist**2!and this is the module of the force between two free particles
			end if
			if((.not.fixed(i)).and.(.not.fixed(j)))then
				do k=1,3
					fvec(k,i)=fvec(k,i)+distvec(k)*force!and these are the components of the force vector for the first particle
					fvec(k,j)=fvec(k,j)-distvec(k)*force!the same for the second particle
				end do
			else if(fixed(j).and.(.not.fixed(i)))then
				do l=1,3
					fvec(k,i)=fvec(k,i)+distvec(k)*force!and these are the components of the force vector for the first particle
				end do
			else if(fixed(i).and.(.not.fixed(j)))then
				do l=1,3
					fvec(k,j)=fvec(k,j)-distvec(k)*force!and these are the components of the force vector for the first particle
				end do
			end if
			if(measure.and.cylinders)then
				!#####HERE MAYBE THERE SHOULD BE A CONTRIBUTION TO THE POTENTIAL ENERGY FOR PARTICLES LOCATED EXCLUSIVELY IN EACH CYLINDER
			else if(measure)then
				epot=epot+4.d0*epsi*(sigm/dist)**6*((sigm/dist)**6-1.d0)+potentialshift!this is the potential energy
			end if
			press=press-force*dist!to calculate the pressure, the sum over all interactions must be carried out
		end do
	end do
	presscount=presscount+1!this counter keeps track of how many timesteps are used to sample the pressure
	return
end subroutine calculate_forces




subroutine measuring()!here the different measurements are performed
	use common_stuff
	implicit none
	integer(8)::i=0!particle counter
	integer(8)::j=0!dimension counter
	real(8)::tmp=0.d0!temporary variable to calculate the distance travelled by the particles
	ekine=0.d0!setting the kinetic energy to zero to calculate the new values
	squaredispl=0.d0!same for the mean square displacement
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
	do i=1,amnttot
		do j=1,3
			tmp=tmp+(rvec(j,i)-rvecini(j,i)-rvecnonper(j,i))**2!calculating the distance travelled by the particles
		end do
		squaredispl=squaredispl+dsqrt(tmp)!computing the mean square displacement needed for the diffusion
		tmp=0.d0!resetting the temporary variable
	end do
	call radial_distribution()!getting a new sample of the radial distribution
	write(30,*) time,epot+ekine,epot,ekine,temp,press,squaredispl/(6.d0*amnttot)!writing data
	press=0.d0!discarding the old value of the pressure so the new one can be computed in the next steps
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
			bintag=nint(dist*dble(maxbin)/l2,8)
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
			cvcy(i)=(3.d0*kboltz)/((4.d0*amntfree*(sumekinecy(i)**2/dble(cvcountcy(i))-sumekinecysqur(i))/(dble(cvcountcy(i))*3.d0*kboltz*tempcy(i)**2))-2.d0)!computing the final value of the heat capacity
			write(30,*) '# The heat capacity in cylinder',i,'is:',cvcy(i)!writing data
		end do
	else
		do i=1,maxbin
			write(40,*) dble(i)*l2/dble(maxbin),dble(distbins(i))/(dble(radcount)*dens)!writing the data of the radial distribution
		end do
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