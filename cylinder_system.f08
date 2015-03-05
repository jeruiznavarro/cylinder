module common_stuff
	implicit none
	save
	real(8)::cellsize(3)=0.d0!size of the cells for each dimension
	real(8)::cv=0.d0!heat capacity of the system
	real(8)::cvcy(10)=0.d0!heat capacity of each cylinder
	real(8)::dens=0.d0!number density of the system
	real(8)::diffconst=0.d0!diffusion constant
	real(8)::ekine=0.d0!kinetic energy of the system
	real(8)::ekinecy(10)=0.d0!kinetic energy of each cylinder in the system
	real(8)::epot=0.d0!potential energy of the system
	real(8)::epsi=0.d0!depth parameter of the Lennard-Jones interaction
	real(8)::kboltz=1.38d-8!Boltzmann's constant
	real(8)::latcon=0.d0!lattice constant
	real(8)::mass=0.d0!generic mass of the particles
	real(8)::measintv=0.d0!interval of time between measurements
	real(8)::press=0.d0!pressure of the system
	real(8)::r(10,10000)=0.d0!positions (1-3,i), speeds (4-6,i), forces (7-9,i) and masses (10,i) for each particle are in this array
	real(8)::intrad=0.d0!interaction radius
	real(8)::radius=0.d0!radius of the non solid zone of the system where particles can flow
	real(8)::sigm=0.d0!size parameter of the Lennard-Jones interaction
	real(8)::siz(3)=0.d0!dimensions of the simulation box
	real(8)::sumekine(2)=0.d0!the sums of the kinetic energy (1) and the kinetic energy squared (2) over time are kept here
	real(8)::sumekinecy(2,10)=0.d0!the sums of the kinetic energy (1) and the kinetic energy squared (2) over time for each cylinder are kept here
	real(8)::time=0.d0!time for which the simulation has been running
	real(8)::temp=0.d0!temperature of the system
	real(8)::tempbath=0.d0!temperature of the heat bath
	real(8)::tempcy(10)=0.d0!temperature of each cylinder in the system
	real(8)::tmeas=0.d0!when time is greater or equal than this measurements will start
	real(8)::trelax=0.d0!relaxation time for the berendsen thermostat
	real(8)::tstep=0.d0!timestep of the simulation
	real(8)::tmax=0.d0!time when the simulation should stop
	real(8)::volum=0.d0!volumen of the system
	integer(8)::amntfix=0!amount of fixed particles in the simulation
	integer(8)::amntfree=0!amount of free particles in the simulation
	integer(8)::amnttot=0!amount of total particles in the simulation
	integer(8)::cells(10000,20,20,20)=0!lists containing the tags of the particles in each cell, the first index are the particles, the second one is the x axis, the third one is y and the last one is z
	integer(8)::cvcount=0!this counter will keep track of how many measurements are used to calculate the heat capacity so the time average can be computed
	integer(8)::cvcountcy(10)=0!this counter will keep track of how many measurements are used to calculate the heat capacity for each cylinder so the time average can be computed
	integer(8)::difftag=0!this variable is the tag of the particle that will be followed to measure diffusion properties
	integer(8)::meascylin=0!this is the number of cylindrical intervals in which there will be independent measurements of variables
	integer(8)::numcells(3)=0!number of cells for each dimension
	integer(8)::occup(20,20,20)=0!number of particles in each cell, occupancy lists
	integer(8)::presscount=0!this counter is used to calculate the average force between particles that is need to obtain the pressure
	logical::andersen=.false.!if true the andersen thermostat will be active
	logical::berendsen=.false.!if true the berendsen thermostat will be active
	logical::cylinders=.false.!if true, independent measurements will be taken in the different cylindrical intervals of the system
	logical::fcc=.false.!if true a face centered cubic structure will be generated
	logical::fixed(10000)=.false.!this is true for fixed particles and false for free ones that can move
	logical::gnuplot=.false.!if true, the program will output data to be visualized in gnuplot
	logical::simpcub=.false.!if true, a simple cubic structure will be generated
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
end module common_stuff

program MD_cylinders_main
	use common_stuff
	implicit none
	call input!reading initial parameters
	call startup!setting up the system
	do
		if(time>tmax)exit!the program ends when it reaches the maximum time
		call kinematics_and_dynamics!moving the system
		if((int(dmod(time,measintv),8)==0).and.(time>=tmeas))then
			if(cylinders)then
				call measuring_with_cylinders!with a certain frequency measurements are performed in each cylindrical section of the fluid cavity of the system
			else
				call measuring!with a certain frequency measurements are performed
			end if
		end if
		time=time+tstep!time advances
	end do
	call output!finishing up the writing of results
end program MD_cylinders_main



subroutine input
	use common_stuff
	implicit none
	open(10,file='input.dat',action='read')
	read(10,*) andersen,berendsen,cylinders,epsi,fcc,gnuplot,intrad,latcon,mass,meascylin,measintv,radius,sigm,simpcub,siz(1),siz(2),siz(3),temp,tempbath,tmeas,trelax,tstep,tmax!reading input data
	close(10)
	tmax=tmax+tstep*1.d-1!to avoid problems with repeated real summations, tmax is slightly modified so the last step is not lost
	return
end subroutine input



subroutine startup!the system to be simulated is set up here
	use common_stuff
	implicit none
	real(8)::dist=0.d0
	real(8)::dran_g!gaussian random number
	real(8)::dran_u!uniform random number
	integer(8)::i1=0!x dimension counter
	integer(8)::i2=0!y dimension counter
	integer(8)::i3=0!z dimension counter
	integer(8)::j=1!particle counter
	integer(8)::k=0!dimension counter
	integer(8)::seed(8)=0!seed for the random number generator
	logical::test=.false.
	call date_and_time(values=seed)!creating a seed using the clock of the computer
	call dran_ini(0)!seed(8))!initializing the random number generator
	open(20,file='position_output.dat',action='write')
	open(30,file='energy_temperature_pressure_output.dat',action='write')
	open(40,file='saved_state.dat',action='write')
	radius=2.5d-1*dmin1(siz(1),siz(2))!radius based on the size of the system
	!radius=2.d1*2.d0**(1.d0/6.d0)*sigm!radius based on 20 equlibrium distances
	volum=siz(1)*siz(2)*siz(3)!volume of the system
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
				if(test)then
					r(1,1)=siz(1)*4.5d-1!x position in the solid lattice of the material, free particles will leave these as they move
					r(2,1)=siz(2)*4.5d-1!y position in the solid lattice of the material, free particles will leave these as they move
					r(3,1)=siz(3)*4.5d-1!z position in the solid lattice of the material, free particles will leave these as they move
					r(1,2)=siz(1)*5.5d-1!x position in the solid lattice of the material, free particles will leave these as they move
					r(2,2)=siz(2)*5.5d-1!y position in the solid lattice of the material, free particles will leave these as they move
					r(3,2)=siz(3)*5.5d-1!z position in the solid lattice of the material, free particles will leave these as they move
					r(4,1)=dsqrt(kboltz*temp/mass)*dran_g()!x position in the solid lattice of the material, free particles will leave these as they move
					r(5,1)=dsqrt(kboltz*temp/mass)*dran_g()!y position in the solid lattice of the material, free particles will leave these as they move
					r(6,1)=dsqrt(kboltz*temp/mass)*dran_g()!z position in the solid lattice of the material, free particles will leave these as they move
					r(4,2)=dsqrt(kboltz*temp/mass)*dran_g()!x position in the solid lattice of the material, free particles will leave these as they move
					r(5,2)=dsqrt(kboltz*temp/mass)*dran_g()!y position in the solid lattice of the material, free particles will leave these as they move
					r(6,2)=dsqrt(kboltz*temp/mass)*dran_g()!z position in the solid lattice of the material, free particles will leave these as they move
					r(10,1)=mass!setting the mass of the first particle
					r(10,2)=mass!setting the mass of the second particle
					if(gnuplot)write(*,*) r(1,1),r(2,1),r(3,1)!data feed to gnuplot
					if(gnuplot)write(*,*) r(1,2),r(2,2),r(3,2)!data feed to gnuplot
					amnttot=2
					go to 1
				else if(simpcub)then!simple cubic structure
					r(1,j)=latcon*dble(i1)!x position in the solid lattice of the material, free particles will leave these as they move
					r(2,j)=latcon*dble(i2)!y position in the solid lattice of the material, free particles will leave these as they move
					r(3,j)=latcon*dble(i3)!z position in the solid lattice of the material, free particles will leave these as they move
					if(gnuplot)write(*,*) r(1,j),r(2,j),r(3,j)!data feed to gnuplot
					if(dsqrt((r(1,j)-5.d-1*siz(1))**2+(r(2,j)-5.d-1*siz(1))**2)>=radius)then
						fixed(j)=.true.!the particles outside of this cylinder are fixed
						amntfix=amntfix+1!counting the number of fixed particles
					else
						do k=1,3
							r(k+3,j)=dsqrt(kboltz*temp/mass)*dran_g()!free particles have a speed
						end do
						amntfree=amntfree+1!counting the number of free particles
					end if
					r(10,j)=mass!setting the mass of the particles
					j=j+1!next particle
				else if(fcc)then!face centered structure
					r(1,j)=latcon*dble(i1)!x position for the first particle of the crystalline cell
					r(2,j)=latcon*dble(i2)!y position for the first particle of the crystalline cell
					r(3,j)=latcon*dble(i3)!z position for the first particle of the crystalline cell
					r(1,j+1)=latcon*(5.d-1+dble(i1))!x position for the second particle of the crystalline cell
					r(2,j+1)=latcon*(5.d-1+dble(i2))!y position for the second particle of the crystalline cell
					r(3,j+1)=latcon*dble(i3)!z position for the second particle of the crystalline cell
					r(1,j+2)=latcon*dble(i1)!x position for the third particle of the crystalline cell
					r(2,j+2)=latcon*(5.d-1+dble(i2))!y position for the third particle of the crystalline cell
					r(3,j+2)=latcon*(5.d-1+dble(i3))!z position for the third particle of the crystalline cell
					r(1,j+3)=latcon*(5.d-1+dble(i1))!x position for the fourth particle of the crystalline cell
					r(2,j+3)=latcon*dble(i2)!y position for the fourth particle of the crystalline cell
					r(3,j+3)=latcon*(5.d-1+dble(i3))!z position for the fourth particle of the crystalline cell
					if(gnuplot)then
						write(*,*) r(1,j),r(2,j),r(3,j)!data feed to gnuplot
						write(*,*) r(1,j+1),r(2,j+1),r(3,j+1)!data feed to gnuplot
						write(*,*) r(1,j+2),r(2,j+2),r(3,j+2)!data feed to gnuplot
						write(*,*) r(1,j+3),r(2,j+3),r(3,j+3)!data feed to gnuplot
					end if
					do k=1,3
						r(k+3,j)=dsqrt(kboltz*temp/mass)*dran_g()!velocity of the first particle
						r(k+3,j+1)=dsqrt(kboltz*temp/mass)*dran_g()!velocity of the second particle
						r(k+3,j+2)=dsqrt(kboltz*temp/mass)*dran_g()!velocity of the third particle
						r(k+3,j+3)=dsqrt(kboltz*temp/mass)*dran_g()!velocity of the fourth particle
					end do
					r(10,j)=mass!setting the mass of the first particle
					r(10,j+1)=mass!setting the mass of the first particle
					r(10,j+2)=mass!setting the mass of the first particle
					r(10,j+3)=mass!setting the mass of the first particle
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
1	if(gnuplot)write(*,*) 'e'!end of gnuplot frame
	if(simpcub)then
		amnttot=j!total amount of particles for the sc structure			#####CAREFUL! YOU MAY BE MISSING THE -1 HERE
	else if(fcc)then
		amnttot=j-1!total amount of particles for the fcc structure
	end if
	call remove_average_momentum
	write(20,*) amnttot!the xyz format needs to have the amount of particles in the first line
	write(20,*) 'This is an optional comment, but it must be here for the file to be readable by the visualizing program.'!and this comment is mandatory
	write(30,*) '#	time	total energy	potential energy	kinetic energy	temperature	pressure'!magnitudes in this file
	if(simpcub)then
		dens=amntfree/(siz(3)*dacos(-1.d0)*radius**2)!this is the number density of the free particles 			#####CAREFUL! IT MAY BE THE TOTAL AMOUNT OF PARTICLES
	else if(fcc)then
		dens=amnttot/(siz(3)*dacos(-1.d0)*radius**2)!this is the number density of the free particles 			#####CAREFUL! IT MAY BE THE TOTAL AMOUNT OF PARTICLES
	end if
	do
		difftag=nint(dble(j)*dran_u(),8)!choosing randomly the particle that will be used to measure diffusion properties
		if(difftag>0)exit!just in case the tag is 0 or negative
	end do
	call create_cell_lists!cells are needed to calculate forces
	call calculate_forces!before the first integration iteration the forces need to be calculated
	return
end subroutine startup



subroutine remove_average_momentum!to avoid the system having a drift, the average momentum must be substracted to all the particles in the system
	use common_stuff
	implicit none
	real(8)::mmntsum(3)=0.d0!the global momentum will be stored here
	integer(8)::i=0!particle counter
	integer(8)::j=0!dimension counter
	do i=1,amnttot
		do j=1,3
			mmntsum(j)=mmntsum(j)+r(j+3,i)*r(10,i)!calculating the global momentum
		end do
	end do
	mmntsum=mmntsum/dble(amnttot)!dividing by the number of particles so that the final result is the average momentum of the system
	do i=1,amnttot
		do j=1,3
			r(j+3,i)=(r(j+3,i)*r(10,i)-mmntsum(j))/r(10,i)!substracting the average momentum to the momentum of each particle and converting it back to velocity
		end do
	end do
	return
end subroutine remove_average_momentum



subroutine create_cell_lists!cell lists are used to speed up the simulation, see the article in wikipedia for more general details
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
			indeces(i)=floor(r(i,j)/cellsize(i),8)+1!particle j is the cell which has these indeces
		end do
		occup(indeces(1),indeces(2),indeces(3))=occup(indeces(1),indeces(2),indeces(3))+1!updating the occupancy list with the new particle
		cells(occup(indeces(1),indeces(2),indeces(3)),indeces(1),indeces(2),indeces(3))=j!adding the new particle to the apropiate cell
	end do
	return
end subroutine create_cell_lists!			#####BE CAREFUL, IF THE SIMULATION BOX IS NOT COMPLETELY ON THE NON NEGATIVE PARTS OF THE AXES SEGMENTATION FAULTS MAY APPEAR



subroutine kinematics_and_dynamics!this is the integrator of the simulation, in this case the Verlet integrator
	use common_stuff
	implicit none
	integer(8)::i=0!x dimension cell counter
	integer(8)::j=0!y dimension cell counter
	integer(8)::k=0!z dimension cell counter
	integer(8)::m=0!dimension counter
	integer(8)::n=0!particle counter in the occupancy lists
	integer(8)::newind(3)=0!indeces of the cell in which a particle is after moving
	if(gnuplot)write(*,*) "splot '-' w p pt 7 lc 0"!new gnuplot frame
	do k=1,numcells(3)
		do j=1,numcells(2)
			do i=1,numcells(1)
				do n=1,occup(i,j,k)
					if(fixed(cells(n,i,j,k)))cycle!fixed particles don't move
					do m=1,3
						r(m+3,cells(n,i,j,k))=r(m+3,cells(n,i,j,k))+r(m+6,cells(n,i,j,k))*tstep*5.d-1/r(10,cells(n,i,j,k))!updating speeds
					end do
				end do
			end do
		end do
	end do
	do k=1,numcells(3)
		do j=1,numcells(2)
			do i=1,numcells(1)
				do n=1,occup(i,j,k)
					if(fixed(cells(n,i,j,k)))cycle!fixed particles don't move
					do m=1,3
						r(m,cells(n,i,j,k))=r(m,cells(n,i,j,k))+r(m+3,cells(n,i,j,k))*tstep!moving the particle
					end do
					do m=1,3
						if(r(m,cells(n,i,j,k))<0.d0)then
							r(m,cells(n,i,j,k))=r(m,cells(n,i,j,k))+siz(m)!if the particle has reached the lower end of one dimension of the simulation box it should be at the other side of it 								(periodic boundary conditions)
						else if(r(m,cells(n,i,j,k))>=siz(m))then
							r(m,cells(n,i,j,k))=r(m,cells(n,i,j,k))-siz(m)!if the particle has reached the higher end of one dimension of the simulation box it should be at the other side of it 								(periodic boundary conditions)
						end if
					end do
					if(gnuplot)write(*,*) r(1,cells(n,i,j,k)),r(2,cells(n,i,j,k)),r(3,cells(n,i,j,k))!data feed to gnuplot
				end do
			end do
		end do
	end do
	if(gnuplot)write(*,*) 'e'!end of gnuplot frame
	call update_particle_in_cells!after the particles have moved, they may have entered a different cell and this needs to be checked and taken care of
	call calculate_forces!the Verlet integrator requires to update the speeds again after moving the particle so the forces need to be calculated again
	do k=1,numcells(3)
		do j=1,numcells(2)
			do i=1,numcells(1)
				do n=1,occup(i,j,k)
					if(fixed(cells(n,i,j,k)))cycle!fixed particles don't move
					do m=1,3
						if(berendsen.and.andersen)then
							stop "Both thermostats can't be active at the the same time, aborting."
						else if(berendsen)then
							r(m+3,cells(n,i,j,k))=dsqrt(1.d0+(tempbath/temp-1.d0)*tstep/trelax)*(r(m+3,cells(n,i,j,k))+r(m+6,cells(n,i,j,k))*tstep*5.d-1/r(10,cells(n,i,j,k)))!updating speeds with 							the berendsen thermostat
						else if(andersen)then
							r(m+3,cells(n,i,j,k))=dsqrt(1.d0+(tempbath/temp-1.d0)*tstep/trelax)*(r(m+3,cells(n,i,j,k))+r(m+6,cells(n,i,j,k))*tstep*5.d-1/r(10,cells(n,i,j,k)))!updating speeds with the 							andersen thermostat 			#####THIS ONE IS MORE COMPLEX, WAIT TILL IT'S EXPLAINED IN CLASS
						else
							r(m+3,cells(n,i,j,k))=r(m+3,cells(n,i,j,k))+r(m+6,cells(n,i,j,k))*tstep*5.d-1/r(10,cells(n,i,j,k))!updating speeds
						end if
					end do
				end do
			end do
		end do
	end do
	return
end subroutine kinematics_and_dynamics



subroutine update_particle_in_cells
	use common_stuff
	implicit none
	integer(8)::i1=0!x dimension, old-cell counter
	integer(8)::j1=0!y dimension, old-cell counter
	integer(8)::k1=0!z dimension, old-cell counter
	integer(8)::i2=0!x dimension, new-cell counter
	integer(8)::j2=0!y dimension, new-cell counter
	integer(8)::k2=0!z dimension, new-cell counter
	integer(8)::n=0!particle counter in the cell lists
	do k1=1,numcells(3)
		do j1=1,numcells(2)
			do i1=1,numcells(1)
				n=1!starting particle counter at the beginning of the list
				do
					if(fixed(cells(n,i1,j1,k1)))cycle!fixed particles don't need updating
					i2=floor(r(1,cells(n,i1,j1,k1))/cellsize(1),8)+1!computing the x index of the cell in which the particle is located after moving
					j2=floor(r(2,cells(n,i1,j1,k1))/cellsize(2),8)+1!computing the y index of the cell in which the particle is located after moving
					k2=floor(r(3,cells(n,i1,j1,k1))/cellsize(3),8)+1!computing the z index of the cell in which the particle is located after moving
					if((i1/=i2).or.(j1/=j2).or.(k1/=k2))then!after the particle has moved, the cell in which it is contained may have changed so it needs to be updated
						occup(i2,j2,k2)=occup(i2,j2,k2)+1!updating the new occupancy list with the particle
						cells(occup(i2,j2,k2),i2,j2,k2)=cells(n,i1,j1,k1)!adding the particle to the corresponding new cell
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



subroutine calculate_forces
	use common_stuff
	implicit none
	real(8)::displ(3)=0.d0!displacement vector to calculate the appropiate distance with the periodic boundary conditions
	real(8)::dist=0.d0!distance between particles
	real(8)::distvec(3)=0.d0!vectorial distance between particles
	real(8)::force=0.d0!modulus of the force vector
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
	logical::interaction(10000,10000)=.false.!this matrix stores which interactions have been calculated already
	logical::measure=.false.!if it is time to make a measurement this variable will be true
	if((nint(dmod(time,measintv),8)==0).and.(time>=tmeas))then
		measure=.true.!to avoid computing the above many times, this logical variable is set
		epot=0.d0!setting the potential energy to zero so it can be calculated
	end if
	interaction=.false.!the interactions must be reset so they are taken into account
	r(7:9,:)=0.d0!setting forces to zero
	do k1=1,numcells(3)
		do j1=1,numcells(2)
			do i1=1,numcells(1)
				do n1=1,occup(i1,j1,k1)
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
									if(((i1==i2).and.(j1==j2)).and.((k1==k2).and.(n1==n2)))cycle!a particle shouldn't interact with itself
									if((interaction(cells(n1,i1,j1,k1),cells(n2,i2,j2,k2))).or.(interaction(cells(n2,i2,j2,k2),cells(n1,i1,j1,k1))))cycle!if this interaction has already been accounted 										for there's no need to go on
									interaction(cells(n1,i1,j1,k1),cells(n2,i2,j2,k2))=.true.!accounting for this interaction
									do l=1,3
										distvec(l)=r(l,cells(n1,i1,j1,k1))-r(l,cells(n2,i2,j2,k2))-displ(l)!these are the components of the vectorial distance between the particles
									end do
									dist=dsqrt(distvec(1)**2+distvec(2)**2+distvec(3)**2)!and this is the distance between the two particles
									if(fixed(cells(n2,i2,j2,k2)).or.fixed(cells(n1,i1,j1,k1)))then!the interaction can be slightly different between fixed and free particles
										force=2.4d1*epsi*(sigm/dist)**6*(2.d0*(sigm/dist)**6-1.d0)/dist**2!this is the module of the force between a free and a fixed particle
									else
										force=2.4d1*epsi*(sigm/dist)**6*(2.d0*(sigm/dist)**6-1.d0)/dist**2!and this is the module of the force between two free particles
									end if
									if((.not.fixed(cells(n1,i1,j1,k1))).and.(.not.fixed(cells(n2,i2,j2,k2))))then
										do l=1,3
											r(l+6,cells(n1,i1,j1,k1))=r(l+6,cells(n1,i1,j1,k1))+distvec(l)*force!and these are the components of the force vector for the first particle
											r(l+6,cells(n2,i2,j2,k2))=r(l+6,cells(n2,i2,j2,k2))-distvec(l)*force!the same for the second particle
										end do
									else if(fixed(cells(n2,i2,j2,k2)).and.(.not.fixed(cells(n1,i1,j1,k1))))then
										do l=1,3
											r(l+6,cells(n1,i1,j1,k1))=r(l+6,cells(n1,i1,j1,k1))+distvec(l)*force!and these are the components of the force vector for the first particle
										end do
									else if(fixed(cells(n1,i1,j1,k1)).and.(.not.fixed(cells(n2,i2,j2,k2))))then
										do l=1,3
											r(l+6,cells(n2,i2,j2,k2))=r(l+6,cells(n2,i2,j2,k2))-distvec(l)*force!and these are the components of the force vector for the first particle
										end do
									end if
									if(measure.and.cylinders)then
										! 			#####HERE MAYBE THERE SHOULD BE A CONTRIBUTION TO THE POTENTIAL ENERGY FOR PARTICLES LOCATED EXCLUSIVELY IN EACH CYLINDER
									else if(measure)then
										epot=epot+4.d0*epsi*(sigm/dist)**6*((sigm/dist)**6-1.d0)+2.d-2!this is the potential energy
									end if
									press=press-force*dist!to calculate the pressure, the sum over all interactions must be carried out
									!press=press+(r(7,cells(n1,i1,j1,k1))*distvec(1)+r(8,cells(n1,i1,j1,k1))*distvec(2)+r(9,cells(n1,i1,j1,k1))*distvec(3))!to calculate the pressure, the force per 										particle is needed 			#####THIS IS THE ACTUAL SCALAR PRODUCT, IF THE PREVIOUS STATEMENT CAUSES PROBLEMS OR DOESN'T WORK USE THIS ONE
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



subroutine measuring
	use common_stuff
	implicit none
	integer(8)::i=0!particle counter
	ekine=0.d0
	do i=1,amnttot
		if(.not.fixed(i))ekine=ekine+r(10,i)*(r(1,i)**2+r(2,i)**2+r(3,i)**2)!measuring the kinetic energy
	end do
	ekine=ekine*5.d-1!don't forget the 1/2
	sumekine(1)=sumekine(1)+ekine!summing the kinetic energy
	sumekine(2)=sumekine(2)+ekine**2!and the squared kinetic energy
	cvcount=cvcount+1!another measurement for the heat capacity
	temp=(2.d0*ekine)/(3.d0*amntfree*kboltz)!and with the kinetic energy the temperature is measured too 			#####CAREFUL WITH THE UNITS OF KB AND DIVIDING BY THE TOTAL AMOUNT OF PARTICLES
	press=dens*temp*kboltz+press/(3.d0*volum*dble(presscount))!final result for the pressure
	presscount=0!setting the counter to zero for the next measurement of the presure
	write(30,*) time, epot+ekine, epot, ekine, temp, press!writing data
	press=0.d0
	call diffusion
	return
end subroutine measuring



subroutine diffusion
	use common_stuff
	implicit none
	
	return
end subroutine diffusion



subroutine measuring_with_cylinders
	use common_stuff
	implicit none
	integer(8)::i=0!particle counter
	integer(8)::j=0!cylinder counter
	ekinecy=0.d0
	do i=1,amnttot
		if(.not.fixed(i))then
			j=floor(1.d1*dsqrt(r(1,i)**2+r(2,i)**2)/radius,8)!finding in which cylinder the particle is located
			if(j>10)j=10!just in case the particle is slightly beyond the radius of the cylinder
			ekinecy(j)=ekinecy(j)+r(10,i)*(r(1,i)**2+r(2,i)**2+r(3,i)**2)!measuring the kinetic energy per cylinder
		end if
	end do
	do j=1,10
		ekinecy(j)=ekinecy(j)*5.d-1!don't forget the 1/2
		sumekinecy(1,j)=sumekinecy(1,j)+ekinecy(j)!summing the kinetic energy per cylinder
		sumekinecy(2,j)=sumekinecy(2,j)+ekinecy(j)**2!and the squared kinetic energy per cylinder
		cvcountcy(j)=cvcountcy(j)+1!another measurement for the heat capacity per cylinder
		tempcy(j)=(2.d0*ekinecy(j))/(3.d0*amntfree*kboltz)!and with the kinetic energy the temperature is measured too per cylinder 			#####CAREFUL WITH THE UNITS OF KB AND DIVIDING BY THE 			TOTAL AMOUNT OF PARTICLES
	end do
	call diffusion_with_cylinders
	return
end subroutine measuring_with_cylinders



subroutine diffusion_with_cylinders
	use common_stuff
	implicit none
	
	return
end subroutine diffusion_with_cylinders



subroutine output
	use common_stuff
	implicit none
	integer(8)::i=0!cylinder counter
	cv=(3.d0*kboltz)/((4.d0*amntfree*(sumekine(1)**2/dble(cvcount)-sumekine(2))/(dble(cvcount)*3.d0*kboltz*temp**2))-2.d0)!computing the final value of the heat capacity
	if(cylinders)then
		do i=1,10
			cvcy(i)=(3.d0*kboltz)/((4.d0*amntfree*(sumekinecy(1,i)**2/dble(cvcountcy(i))-sumekinecy(2,i))/(dble(cvcountcy(i))*3.d0*kboltz*tempcy(i)**2))-2.d0)!computing the final value of the heat 				capacity
		end do
	end if
	write(30,*) '# The heat capacity is: ',cvcy(1),cvcy(2),cvcy(3),cvcy(4),cvcy(5),cvcy(6),cvcy(7),cvcy(8),cvcy(9),cvcy(10)!writing data
	close(20)
	close(30)
	close(40)
	return
end subroutine output
