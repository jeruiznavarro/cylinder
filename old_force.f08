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
		do l=1,3
			fvec(l,n1)=0.d0!setting forces to zero
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
