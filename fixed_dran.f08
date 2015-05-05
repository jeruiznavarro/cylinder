subroutine dran_ini(iseed0)
implicit double precision(a-h,o-z)
parameter(ip=1279)
parameter(np=14)
parameter(nbit=31)
parameter(m=2**np,np1=nbit-np,nn=2**np1-1,nn1=nn+1)
integer ix(ip)
dimension g(0:m)

data c0,c1,c2/2.515517,0.802853,0.010328/
data d1,d2,d3/1.432788,0.189269,0.001308/
data pi/3.141592653589793d0/

common /ixx/ ix
common /icc/ ic
common /gg/ g

dseed=iseed0
do i=1,ip
 ix(i)=0
 do j=0,nbit-1
    if(rand_xx(dseed).lt.0.5) ix(i)=ibset(ix(i),j)
 enddo
enddo
ic=0

do i=m/2,m
	p=1.0-dble(i+1)/(m+2)
	t=sqrt(-2.0*log(p))
	x=t-(c0+t*(c1+c2*t))/(1.0+t*(d1+t*(d2+t*d3)))
	g(i)=x
	g(m-i)=-x
enddo

u2th=1.0-dble(m+2)/m*sqrt(2.0/pi)*g(m)*exp(-g(m)*g(m)/2)
u2th=nn1*sqrt(u2th)
g=g/u2th

return
end



function i_dran(n)
implicit double precision(a-h,o-z)
parameter(ip=1279)
parameter(iq=418)
parameter(is=ip-iq)
integer ix(ip)
common /ixx/ ix
common /icc/ ic
ic=ic+1
if (ic > ip) ic=1
if (ic > iq) then
ix(ic)=ieor(ix(ic),ix(ic-iq))
else
ix(ic)=ieor(ix(ic),ix(ic+is))
endif
i_dran=ix(ic)
if (n > 0) i_dran=mod(i_dran,n)+1
return
end




function dran_u()
implicit double precision(a-h,o-z)
parameter(ip=1279)
parameter(iq=418)
parameter(is=ip-iq)
parameter (rmax=2147483648.0d0)
integer ix(ip)
common /ixx/ ix
common /icc/ ic
ic=ic+1
if (ic > ip) ic=1
if (ic > iq) then
ix(ic)=ieor(ix(ic),ix(ic-iq))
else
ix(ic)=ieor(ix(ic),ix(ic+is))
endif
dran_u=(dble(ix(ic))+0.5d0)/rmax
return
end



function dran_g()
implicit double precision(a-h,o-z)
parameter(ip=1279)
parameter(iq=418)
parameter(np=14)
parameter(nbit=31)
parameter(m=2**np,np1=nbit-np,nn=2**np1-1,nn1=nn+1)
parameter(is=ip-iq)

integer ix(ip)
dimension g(0:m)

common /ixx/ ix
common /icc/ ic
common /gg/ g

ic=ic+1
if(ic > ip) ic=1
if(ic > iq) then
ix(ic)=ieor(ix(ic),ix(ic-iq))
else
ix(ic)=ieor(ix(ic),ix(ic+is))
endif
i=ishft(ix(ic),-np1)
i2=iand(ix(ic),nn)
dran_g=i2*g(i+1)+(nn1-i2)*g(i)
return
end



function dran_gbmw()
implicit double precision(a-h,o-z)
parameter(ip=1279)
parameter(iq=418)
parameter(is=ip-iq)
parameter (rmax=2147483648.0d0)
integer ix(ip)
integer, save :: icount=1
double precision, save :: u,v
common /ixx/ ix
common /icc/ ic
data pi2 /6.283185307179586d0/
if (icount == 1) then
	ic=ic+1
	if (ic > ip) ic=1
	if (ic > iq) then
	ix(ic)=ieor(ix(ic),ix(ic-iq))
	else
	ix(ic)=ieor(ix(ic),ix(ic+is))
	endif
	u=pi2*dble(ix(ic)+0.5d0)/rmax
	ic=ic+1
	if(ic > ip) ic=1
	if(ic > iq) then
	ix(ic)=ieor(ix(ic),ix(ic-iq))
	else
	ix(ic)=ieor(ix(ic),ix(ic+is))
	endif
	v=(dble(ix(ic))+0.5d0)/rmax
	v=dsqrt(-2.0d0*log(v))
	dran_gbmw=dcos(u)*v
	icount=2
else
  dran_gbmw=dsin(u)*v
  icount=1
endif
return
end



function rand_xx(dseed)
double precision a,c,xm,rm,dseed,rand_xx
parameter (xm=2.d0**32,rm=1.d0/xm,a=69069.d0,c=1.d0)
dseed=mod(dseed*a+c,xm)
rand_xx=dseed*rm
return
end