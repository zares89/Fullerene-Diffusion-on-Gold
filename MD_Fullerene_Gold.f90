! This code runs a molecular dynamic simulation in NVT ensemble for the diffusion of a fullerene particle on a gold substrate. Sutton-Chen and Tersoff force fields are used. 
Program MDgold_Fullerine
use DFPORT
use LIBM
use DFLIB
implicit none

real*8 :: rx,ry,rz,eps,c,a,parsc,parter,parlj,scfx,scfy,scfz,acell,phi,t1,t2,timecpu,rcut,c11,boxx,boxy,boxz,Esc,Et,fterx,ftery,fterz
real*8 :: rxf,ryf,rzf,gx,gy,gz
real*8 :: temp,dt
real*8 :: B,lambda1,lambda2,beta,alpha,nt,ct,ds,h,lambda3,R,D
real*8 :: sigma
integer :: nmax,m,n,nall,i,natoms,nmd,nallt,nallf,nalltf
parameter (Nmax=5000)
dimension rx(nmax),ry(nmax),rz(nmax),parsc(5),parter(13),parlj(2),scfx(nmax),scfy(nmax),scfz(nmax),rxf(nmax),ryf(nmax),rzf(nmax),fterx(nmax),ftery(nmax),fterz(nmax)
call cpu_time(t1)
! Sutton-Chen paramesuts for Au
eps=1.2793e-2;m=8;n=10;c=34.408;a=4.08
parsc(1)=eps;parsc(2)=m;parsc(3)=n;parsc(4)=c;parsc(5)=a
acell=4.08;rcut=2*acell
! Tersoff parameters for Carbon (Tersoff 1989)
A=1.3936e3;B=3.467e2;lambda1=3.4879;lambda2=2.2119;beta=1.5724e-7;
alpha=0.;nt=7.2751e-1;ct=3.8049e4;ds=4.384;h=-5.7058e-1;lambda3=0.;R=1.95;D=0.15;
!par=[A B lambda1 lambda2 lambda3 alpha beta n c ds h R D];
parter(1)=A;parter(2)=B;parter(3)=lambda1;parter(4)=lambda2;parter(5)=lambda3;parter(6)=alpha;
parter(7)=beta;parter(8)=nt;parter(9)=ct;parter(10)=ds;parter(11)=h;parter(12)=R;parter(13)=D
! LJ parameters for C-Au
eps=0.01273;sigma=2.9943
parlj(1)=eps;parlj(2)=sigma
!------------------------------------------------------------------
! MD setup
gx=6;gy=6;gz=2
nmd=3000
temp=300
dt=1.e-15
open(22,file='xyzgold.dat',status='old')
call incrd(nall,nallt,nmax,rx,ry,rz)
close(22)

!call periodic(rx1,ry1,rz1,rx,ry,rz,natoms,nall,acell,nmax,boxx,boxy,boxz)
!nall=natoms
rx=rx*acell;ry=ry*acell;rz=rz*acell
boxx=maxval(rx)-minval(rx)+acell
boxy=maxval(ry)-minval(ry)+acell
boxz=(maxval(rz)-minval(rz))*10
open(22,file='xyzc60.dat',status='old')
! Call subroutine to initialze the system with random velocities according to Maxwell-Boltzmann distribution
call incrd(nallf,nalltf,nmax,rxf,ryf,rzf)
close(22)
! Run MD simulation
call  MD(rx,ry,rz,rxf,ryf,rzf,nall,nallt,nallf,nalltf,temp,parsc,parter,parlj,acell,rcut,nmd,boxx,boxy,boxz)

call cpu_time(t2)
timecpu=t2-t1
open(44,file='elapsed_time.dat',status='unknown')
write(44,*) timecpu
close(44)

contains

subroutine incrd(nall,nallt,nmax,gcoordmx,gcoordmy,gcoordmz)

implicit none
integer :: nall,nmax,i,nallt
PARAMETER (Nmax=5000)
real*8 :: gcoordmx,gcoordmy,gcoordmz
dimension gcoordmx(nmax),gcoordmy(nmax),gcoordmz(nmax)

read(22, *) nall
read(22, *) nallt
!------------
read(22, *) (gcoordmx(i),gcoordmy(i),gcoordmz(i),i=1,nall)
end
!----The subroutine that sets up the initial random velocities
!----according to Maxwell-Boltzmann distribution-------------
subroutine strtup(vx,vy,vz,nallt,temp,massp)

implicit none

real*8 :: rtemp,dummy,gauss,sumx,sumy,sumz,bol
real*8 :: vx,vy,vz,temp
real*8 :: ZBQLNOR
real*8 :: massp,mass,rmass,boltzct,rbolt
real*8 :: rmaxy
integer :: nall,i,ip,count,nallt
integer :: pp,qq,intsiz,idelay,mm(100),max,ii,k1
integer :: j,k,n,seed
dimension vx(nallt),vy(nallt),vz(nallt)

boltzct=8.617332478e-5
rbolt=sqrt(boltzct)
rtemp=sqrt(temp)
mass=massp
rmass=sqrt(mass)
do i=1,nallt
	vx(i)=rtemp*rbolt*sqrt(16.02*3)*ZBQLNOR(0.,1.)/(1.5*rmass)
	vy(i)=rtemp*rbolt*sqrt(16.02*3)*ZBQLNOR(0.,1.)/(1.5*rmass)
	vz(i)=rtemp*rbolt*sqrt(16.02*3)*ZBQLNOR(0.,1.)/(1.5*rmass)
end do

sumx=0.
sumy=0.
sumz=0.
do i=1,nallt
	sumx=sumx+mass*vx(i)
	sumy=sumy+mass*vy(i)
	sumz=sumz+mass*vz(i)
end do
sumx=sumx/dble(nallt)
sumy=sumy/dble(nallt)
sumz=sumz/dble(nallt)
open(57,file='IniVel.dat',status='unknown')
do i=1,nallt
	vx(i)=vx(i)-sumx/mass
	vy(i)=vy(i)-sumy/mass
	vz(i)=vz(i)-sumz/mass
	write(57,*) i,vx(i)
end do
close(57)
end
!****************************
function gauss(dummy)

implicit none

real*8 :: dummy,a1,a3,a5,a7,a9,sum,r,r2,gauss,y
real*8 :: x
integer :: i,count,seed
parameter (a1=3.949846138,a3=0.252408784,a5=0.076542912)
parameter (a7=0.008355968,a9=0.029899776)

call random_seed(seed) 
sum=0.
do i=1,12
	call random_number(r)
	sum=sum+r
end do
!
r=(sum-6.0)/4.0
r2=r*r
gauss=((((a9*r2+a7)*r2+a5)*r2+a3)*r2+a1)*r
end
!--------------------------------------------------------
! Subroutine to evaluate Force and Potential energy
subroutine SCforce(nall,boxx,boxy,boxz,rcut,rx,ry,rz,par,scfx,scfy,scfz,Esc)

implicit none

real*8 :: scfx,scfy,scfz,rx,ry,rz,par,phi,rouk,rouj,vr,eps,d,dx,dy,dz,c,rcut,boxx,boxy,boxz,Esc
dimension scfx(nall),scfy(nall),scfz(nall),rx(nall),ry(nall),rz(nall),par(5)
integer :: m,n
integer :: k,j,nall
integer :: minim

eps=par(1);m=par(2);n=par(3);c=par(4)
scfx(:)=0.;scfy(:)=0.;scfz(:)=0.

Esc=0.
do k=1,nall
	call roum(nall,boxx,boxy,boxz,rcut,k,rx,ry,rz,par,rouk)
	do j=1,nall
		dx=rx(k)-rx(j)
		dy=ry(k)-ry(j)
		dz=rz(k)-rz(j)
		dx=dx-boxx*anint(dx/boxx)
		dy=dy-boxy*anint(dy/boxy)
		d=dx*dx+dy*dy+dz*dz
		if (d<=rcut*rcut .and. d>0.) then
		call vpair(par,vr,d)
		call phim(par,phi,d)
		call roum(nall,boxx,boxy,boxz,rcut,j,rx,ry,rz,par,rouj)
		scfx(k)=scfx(k)+eps*(n*vr/d*dx-c*+1./2.*m*phi/d*dx* &
		(1/sqrt(rouk)+1/sqrt(rouj)))
		scfy(k)=scfy(k)+eps*(n*vr/d*dy-c*+1./2.*m*phi/d*dy* &
		(1/sqrt(rouk)+1/sqrt(rouj)))
		scfz(k)=scfz(k)+eps*(n*vr/d*dz-c*+1./2.*m*phi/d*dz* &
		(1/sqrt(rouk)+1/sqrt(rouj)))
		Esc=eps*1./2.*vr+Esc
		end if
	enddo
	Esc=Esc-eps*c*sqrt(rouk)
enddo
write(*,*) Esc
end
!-------------------------------------------------------------
subroutine vpair(par,vr,d)

implicit none

real*8 :: vr,par,a,d
integer :: n
dimension par(5)
a=par(5);n=par(3)
vr=(a/sqrt(d))**n

end
!-----------------------------------------------------------------------
subroutine phim(par,phi,d)

implicit none

real*8 :: par,phi,d,a
integer :: m
dimension par(5)

a=par(5);m=par(2)

phi=(a/sqrt(d))**m

end
!-------------------------------------------------------------------------
subroutine roum(nall,boxx,boxy,boxz,rcut,j,rx,ry,rz,par,rou)

implicit none

real*8 :: rx,ry,rz,par,rou,d,rcut,dx,dy,dz,boxx,boxy,boxz
dimension rx(nall),ry(nall),rz(nall),par(5)
integer :: i,j,nall
rou=0.
do i=1,nall
		dx=rx(i)-rx(j)
		dy=ry(i)-ry(j)
		dz=rz(i)-rz(j)
		dx=dx-boxx*anint(dx/boxx)
		dy=dy-boxy*anint(dy/boxy)
		d=dx*dx+dy*dy+dz*dz
		if (d<=rcut*rcut .and. d>0.) then
			call phim(par,phi,d)
			rou=rou+phi
		end if
enddo
end 
!------------------------------------------------------------------------------------	
subroutine verletgold(nall,nallt,par,massp,rx,ry,rz,vx,vy,vz,dt,boxx,boxy,boxz,temp,curtemp,ksi,Esc)

implicit none

real*8 :: rx,ry,rz,vx,vy,vz,ax,ay,az,ax1,ay1,az1
real*8 :: dt
real*8 :: boxx,boxy,boxz
real*8 :: ksi,massp,curtemp,ksi1,temp,Esc
real*8 :: Q,boltzct,tau,sumk,sumk1,q1
real*8 :: par
dimension par(5)
integer :: minim,j,i
integer :: nall,nallt
dimension rx(nall),ry(nall),rz(nall),vx(nallt),vy(nallt),vz(nallt),ax(nall),ay(nall),az(nall),ax1(nall),ay1(nall),az1(nall)
minim=0
dt=1.e-15
tau=100*dt
boltzct=8.617332478e-5
call SCforce(nall,boxx,boxy,boxz,rcut,rx,ry,rz,par,scfx,scfy,scfz,Esc)
call acc(nall,scfx,scfy,scfz,ax1,ay1,az1,massp)
Q=3*(nallt-1)*boltzct*temp*tau**2
do i=1,nallt
	rx(i)=rx(i)+vx(i)*dt+0.5*ax1(i)*dt**2-0.5*dt**2*ksi*vx(i)
	ry(i)=ry(i)+vy(i)*dt+0.5*ay1(i)*dt**2-0.5*dt**2*ksi*vy(i)
	rz(i)=rz(i)+vz(i)*dt+0.5*az1(i)*dt**2-0.5*dt**2*ksi*vz(i)
end do
ksi1=ksi
sumk=0.
sumk1=0.
do i=1,nallt
	sumk=sumk+massp*(vx(i)+dt*0.5*(ax1(i)-ksi1*vx(i)))**2+massp*(vy(i)+dt*0.5*(ay1(i)-ksi1*vy(i)))**2+ &
	massp*(vz(i)+dt*0.5*(az1(i)-ksi1*vz(i)))**2
	sumk1=sumk1+massp*(vx(i)**2+vy(i)**2+vz(i)**2)
end do
sumk=sumk*6.24181*0.01
sumk1=sumk1*6.24181*0.01
ksi=ksi+dt/(2*Q)*(sumk1-3*(nallt-1)*boltzct*temp)+dt/(2*Q)*(sumk-3*(nallt-1)*boltzct*temp)
do j=1,nall
		rx(j)=rx(j)-boxx*anint(rx(j)/boxx)
		ry(j)=ry(j)-boxy*anint(ry(j)/boxy)
		rz(j)=rz(j)-boxz*anint(rz(j)/boxz)
end do
call SCforce(nall,boxx,boxy,boxz,rcut,rx,ry,rz,par,scfx,scfy,scfz,Esc)
call acc(nallt,scfx,scfy,scfz,ax,ay,az,massp)
do i=1,nallt
	vx(i)=2./(2+ksi*dt)*(vx(i)+0.5*dt*(ax1(i)-ksi1*vx(i))+0.5*dt*ax(i))
	vy(i)=2./(2+ksi*dt)*(vy(i)+0.5*dt*(ay1(i)-ksi1*vy(i))+0.5*dt*ay(i))
	vz(i)=2./(2+ksi*dt)*(vz(i)+0.5*dt*(az1(i)-ksi1*vz(i))+0.5*dt*az(i))
end do
q1=ax(1)/ax1(1)
curtemp=sumk1/(boltzct*dble((3*nallt)))
end

subroutine acc(nallt,scfx,scfy,scfz,ax,ay,az,massp)

real*8 :: scfx,scfy,scfz,ax,ay,az
real*8 :: massp
dimension scfx(nall),scfy(nall),scfz(nall)
dimension ax(nall),ay(nall),az(nall)
integer :: nallt

do i=1,nallt
ax(i)=16.0210*scfx(i)/massp;ay(i)=16.0210*scfy(i)/massp;az(i)=16.0210*scfz(i)/massp
end do

end
!-------------------------------------------------------------------------------
subroutine tersoff(Et,Fterx,Ftery,Fterz,rx,ry,rz,par,nall)

implicit none

real*8 :: Et,rx,ry,rz,Fterx,Ftery,Fterz
real*8 :: beta,lambda3,n,R,D,c,ds,h,rs,vr,fc,fr,fa,bs
real*8 :: dx,dy,dz,ksij,dij,cont,sumb,cont1,sumb1,kesi,kesi1,cont2
integer :: i,j,nall,k
dimension rx(nall),ry(nall),rz(nall),Fterx(nall),Ftery(nall),Fterz(nall)
real*8, dimension(13) :: par
real*8, dimension(3) :: partfr,partfa,partr,partrn,partfci,partfcj,partb,partfc
!par=[A B lambda1 lambda2 lambda3 alpha beta n c ds h R D];
beta=par(7)
lambda3=par(5)
n=par(8)
R=par(12)
D=par(13)
c=par(9)
ds=par(10)
h=par(11)
lambda1=par(3);lambda2=par(4)
Et=0.
!parsay1=parsay(1)
sumb1=sumb
kesi1=kesi
kesi=0
cont=fterx(1)
sumb=0
Fterx(:)=0.;Ftery(:)=0.;Fterz(:)=0.;

do i=1,nall

	do j=1,nall
		if (i/=j) then
		dx=rx(i)-rx(j)
		dy=ry(i)-ry(j)
		dz=rz(i)-rz(j)
		dx=dx-boxx*anint(dx/boxx)
		dy=dy-boxy*anint(dy/boxy)
		dz=dz-boxz*anint(dz/boxz)
		dij=dx*dx+dy*dy+dz*dz
		rs=SQRT(dij)
		call fcut(rs,R,D,fc)
		!fc=fcut(r,R,D)
		if (fc==0.) then
			vr=0.
		else
			!Calling manybody potential term
			call manybody(i,j,dx,dy,dz,rs,rx,ry,rz,beta,lambda3,n,R,D,c,ds,h,bs,nall,ksij)
		
			call frep(par(1),rs,par(3),fr)
			call fatt(rs,par(2),par(4),fa)
			vr=fc*(fr+bs*fa)
			partr(1)=1/rs*dx
			partr(2)=1/rs*dy
			partr(3)=1/rs*dz
			partrn=-partr
			do k=1,nall
				if (k==i)  then
					call partialfc(rs,R,D,partr,partfc)
					call partmanybody(i,j,k,dx,dy,dz,rs,rx,ry,rz,partr,beta,lambda3,n,R,D,c,ds,h,partb,nall,ksij)
					partfr(1)=-lambda1*fr*partr(1);partfr(2)=-lambda1*fr*partr(2);partfr(3)=-lambda1*fr*partr(3)
					partfa(1)=-lambda2*fa*partr(1);partfa(2)=-lambda2*fa*partr(2);partfa(3)=-lambda2*fa*partr(3)
					Fterx(k)=Fterx(k)-0.5*(partfc(1)*(fr+bs*fa)+fc*(partfr(1)+partb(1)*fa+bs*partfa(1)))
					Ftery(k)=Ftery(k)-0.5*(partfc(2)*(fr+bs*fa)+fc*(partfr(2)+partb(2)*fa+bs*partfa(2)))
					Fterz(k)=Fterz(k)-0.5*(partfc(3)*(fr+bs*fa)+fc*(partfr(3)+partb(3)*fa+bs*partfa(3)))
				else if (k==j) then
					call partialfc(rs,R,D,partrn,partfc)
					call partmanybody(i,j,k,dx,dy,dz,rs,rx,ry,rz,partr,beta,lambda3,n,R,D,c,ds,h,partb,nall,ksij)
					partfr(1)=-lambda1*fr*partrn(1);partfr(2)=-lambda1*fr*partrn(2);partfr(3)=-lambda1*fr*partrn(3)
					partfa(1)=-lambda2*fa*partrn(1);partfa(2)=-lambda2*fa*partrn(2);partfa(3)=-lambda2*fa*partrn(3)
					Fterx(k)=Fterx(k)-0.5*(partfc(1)*(fr+bs*fa)+fc*(partfr(1)+partb(1)*fa+bs*partfa(1)))
					Ftery(k)=Ftery(k)-0.5*(partfc(2)*(fr+bs*fa)+fc*(partfr(2)+partb(2)*fa+bs*partfa(2)))
					Fterz(k)=Fterz(k)-0.5*(partfc(3)*(fr+bs*fa)+fc*(partfr(3)+partb(3)*fa+bs*partfa(3)))
				else
					call partmanybody(i,j,k,dx,dy,dz,rs,rx,ry,rz,partr,beta,lambda3,n,R,D,c,ds,h,partb,nall,ksij)
					Fterx(k)=Fterx(k)-0.5*(fc*partb(1)*fa)
					Ftery(k)=Ftery(k)-0.5*(fc*partb(2)*fa)
					Fterz(k)=Fterz(k)-0.5*(fc*partb(3)*fa)
				endif
			
			kesi=ksij**2+kesi
			end do
		end if
		Et=Et+0.5*vr
		
		end if
	enddo
	sumb=Fterx(1)
enddo
cont=fterx(1)/cont
cont1=sumb/sumb1
!cont2=parsay(1)/parsay1
write(*,*) Et
end
!--------------------------------------------
subroutine fcut(rs,R,D,fc)

implicit none
real*8 :: rs,R,D,fc,rpd,rmd
real*8, parameter :: pi1 = 3.141592653589793239
rpd=R+D
rmd=R-D
if (rs<rmd) then
	fc=1.0
else if (rmd<rs .and. rs<rpd) then
	fc=1.0/2.0-1.0/2.0*sin(pi1/2.*(rs-R)/D)
else
	fc=0.
end if
end
!--------------------------------------------
!--------------------------------
! Repulsive potential subroutine
subroutine frep(A,rs,lambda1,fr)
implicit none 
real*8 :: A,rs,lambda1,fr
fr=A*EXP(-lambda1*rs)
end
!---------------------------
! Attractive potential subroutine
subroutine fatt(rs,Bs,lambda2,fa)
implicit none
real*8 :: rs,Bs,lambda2,fa
fa=-Bs*EXP(-lambda2*rs)
end
!-------------------------------------

!------------------------
subroutine partialfc(rs,R,D,partr,partfc)

implicit none

real*8 :: rs,R,D
real*8, dimension(3) :: partr,partfc
real*8, parameter :: pi1 = 3.141592653589793239

if (rs < R-D .OR. rs > R+D) then
	partfc(1)=0.
	partfc(2)=0.
	partfc(3)=0.
else
	partfc(1)=-pi1/(4.*D)*partr(1)*cos(pi1/2.*(rs-R)/D)
	partfc(2)=-pi1/(4.*D)*partr(2)*cos(pi1/2.*(rs-R)/D)
	partfc(3)=-pi1/(4.*D)*partr(3)*cos(pi1/2.*(rs-R)/D)
end if
end

!------------------------------------------------------------------------------------------
subroutine partmanybody(i,j,k,dxij,dyij,dzij,rij,rx,ry,rz,partr,beta,lambda3,n,R,D,c,ds,h,partb,nall,ksij)

implicit none

real*8 :: beta,lambda3,n,R,D,ds,h,ksij,cteta,gtet,steta,fc,rij,rik,rim,c,vcn
real*8, dimension(nall) :: rx,ry,rz
real*8 :: dxim,dyim,dzim,dim,dxij,dyij,dzij,dxjm,dyjm,dzjm,dxg,dyg,dzg
real*8 :: partfci,partfcj
real*8, dimension(3) :: partb,parsay,parsay1,parsay2,partrm,va,vb,vc,partrijk,partrimk,partetm,pargk,partfc,partr
integer :: i,j,k,m,nall
parsay(:)=0
partb(:)=0
do m=1,nall
	if (m/=i .and. m/=j) then
		if (k==i .or. k==j .or. k==m) then
			dxim=rx(i)-rx(m)
			dyim=ry(i)-ry(m)
			dzim=rz(i)-rz(m)
			dxim=dxim-boxx*anint(dxim/boxx)
			dyim=dyim-boxy*anint(dyim/boxy)
			dim=dxim*dxim+dyim*dyim+dzim*dzim
			rim=SQRT(dim)
			call fcut(rim,R,D,fc)
			if (fc/=0) then
			va(1)=dxim
			va(2)=dyim
			va(3)=dzim
			vb(1)=dxij
			vb(2)=dyij
			vb(3)=dzij
			cteta=DOT_PRODUCT(va,vb)/(rij*rim)
			gtet=1.+c**2/ds**2-c**2/(ds**2+(h-cteta)**2)
			if (k==i) then
				dxg=(dxim+dxij)
				dyg=(dyim+dyij)
				dzg=(dzim+dzij)
				partrimk(1)=1/rim*dxim
				partrimk(2)=1/rim*dyim
				partrimk(3)=1/rim*dzim
				partrijk(1)=partr(1)
				partrijk(2)=partr(2)
				partrijk(3)=partr(3)
				else if (k==j) then
				dxg=-dxim
				dyg=-dyim
				dzg=-dzim
				partrimk(1)=0
				partrimk(2)=0
				partrimk(3)=0
				partrijk(1)=-partr(1)
				partrijk(2)=-partr(2)
				partrijk(3)=-partr(3)
				else
				dxg=-dxij
				dyg=-dyij
				dzg=-dzij
				partrimk(1)=-1/rim*dxim
				partrimk(2)=-1/rim*dyim
				partrimk(3)=-1/rim*dzim
				partrijk(1)=0
				partrijk(2)=0
				partrijk(3)=0
			end if
			pargk(1)=-(h-cteta)*2.*c**2/(ds**2+(h-cteta)**2)**2*(rij*rim*dxg-(dxim*dxij+dyim*dyij+dzim*dzij)* &
			(partrijk(1)*rim+partrimk(1)*rij))/(rij*rim)**2
			pargk(2)=-(h-cteta)*2.*c**2/(ds**2+(h-cteta)**2)**2*(rij*rim*dyg-(dxim*dxij+dyim*dyij+dzim*dzij)* &
			(partrijk(2)*rim+partrimk(2)*rij))/(rij*rim)**2
			pargk(3)=-(h-cteta)*2.*c**2/(ds**2+(h-cteta)**2)**2*(rij*rim*dzg-(dxim*dxij+dyim*dyij+dzim*dzij)* &
			(partrijk(3)*rim+partrimk(3)*rij))/(rij*rim)**2
			
			call partialfc(rim,R,D,partrimk,partfc)

			parsay(1)=parsay(1)+(partfc(1)*gtet+fc*pargk(1))
			parsay(2)=parsay(2)+(partfc(2)*gtet+fc*pargk(2))
			parsay(3)=parsay(3)+(partfc(3)*gtet+fc*pargk(3))
			end if
		end if
	end if
end do
if (ksij/=0) then
partb(1)=-1./2.*(1+beta**n*ksij**n)**(-1./dble(2*n)-1)*(beta**n*ksij**(n-1))*parsay(1)
partb(2)=-1./2.*(1+beta**n*ksij**n)**(-1./dble(2*n)-1)*(beta**n*ksij**(n-1))*parsay(2)
partb(3)=-1./2.*(1+beta**n*ksij**n)**(-1./dble(2*n)-1)*(beta**n*ksij**(n-1))*parsay(3)
end if
end

!---------------------------------------------------------------------------------------------

subroutine manybody(i,j,dxij,dyij,dzij,rij,rx,ry,rz,beta,lambda3,n,R,D,c,ds,h,bs,nall,ksij)
implicit none
integer :: i,j,k,nall
real*8 :: dx,dy,dz,dik,dxij,dyij,dzij
real*8 :: rx,ry,rz,beta,lambda3,n,R,D,c,ds,h,bs,ksij,cteta,gtet,rij,rik,fc
real*8,dimension(3) :: va,vb
dimension rx(nall),ry(nall),rz(nall)

ksij=0.
do k=1,nall
	if (k/=i .and. k/=j) then
		dx=rx(i)-rx(k)
		dy=ry(i)-ry(k)
		dz=rz(i)-rz(k)
		dx=dx-boxx*anint(dx/boxx)
		dy=dy-boxy*anint(dy/boxy)
		dik=dx*dx+dy*dy+dz*dz
		rik=sqrt(dik)
		call fcut(rik,R,D,fc)
		if (fc/=0) then
		va(1)=dx
		va(2)=dy
		va(3)=dz
		vb(1)=dxij
		vb(2)=dyij
		vb(3)=dzij
		cteta=DOT_PRODUCT(va,vb)/(rij*rik)
		gtet=1.+c**2/ds**2-c**2/(ds**2+(h-cteta)**2)
		ksij=ksij+fc*gtet*EXP(lambda3**3*(rij-rik)**3)
		end if
	end if
end do
bs=(1.+beta**n*ksij**n)**(-1.0/(2.0*n))
end
!----------------------------------------------------------------------------------------
subroutine verletful(nall,nallt,par,massp,rx,ry,rz,vx,vy,vz,dt,boxx,boxy,boxz,temp,curtemp,ksi,Et,sumb)

implicit none

real*8 :: rx,ry,rz,vx,vy,vz,ax,ay,az,ax1,ay1,az1,Fterx,Ftery,Fterz
real*8 :: dt
real*8 :: boxx,boxy,boxz
real*8 :: ksi,massp,curtemp,ksi1,temp,Et
real*8 :: Q,boltzct,tau,sumk,sumk1
real*8 :: par,rat,ratk,q1,sumb,kesi
dimension par(5)
integer :: minim,j,i
integer :: nall,nallt
dimension rx(nall),ry(nall),rz(nall),vx(nallt),vy(nallt),vz(nallt),ax(nall),ay(nall),az(nall),ax1(nall),ay1(nall),az1(nall)
dimension Fterx(nall),Ftery(nall),Fterz(nall)
minim=0
dt=1.e-15
tau=5*dt
boltzct=8.617332478e-5
call tersoff(Et,Fterx,Ftery,Fterz,rx,ry,rz,par,nall)
call acc(nall,Fterx,Ftery,Fterz,ax1,ay1,az1,massp)
Q=3*(nallt-1)*boltzct*temp*tau**2
do i=1,nallt
	rx(i)=rx(i)+vx(i)*dt+0.5*ax1(i)*dt**2-0.5*dt**2*ksi*vx(i)
	ry(i)=ry(i)+vy(i)*dt+0.5*ay1(i)*dt**2-0.5*dt**2*ksi*vy(i)
	rz(i)=rz(i)+vz(i)*dt+0.5*az1(i)*dt**2-0.5*dt**2*ksi*vz(i)
end do
ksi1=ksi
sumk=0.
sumk1=0.
do i=1,nallt
	sumk=sumk+massp*(vx(i)+dt*0.5*(ax1(i)-ksi1*vx(i)))**2+massp*(vy(i)+dt*0.5*(ay1(i)-ksi1*vy(i)))**2+ &
	massp*(vz(i)+dt*0.5*(az1(i)-ksi1*vz(i)))**2
	sumk1=sumk1+massp*(vx(i)**2+vy(i)**2+vz(i)**2)
end do
sumk=sumk*6.24181*0.01
sumk1=sumk1*6.24181*0.01
ksi=ksi+dt/(2*Q)*(sumk1-3*(nallt-1)*boltzct*temp)+dt/(2*Q)*(sumk-3*(nallt-1)*boltzct*temp)
do j=1,nall
		rx(j)=rx(j)-boxx*anint(rx(j)/boxx)
		ry(j)=ry(j)-boxy*anint(ry(j)/boxy)
		rz(j)=rz(j)-boxz*anint(rz(j)/boxz)
end do
call tersoff(Et,Fterx,Ftery,Fterz,rx,ry,rz,par,nall)
call acc(nallt,Fterx,Ftery,Fterz,ax,ay,az,massp)
open(59,file='Vel.dat',status='unknown')
do i=1,nallt
	vx(i)=2./(2+ksi*dt)*(vx(i)+0.5*dt*(ax1(i)-ksi1*vx(i))+0.5*dt*ax(i))
	vy(i)=2./(2+ksi*dt)*(vy(i)+0.5*dt*(ay1(i)-ksi1*vy(i))+0.5*dt*ay(i))
	vz(i)=2./(2+ksi*dt)*(vz(i)+0.5*dt*(az1(i)-ksi1*vz(i))+0.5*dt*az(i))
	write(59,*) i,vx(i)
end do
close(59)
q1=ax(1)/ax1(1)
curtemp=sumk1/(boltzct*dble((3*nall)))
rat=ksi/ksi1
ratk=sumk/sumk1
end
!---------------------------------------------
subroutine verletlj(nall,nallt,nallf,nalltf,parlj,parter,parsc,rx,ry,rz,rxf,ryf,rzf,vx,vy,vz,vxf,vyf,vzf,boxx,boxy,boxz,temp,curtemp,ksi,Elj,rf)
implicit none
real*8 :: rx,ry,rz,vx,vy,vz,rxf,ryf,rzf,vxf,vyf,vzf,axf,ayf,azf,axf1,ayf1,azf1,ax,ay,az,ax1,ay1,az1
real*8,dimension(2) :: parlj
real*8 :: temp,curtemp,ksi,Elj,ksi1,Et
real*8 :: boxx,boxy,boxz,dt,parter,parsc
real*8 :: Q,boltzct,tau,sumk,sumk1
real*8 :: masspc,masspau,rcutsc,rat,ratk,sumrf,rf,sumx,sumy,sumz
integer :: j,i
integer :: nall,nallt,nallf,nalltf
dimension rx(nall),ry(nall),rz(nall),rxf(nallf),ryf(nallf),rzf(nallf),vx(nallt),vy(nallt),vz(nallt),vxf(nalltf),vyf(nalltf),vzf(nalltf)
dimension ax(nall),ay(nall),az(nall),ax1(nall),ay1(nall),az1(nall),axf(nallf),ayf(nallf),azf(nallf),axf1(nallf),ayf1(nallf),azf1(nallf)
real*8 :: Fljx,Fljy,Fljz,Fljfx,Fljfy,Fljfz,Fterx,Ftery,Fterz,Fscx,Fscy,Fscz
dimension Fljx(nall),Fljy(nall),Fljz(nall),Fljfx(nallf),Fljfy(nallf),Fljfz(nallf),Fterx(nallf),Ftery(nallf),Fterz(nallf),Fscx(nall),Fscy(nall),Fscz(nall)
dimension parter(13),parsc(5),sumrf(3),rf(3)
masspc=12.0107/(6.022e26)
masspau=196.97/(6.022e26)
dt=1.e-15
tau=10*dt
boltzct=8.617332478e-5
rcutsc=2*4.08
call lennard(Elj,Fljx,Fljy,Fljz,Fljfx,Fljfy,Fljfz,rx,ry,rz,rxf,ryf,rzf,parlj,nall,nallf)
call tersoff(Et,Fterx,Ftery,Fterz,rxf,ryf,rzf,parter,nallf)
call SCforce(nall,boxx,boxy,boxz,rcutsc,rx,ry,rz,parsc,Fscx,Fscy,Fscz,Esc)
call acca(nallt,nalltf,Fterx,Ftery,Fterz,Fljx,Fljy,Fljz,Fljfx,Fljfy,Fljfz,Fscx,Fscy,Fscz,ax1,ay1,az1,axf1,ayf1,azf1,masspc,masspau)
Q=3*(nallt+nalltf-1)*boltzct*temp*tau**2
do i=1,nallt
	rx(i)=rx(i)+vx(i)*dt+0.5*ax1(i)*dt**2-0.5*dt**2*ksi*vx(i)
	ry(i)=ry(i)+vy(i)*dt+0.5*ay1(i)*dt**2-0.5*dt**2*ksi*vy(i)
	rz(i)=rz(i)+vz(i)*dt+0.5*az1(i)*dt**2-0.5*dt**2*ksi*vz(i)
end do
sumrf(:)=0.
do i=1,nalltf
	rxf(i)=rxf(i)+vxf(i)*dt+0.5*axf1(i)*dt**2-0.5*dt**2*ksi*vxf(i)
	ryf(i)=ryf(i)+vyf(i)*dt+0.5*ayf1(i)*dt**2-0.5*dt**2*ksi*vyf(i)
	rzf(i)=rzf(i)+vzf(i)*dt+0.5*azf1(i)*dt**2-0.5*dt**2*ksi*vzf(i)
	sumrf(1)=rxf(i)+sumrf(1)
	sumrf(2)=ryf(i)+sumrf(2)
	sumrf(3)=rzf(i)+sumrf(3)
end do
do i=1,3
rf(i)=1./dble(nalltf)*sumrf(i)
end do
ksi1=ksi
sumk=0.
sumk1=0.
do i=1,nallt
	sumk=sumk+masspau*(vx(i)+dt*0.5*(ax1(i)-ksi1*vx(i)))**2+masspau*(vy(i)+dt*0.5*(ay1(i)-ksi1*vy(i)))**2+ &
	masspau*(vz(i)+dt*0.5*(az1(i)-ksi1*vz(i)))**2
	sumk1=sumk1+masspau*(vx(i)**2+vy(i)**2+vz(i)**2)
end do
do i=1,nalltf
	sumk=sumk+masspc*(vxf(i)+dt*0.5*(axf1(i)-ksi1*vxf(i)))**2+masspc*(vyf(i)+dt*0.5*(ayf1(i)-ksi1*vyf(i)))**2+ &
	masspc*(vzf(i)+dt*0.5*(azf1(i)-ksi1*vzf(i)))**2
	sumk1=sumk1+masspc*(vxf(i)**2+vyf(i)**2+vzf(i)**2)
end do
sumk=sumk*6.24181*0.01
sumk1=sumk1*6.24181*0.01
ksi=ksi+dt/(2*Q)*(sumk1-3*(nallt+nalltf-1)*boltzct*temp)+dt/(2*Q)*(sumk-3*(nallt+nalltf-1)*boltzct*temp)
do j=1,nall
		rx(j)=rx(j)-boxx*anint(rx(j)/boxx)
		ry(j)=ry(j)-boxy*anint(ry(j)/boxy)
end do
do j=1,nallf
		rxf(j)=rxf(j)-boxx*anint(rxf(j)/boxx)
		ryf(j)=ryf(j)-boxy*anint(ryf(j)/boxy)
		rzf(j)=rzf(j)-boxz*anint(rzf(j)/boxz)
end do
call lennard(Elj,Fljx,Fljy,Fljz,Fljfx,Fljfy,Fljfz,rx,ry,rz,rxf,ryf,rzf,parlj,nall,nallf)
call tersoff(Et,Fterx,Ftery,Fterz,rxf,ryf,rzf,parter,nallf)
call SCforce(nall,boxx,boxy,boxz,rcutsc,rx,ry,rz,parsc,Fscx,Fscy,Fscz,Esc)
call acca(nallt,nalltf,Fterx,Ftery,Fterz,Fljx,Fljy,Fljz,Fljfx,Fljfy,Fljfz,Fscx,Fscy,Fscz,ax,ay,az,axf,ayf,azf,masspc,masspau)
do i=1,nallt
	vx(i)=2./(2+ksi*dt)*(vx(i)+0.5*dt*(ax1(i)-ksi1*vx(i))+0.5*dt*ax(i))
	vy(i)=2./(2+ksi*dt)*(vy(i)+0.5*dt*(ay1(i)-ksi1*vy(i))+0.5*dt*ay(i))
	vz(i)=2./(2+ksi*dt)*(vz(i)+0.5*dt*(az1(i)-ksi1*vz(i))+0.5*dt*az(i))
end do
do i=1,nalltf
	vxf(i)=2./(2+ksi*dt)*(vxf(i)+0.5*dt*(axf1(i)-ksi1*vxf(i))+0.5*dt*axf(i))
	vyf(i)=2./(2+ksi*dt)*(vyf(i)+0.5*dt*(ayf1(i)-ksi1*vyf(i))+0.5*dt*ayf(i))
	vzf(i)=2./(2+ksi*dt)*(vzf(i)+0.5*dt*(azf1(i)-ksi1*vzf(i))+0.5*dt*azf(i))
end do

q=ax(1)/ax1(1)
curtemp=sumk1/(boltzct*dble(3*(nallt+nalltf)))
rat=ksi/ksi1
ratk=sumk/sumk1
end

!--------------------------------
subroutine lennard(Elj,Fljx,Fljy,Fljz,Fljfx,Fljfy,Fljfz,rx,ry,rz,rxf,ryf,rzf,parlj,nall,nallf)
implicit none
real*8 :: rx,ry,rz,rxf,ryf,rzf
real*8 :: Fljx,Fljy,Fljz,Fljfx,Fljfy,Fljfz
real*8 :: Elj,parlj,dx,dy,dz,rs,rcut
integer :: nall,nallf,i,j
dimension rx(nall),ry(nall),rz(nall),rxf(nallf),ryf(nallf),rzf(nallf)
dimension Fljx(nall),fljy(nall),fljz(nall),fljfx(nallf),fljfy(nallf),fljfz(nallf),parlj(2)

eps=parlj(1);sigma=parlj(2)
rcut=15.
Elj=0;Fljx(:)=0;Fljy(:)=0;Fljz(:)=0;Fljfx(:)=0;Fljfy(:)=0;Fljfz(:)=0
Do i=1,nall
	Do j=1,nallf
		dx=rx(i)-rxf(j)
		dy=ry(i)-ryf(j)
		dz=rz(i)-rzf(j)
        dx=dx-boxx*anint(dx/boxx)
		dy=dy-boxy*anint(dy/boxy)
		rs=sqrt(dx*dx+dy*dy+dz*dz)
		if (rs<=rcut) then
		Elj=Elj-eps*((sigma/rs)**12-(sigma/rs)**6)
		Fljx(i)=Fljx(i)+4.*eps/rs**2*dx*(12*(sigma/rs)**12-6*(sigma/rs)**6)
		Fljy(i)=Fljy(i)+4.*eps/rs**2*dy*(12*(sigma/rs)**12-6*(sigma/rs)**6)
		Fljz(i)=Fljz(i)+4.*eps/rs**2*dz*(12*(sigma/rs)**12-6*(sigma/rs)**6)
		Fljfx(j)=Fljfx(j)-4.*eps/rs**2*dx*(12*(sigma/rs)**12-6*(sigma/rs)**6)
		Fljfy(j)=Fljfy(j)-4.*eps/rs**2*dy*(12*(sigma/rs)**12-6*(sigma/rs)**6)
		Fljfz(j)=Fljfz(j)-4.*eps/rs**2*dz*(12*(sigma/rs)**12-6*(sigma/rs)**6)
		end if
	end do
end do
end
!---------------------------------------------
subroutine acca(nallt,nalltf,Fterx,Ftery,Fterz,Fljx,Fljy,Fljz,Fljfx,Fljfy,Fljfz,Fscx,Fscy,Fscz,ax,ay,az,axf,ayf,azf,masspc,masspau)
implicit none
real*8 :: Fterx,Ftery,Fterz,Fljx,Fljy,Fljz,Fljfx,Fljfy,Fljfz,Fscx,Fscy,Fscz,ax,ay,az,axf,ayf,azf
real*8 :: masspc,masspau
dimension Fterx(nalltf),Ftery(nalltf),Fterz(nalltf),Fljx(nallt),Fljy(nallt),Fljz(nallt),Fscx(nallt),Fscy(nallt),Fscz(nallt),Fljfx(nalltf),Fljfy(nalltf),Fljfz(nalltf)
dimension ax(nallt),ay(nallt),az(nallt),axf(nalltf),ayf(nalltf),azf(nalltf)
integer :: nallt,nalltf
ax(:)=0;ay(:)=0;az(:)=0;axf(:)=0;ayf(:)=0;azf(:)=0
do i=1,nallt
ax(i)=ax(i)+16.0210*Fscx(i)/masspau+16.0210*Fljx(i)/masspau
ay(i)=ay(i)+16.0210*Fscy(i)/masspau+16.0210*Fljy(i)/masspau
az(i)=az(i)+16.0210*Fscz(i)/masspau+16.0210*Fljz(i)/masspau
end do
do i=1,nalltf
axf(i)=axf(i)+16.0210*Fterx(i)/masspc+16.0210*Fljfx(i)/masspc
ayf(i)=ayf(i)+16.0210*Ftery(i)/masspc+16.0210*Fljfy(i)/masspc
azf(i)=azf(i)+16.0210*Fterz(i)/masspc+16.0210*Fljfz(i)/masspc
end do
end

!----------------------------
subroutine MD(rx,ry,rz,rxf,ryf,rzf,nall,nallt,nallf,nalltf,temp,parsc,parter,parlj,acell,rcut,nmd,boxx,boxy,boxz)
implicit none
real*8 :: rx,ry,rz,vx,vy,vz,rxf,ryf,rzf,vxf,vyf,vzf
real*8 :: temp,parsc,parter,parlj,acell,rcut,curtemp,massp,ksi,Esc,Et,Elj,sumb,sumx,sumy,sumz,masspau,masspc
real*8 :: boxx,boxy,boxz,rf
real*8 :: dt,qa
integer :: nall,nallt,i,nmd,nallf,nalltf,j
dimension rx(nall),ry(nall),rz(nall),parsc(5),parter(13),parlj(2),vx(nallt),vy(nallt),vz(nallt)
dimension rxf(nallf),ryf(nallf),rzf(nallf),vxf(nallf),vyf(nallf),vzf(nallf),rf(3)

!Relaxation of Gold
massp=196.97/(6.022e26) !Mass in Kg
masspau=massp
ksi=0.
call strtup(vx,vy,vz,nallt,temp,massp)
open(55,file='Energy_step_gold.dat',status='unknown')
open(110,file='Temp_step_gold.dat',status='unknown')
do i=1,5000

call verletgold(nall,nallt,parsc,massp,rx,ry,rz,vx,vy,vz,dt,boxx,boxy,boxz,temp,curtemp,ksi,Esc)
write(*,*) "Substrate temperature is", curtemp
write(55,*) i,Esc
write(110,*) i,curtemp
end do
close(55)
close(110)
!!Relaxation of Fullerene
massp=12.0107/(6.022e26) !Mass in Kg
masspc=massp
ksi=0.
call strtup(vxf,vyf,vzf,nalltf,temp,massp)
open(55,file='Energy_step_f.dat',status='unknown')
open(110,file='Temp_step_f.dat',status='unknown')
do i=1,10000
call verletful(nallf,nalltf,parter,massp,rxf,ryf,rzf,vxf,vyf,vzf,dt,boxx,boxy,boxz,temp,curtemp,ksi,Et,sumb)
write(*,*) "Fullerene temperature is", curtemp
write(55,*) i,Et
write(110,*) i,curtemp
end do
close(55)
close(110)
qa=minval(rzf)-maxval(rz)
do i=1,nallf
rzf(i)=rzf(i)-qa+3.
end do
qa=minval(rzf)-maxval(rz)
qa=minval(rxf)-minval(rx)
qa=minval(ryf)-minval(ry)
!-----------------------------------------------------
open(77,file='Full_traj.dat',status='unknown')
open(55,file='Energy_step.dat',status='unknown')
open(110,file='Temp_step.dat',status='unknown')
do i=1,20000
call verletlj(nall,nallt,nallf,nalltf,parlj,parter,parsc,rx,ry,rz,rxf,ryf,rzf,vx,vy,vz,vxf,vyf,vzf,boxx,boxy,boxz,temp,curtemp,ksi,Elj,rf)
write(*,*) "System temperature is", curtemp
write(77,*) i,rf(1),rf(2),rf(3)
write(55,*) i,Elj
write(110,*) i,curtemp
end do
close(55)
close(110)
close(77)
end
!------------------------------------

subroutine potin(massp,munit,eunit,lunit,dt,temp,nall,tmp)

implicit none

real*8 :: massp,munit,eunit,lunit,dt,temp,tmp
real*8 :: time
integer :: nall

munit=196.96657
lunit=4.08
eunit=1.

massp=196.96657/munit
time=munit*lunit*lunit/eunit
time=sqrt(time*0.1/6.022/1.602)*1.e-13
dt=dt/time

tmp=eunit*1.602e-19/1.381e-23
temp=temp/tmp

end

      FUNCTION ZBQLNOR(MU,SIGMA)
!
!       Returns a random number Normally distributed with mean
!       MU and standard deviation |SIGMA|, using the Box-Muller
!       algorithm
!
!OBTAINED FROM http://www.ucl.ac.uk/~ucakarc/work/software/randgen.f
! Andreas Kloeckner, 7/5/2006
! FROM http://www.ucl.ac.uk/~ucakarc/work/software/randgen.txt
      DOUBLE PRECISION THETA,R,ZBQLNOR,ZBQLU01,PI,MU,SIGMA
      DOUBLE PRECISION SPARE
      INTEGER STATUS
      SAVE STATUS,SPARE,PI
      DATA STATUS /-1/
	  SIGMA=1.
	  MU=0.
      IF (STATUS.EQ.-1) PI = 4.0D0*DATAN(1.0D0)

      IF (STATUS.LE.0) THEN
       THETA = 2.0D0*PI*ZBQLU01(0.0D0)
       R = DSQRT( -2.0D0*DLOG(ZBQLU01(0.0D0)) )
       ZBQLNOR = (R*DCOS(THETA))
       SPARE = (R*DSIN(THETA))
       STATUS = 1
      ELSE
       ZBQLNOR = SPARE
       STATUS = 0
      ENDIF
      
      ZBQLNOR = MU + (SIGMA*ZBQLNOR)

      END
	  
	   FUNCTION ZBQLU01(DUMMY)
!
!       Returns a uniform random number between 0 & 1, using
!       a Marsaglia-Zaman type subtract-with-borrow generator.
!       Uses double precision, rather than integer, arithmetic 
!       throughout because MZ's integer constants overflow
!       32-bit integer storage (which goes from -2^31 to 2^31).
!       Ideally, we would explicitly truncate all integer 
!       quantities at each stage to ensure that the double
!       precision representations do not accumulate approximation
!       error; however, on some machines the use of DNINT to
!       accomplish this is *seriously* slow (run-time increased
!       by a factor of about 3). This double precision version 
!       has been tested against an integer implementation that
!       uses long integers (non-standard and, again, slow) -
!       the output was identical up to the 16th decimal place
!       after 10^10 calls, so we're probably OK ...
!
!OBTAINED FROM http://www.ucl.ac.uk/~ucakarc/work/software/randgen.f
! Andreas Kloeckner, 7/5/2006
! FROM http://www.ucl.ac.uk/~ucakarc/work/software/randgen.txt

      DOUBLE PRECISION ZBQLU01,DUMMY,B,C,ZBQLIX(43),X,B2,BINV
      INTEGER CURPOS,ID22,ID43

      COMMON /ZBQL0001/ ZBQLIX,B,C
      SAVE /ZBQL0001/
      SAVE CURPOS,ID22,ID43
      DATA CURPOS,ID22,ID43 /1,22,43/

      B2 = B
      BINV = 1.0D0/B
 5    X = ZBQLIX(ID22) - ZBQLIX(ID43) - C
      IF (X.LT.0.0D0) THEN
       X = X + B
       C = 1.0D0
      ELSE
       C = 0.0D0
      ENDIF
      ZBQLIX(ID43) = X
!
!     Update array pointers. Do explicit check for bounds of each to
!     avoid expense of modular arithmetic. If one of them is 0 the others
!     won't be
!
      CURPOS = CURPOS - 1
      ID22 = ID22 - 1
      ID43 = ID43 - 1
      IF (CURPOS.EQ.0) THEN
       CURPOS=43
      ELSEIF (ID22.EQ.0) THEN
       ID22 = 43
      ELSEIF (ID43.EQ.0) THEN
       ID43 = 43
      ENDIF
!
!     The integer arithmetic there can yield X=0, which can cause 
!     problems in subsequent routines (e.g. ZBQLEXP). The problem
!     is simply that X is discrete whereas U is supposed to 
!     be continuous - hence if X is 0, go back and generate another
!     X and return X/B^2 (etc.), which will be uniform on (0,1/B). 
!
      IF (X.LT.BINV) THEN
       B2 = B2*B
       GOTO 5
      ENDIF

      ZBQLU01 = X/B2

      END
 
end program
!----------------------------------------
