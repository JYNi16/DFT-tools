Program main
implicit none
integer, parameter :: N=100
real, parameter :: pi=3.1415926
integer :: i, j, k(4), nx, ny, nz, num(2), e(4), r1, r2, nt, n1, n2, x, y, z, k1, k2, k3, &
s1, s2, s3, num1, num2
real(8) :: a0, ap1, bp1, cp1, ap2, bp2, cp2, ap, bp, cp, fr(N,3), fr1(N,3), &
fr2(N,3), fr3(N,3), lp(3,3), lp1(3,3), lp2(3,3), lp3(3,3), car(10,N,3), car1(10,N,3), &
car2(10,N,3), car3(10,N,3), ang1, ang2, ang3, dist2(3), dtm, do1, Cri, m(3,3), mr, ang(3), hei(2), &
carp(5,3), carp1(10000,3), sc(6,2000,3), frt(10000,3)

write(*,*) "Please input the pseudocubic lattic parameter a0!"
!read(*,*) a0
a0=3.8
write(*,*) "Please input the supercell that will be created!"
!do i=1, 3
!read(*,*) m(i,1), m(i,2), m(i,3)
!end do
m(1,1)=2; m(1,2)=0; m(1,3)=0
m(2,1)=0; m(2,2)=2; m(2,3)=0
m(3,1)=0; m(3,2)=0; m(3,3)=2

write(*,*) "Please input the rotational mode for the three directions (x,y,z): 1 for M and -1 for R!"
read(*,*) nx, ny, nz

write(*,*) "Please input the ang1 by which the Oct. rotates around x axis, &
the ang2 for y axis, ang3 for z axis, dist2 (three components) for the &
movements of B ion!"
read(*,*) ang(1), ang(2), ang(3), dist2(1), dist2(2), dist2(3)

write(*,*) "Please input the ratio of two compounds!"
!read(*,*) r1, r2
r1=1;r2=1
write(*,*) "Please input the random numbers indicating the different elements at A and B positions!"
!read(*,*) e(1), e(2), e(3), e(4)
e(1)=1; e(2)=1; e(3)=2; e(4)=2
!Creating supercell
fr(1,1)=0.0
fr(1,2)=0.0
fr(1,3)=0.0
fr(2,1)=0.5
fr(2,2)=0.5
fr(2,3)=0.5
fr(3,1)=0.0
fr(3,2)=0.5
fr(3,3)=0.5
fr(4,1)=0.5
fr(4,2)=0.0
fr(4,3)=0.5
fr(5,1)=0.5
fr(5,2)=0.5
fr(5,3)=0.0

do i=1, 5
carp(i,1)=fr(i,1)*a0
carp(i,2)=fr(i,2)*a0
carp(i,3)=fr(i,3)*a0
end do

s1=9
s2=9
s3=20

do k1=1, 5
num2=0
do z=1, s3
do y=1, s2
do x=1, s1
num2=num2+1
if(MOD(z,2).NE.0) then
if(MOD(y,2).NE.0) then
sc(k1,num2,1)=carp(k1,1)+(x-1)*a0
sc(k1,num2,2)=carp(k1,2)+(y-1)*a0
sc(k1,num2,3)=carp(k1,3)+(z-1)*a0
else if(MOD(y,2)==0) then
k2=s1-x+1
sc(k1,num2,1)=carp(k1,1)+(k2-1)*a0
sc(k1,num2,2)=carp(k1,2)+(y-1)*a0
sc(k1,num2,3)=carp(k1,3)+(z-1)*a0
end if
else if(MOD(z,2)==0) then
k3=s2-y+1
if(MOD(y,2).NE.0) then
k2=s1-x+1
sc(k1,num2,1)=carp(k1,1)+(k2-1)*a0
sc(k1,num2,2)=carp(k1,2)+(k3-1)*a0
sc(k1,num2,3)=carp(k1,3)+(z-1)*a0
else if(MOD(y,2)==0) then
sc(k1,num2,1)=carp(k1,1)+(x-1)*a0
sc(k1,num2,2)=carp(k1,2)+(y-1)*a0
sc(k1,num2,3)=carp(k1,3)+(z-1)*a0
end if
end if
end do
end do
end do
end do

num2=0
do k1=1, 5
do i=1, s1*s2*s3
num2=num2+1
carp1(num2,1)=sc(k1,i,1)-sc(1,345,1)
carp1(num2,2)=sc(k1,i,2)-sc(1,345,2)
carp1(num2,3)=sc(k1,i,3)-sc(1,345,3)
end do
end do

do i=1, 3
lp1(i,1)=m(i,1)*a0
lp1(i,2)=m(i,2)*a0
lp1(i,3)=m(i,3)*a0
end do

mr=lp1(1,1)*lp1(2,2)*lp1(3,3)+lp1(1,2)*lp1(2,3)*lp1(3,1)+lp1(1,3)*lp1(3,2)*lp1(2,1)&
-lp1(1,3)*lp1(2,2)*lp1(3,1)-lp1(1,1)*lp1(3,2)*lp1(2,3)-lp1(1,2)*lp1(2,1)*lp1(3,3)

do i=1, num2

frt(i,1)=(carp1(i,1)*(lp1(2,2)*lp1(3,3)-lp1(2,3)*lp1(3,2))+carp1(i,2)*(-lp1(2,1)*&
lp1(3,3)+lp1(2,3)*lp1(3,1))+carp1(i,3)*(lp1(2,1)*lp1(3,2)-lp1(2,2)*lp1(3,1)))/mr
frt(i,2)=(carp1(i,1)*(-lp1(1,2)*lp1(3,3)+lp1(1,3)*lp1(3,2))+carp1(i,2)*(lp1(1,1)*&
lp1(3,3)-lp1(1,3)*lp1(3,1))+carp1(i,3)*(-lp1(1,1)*lp1(3,2)+lp1(1,2)*lp1(3,1)))/mr
frt(i,3)=(carp1(i,1)*(lp1(1,2)*lp1(2,3)-lp1(1,3)*lp1(2,2))+carp1(i,2)*(-lp1(1,1)*&
lp1(2,3)+lp1(1,3)*lp1(2,1))+carp1(i,3)*(lp1(1,1)*lp1(2,2)-lp1(1,2)*lp1(2,1)))/mr

end do


num1=0
do i=1, num2
if(frt(i,1)>-0.001.and.frt(i,1)<0.99) then
if(frt(i,2)>-0.001.and.frt(i,2)<0.99) then
if(frt(i,3)>-0.001.and.frt(i,3)<0.99) then
num1=num1+1
fr(num1,1)=frt(i,1)
fr(num1,2)=frt(i,2)
fr(num1,3)=frt(i,3)
end if
end if
end if
end do
!End

nt=num1
n1=sqrt(m(1,1)**2+m(1,2)**2+m(1,3)**2)*sqrt(m(2,1)**2+m(2,2)**2+m(2,3)**2)*sqrt(m(3,1)**2+m(3,2)**2+m(3,3)**2)+1
n2=(n1-1)*2

dtm=a0 !tm means 'transition metal'.
do1=a0/2.0 !do1 means 'tm-oxygen'.
Cri=0.21 ! The criterion to judge the distance between tm pair.

ap=a0*cos(2*pi*ang(2)/360)*cos(2*pi*ang(3)/360)
bp=a0*cos(2*pi*ang(1)/360)*cos(2*pi*ang(3)/360)
cp=a0*cos(2*pi*ang(1)/360)*cos(2*pi*ang(2)/360)

ang1=do1*tan(2*pi*ang(1)/360)
ang2=do1*tan(2*pi*ang(2)/360)
ang3=do1*tan(2*pi*ang(3)/360)

do i=1, nt
car(1,i,1)=fr(i,1)*lp1(1,1)+fr(i,2)*lp1(2,1)+fr(i,3)*lp1(3,1)
car(1,i,2)=fr(i,1)*lp1(1,2)+fr(i,2)*lp1(2,2)+fr(i,3)*lp1(3,2)
car(1,i,3)=fr(i,1)*lp1(1,3)+fr(i,2)*lp1(2,3)+fr(i,3)*lp1(3,3)
end do

do i=n2+1, nt
car(2,i,1)=(fr(i,1)+1)*lp1(1,1)+fr(i,2)*lp1(2,1)+fr(i,3)*lp1(3,1)
car(2,i,2)=(fr(i,1)+1)*lp1(1,2)+fr(i,2)*lp1(2,2)+fr(i,3)*lp1(3,2)
car(2,i,3)=(fr(i,1)+1)*lp1(1,3)+fr(i,2)*lp1(2,3)+fr(i,3)*lp1(3,3)
end do

do i=n2+1, nt
car(3,i,1)=fr(i,1)*lp1(1,1)+(fr(i,2)+1)*lp1(2,1)+fr(i,3)*lp1(3,1)
car(3,i,2)=fr(i,1)*lp1(1,2)+(fr(i,2)+1)*lp1(2,2)+fr(i,3)*lp1(3,2)
car(3,i,3)=fr(i,1)*lp1(1,3)+(fr(i,2)+1)*lp1(2,3)+fr(i,3)*lp1(3,3)
end do

do i=n2+1, nt
car(4,i,1)=fr(i,1)*lp1(1,1)+fr(i,2)*lp1(2,1)+(fr(i,3)+1)*lp1(3,1)
car(4,i,2)=fr(i,1)*lp1(1,2)+fr(i,2)*lp1(2,2)+(fr(i,3)+1)*lp1(3,2)
car(4,i,3)=fr(i,1)*lp1(1,3)+fr(i,2)*lp1(2,3)+(fr(i,3)+1)*lp1(3,3)
end do

do i=n2+1, nt
car(5,i,1)=(fr(i,1)-1)*lp1(1,1)+fr(i,2)*lp1(2,1)+fr(i,3)*lp1(3,1)
car(5,i,2)=(fr(i,1)-1)*lp1(1,2)+fr(i,2)*lp1(2,2)+fr(i,3)*lp1(3,2)
car(5,i,3)=(fr(i,1)-1)*lp1(1,3)+fr(i,2)*lp1(2,3)+fr(i,3)*lp1(3,3)
end do

do i=n2+1, nt
car(6,i,1)=fr(i,1)*lp1(1,1)+(fr(i,2)-1)*lp1(2,1)+fr(i,3)*lp1(3,1)
car(6,i,2)=fr(i,1)*lp1(1,2)+(fr(i,2)-1)*lp1(2,2)+fr(i,3)*lp1(3,2)
car(6,i,3)=fr(i,1)*lp1(1,3)+(fr(i,2)-1)*lp1(2,3)+fr(i,3)*lp1(3,3)
end do

do i=n2+1, nt
car(7,i,1)=fr(i,1)*lp1(1,1)+fr(i,2)*lp1(2,1)+(fr(i,3)-1)*lp1(3,1)
car(7,i,2)=fr(i,1)*lp1(1,2)+fr(i,2)*lp1(2,2)+(fr(i,3)-1)*lp1(3,2)
car(7,i,3)=fr(i,1)*lp1(1,3)+fr(i,2)*lp1(2,3)+(fr(i,3)-1)*lp1(3,3)
end do

!Rotation around the x axis
do i=n1, n2

if(i/=n1) then
if(nx==1) then
if(abs(sqrt((car(1,i,1)-car(1,i-1,1))**2+(car(1,i,2)-car(1,i-1,2))**2+(car(1,i,3)-car(1,i-1,3))**2)-dtm)<=Cri) then
if(abs(car(1,i,1)-car(1,i-1,1))<=Cri) then

ang1=-ang1

end if
end if

else if(nx==-1) then
if(abs(sqrt((car(1,i,1)-car(1,i-1,1))**2+(car(1,i,2)-car(1,i-1,2))**2+(car(1,i,3)-car(1,i-1,3))**2)-dtm)<=Cri) then

ang1=-ang1

end if

end if

end if

do j=n2+1, nt

if(abs(sqrt((car(1,j,1)-car(1,i,1))**2+(car(1,j,2)-car(1,i,2))**2+(car(1,j,3)-car(1,i,3))**2)-do1)<=Cri) then

if((car(1,j,3)-car(1,i,3))<-1.0) then
car1(1,j,1)=car(1,j,1)
car1(1,j,2)=car(1,j,2)+ang1
car1(1,j,3)=car(1,j,3)
else if((car(1,j,3)-car(1,i,3))>1.0) then
car1(1,j,1)=car(1,j,1)
car1(1,j,2)=car(1,j,2)-ang1
car1(1,j,3)=car(1,j,3)
else if((car(1,j,1)-car(1,i,1))>1.0) then
car1(1,j,1)=car(1,j,1)
car1(1,j,2)=car(1,j,2)
car1(1,j,3)=car(1,j,3)
else if((car(1,j,1)-car(1,i,1))<-1.0) then
car1(1,j,1)=car(1,j,1)
car1(1,j,2)=car(1,j,2)
car1(1,j,3)=car(1,j,3)
else if((car(1,j,2)-car(1,i,2))<-1.0) then
car1(1,j,1)=car(1,j,1)
car1(1,j,2)=car(1,j,2)
car1(1,j,3)=car(1,j,3)-ang1
else if((car(1,j,2)-car(1,i,2))>1.0) then
car1(1,j,1)=car(1,j,1)
car1(1,j,2)=car(1,j,2)
car1(1,j,3)=car(1,j,3)+ang1
end if

end if

if(abs(sqrt((car(2,j,1)-car(1,i,1))**2+(car(2,j,2)-car(1,i,2))**2+(car(2,j,3)-car(1,i,3))**2)-do1)<=Cri) then

if((car(2,j,3)-car(1,i,3))<-1.0) then
car1(1,j,1)=car(1,j,1)
car1(1,j,2)=car(1,j,2)+ang1
car1(1,j,3)=car(1,j,3)
else if((car(2,j,3)-car(1,i,3))>1.0) then
car1(1,j,1)=car(1,j,1)
car1(1,j,2)=car(1,j,2)-ang1
car1(1,j,3)=car(1,j,3)
else if((car(2,j,1)-car(1,i,1))>1.0) then
car1(1,j,1)=car(1,j,1)
car1(1,j,2)=car(1,j,2)
car1(1,j,3)=car(1,j,3)
else if((car(2,j,1)-car(1,i,1))<-1.0) then
car1(1,j,1)=car(1,j,1)
car1(1,j,2)=car(1,j,2)
car1(1,j,3)=car(1,j,3)
else if((car(2,j,2)-car(1,i,2))<-1.0) then
car1(1,j,1)=car(1,j,1)
car1(1,j,2)=car(1,j,2)
car1(1,j,3)=car(1,j,3)-ang1
else if((car(2,j,2)-car(1,i,2))>1.0) then
car1(1,j,1)=car(1,j,1)
car1(1,j,2)=car(1,j,2)
car1(1,j,3)=car(1,j,3)+ang1
end if

end if

if(abs(sqrt((car(3,j,1)-car(1,i,1))**2+(car(3,j,2)-car(1,i,2))**2+(car(3,j,3)-car(1,i,3))**2)-do1)<=Cri) then

if((car(3,j,3)-car(1,i,3))<-1.0) then
car1(1,j,1)=car(1,j,1)
car1(1,j,2)=car(1,j,2)+ang1
car1(1,j,3)=car(1,j,3)
else if((car(3,j,3)-car(1,i,3))>1.0) then
car1(1,j,1)=car(1,j,1)
car1(1,j,2)=car(1,j,2)-ang1
car1(1,j,3)=car(1,j,3)
else if((car(3,j,1)-car(1,i,1))>1.0) then
car1(1,j,1)=car(1,j,1)
car1(1,j,2)=car(1,j,2)
car1(1,j,3)=car(1,j,3)
else if((car(3,j,1)-car(1,i,1))<-1.0) then
car1(1,j,1)=car(1,j,1)
car1(1,j,2)=car(1,j,2)
car1(1,j,3)=car(1,j,3)
else if((car(3,j,2)-car(1,i,2))<-1.0) then
car1(1,j,1)=car(1,j,1)
car1(1,j,2)=car(1,j,2)
car1(1,j,3)=car(1,j,3)-ang1
else if((car(3,j,2)-car(1,i,2))>1.0) then
car1(1,j,1)=car(1,j,1)
car1(1,j,2)=car(1,j,2)
car1(1,j,3)=car(1,j,3)+ang1
end if

end if

if(abs(sqrt((car(4,j,1)-car(1,i,1))**2+(car(4,j,2)-car(1,i,2))**2+(car(4,j,3)-car(1,i,3))**2)-do1)<=Cri) then

if((car(4,j,3)-car(1,i,3))<-1.0) then
car1(1,j,1)=car(1,j,1)
car1(1,j,2)=car(1,j,2)+ang1
car1(1,j,3)=car(1,j,3)
else if((car(4,j,3)-car(1,i,3))>1.0) then
car1(1,j,1)=car(1,j,1)
car1(1,j,2)=car(1,j,2)-ang1
car1(1,j,3)=car(1,j,3)
else if((car(4,j,1)-car(1,i,1))>1.0) then
car1(1,j,1)=car(1,j,1)
car1(1,j,2)=car(1,j,2)
car1(1,j,3)=car(1,j,3)
else if((car(4,j,1)-car(1,i,1))<-1.0) then
car1(1,j,1)=car(1,j,1)
car1(1,j,2)=car(1,j,2)
car1(1,j,3)=car(1,j,3)
else if((car(4,j,2)-car(1,i,2))<-1.0) then
car1(1,j,1)=car(1,j,1)
car1(1,j,2)=car(1,j,2)
car1(1,j,3)=car(1,j,3)-ang1
else if((car(4,j,2)-car(1,i,2))>1.0) then
car1(1,j,1)=car(1,j,1)
car1(1,j,2)=car(1,j,2)
car1(1,j,3)=car(1,j,3)+ang1
end if

end if

if(abs(sqrt((car(5,j,1)-car(1,i,1))**2+(car(5,j,2)-car(1,i,2))**2+(car(5,j,3)-car(1,i,3))**2)-do1)<=Cri) then

if((car(5,j,3)-car(1,i,3))<-1.0) then
car1(1,j,1)=car(1,j,1)
car1(1,j,2)=car(1,j,2)+ang1
car1(1,j,3)=car(1,j,3)
else if((car(5,j,3)-car(1,i,3))>1.0) then
car1(1,j,1)=car(1,j,1)
car1(1,j,2)=car(1,j,2)-ang1
car1(1,j,3)=car(1,j,3)
else if((car(5,j,1)-car(1,i,1))>1.0) then
car1(1,j,1)=car(1,j,1)
car1(1,j,2)=car(1,j,2)
car1(1,j,3)=car(1,j,3)
else if((car(5,j,1)-car(1,i,1))<-1.0) then
car1(1,j,1)=car(1,j,1)
car1(1,j,2)=car(1,j,2)
car1(1,j,3)=car(1,j,3)
else if((car(5,j,2)-car(1,i,2))<-1.0) then
car1(1,j,1)=car(1,j,1)
car1(1,j,2)=car(1,j,2)
car1(1,j,3)=car(1,j,3)-ang1
else if((car(5,j,2)-car(1,i,2))>1.0) then
car1(1,j,1)=car(1,j,1)
car1(1,j,2)=car(1,j,2)
car1(1,j,3)=car(1,j,3)+ang1
end if

end if

if(abs(sqrt((car(6,j,1)-car(1,i,1))**2+(car(6,j,2)-car(1,i,2))**2+(car(6,j,3)-car(1,i,3))**2)-do1)<=Cri) then

if((car(6,j,3)-car(1,i,3))<-1.0) then
car1(1,j,1)=car(1,j,1)
car1(1,j,2)=car(1,j,2)+ang1
car1(1,j,3)=car(1,j,3)
else if((car(6,j,3)-car(1,i,3))>1.0) then
car1(1,j,1)=car(1,j,1)
car1(1,j,2)=car(1,j,2)-ang1
car1(1,j,3)=car(1,j,3)
else if((car(6,j,1)-car(1,i,1))>1.0) then
car1(1,j,1)=car(1,j,1)
car1(1,j,2)=car(1,j,2)
car1(1,j,3)=car(1,j,3)
else if((car(6,j,1)-car(1,i,1))<-1.0) then
car1(1,j,1)=car(1,j,1)
car1(1,j,2)=car(1,j,2)
car1(1,j,3)=car(1,j,3)
else if((car(6,j,2)-car(1,i,2))<-1.0) then
car1(1,j,1)=car(1,j,1)
car1(1,j,2)=car(1,j,2)
car1(1,j,3)=car(1,j,3)-ang1
else if((car(6,j,2)-car(1,i,2))>1.0) then
car1(1,j,1)=car(1,j,1)
car1(1,j,2)=car(1,j,2)
car1(1,j,3)=car(1,j,3)+ang1
end if

end if

if(abs(sqrt((car(7,j,1)-car(1,i,1))**2+(car(7,j,2)-car(1,i,2))**2+(car(7,j,3)-car(1,i,3))**2)-do1)<=Cri) then

if((car(7,j,3)-car(1,i,3))<-1.0) then
car1(1,j,1)=car(1,j,1)
car1(1,j,2)=car(1,j,2)+ang1
car1(1,j,3)=car(1,j,3)
else if((car(7,j,3)-car(1,i,3))>1.0) then
car1(1,j,1)=car(1,j,1)
car1(1,j,2)=car(1,j,2)-ang1
car1(1,j,3)=car(1,j,3)
else if((car(7,j,1)-car(1,i,1))>1.0) then
car1(1,j,1)=car(1,j,1)
car1(1,j,2)=car(1,j,2)
car1(1,j,3)=car(1,j,3)
else if((car(7,j,1)-car(1,i,1))<-1.0) then
car1(1,j,1)=car(1,j,1)
car1(1,j,2)=car(1,j,2)
car1(1,j,3)=car(1,j,3)
else if((car(7,j,2)-car(1,i,2))<-1.0) then
car1(1,j,1)=car(1,j,1)
car1(1,j,2)=car(1,j,2)
car1(1,j,3)=car(1,j,3)-ang1
else if((car(7,j,2)-car(1,i,2))>1.0) then
car1(1,j,1)=car(1,j,1)
car1(1,j,2)=car(1,j,2)
car1(1,j,3)=car(1,j,3)+ang1
end if

end if

end do
end do
!End

!The transformation from car to dir!
mr=lp1(1,1)*lp1(2,2)*lp1(3,3)+lp1(1,2)*lp1(2,3)*lp1(3,1)+lp1(1,3)*lp1(3,2)*lp1(2,1)&
-lp1(1,3)*lp1(2,2)*lp1(3,1)-lp1(1,1)*lp1(3,2)*lp1(2,3)-lp1(1,2)*lp1(2,1)*lp1(3,3)

do i=n2+1, nt

fr1(i,1)=(car1(1,i,1)*(lp1(2,2)*lp1(3,3)-lp1(2,3)*lp1(3,2))+car1(1,i,2)*(-lp1(2,1)*&
lp1(3,3)+lp1(2,3)*lp1(3,1))+car1(1,i,3)*(lp1(2,1)*lp1(3,2)-lp1(2,2)*lp1(3,1)))/mr
fr1(i,2)=(car1(1,i,1)*(-lp1(1,2)*lp1(3,3)+lp1(1,3)*lp1(3,2))+car1(1,i,2)*(lp1(1,1)*&
lp1(3,3)-lp1(1,3)*lp1(3,1))+car1(1,i,3)*(-lp1(1,1)*lp1(3,2)+lp1(1,2)*lp1(3,1)))/mr
fr1(i,3)=(car1(1,i,1)*(lp1(1,2)*lp1(2,3)-lp1(1,3)*lp1(2,2))+car1(1,i,2)*(-lp1(1,1)*&
lp1(2,3)+lp1(1,3)*lp1(2,1))+car1(1,i,3)*(lp1(1,1)*lp1(2,2)-lp1(1,2)*lp1(2,1)))/mr

end do
!End

do i=n2+1, nt
car1(2,i,1)=(fr1(i,1)+1)*lp1(1,1)+fr1(i,2)*lp1(2,1)+fr1(i,3)*lp1(3,1)
car1(2,i,2)=(fr1(i,1)+1)*lp1(1,2)+fr1(i,2)*lp1(2,2)+fr1(i,3)*lp1(3,2)
car1(2,i,3)=(fr1(i,1)+1)*lp1(1,3)+fr1(i,2)*lp1(2,3)+fr1(i,3)*lp1(3,3)
end do

do i=n2+1, nt
car1(3,i,1)=fr1(i,1)*lp1(1,1)+(fr1(i,2)+1)*lp1(2,1)+fr1(i,3)*lp1(3,1)
car1(3,i,2)=fr1(i,1)*lp1(1,2)+(fr1(i,2)+1)*lp1(2,2)+fr1(i,3)*lp1(3,2)
car1(3,i,3)=fr1(i,1)*lp1(1,3)+(fr1(i,2)+1)*lp1(2,3)+fr1(i,3)*lp1(3,3)
end do

do i=n2+1, nt
car1(4,i,1)=fr1(i,1)*lp1(1,1)+fr1(i,2)*lp1(2,1)+(fr1(i,3)+1)*lp1(3,1)
car1(4,i,2)=fr1(i,1)*lp1(1,2)+fr1(i,2)*lp1(2,2)+(fr1(i,3)+1)*lp1(3,2)
car1(4,i,3)=fr1(i,1)*lp1(1,3)+fr1(i,2)*lp1(2,3)+(fr1(i,3)+1)*lp1(3,3)
end do

do i=n2+1, nt
car1(5,i,1)=(fr1(i,1)-1)*lp1(1,1)+fr1(i,2)*lp1(2,1)+fr1(i,3)*lp1(3,1)
car1(5,i,2)=(fr1(i,1)-1)*lp1(1,2)+fr1(i,2)*lp1(2,2)+fr1(i,3)*lp1(3,2)
car1(5,i,3)=(fr1(i,1)-1)*lp1(1,3)+fr1(i,2)*lp1(2,3)+fr1(i,3)*lp1(3,3)
end do

do i=n2+1, nt
car1(6,i,1)=fr1(i,1)*lp1(1,1)+(fr1(i,2)-1)*lp1(2,1)+fr1(i,3)*lp1(3,1)
car1(6,i,2)=fr1(i,1)*lp1(1,2)+(fr1(i,2)-1)*lp1(2,2)+fr1(i,3)*lp1(3,2)
car1(6,i,3)=fr1(i,1)*lp1(1,3)+(fr1(i,2)-1)*lp1(2,3)+fr1(i,3)*lp1(3,3)
end do

do i=n2+1, nt
car1(7,i,1)=fr1(i,1)*lp1(1,1)+fr1(i,2)*lp1(2,1)+(fr1(i,3)-1)*lp1(3,1)
car1(7,i,2)=fr1(i,1)*lp1(1,2)+fr1(i,2)*lp1(2,2)+(fr1(i,3)-1)*lp1(3,2)
car1(7,i,3)=fr1(i,1)*lp1(1,3)+fr1(i,2)*lp1(2,3)+(fr1(i,3)-1)*lp1(3,3)
end do

!Rotation around y axis
do i=n1, n2

if(i/=n1) then
if(ny==1) then
if(abs(sqrt((car(1,i,1)-car(1,i-1,1))**2+(car(1,i,2)-car(1,i-1,2))**2+(car(1,i,3)-car(1,i-1,3))**2)-dtm)<=Cri) then
if(abs(car(1,i,2)-car(1,i-1,2))<=Cri) then

ang2=-ang2

end if
end if

else if(ny==-1) then
if(abs(sqrt((car(1,i,1)-car(1,i-1,1))**2+(car(1,i,2)-car(1,i-1,2))**2+(car(1,i,3)-car(1,i-1,3))**2)-dtm)<=Cri) then

ang2=-ang2

end if

end if

end if

do j=n2+1, nt

if(abs(sqrt((car1(1,j,1)-car(1,i,1))**2+(car1(1,j,2)-car(1,i,2))**2+(car1(1,j,3)-car(1,i,3))**2)-do1)<=Cri) then

if((car1(1,j,3)-car(1,i,3))<-1.0) then
car2(1,j,1)=car1(1,j,1)-ang2
car2(1,j,2)=car1(1,j,2)
car2(1,j,3)=car1(1,j,3)
else if((car1(1,j,3)-car(1,i,3))>1.0) then
car2(1,j,1)=car1(1,j,1)+ang2
car2(1,j,2)=car1(1,j,2)
car2(1,j,3)=car1(1,j,3)
else if((car1(1,j,1)-car(1,i,1))>1.0) then
car2(1,j,1)=car1(1,j,1)
car2(1,j,2)=car1(1,j,2)
car2(1,j,3)=car1(1,j,3)-ang2
else if((car1(1,j,1)-car(1,i,1))<-1.0) then
car2(1,j,1)=car1(1,j,1)
car2(1,j,2)=car1(1,j,2)
car2(1,j,3)=car1(1,j,3)+ang2
else if((car1(1,j,2)-car(1,i,2))<-1.0) then
car2(1,j,1)=car1(1,j,1)
car2(1,j,2)=car1(1,j,2)
car2(1,j,3)=car1(1,j,3)
else if((car1(1,j,2)-car(1,i,2))>1.0) then
car2(1,j,1)=car1(1,j,1)
car2(1,j,2)=car1(1,j,2)
car2(1,j,3)=car1(1,j,3)
end if

end if

if(abs(sqrt((car1(2,j,1)-car(1,i,1))**2+(car1(2,j,2)-car(1,i,2))**2+(car1(2,j,3)-car(1,i,3))**2)-do1)<=Cri) then

if((car1(2,j,3)-car(1,i,3))<-1.0) then
car2(1,j,1)=car1(1,j,1)-ang2
car2(1,j,2)=car1(1,j,2)
car2(1,j,3)=car1(1,j,3)
else if((car1(2,j,3)-car(1,i,3))>1.0) then
car2(1,j,1)=car1(1,j,1)+ang2
car2(1,j,2)=car1(1,j,2)
car2(1,j,3)=car1(1,j,3)
else if((car1(2,j,1)-car(1,i,1))>1.0) then
car2(1,j,1)=car1(1,j,1)
car2(1,j,2)=car1(1,j,2)
car2(1,j,3)=car1(1,j,3)-ang2
else if((car1(2,j,1)-car(1,i,1))<-1.0) then
car2(1,j,1)=car1(1,j,1)
car2(1,j,2)=car1(1,j,2)
car2(1,j,3)=car1(1,j,3)+ang2
else if((car1(2,j,2)-car(1,i,2))<-1.0) then
car2(1,j,1)=car1(1,j,1)
car2(1,j,2)=car1(1,j,2)
car2(1,j,3)=car1(1,j,3)
else if((car1(2,j,2)-car(1,i,2))>1.0) then
car2(1,j,1)=car1(1,j,1)
car2(1,j,2)=car1(1,j,2)
car2(1,j,3)=car1(1,j,3)
end if

end if

if(abs(sqrt((car1(3,j,1)-car(1,i,1))**2+(car1(3,j,2)-car(1,i,2))**2+(car1(3,j,3)-car(1,i,3))**2)-do1)<=Cri) then

if((car1(3,j,3)-car(1,i,3))<-1.0) then
car2(1,j,1)=car1(1,j,1)-ang2
car2(1,j,2)=car1(1,j,2)
car2(1,j,3)=car1(1,j,3)
else if((car1(3,j,3)-car(1,i,3))>1.0) then
car2(1,j,1)=car1(1,j,1)+ang2
car2(1,j,2)=car1(1,j,2)
car2(1,j,3)=car1(1,j,3)
else if((car1(3,j,1)-car(1,i,1))>1.0) then
car2(1,j,1)=car1(1,j,1)
car2(1,j,2)=car1(1,j,2)
car2(1,j,3)=car1(1,j,3)-ang2
else if((car1(3,j,1)-car(1,i,1))<-1.0) then
car2(1,j,1)=car1(1,j,1)
car2(1,j,2)=car1(1,j,2)
car2(1,j,3)=car1(1,j,3)+ang2
else if((car1(3,j,2)-car(1,i,2))<-1.0) then
car2(1,j,1)=car1(1,j,1)
car2(1,j,2)=car1(1,j,2)
car2(1,j,3)=car1(1,j,3)
else if((car1(3,j,2)-car(1,i,2))>1.0) then
car2(1,j,1)=car1(1,j,1)
car2(1,j,2)=car1(1,j,2)
car2(1,j,3)=car1(1,j,3)
end if

end if

if(abs(sqrt((car1(4,j,1)-car(1,i,1))**2+(car1(4,j,2)-car(1,i,2))**2+(car1(4,j,3)-car(1,i,3))**2)-do1)<=Cri) then

if((car1(4,j,3)-car(1,i,3))<-1.0) then
car2(1,j,1)=car1(1,j,1)-ang2
car2(1,j,2)=car1(1,j,2)
car2(1,j,3)=car1(1,j,3)
else if((car1(4,j,3)-car(1,i,3))>1.0) then
car2(1,j,1)=car1(1,j,1)+ang2
car2(1,j,2)=car1(1,j,2)
car2(1,j,3)=car1(1,j,3)
else if((car1(4,j,1)-car(1,i,1))>1.0) then
car2(1,j,1)=car1(1,j,1)
car2(1,j,2)=car1(1,j,2)
car2(1,j,3)=car1(1,j,3)-ang2
else if((car1(4,j,1)-car(1,i,1))<-1.0) then
car2(1,j,1)=car1(1,j,1)
car2(1,j,2)=car1(1,j,2)
car2(1,j,3)=car1(1,j,3)+ang2
else if((car1(4,j,2)-car(1,i,2))<-1.0) then
car2(1,j,1)=car1(1,j,1)
car2(1,j,2)=car1(1,j,2)
car2(1,j,3)=car1(1,j,3)
else if((car1(4,j,2)-car(1,i,2))>1.0) then
car2(1,j,1)=car1(1,j,1)
car2(1,j,2)=car1(1,j,2)
car2(1,j,3)=car1(1,j,3)
end if

end if

if(abs(sqrt((car1(5,j,1)-car(1,i,1))**2+(car1(5,j,2)-car(1,i,2))**2+(car1(5,j,3)-car(1,i,3))**2)-do1)<=Cri) then

if((car1(5,j,3)-car(1,i,3))<-1.0) then
car2(1,j,1)=car1(1,j,1)-ang2
car2(1,j,2)=car1(1,j,2)
car2(1,j,3)=car1(1,j,3)
else if((car1(5,j,3)-car(1,i,3))>1.0) then
car2(1,j,1)=car1(1,j,1)+ang2
car2(1,j,2)=car1(1,j,2)
car2(1,j,3)=car1(1,j,3)
else if((car1(5,j,1)-car(1,i,1))>1.0) then
car2(1,j,1)=car1(1,j,1)
car2(1,j,2)=car1(1,j,2)
car2(1,j,3)=car1(1,j,3)-ang2
else if((car1(5,j,1)-car(1,i,1))<-1.0) then
car2(1,j,1)=car1(1,j,1)
car2(1,j,2)=car1(1,j,2)
car2(1,j,3)=car1(1,j,3)+ang2
else if((car1(5,j,2)-car(1,i,2))<-1.0) then
car2(1,j,1)=car1(1,j,1)
car2(1,j,2)=car1(1,j,2)
car2(1,j,3)=car1(1,j,3)
else if((car1(5,j,2)-car(1,i,2))>1.0) then
car2(1,j,1)=car1(1,j,1)
car2(1,j,2)=car1(1,j,2)
car2(1,j,3)=car1(1,j,3)
end if

end if

if(abs(sqrt((car1(6,j,1)-car(1,i,1))**2+(car1(6,j,2)-car(1,i,2))**2+(car1(6,j,3)-car(1,i,3))**2)-do1)<=Cri) then

if((car1(6,j,3)-car(1,i,3))<-1.0) then
car2(1,j,1)=car1(1,j,1)-ang2
car2(1,j,2)=car1(1,j,2)
car2(1,j,3)=car1(1,j,3)
else if((car1(6,j,3)-car(1,i,3))>1.0) then
car2(1,j,1)=car1(1,j,1)+ang2
car2(1,j,2)=car1(1,j,2)
car2(1,j,3)=car1(1,j,3)
else if((car1(6,j,1)-car(1,i,1))>1.0) then
car2(1,j,1)=car1(1,j,1)
car2(1,j,2)=car1(1,j,2)
car2(1,j,3)=car1(1,j,3)-ang2
else if((car1(6,j,1)-car(1,i,1))<-1.0) then
car2(1,j,1)=car1(1,j,1)
car2(1,j,2)=car1(1,j,2)
car2(1,j,3)=car1(1,j,3)+ang2
else if((car1(6,j,2)-car(1,i,2))<-1.0) then
car2(1,j,1)=car1(1,j,1)
car2(1,j,2)=car1(1,j,2)
car2(1,j,3)=car1(1,j,3)
else if((car1(6,j,2)-car(1,i,2))>1.0) then
car2(1,j,1)=car1(1,j,1)
car2(1,j,2)=car1(1,j,2)
car2(1,j,3)=car1(1,j,3)
end if

end if

if(abs(sqrt((car1(7,j,1)-car(1,i,1))**2+(car1(7,j,2)-car(1,i,2))**2+(car1(7,j,3)-car(1,i,3))**2)-do1)<=Cri) then

if((car1(7,j,3)-car(1,i,3))<-1.0) then
car2(1,j,1)=car1(1,j,1)-ang2
car2(1,j,2)=car1(1,j,2)
car2(1,j,3)=car1(1,j,3)
else if((car1(7,j,3)-car(1,i,3))>1.0) then
car2(1,j,1)=car1(1,j,1)+ang2
car2(1,j,2)=car1(1,j,2)
car2(1,j,3)=car1(1,j,3)
else if((car1(7,j,1)-car(1,i,1))>1.0) then
car2(1,j,1)=car1(1,j,1)
car2(1,j,2)=car1(1,j,2)
car2(1,j,3)=car1(1,j,3)-ang2
else if((car1(7,j,1)-car(1,i,1))<-1.0) then
car2(1,j,1)=car1(1,j,1)
car2(1,j,2)=car1(1,j,2)
car2(1,j,3)=car1(1,j,3)+ang2
else if((car1(7,j,2)-car(1,i,2))<-1.0) then
car2(1,j,1)=car1(1,j,1)
car2(1,j,2)=car1(1,j,2)
car2(1,j,3)=car1(1,j,3)
else if((car1(7,j,2)-car(1,i,2))>1.0) then
car2(1,j,1)=car1(1,j,1)
car2(1,j,2)=car1(1,j,2)
car2(1,j,3)=car1(1,j,3)
end if

end if

end do
end do
!End

!The transformation from car to dir!
mr=lp1(1,1)*lp1(2,2)*lp1(3,3)+lp1(1,2)*lp1(2,3)*lp1(3,1)+lp1(1,3)*lp1(3,2)*lp1(2,1)&
-lp1(1,3)*lp1(2,2)*lp1(3,1)-lp1(1,1)*lp1(3,2)*lp1(2,3)-lp1(1,2)*lp1(2,1)*lp1(3,3)

do i=n2+1, nt

fr2(i,1)=(car2(1,i,1)*(lp1(2,2)*lp1(3,3)-lp1(2,3)*lp1(3,2))+car2(1,i,2)*(-lp1(2,1)*&
lp1(3,3)+lp1(2,3)*lp1(3,1))+car2(1,i,3)*(lp1(2,1)*lp1(3,2)-lp1(2,2)*lp1(3,1)))/mr
fr2(i,2)=(car2(1,i,1)*(-lp1(1,2)*lp1(3,3)+lp1(1,3)*lp1(3,2))+car2(1,i,2)*(lp1(1,1)*&
lp1(3,3)-lp1(1,3)*lp1(3,1))+car2(1,i,3)*(-lp1(1,1)*lp1(3,2)+lp1(1,2)*lp1(3,1)))/mr
fr2(i,3)=(car2(1,i,1)*(lp1(1,2)*lp1(2,3)-lp1(1,3)*lp1(2,2))+car2(1,i,2)*(-lp1(1,1)*&
lp1(2,3)+lp1(1,3)*lp1(2,1))+car2(1,i,3)*(lp1(1,1)*lp1(2,2)-lp1(1,2)*lp1(2,1)))/mr

end do
!End

do i=n2+1, nt
car2(2,i,1)=(fr2(i,1)+1)*lp1(1,1)+fr2(i,2)*lp1(2,1)+fr2(i,3)*lp1(3,1)
car2(2,i,2)=(fr2(i,1)+1)*lp1(1,2)+fr2(i,2)*lp1(2,2)+fr2(i,3)*lp1(3,2)
car2(2,i,3)=(fr2(i,1)+1)*lp1(1,3)+fr2(i,2)*lp1(2,3)+fr2(i,3)*lp1(3,3)
end do

do i=n2+1, nt
car2(3,i,1)=fr2(i,1)*lp1(1,1)+(fr2(i,2)+1)*lp1(2,1)+fr2(i,3)*lp1(3,1)
car2(3,i,2)=fr2(i,1)*lp1(1,2)+(fr2(i,2)+1)*lp1(2,2)+fr2(i,3)*lp1(3,2)
car2(3,i,3)=fr2(i,1)*lp1(1,3)+(fr2(i,2)+1)*lp1(2,3)+fr2(i,3)*lp1(3,3)
end do

do i=n2+1, nt
car2(4,i,1)=fr2(i,1)*lp1(1,1)+fr2(i,2)*lp1(2,1)+(fr2(i,3)+1)*lp1(3,1)
car2(4,i,2)=fr2(i,1)*lp1(1,2)+fr2(i,2)*lp1(2,2)+(fr2(i,3)+1)*lp1(3,2)
car2(4,i,3)=fr2(i,1)*lp1(1,3)+fr2(i,2)*lp1(2,3)+(fr2(i,3)+1)*lp1(3,3)
end do

do i=n2+1, nt
car2(5,i,1)=(fr2(i,1)-1)*lp1(1,1)+fr2(i,2)*lp1(2,1)+fr2(i,3)*lp1(3,1)
car2(5,i,2)=(fr2(i,1)-1)*lp1(1,2)+fr2(i,2)*lp1(2,2)+fr2(i,3)*lp1(3,2)
car2(5,i,3)=(fr2(i,1)-1)*lp1(1,3)+fr2(i,2)*lp1(2,3)+fr2(i,3)*lp1(3,3)
end do

do i=n2+1, nt
car2(6,i,1)=fr2(i,1)*lp1(1,1)+(fr2(i,2)-1)*lp1(2,1)+fr2(i,3)*lp1(3,1)
car2(6,i,2)=fr2(i,1)*lp1(1,2)+(fr2(i,2)-1)*lp1(2,2)+fr2(i,3)*lp1(3,2)
car2(6,i,3)=fr2(i,1)*lp1(1,3)+(fr2(i,2)-1)*lp1(2,3)+fr2(i,3)*lp1(3,3)
end do

do i=n2+1, nt
car2(7,i,1)=fr2(i,1)*lp1(1,1)+fr2(i,2)*lp1(2,1)+(fr2(i,3)-1)*lp1(3,1)
car2(7,i,2)=fr2(i,1)*lp1(1,2)+fr2(i,2)*lp1(2,2)+(fr2(i,3)-1)*lp1(3,2)
car2(7,i,3)=fr2(i,1)*lp1(1,3)+fr2(i,2)*lp1(2,3)+(fr2(i,3)-1)*lp1(3,3)
end do

!Rotation around z axis
do i=n1, n2

if(i/=n1) then
if(nz==1) then
if(abs(sqrt((car(1,i,1)-car(1,i-1,1))**2+(car(1,i,2)-car(1,i-1,2))**2+(car(1,i,3)-car(1,i-1,3))**2)-dtm)<=Cri) then
if(abs(car(1,i,3)-car(1,i-1,3))<=Cri) then

ang3=-ang3

end if
end if

else if(nz==-1) then
if(abs(sqrt((car(1,i,1)-car(1,i-1,1))**2+(car(1,i,2)-car(1,i-1,2))**2+(car(1,i,3)-car(1,i-1,3))**2)-dtm)<=Cri) then

ang3=-ang3

end if

end if

end if

do j=n2+1, nt

if(abs(sqrt((car2(1,j,1)-car(1,i,1))**2+(car2(1,j,2)-car(1,i,2))**2+(car2(1,j,3)-car(1,i,3))**2)-do1)<=Cri) then

if((car2(1,j,3)-car(1,i,3))<-1.0) then
car3(1,j,1)=car2(1,j,1)
car3(1,j,2)=car2(1,j,2)
car3(1,j,3)=car2(1,j,3)
else if((car2(1,j,3)-car(1,i,3))>1.0) then
car3(1,j,1)=car2(1,j,1)
car3(1,j,2)=car2(1,j,2)
car3(1,j,3)=car2(1,j,3)
else if((car2(1,j,1)-car(1,i,1))>1.0) then
car3(1,j,1)=car2(1,j,1)
car3(1,j,2)=car2(1,j,2)+ang3
car3(1,j,3)=car2(1,j,3)
else if((car2(1,j,1)-car(1,i,1))<-1.0) then
car3(1,j,1)=car2(1,j,1)
car3(1,j,2)=car2(1,j,2)-ang3
car3(1,j,3)=car2(1,j,3)
else if((car2(1,j,2)-car(1,i,2))<-1.0) then
car3(1,j,1)=car2(1,j,1)+ang3
car3(1,j,2)=car2(1,j,2)
car3(1,j,3)=car2(1,j,3)
else if((car2(1,j,2)-car(1,i,2))>1.0) then
car3(1,j,1)=car2(1,j,1)-ang3
car3(1,j,2)=car2(1,j,2)
car3(1,j,3)=car2(1,j,3)
end if

end if

if(abs(sqrt((car2(2,j,1)-car(1,i,1))**2+(car2(2,j,2)-car(1,i,2))**2+(car2(2,j,3)-car(1,i,3))**2)-do1)<=Cri) then

if((car2(2,j,3)-car(1,i,3))<-1.0) then
car3(1,j,1)=car2(1,j,1)
car3(1,j,2)=car2(1,j,2)
car3(1,j,3)=car2(1,j,3)
else if((car2(2,j,3)-car(1,i,3))>1.0) then
car3(1,j,1)=car2(1,j,1)
car3(1,j,2)=car2(1,j,2)
car3(1,j,3)=car2(1,j,3)
else if((car2(2,j,1)-car(1,i,1))>1.0) then
car3(1,j,1)=car2(1,j,1)
car3(1,j,2)=car2(1,j,2)+ang3
car3(1,j,3)=car2(1,j,3)
else if((car2(2,j,1)-car(1,i,1))<-1.0) then
car3(1,j,1)=car2(1,j,1)
car3(1,j,2)=car2(1,j,2)-ang3
car3(1,j,3)=car2(1,j,3)
else if((car2(2,j,2)-car(1,i,2))<-1.0) then
car3(1,j,1)=car2(1,j,1)+ang3
car3(1,j,2)=car2(1,j,2)
car3(1,j,3)=car2(1,j,3)
else if((car2(2,j,2)-car(1,i,2))>1.0) then
car3(1,j,1)=car2(1,j,1)-ang3
car3(1,j,2)=car2(1,j,2)
car3(1,j,3)=car2(1,j,3)
end if

end if

if(abs(sqrt((car2(3,j,1)-car(1,i,1))**2+(car2(3,j,2)-car(1,i,2))**2+(car2(3,j,3)-car(1,i,3))**2)-do1)<=Cri) then

if((car2(3,j,3)-car(1,i,3))<-1.0) then
car3(1,j,1)=car2(1,j,1)
car3(1,j,2)=car2(1,j,2)
car3(1,j,3)=car2(1,j,3)
else if((car2(3,j,3)-car(1,i,3))>1.0) then
car3(1,j,1)=car2(1,j,1)
car3(1,j,2)=car2(1,j,2)
car3(1,j,3)=car2(1,j,3)
else if((car2(3,j,1)-car(1,i,1))>1.0) then
car3(1,j,1)=car2(1,j,1)
car3(1,j,2)=car2(1,j,2)+ang3
car3(1,j,3)=car2(1,j,3)
else if((car2(3,j,1)-car(1,i,1))<-1.0) then
car3(1,j,1)=car2(1,j,1)
car3(1,j,2)=car2(1,j,2)-ang3
car3(1,j,3)=car2(1,j,3)
else if((car2(3,j,2)-car(1,i,2))<-1.0) then
car3(1,j,1)=car2(1,j,1)+ang3
car3(1,j,2)=car2(1,j,2)
car3(1,j,3)=car2(1,j,3)
else if((car2(3,j,2)-car(1,i,2))>1.0) then
car3(1,j,1)=car2(1,j,1)-ang3
car3(1,j,2)=car2(1,j,2)
car3(1,j,3)=car2(1,j,3)
end if

end if

if(abs(sqrt((car2(4,j,1)-car(1,i,1))**2+(car2(4,j,2)-car(1,i,2))**2+(car2(4,j,3)-car(1,i,3))**2)-do1)<=Cri) then

if((car2(4,j,3)-car(1,i,3))<-1.0) then
car3(1,j,1)=car2(1,j,1)
car3(1,j,2)=car2(1,j,2)
car3(1,j,3)=car2(1,j,3)
else if((car2(4,j,3)-car(1,i,3))>1.0) then
car3(1,j,1)=car2(1,j,1)
car3(1,j,2)=car2(1,j,2)
car3(1,j,3)=car2(1,j,3)
else if((car2(4,j,1)-car(1,i,1))>1.0) then
car3(1,j,1)=car2(1,j,1)
car3(1,j,2)=car2(1,j,2)+ang3
car3(1,j,3)=car2(1,j,3)
else if((car2(4,j,1)-car(1,i,1))<-1.0) then
car3(1,j,1)=car2(1,j,1)
car3(1,j,2)=car2(1,j,2)-ang3
car3(1,j,3)=car2(1,j,3)
else if((car2(4,j,2)-car(1,i,2))<-1.0) then
car3(1,j,1)=car2(1,j,1)+ang3
car3(1,j,2)=car2(1,j,2)
car3(1,j,3)=car2(1,j,3)
else if((car2(4,j,2)-car(1,i,2))>1.0) then
car3(1,j,1)=car2(1,j,1)-ang3
car3(1,j,2)=car2(1,j,2)
car3(1,j,3)=car2(1,j,3)
end if

end if

if(abs(sqrt((car2(5,j,1)-car(1,i,1))**2+(car2(5,j,2)-car(1,i,2))**2+(car2(5,j,3)-car(1,i,3))**2)-do1)<=Cri) then

if((car2(5,j,3)-car(1,i,3))<-1.0) then
car3(1,j,1)=car2(1,j,1)
car3(1,j,2)=car2(1,j,2)
car3(1,j,3)=car2(1,j,3)
else if((car2(5,j,3)-car(1,i,3))>1.0) then
car3(1,j,1)=car2(1,j,1)
car3(1,j,2)=car2(1,j,2)
car3(1,j,3)=car2(1,j,3)
else if((car2(5,j,1)-car(1,i,1))>1.0) then
car3(1,j,1)=car2(1,j,1)
car3(1,j,2)=car2(1,j,2)+ang3
car3(1,j,3)=car2(1,j,3)
else if((car2(5,j,1)-car(1,i,1))<-1.0) then
car3(1,j,1)=car2(1,j,1)
car3(1,j,2)=car2(1,j,2)-ang3
car3(1,j,3)=car2(1,j,3)
else if((car2(5,j,2)-car(1,i,2))<-1.0) then
car3(1,j,1)=car2(1,j,1)+ang3
car3(1,j,2)=car2(1,j,2)
car3(1,j,3)=car2(1,j,3)
else if((car2(5,j,2)-car(1,i,2))>1.0) then
car3(1,j,1)=car2(1,j,1)-ang3
car3(1,j,2)=car2(1,j,2)
car3(1,j,3)=car2(1,j,3)
end if

end if

if(abs(sqrt((car2(6,j,1)-car(1,i,1))**2+(car2(6,j,2)-car(1,i,2))**2+(car2(6,j,3)-car(1,i,3))**2)-do1)<=Cri) then

if((car2(6,j,3)-car(1,i,3))<-1.0) then
car3(1,j,1)=car2(1,j,1)
car3(1,j,2)=car2(1,j,2)
car3(1,j,3)=car2(1,j,3)
else if((car2(6,j,3)-car(1,i,3))>1.0) then
car3(1,j,1)=car2(1,j,1)
car3(1,j,2)=car2(1,j,2)
car3(1,j,3)=car2(1,j,3)
else if((car2(6,j,1)-car(1,i,1))>1.0) then
car3(1,j,1)=car2(1,j,1)
car3(1,j,2)=car2(1,j,2)+ang3
car3(1,j,3)=car2(1,j,3)
else if((car2(6,j,1)-car(1,i,1))<-1.0) then
car3(1,j,1)=car2(1,j,1)
car3(1,j,2)=car2(1,j,2)-ang3
car3(1,j,3)=car2(1,j,3)
else if((car2(6,j,2)-car(1,i,2))<-1.0) then
car3(1,j,1)=car2(1,j,1)+ang3
car3(1,j,2)=car2(1,j,2)
car3(1,j,3)=car2(1,j,3)
else if((car2(6,j,2)-car(1,i,2))>1.0) then
car3(1,j,1)=car2(1,j,1)-ang3
car3(1,j,2)=car2(1,j,2)
car3(1,j,3)=car2(1,j,3)
end if

end if

if(abs(sqrt((car2(7,j,1)-car(1,i,1))**2+(car2(7,j,2)-car(1,i,2))**2+(car2(7,j,3)-car(1,i,3))**2)-do1)<=Cri) then

if((car2(7,j,3)-car(1,i,3))<-1.0) then
car3(1,j,1)=car2(1,j,1)
car3(1,j,2)=car2(1,j,2)
car3(1,j,3)=car2(1,j,3)
else if((car2(7,j,3)-car(1,i,3))>1.0) then
car3(1,j,1)=car2(1,j,1)
car3(1,j,2)=car2(1,j,2)
car3(1,j,3)=car2(1,j,3)
else if((car2(7,j,1)-car(1,i,1))>1.0) then
car3(1,j,1)=car2(1,j,1)
car3(1,j,2)=car2(1,j,2)+ang3
car3(1,j,3)=car2(1,j,3)
else if((car2(7,j,1)-car(1,i,1))<-1.0) then
car3(1,j,1)=car2(1,j,1)
car3(1,j,2)=car2(1,j,2)-ang3
car3(1,j,3)=car2(1,j,3)
else if((car2(7,j,2)-car(1,i,2))<-1.0) then
car3(1,j,1)=car2(1,j,1)+ang3
car3(1,j,2)=car2(1,j,2)
car3(1,j,3)=car2(1,j,3)
else if((car2(7,j,2)-car(1,i,2))>1.0) then
car3(1,j,1)=car2(1,j,1)-ang3
car3(1,j,2)=car2(1,j,2)
car3(1,j,3)=car2(1,j,3)
end if

end if

end do
end do
!End

do i=1, 3
lp(i,1)=m(i,1)*ap
lp(i,2)=m(i,2)*bp
lp(i,3)=m(i,3)*cp
end do

hei(1)=sqrt(lp1(3,1)**2+lp1(3,2)**2+lp1(3,3)**2)

num(1)=r1*(n1-1)/(r1+r2)
i=1
k(1)=1
k(2)=1

do j=1, n1-1

if(car(1,i,3)-r2*hei(1)/(r1+r2)>-Cri) then
car3(1,k(1),1)=car(1,i,1)
car3(1,k(1),2)=car(1,i,2)
car3(1,k(1),3)=car(1,i,3)
k(1)=k(1)+1
else
car3(1,k(2)+num(1),1)=car(1,i,1)
car3(1,k(2)+num(1),2)=car(1,i,2)
car3(1,k(2)+num(1),3)=car(1,i,3)
k(2)=k(2)+1
end if
i=i+1
end do

!Moving the ions at B positions
num(2)=r1*(n2-n1+1)/(r1+r2)-1
i=n1
k(3)=n1
k(4)=n1

do j=n1, n2

if(car(1,i,3)-r2*hei(1)/(r1+r2)>-Cri) then
car3(1,k(3),1)=car(1,i,1)+dist2(1)
car3(1,k(3),2)=car(1,i,2)+dist2(2)
car3(1,k(3),3)=car(1,i,3)+dist2(3)
k(3)=k(3)+1
else
car3(1,k(4)+num(2)+1,1)=car(1,i,1)+dist2(1)
car3(1,k(4)+num(2)+1,2)=car(1,i,2)+dist2(2)
car3(1,k(4)+num(2)+1,3)=car(1,i,3)+dist2(3)
k(4)=k(4)+1
end if
i=i+1

end do
!End

mr=lp1(1,1)*lp1(2,2)*lp1(3,3)+lp1(1,2)*lp1(2,3)*lp1(3,1)+lp1(1,3)*lp1(3,2)*lp1(2,1)&
-lp1(1,3)*lp1(2,2)*lp1(3,1)-lp1(1,1)*lp1(3,2)*lp1(2,3)-lp1(1,2)*lp1(2,1)*lp1(3,3)

do i=1, nt

fr3(i,1)=(car3(1,i,1)*(lp1(2,2)*lp1(3,3)-lp1(2,3)*lp1(3,2))+car3(1,i,2)*(-lp1(2,1)*&
lp1(3,3)+lp1(2,3)*lp1(3,1))+car3(1,i,3)*(lp1(2,1)*lp1(3,2)-lp1(2,2)*lp1(3,1)))/mr
fr3(i,2)=(car3(1,i,1)*(-lp1(1,2)*lp1(3,3)+lp1(1,3)*lp1(3,2))+car3(1,i,2)*(lp1(1,1)*&
lp1(3,3)-lp1(1,3)*lp1(3,1))+car3(1,i,3)*(-lp1(1,1)*lp1(3,2)+lp1(1,2)*lp1(3,1)))/mr
fr3(i,3)=(car3(1,i,1)*(lp1(1,2)*lp1(2,3)-lp1(1,3)*lp1(2,2))+car3(1,i,2)*(-lp1(1,1)*&
lp1(2,3)+lp1(1,3)*lp1(2,1))+car3(1,i,3)*(lp1(1,1)*lp1(2,2)-lp1(1,2)*lp1(2,1)))/mr

end do

do i=1, nt
if(fr3(i,1)<-1E-5) then
fr3(i,1)=fr3(i,1)+1
end if
if(fr3(i,2)<-1E-5) then
fr3(i,2)=fr3(i,2)+1
end if
if(fr3(i,3)<-1E-5) then
fr3(i,3)=fr3(i,3)+1
end if
end do

do i=1, nt
if(fr3(i,1)>=1.0) then
fr3(i,1)=fr3(i,1)-1
end if
if(fr3(i,2)>1.0) then
fr3(i,2)=fr3(i,2)-1
end if
if(fr3(i,3)>1.0) then
fr3(i,3)=fr3(i,3)-1
end if
end do

open(10, file="POSCAR")
write(10,*) "Perovskite system"
write(10,*) 1.000
do i=1, 3
write(10,"(3(1X,F16.11))") lp(i,1), lp(i,2), lp(i,3)
end do
if(e(1)==e(2).and.e(3)/=e(4)) then
write(10,*) n1-1, num(2)+1, n2-n1-num(2), nt-n2
else if(e(1)/=e(2).and.e(3)==e(4)) then
write(10,*) num(1), n1-1-num(1), n2-n1+1, nt-n2
else if(e(1)/=e(2).and.e(3)/=e(4)) then
write(10,*) num(1), n1-1-num(1), num(2)+1, n2-n1-num(2), nt-n2
else if(e(1)==e(2).and.e(3)==e(4)) then
write(10,*) n1-1, n2-n1+1, nt-n2
end if
write(10,*) "Direct"
do i=1, nt
write(10,"(3(1X,F16.11))") fr3(i,1), fr3(i,2), fr3(i,3)
end do
close(10)

stop
end
