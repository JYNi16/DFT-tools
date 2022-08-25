program writeband
implicit none

integer i,j,k,mistake,ispin
integer nk,nvolue,nkd,nband,num,kpoints
real lattice,volum,pi,ax,ay,az,z,y,q,fermi,a
real vector(3,3),b(3,3),c(3,3),d(3),e(3),x(3),dis1(3,3)
real,allocatable::position(:,:),eigenval(:,:),weight(:),distance(:),dis2(:),eigenup(:,:),eigendn(:,:)

pi=2.0*asin(1.0)
!write(*,*) pi
write(*,*)'*****  please prepare the POSCAR DOSCAR EIGENVAL KPOINTS well  ******'
write(*,*)'test 1'

open(unit=10,file='EIGENVAL',status='old',iostat=mistake)
if (mistake.ne.0)then
  write(*,*)"*********   can't search EIGENVAL,please check it   *********"
  stop
end if

open(unit=9,file='POSCAR',status='old',iostat=mistake)
if (mistake.ne.0)then
  write(*,*) "*********   can't search POSCAR,please check it   *********"
  stop
end if

open(unit=7,file='DOSCAR',status='old',iostat=mistake)
if (mistake.ne.0)then
  write(*,*) "*********   can't search DOSCAR,please check it   *********"
  stop
end if

open(unit=8,file='KPOINTS',status='old',iostat=mistake)
if (mistake.ne.0)then
 ! write(*,*)'******    can't search KPOINTS,please check it   *********'
  stop
end if
read(8,*)
read(8,*) kpoints
write(*,*)kpoints
close(8)

write(*,*)'test 2'

write(*,*)'*******   please input the number of ispin   *******'
read (*,*) ispin
if (ispin==1) then
  open(unit=11,file='band.dat',status='replace')
  else if(ispin==2) then
    open(unit=11,file='bandup.dat',status='replace')
	open(unit=12,file='banddn.dat',status='replace')

	else 
	write(*,*)'********   please check the ispin=1 or 2   ********'
end if 


write (*,*) 'test 3'

do i=1,5
  read (10,*)
  read (7,*)
enddo
read (7,*)z,y,q,fermi
read (10,*)nvolue,nk,nband
close(7)
!write(*,*)fermi

allocate(position(nk,3))
allocate(weight(nk))
allocate(eigenval(nk,nband))
allocate(eigenup(nk,nband))
allocate(eigendn(nk,nband))
allocate(distance(nk))
allocate(dis2(nk))

do i=1,nk
  read(10,*)
  read(10,*)(position(i,k),k=1,3),weight(i)
  do j=1,nband
    if(ispin==1)then
	  read(10,*)num,eigenval(i,j)
	else
	  read(10,*)num,eigenup(i,j),eigendn(i,j)
	  end if
  end do
end do
close(10)

!100 format(3E14.7,E12.7)
!101 format(I3,f8.6)
!102 format(I3,2f10.6)

write(*,*) 'test 4'

read (9,*)
read(9,*) lattice
do i=1,3
  read(9,*)(vector(i,j),j=1,3)
end do
close(9)

write(*,*) 'test 5'

do i=1,3
  do j=1,3
    c(i,j)=lattice*vector(i,j)
  end do
end do

d(1)=c(2,2)*c(3,3)-c(2,3)*c(3,2)
d(2)=c(2,3)*c(3,1)-c(2,1)*c(3,3)
d(3)=c(2,1)*c(3,2)-c(2,2)*c(3,1)

volum=c(1,1)*d(1)+c(1,2)*d(2)+c(1,3)*d(3)   !混合积

!write(*,*) volum

write(*,*) 'test 6'

!b(1,1)=2*pi*d(2)/volum
!b(1,2)=2*pi*d(3)/volum
!b(1,3)=2*pi*d(1)/volum

do i=1,3
  if(i==1) then
    j=2
	k=3
  else if(i==2) then
	  j=3
	  k=1
    else
    	j=1
    	k=2
    end if
 
  e(1)=c(j,2)*c(k,3)-c(j,3)*c(k,2)
  e(2)=c(j,3)*c(k,1)-c(j,1)*c(k,3)
  e(3)=c(j,1)*c(k,2)-c(j,2)*c(k,1)
  do j=1,3
    b(i,j)=2*pi*e(j)/volum 
	!write(*,*)b(i,j)
  end do   !倒格子矢量
end do

write(*,*) 'test 7'

distance(1)=0.00
!do i=2,nk
!  do j=1,3   !坐标轴信息
 !   do k=1,3  !晶矢信息
!      dis1(j,k)=(position(i,k)-position(i-1,k))*b(k,j)
!	  x(j)=x(j)+dis1(j,k)
!	end do     
!  end do 
!  dis2(i)=sqrt(x(1)**2+x(2)**2+x(3)**2)
!  distance(i)=distance(i-1)+dis2(i)
!end do

do i=2,nk
  ax=(position(i,1)-position(i-1,1))*b(1,1)+(position(i,2)-position(i-1,2))*b(2,1)+(position(i,3)-position(i-1,3))*b(3,1)
  ay=(position(i,1)-position(i-1,1))*b(1,2)+(position(i,2)-position(i-1,2))*b(2,2)+(position(i,3)-position(i-1,3))*b(3,2)
  az=(position(i,1)-position(i-1,1))*b(1,3)+(position(i,2)-position(i-1,2))*b(2,3)+(position(i,3)-position(i-1,3))*b(3,3)
  a=sqrt(ax**2+ay**2+az**2)
  distance(i)=distance(I-1)+a
end do

write(*,*) 'test 8'

if(ispin==1)then
  do j=1,nband
    do i=1,nk
      write(11,103)distance(i),eigenval(i,j)-fermi
	  end do
	write(11,*)
  end do
do i=1,nk
  if(mod(i,kpoints)==0)then
    write(11,*)distance(i)
  end if
end do


!if( i==56)then
!write(11,*)distance(i)
!endif
!if( i==76)then
!write(11,*)distance(i)
!endif
!if( i==96)then
!!write(11,*)distance(i)
!endif
!if( i==116)then
!write(11,*)distance(i)
!endif
!if( i==136)then
!write(11,*)distance(i)
!endif

!enddo

  close(11)
else if(ispin==2)then
  do j=1,nband
    do i=1,nk
	  write(11,103)distance(i),eigenup(i,j)-fermi
	  write(12,103)distance(i),eigendn(i,j)-fermi
	enddo
  write(11,*)
  write(12,*)
  end do
do i=1,nk
  if(mod(i,kpoints)==0)then
    write(11,*)distance(i)
	write(12,*)distance(i)
  end if
end do
  close(11)
  close(12)
end if

deallocate(position)
deallocate(weight)
deallocate(eigenval)
deallocate(eigenup)
deallocate(eigendn)
deallocate(distance)
deallocate(dis2)


103 format(f10.6,f12.6)

end







