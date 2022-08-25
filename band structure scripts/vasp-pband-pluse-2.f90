8program writepband
implicit none
      !parameter (nbtot = 200) 
      !parameter (nktot = 500) 
      !parameter (nxd = 300) 
      !parameter (niontot = 20) 
      real*8 emin,emax,ef,quanzhong,volume,del,dkx,dky,dkz,pi,aa,c(3),aaa(3,3),b(3,3),w,zero
      integer ispin,nkpt,nband,nion,nion_plot,i,j,k,nk,nb,ni,nn,ndiv,kp,kk,n,natom,x
	  integer,allocatable ::atom(:)
        
      !real*8  xx(nxd)
      real*8  dump(20)
real*8,allocatable::eig(:,:,:),kpt(:,:),wei(:,:,:,:,:),xx(:),ee(:,:,:),twei(:,:,:)
!character(30) null
pi=3.1415926

open(22,file='KPOINTS',status='old')
read(22,*)
read(22,*)ndiv

write(*,*)
write(*,*)"================= Spin options ======================"
write(*,*)"1: No Spin-polarized calculation                     "
write(*,*)"2: Spin-polarized calculation                        "
write(*,*)
write(*,*) "------------->>"
!read(*,*)ispin
ispin=1
write(*,*)
!write(*,*)"input energy range"
!write(*,*) "------------->>"
!read(*,*)emin,emax
write(*,*) "------------->>"
!write(*,*)"input how many kpoints in calculation:"
write(*,*) "------------->>"
!read(*,*)nkpt

!write(*,*)"input amount of bands:"
write(*,*) "------------->>"
!read(*,*)nband

!write(*,*)"input how many inos in calculation:"
write(*,*) "------------->>"
!read(*,*)nion

write(*,*) "------------->>"
write(*,*)"how many atom you want to plot"
read (*,*)natom
write(*,*)natom

call readfermi(ef)
write(*,*)ef
write(*,*) 'Which atom you want in your projection band ' 
write(*,*)

allocate (atom(natom))

do i =1,natom
read(*,*)atom(i)
enddo

write(*,*)
write(*,*) 'Enter the scaling factor in your projection band ' 
write(*,*)
read(*,*) quanzhong
write(*,*)

!if (ispin.eq.1) then 
!open(11,file='band-p.dat') 
!elseif (ispin.eq.2) then 


!endif 
!define dimension
!allocate(aaa(3,3))
!allocate(b(3,3))
!read lattice constant from POSCAR
open (77,file='POSCAR',status='old')
read(77,*)
read(77,*)aa
do i=1,3
read(77,*) (aaa(i,j),j=1,3)
enddo
do i=1,3
 do j=1,3
 aaa(i,j)=aa*aaa(i,j)
 enddo
enddo
!500  format (3f12.8)
!calcalate volume
volume=aaa(1,1)*aaa(2,2)*aaa(3,3)+aaa(1,2)*aaa(2,3)*aaa(3,1)+aaa(1,3)*aaa(2,1)*aaa(3,2)-aaa(1,1)*aaa(2,3)*aaa(3,2)-aaa(1,2)*aaa(2,1)*aaa(1,3)-aaa(1,3)*aaa(2,2)*aaa(3,1)
do i=1,3
 if (i .eq. 1) then
    j=2
	k=3
 else if (i .eq. 2) then
    j=3
	k=1
 else
    j=1
	k=2
 endif
c(1)=aaa(j,2)*aaa(k,3)-aaa(j,3)*aaa(k,2)
c(2)=aaa(j,3)*aaa(k,1)-aaa(j,1)*aaa(k,3)
c(3)=aaa(j,1)*aaa(k,2)-aaa(j,2)*aaa(k,1)
do j=1,3
  b(i,j)=2*pi*c(j)/volume
  !write(*,*)b(i,j)
enddo
enddo
!read eigenval and kpoints coordinates
open(7,file='PROCAR',form='FORMATTED',status='OLD')
read(7,*)
read(7,104)nkpt,nband,nion
write(*,*)nkpt,nband,nion
allocate(eig(ispin,nkpt,nband),kpt(nkpt,3),wei(ispin,nkpt,nband,nion,10),xx(nkpt),ee(ispin,nkpt,nband),twei(ispin,nkpt,nband))
backspace(7)
!read(7,*)
do i=1,ispin
read(7,*)
!write(*,*)nkpt
   do k=1,nkpt
   read(7,*)
   !write(*,*)k
   read(7,105)kp,(kpt(k,j),j=1,3),w
   !write(*,*)kpt(k,1)
   !write(6,fmt='(f12.8)')(kpt(k,1))
   read(7,*)
     do nb=1,nband
	 read(7,106)eig(i,k,nb)
	 !write(*,*)eig(i,k,nb)
	 read(7,*)
	 read(7,*)
	   do ni=1,nion
	   read(7,107)wei(i,k,nb,ni,1:10)
	   enddo
	 read(7,*)
	 read(7,*)
	 enddo
   enddo
   !write(*,*)kpt(2,1:3)
   !write(*,*)kpt(2,3)-kpt(2,1)
   !continue
!   *** find reciprocal lattice vector ***
   xx(1) = 0.0
   nn = 1
   do k=1,nkpt-1
         dkx=(kpt(k+1,1)-kpt(k,1))*b(1,1) + (kpt(k+1,2)-kpt(k,2))*b(2,1)       &
     &   + (kpt(k+1,3)- kpt(k,3))*b(3,1)                                  
      dky=(kpt(k+1,1)-kpt(k,1))*b(1,2) + (kpt(k+1,2)-kpt(k,2))*b(2,2)       &
     &   + (kpt(k+1,3)- kpt(k,3))*b(3,2)                                  
      dkz=(kpt(k+1,1)-kpt(k,1))*b(1,3) + (kpt(k+1,2)-kpt(k,2))*b(2,3)       &
     &   + (kpt(k+1,3)- kpt(k,3))*b(3,3)
	 !write(*,*)k                                  
      del =  sqrt ( dkx**2 + dky**2 + dkz**2 ) 
	 ! write(*,*)del
      nn = nn +1 
      xx(nn) = xx(nn-1) + del 
	  !write(*,*)xx(nn)
   enddo
   
   if(i==1)then
   open(50,file="atom.dat",status='replace')
  
    do nb=1,nband
     do k=1,nkpt
	 ee(i,k,nb) = eig(i,k,nb)-ef
!write(*,*)natom
	   do x=1,natom
	     twei(i,k,nb)=twei(i,k,nb)+wei(i,k,nb,atom(x),10)
	   enddo
!	 write(*,*)twei(i,k,nb)
	 !if (ee(i,k,nb) .gt. emax) ee(i,k,nb)=emax
	 !if (ee(i,k,nb) .lt. emin) ee(i,k,nb)=emin
     if(twei(i,k,nb)>quanzhong)write(50,300)xx(k),ee(i,k,nb),twei(i,k,nb)
    
	 enddo
	 write(50,*)
	
    enddo

   endif

 if(i==2)then
  open(51,file="atom.2.dat",status='replace')
    do nb=1,nband
     do k=1,nkpt
	 ee(i,k,nb) = eig(i,k,nb)-ef
!write(*,*)natom
	   do x=1,natom
	     twei(i,k,nb)=twei(i,k,nb)+wei(i,k,nb,atom(x),10)
	   enddo
!	 write(*,*)twei(i,k,nb)
	 !if (ee(i,k,nb) .gt. emax) ee(i,k,nb)=emax
	 !if (ee(i,k,nb) .lt. emin) ee(i,k,nb)=emin
     if(twei(i,k,nb)>quanzhong)write(51,300)xx(k),ee(i,k,nb),twei(i,k,nb)
    
	 enddo
	 write(51,*)
	
    enddo

   endif

800 format (f10.6,f11.6)
300 format (f11.6,f12.6,f13.5)
!800 format (21x)
enddo
!write(*,*)wei(1,1,4,1:9)
open(15,file='line.dat',status='replace')
!write(15,*)
nk=nkpt/ndiv
do n=1,nk
kk=n*ndiv
write (15,400)xx(kk),eig(1,kk,1)-ef
write (15,400)xx(kk),eig(1,kk,nband)-ef
write(15,*)
enddo
write (15,400)xx(1),eig(1,1,1)-ef
write (15,400)xx(1),eig(1,1,nband)-ef
write(15,*)
zero=0.0
write (15,400) xx(1),zero
write (15,400) xx(nkpt),zero

400 format (f10.6,f12.6)
!find reciprocal lattice vector!!!


101 format(10x,f9.5)
102 format(f10.5)
103 format(20a4)
104 format(16x,i3,20x,i5,19x,i4)
!105 format(18x,3f11.8,24x)
105 format(10x,i3,5x,3f11.8,13x,f11.8)
106 format(17x,f14.8,19x)
107 format(3x,10f7.3,7x)
deallocate(eig,kpt,xx,ee)
close(50)
close(51)
close(52)
close(53)
close(54)
close(55)
close(60)
close(61)
close(62)
close(63)
close(64)
close(65)

end

!--------subroutine redafermi-----------------
subroutine readfermi(fermi)
! Copyright (C) 2009 W. Wang
! This file is distributed under the terms of the GNU General Public License.
! See(i,k,nb) the file COPYING for license details.
implicit none
real*8 :: fermi,normalization,enmax,enmin,step
integer :: i,rows,columns,nth,iserror
open(unit=11,file='DOSCAR',form='formatted',status='old',iostat=iserror)
if (iserror>0) then
    write(*,*)"The DOSCAR file does not exist"
    goto 200
end if

!kptip iformation of start
do i=1,5
   read (11,*)
end do

write(*,*)
!read fermi energy and rows
read (11,*)enmax,enmin,rows,fermi,normalization
step=(enmax-enmin)/(rows-1)
close(11)
return
200 end subroutine
