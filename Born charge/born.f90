program borncharge
   implicit none
   real(8) :: b(100,3), a(100,3), diff(100,3), P(3), lattice(3,3), diff_xyz(100,3), born(300,3) 
   integer :: i, j, k, num, m, n 
   integer :: Atom=18
   open(10,file='POSCAR', status='old')
   do i= 1,7   
   read(10,*)
   end do
  
   do i= 1, Atom
   read(10,*) b(i,1), b(i,2), b(i,3)
   end do
   close(10) 
   400 format( 1X,f19.16,1X,f19.16,1X,f19.16)
   
   open(10,file='tmp') 
   do i= 1, Atom
   write(10,400) b(i,1), b(i,2), b(i,3)
   end do
   close(10)
  open(10,file='POSCAR.sym')
   do i= 1, Atom
   read(10,*) a(i,1), a(i,2), a(i,3)
   end do
   close(10)

!   open(10,file='POSCAR.sym') 
!   do i=1, Atom
!      do j = 1,3
!        if (abs(b(i,j)-0.0d0)< 0.2) then
!           a(i,j) = 0.0d0
!        else if (abs(b(i,j)-0.25d0)< 0.2) then
!           a(i,j) = 0.25d0
!        else if (abs(b(i,j)-0.50d0)< 0.2) then
!           a(i,j) = 0.50d0
!        else if (abs(b(i,j)-0.75d0)< 0.2) then
!           a(i,j) = 0.75d0
!        else if (abs(b(i,j)-1.0d0)< 0.2) then
!           a(i,j) = 1.0d0
!       end if
!      end do
!   write(10,400) a(i,1), a(i,2), a(i,3)
!   end do
! a  close(10)

   open(10,file="POSCAR.diff")
   do i=1, Atom
     do j =1,3
     diff(i,j)=b(i,j)-a(i,j)
     end do
   write(10,400) diff(i,1), diff(i,2), diff(i,3)
   end do
   close(10)
  
   open(10,file="POSCAR",status='old')
   read(10,*)
   read(10,*)
      do i = 1,3
        read(10,*) lattice(i,1), lattice(i,2), lattice(i,3)  
        write(*,400) lattice(i,1), lattice(i,2), lattice(i,3)  
      end do
   close(10)

  open(10,file="POSCAR.diff.xyz")
       do i=1, Atom
          do j =1,3
          diff_xyz(i,j)= diff(i,1)* lattice(1,j) + diff(i,2)*lattice(2,j) + diff(i,3)* lattice(3,j) 
          end do
       write(10,400) diff_xyz(i,1),diff_xyz(i,2),diff_xyz(i,3)
       end do 
  close(10)       

 open(10,file="Zeff",status="old")
 k=1
 do j=1, Atom
    read(10,*)
    do i=1,3
    read(10,*) num, born(k,1), born(k,2), born(k,3)
    k=k+1 
    end do
 end do
 close(10)
 P=0.0d0 
 n=1 
 do i=1, Atom
 P(1)= P(1) + diff_xyz(i,1)*born(n,1) +   diff_xyz(i,2)*born(n,2) +   diff_xyz(i,3)*born(n,3)
 P(2)= P(2) + diff_xyz(i,1)*born(n+1,1) + diff_xyz(i,2)*born(n+1,2) + diff_xyz(i,3)*born(n+1,3)
 P(3)= P(3) + diff_xyz(i,1)*born(n+2,1) + diff_xyz(i,2)*born(n+2,2) + diff_xyz(i,3)*born(n+2,3)
 n=n+3 
 end do 

  write (*,*) P 
  write (*,*) P*1600/343.7147
end    
