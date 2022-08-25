!
!  Notice === A.B =transpose(matmul(transpose(A),transpose(B))  = matmul(B,A)
!
!
!
! This codes will realize:
! (1) the Rotation of coordination:
!   R1 ===L===R2 in which L means ligand atom, usually IIV(Oxegen) and IIIV(Cl) types atoms
!   R1 and R2  are transition metal atoms. Usually, the angle among the
!   cluster is not pi/2, but the bond R1L  is nearly perpendicular to the bond
!   LR2,   the configuration is shown as :
!                          R1 ====== L
!                                    ||
!                                    ||
!                                    ||
!                                    R2
!
! (2) generation of the effective hopping matrix along different cluster (path) :
!    
! T1ij
! T2ij
! T3ij
!
! (3) generation of the effective on-site energy matrix
!
! eetg1 =
! eeg1=
! eetg2 =
! eeg2 =
!
!
! (4) generation of the explicit.fdf file for lat_mod.x
!  including  hopping matrix and on-site energy
!
!!
!! two input files  POSCAR.0  and para.dat you should prepare
!! POSCAR must be in fractional form 
!! two outfile  explicit.fdf, newposar.0 
!! explicit.fdf 
!!

!!
!! we will deal with four poscars, i.e., CrI3, WI3, CrWI6, CrGeTe6
!! the former three poscars have the same structure   I6Cr2, I6W2, I6CrW
!! we prepare the poscar should as above structure..
!! the forth poscar is different from the former structure, which includes ten
!! atoms  Te6 Cr2 Ge2
!!
!! note that the transition metal atoms must be in 7 and 8 sites of poscar...
!!

program main
    implicit none
    integer :: id, iu ! unit id  and iu
    real(8) :: hd(3,3), hd_new(3,3) ! head_poscar(3,3)
    real(8) :: unit_mat(3,3) ! identity matrix
    real(8) :: iv(3),jv(3),kv(3),ivnew(3),jvnew(3), kvnew(3) !
    real(8) :: av(3),bv(3),cv(3), t1(3), t2(3),t3(3)!
	real(8) :: z1(3),z2(3),z3(3),x1(3),x2(3),x3(3),Rz(3:3),Rx(3:3),z(3)=(0.0,0.0,1.0), x(3) = (1.0,0.0,0.0),pos_car_z(3,8),pos_car_x(3,8),z3_n(3),x3_n(3)!njy_Roration_oct
    !real(8) :: pos(3,8), pos_car(3,8)
    real(8), allocatable :: pos(:,:), pos_car(:,:)
    real(8), parameter :: autoa = 0.529177249d0
    real(8), parameter :: rad2deg = 180.0d0 / 3.1415926d0 
    real(8), parameter :: deg2rad = 3.1415926d0 / 180.0d0
    real(8), parameter ::   = 0.1d0
    character(120) :: ctmp
    integer :: iline,iloop,jloop! dummy index
    integer :: itmp
    integer, allocatable :: indx_nest(:)
    real(8) :: modulus
    logical :: lfrac  = .TRUE.
    real(8) :: rtmp(3)
    real(8) :: angle_vec
    integer :: lvert
    real(8) :: R2R1(3) ! R2 - R1
    real(8) :: taxis(3), alpha
    real(8) :: rotmat(3,3), Rmat(3,3)
    real(8) :: rescale1, rescale2
    real(8) :: matrix1(3,3),matrix2(3,3), matrix3(3,3)
    logical :: ltest
    !
    real(8), allocatable :: hmat1(:,:), hmat2(:,:), hmat3(:,:)
    real(8) :: vpdp_s, vpds_s, vddd, vddp,vdds, vpdp,vpds
    real(8),allocatable :: et2g(:), eeg(:)
    integer :: lsimp = 0
    real(8) :: cr3(3),cr4(3),cr3_cart(3), cr4_cart(3)
    
    real(8) :: deltaE = 1.0d0
    integer :: npara
    !
    integer :: is(2), norb
    real(8) :: diff_vec(3,5)
    real(8) :: Hcr1site(5,5),Hcr2site(5,5)
    real(8) :: tzero
    character(20) :: file_in, file_out, file_tmp
    logical :: alive
    integer :: num_atom
    character(2) :: sym_atom(2) ! which describe the transition metal atoms
    character(2) :: ligand
    real(8), allocatable :: new_pos_car(:,:)
    real(8), allocatable :: new_pos_frac(:,:)
    
    !
    ! and read from para.dat
    !
    
    ! iu = 10
    !file_in = 'para.dat'
    !alive = .FALSE.
    !inquire(file = file_in, exist = alive)
    !if(alive) then
    !    open(iu,file = file_in)
    !    print*, 'open ', file_in, 'successfully!'
    !else
    !   stop 'do not exist input file para.dat! please prepare the file!'
    !end if
    !read(iu, *) file_tmp
    !print*, 'we will deal with ', file_tmp  
    !read(iu, *) sym_atom(:)
    !read(iu, *) ligand
    !read(iu, *) deltaE
    !
    ! para.dat
    ! ICr or IW or ICrW or TeCrGe
    !
    if( trim(file_tmp) == 'TeCrGe') then
        num_atom = 10
        allocate(pos(3,num_atom))

        allocate(new_pos_car(3,num_atom))
        allocate(new_pos_frac(3,num_atom))
  
        ! position array reading from POSCAR in fractional coordination

        allocate(pos_car(3, num_atom)) 

        ! position array in cartesian form
    else
        num_atom = 8
        allocate(pos(3,num_atom)) 

        ! position array reading from POSCAR in fractional coordination
        
        allocate(pos_car(3,num_atom))

        allocate(new_pos_car(3,num_atom))
        allocate(new_pos_frac(3,num_atom))

        ! position array in cartesian form
    
    end if
    
    pos_car = 0.0d0
    unit_mat= 0.0d0
    do iline = 1, 3
        unit_mat(iline,iline) = 1.0d0
    end do

    id = 20
    alive = .FALSE.
    file_in = 'POSCAR.0'

    inquire(file=file_in, exist = alive)
    if(alive) then
        open(id,file = file_in)
        print*, 'open ', file_in, 'successfully!'
    else
       stop 'do not exist input file POSCAR.0! please prepare the file!'
    end if

    read(id,*) ctmp
    read(id,*)
    print*, ctmp
    do iline = 1,3
        read(id,*) hd(:,iline)
        write(*,'(6x, 3f18.12)') hd(:,iline)
    end do
    read(id,*)
    read(id,*)
    read(id,*)
    do iline =1, num_atom
        read(id,*) pos(:,iline)
        write(*,'(6x, 3f18.12)') pos(:,iline)
    end do
    close(id)
    !
    ! read poscar over
    !
    print*, ''
    print*, '=========  information about lattice ================ ' 
    print*, ''
    call cell_info(hd)
    print*, ''
    print*, '=========  information about lattice ================ ' 
    print*, ''
    !
    ! print the unit cell informations
    !
    iv(:) = unit_mat(:,1)  ! x direction unit vector in row form
    jv(:) = unit_mat(:,2)  ! y direction unit vector in row form 
    kv(:) = unit_mat(:,3)  ! z direction unit vectro in row form

    print*, 'Header of poscar'
    write(*, '(6x, 3f18.12)') hd(1,1)*iv + hd(2,1)*jv + hd(3,1)*kv
    write(*, '(6x, 3f18.12)') hd(1,2)*iv + hd(2,2)*jv + hd(3,2)*kv
    write(*, '(6x, 3f18.12)') hd(1,3)*iv + hd(2,3)*jv + hd(3,3)*kv

    ! we should transfer the fractional coordination into cartesian coordination to identify the orientation
    ! frac to cart frac2cart()
    ! cart to frac cart2cart()
    ! 
    print*, '     ======== cartesian coordinate ============'
    do iline = 1, num_atom
        call frac2cart( hd, pos(:,iline), pos_car(:,iline))
        write(*,'(6x,3f18.10)') pos_car(:,iline)
    end do
    ! 
    lfrac = .TRUE.
    if(lfrac ) then
        pos = 0.0d0
        print*, '    ======= fractional coordinate ============ '
        do iline = 1, num_atom
            call cart2frac( hd, pos_car(:,iline), pos(:,iline))
            write(*,'(6x,3f18.12)') pos(:,iline)
        end do
    end if

    !  in all poscar, the transition metal atoms Cr and W is in 7,8 position.
    !  find the nearest neighbour of Cr1
    print*,''
    write(*,'(6x,2a8)') "R1 === > ", pos_car(:,7)
    do iline = 1, 6
        print*, modulus(pos_car(:,iline) - pos_car(:,7),3)
    end do
    !
    !  pos(:,7) and pos(:,8) describe Cr atoms
    !  pos(:,1) to pos(:,6) denote I atoms
    !
    print*,''
    write(*,'(6x,2a8)')'R2 === > ', pos_car(:,8)
    do iline = 1, 6
        print*, modulus(pos_car(:,iline) - pos_car(:,8),3)
    end do
 
1000 format(6x,3a4,i2,4a4,i2,a25)

     print*, 'Rotation z-axis start'
     z1(:) = pos_car(:,1) - pos_car(:,7)
     z2(:) = pos_car(:,4) - pos_car(:,7)
	 z3(:) = z1(:) - z2(:)!z3 present the axis along z1+z2
	 call normalize(z3,z3_n,3)
	 alpha_z = angle_vec(z3_n,z,3)*rad2deg
	 write(*,'(6x,a15,3x, f12.6,a10)')'angle = ', angle_z,'degree'
	 call cross_product(z,z3_n,rtmp)
     call normalize(rtmp,kvnew,3)
	 write(*,'(6x,a20,6x,3f15.8)') 'Rotation axis is : ',kvnew 
	 call rot_mat(alpha_z,kvnew,rotmat)
	 do i=1,3
	    write(*,'(3f12.6)') rotmat(:,i)
     end do 
	 do i=1,8 
	    pos_car_z(:,i) =  matmul(pos_car(:,i),rotmat)
     end do 
	 do i=1,8
	   write(*,'(6x, 3f18.12)') pos_car_z(:,i)
	 end do 
	 print*, 'Rotation z-axis over'
	 
	 print*,'Rotation x-axis start'
	 x1(:) = pos_car(:,6) - pos_car(:,7)
     x2(:) = pos_car(:,3) - pos_car(:,7)
	 x3(:) = x1(:) - x2(:) ! x3 present the axis along x1+x2
	 call normalize (x3,x3_n,3)
	 write(*,'(6x,a15,3x, f12.6,a10)')'angle = ', angle_vec(x,x3_n,3) * rad2deg,'degree'
	 alpha_x = angle_vec(x,x3_n,3)*rad2deg
	 write(*,'(6x,a15,3x, f12.6,a10)')'angle = ', angle_x,'degree'
	 !call cross_product(x,x3_n,rtmp)
     !call normalize(rtmp,avnew,3)
	 !write(*,'(6x,a15,3x,f15.8)') 'Rotation axis is :', avnew
	 !call cross_product(x,z3,3)
     !call normalize(rtmp,kvnew,3)
	 call rot_mat(alpha_x,kvnew,rotmat) ! Rotation along z axis
	 do i=1,3
	    write(*,'(3f12.6)') rotmat(:,i)
     end do 
	 do i=1,8 
	    pos_car_x(:,i) =  matmul(pos_car_z(:,i),rotmat)
     end do 
	 do i=1,8
	   write(*,'(6x, 3f18.12)') pos_car_x(:,i)
	 end do 
	 print*, 'Rotation x-axis over'
	

end program


function modulus( x, n )
    implicit none
    integer, intent(in) :: n
    real(8),intent(in) ::x(n)
    real(8) :: modulus
    integer :: i
    real(8) :: rtmp
    modulus = 0.0
    rtmp = 0.0
    do i = 1, n
        rtmp = rtmp + x(i)**2
    end do
    modulus = dsqrt(rtmp)
    return
end function

function angle_vec(x1,x2,n)! assuming n =2 or 3,  return angle_vec in unit of rad
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: x1(n), x2(n) ! x1 and x2 must have the same shape
    real(8) :: angle_vec
    real(8) :: modulus
    angle_vec = acos(dot_product(x1,x2)/(modulus(x1,n)*modulus(x2,n)))
    return
end function

subroutine frac2cart(abc, posin, posout)
    implicit none
    real(8), intent(in) :: abc(3,3)
    real(8), intent(in) :: posin(3)
    real(8), intent(out) :: posout(3)
    !
    !
        posout(:) = abc(:,1) * posin(1) + abc(:,2) * posin(2) + abc(:,3) * posin(3)

end subroutine
subroutine cross_product(x1_in, x2_in, x_out)! juse suitable to the three dimensional vectors.
    implicit none
    real(8), intent(in) :: x1_in(3), x2_in(3)
    real(8), intent(out) :: x_out(3)
    x_out(1) = x1_in(2) * x2_in(3) - x1_in(3) * x2_in(2)
    x_out(2) = -(x1_in(1) * x2_in(3)- x1_in(3) *x2_in(1))
    x_out(3) = x1_in(1)*x2_in(2) - x1_in(2) * x2_in(1)

end subroutine


subroutine cart2frac(abc, posin, posout)
    implicit none
    real(8), intent(in) :: abc(3,3)
    real(8), intent(in) :: posin(3)
    real(8), intent(out) :: posout(3)
    real(8) :: vol, u,v,w
    real(8) :: vtmp(3)! temperory vector
    !
    !  abc is an oblique coordination:
    !  abc(1,1) abc(2,1) abc(3,1)
    !  abc(1,2) abc(2,2) abc(3,2)
    !  abc(1,3) abc(2,3) abc(3,3)
    !
    call cross_product(abc(:,2),abc(:,3), vtmp(:)) ! sigma_a = b x c
    vol = dot_product(abc(:,1), vtmp(:))

    if( vol <= 0.0) then
       stop "Error, left hand coordination is adopted, please check POSCAR"
    end if
    u = dot_product(posin,vtmp(:))/vol

    call cross_product(abc(:,3),abc(:,1), vtmp(:)) ! sigma_b = c x a
    v = dot_product(posin,vtmp(:))/vol

    call cross_product(abc(:,1),abc(:,2), vtmp(:)) ! sigma_c = a x b
    w = dot_product(posin,vtmp(:))/vol

    posout(1) = u
    posout(2) = v
    posout(3) = w
end subroutine

subroutine cell_info(head)
     implicit none
	real(8), intent(in) :: head(3,3)
	real(8) ::  a,b,c, vol
	real(8) ::  sigma(3)
	real(8) ::  modulus,angle_vec
        real(8) ::  rad2deg_
        
        rad2deg_ = 180d0 / 3.1415926d0 
    
        a = modulus(head(:,1),3)
        b = modulus(head(:,2),3)
        c = modulus(head(:,3),3)
        write(*,'(a10,6x,f15.6)') 'a = ', a
        write(*,'(a10,6x,f15.6)') 'b = ', b
        write(*,'(a10,6x,f15.6)') 'c = ', c

        write(*,'(a40,6x,f12.4)') ' the angle between a and b vector =',angle_vec(head(:,1),head(:,2),3) * rad2deg_
        write(*,'(a40,6x,f12.4)') ' the angle between b and c vector =',angle_vec(head(:,2),head(:,3),3) * rad2deg_
        write(*,'(a40,6x,f12.4)') ' the angle between c and a vector =',angle_vec(head(:,3),head(:,1),3) * rad2deg_
        call cross_product(head(:,2),head(:,3), sigma(:))
        write(*,'(a30,6x,f12.6)')  ' Volume of unit cell =', dot_product( head(:,1),sigma(:))

        return
end subroutine
subroutine normalize(x_in,x_out,n)
        implicit none
        integer, intent(in) :: n
	real(8), intent(in) :: x_in(n)
	real(8), intent(out) :: x_out(n)
	real(8) :: modulus
        integer :: i

        do i = 1, n
            x_out(i) = x_in(i)/modulus(x_in, n)
        end do

end subroutine

!subroutine rot_vec(x_in, x_out, axis)! 3D rotation with respect to axis
!    implicit none
!    real(8), intent(in), x_in(3)
!    real(8), intent(in), axis(3)

!end subroutine
        subroutine rot_mat(theta, axis,mat) 
        ! note that axis is the unit vector, theta/180 is in unit of rad
        implicit none
        real(8), intent(in) :: theta ! [-pi, pi]
        real(8), intent(in) :: axis(3)
        real(8), intent(out) :: mat(3,3)
        real(8) :: theta1
        real(8) :: modulus
        real(8) :: axis1(3)
        real(8) :: a, b, c

        theta1  = theta*3.1415926d0 / 180.0d0
        if(abs(modulus(axis,3)-1.0)<0.0000001d0) then
                print*, 'input axis is unit vector'
                axis1(:) = axis(:)
        else 
                print*, 'input axis is not unit vector, axis will be normalized '
                call  normalize(axis,axis1,3)
        end if
        a = axis1(1)
        b = axis1(2)
        c = axis1(3)
        mat = 0.0d0
        if( (a**2 + b**2) == 0.0d0 ) then
                mat(1,1) = a**2 + dcos(theta1)
        else
                mat(1,1) =  a**2 + (b**2 + (a**2) * (c**2)) * dcos(theta1)/(a**2 + b**2)
        end if
        mat(2,1) =  a * b - a * b* dcos(theta1) - c* dsin(theta1) 
        mat(3,1) =  a * c - a * c* dcos(theta1) + b* dsin(theta1)

        mat(1,2) = a*b - a * b *dcos(theta1) + c* dsin(theta1)
        if( (a**2 + b**2) == 0.0d0 ) then
                mat(2,2) = b**2 + dcos(theta1)
        else
                mat(2,2) = b**2 + ((a**2) + (b**2) * (c**2))*dcos(theta1)/(a**2+ b**2)
        end if
        mat(3,2) = b*c - b*c* dcos(theta1)- a * dsin(theta1) 

        mat(1,3) = a*c - a * c* dcos(theta1) - b*dsin(theta1)
        mat(2,3) = b*c - b * c * dcos(theta1) + a*dsin(theta1)
        mat(3,3) = - a**2 - b**2 + (a**2 + b**2)*dcos(theta1) + 1.0d0

        end subroutine

