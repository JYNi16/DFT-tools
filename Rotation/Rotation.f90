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
    
    iu = 10
    file_in = 'para.dat'
    alive = .FALSE.
    inquire(file = file_in, exist = alive)
    if(alive) then
        open(iu,file = file_in)
        print*, 'open ', file_in, 'successfully!'
    else
       stop 'do not exist input file para.dat! please prepare the file!'
    end if
    read(iu, *) file_tmp
    print*, 'we will deal with ', file_tmp  
    read(iu, *) sym_atom(:)
    read(iu, *) ligand
    read(iu, *) deltaE
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
    write(*,'(6x,2a8)') "R1 === > ", sym_atom(1)
    do iline = 1, 6
        print*, modulus(pos_car(:,iline) - pos_car(:,7),3)
    end do
    !
    !  pos(:,7) and pos(:,8) describe Cr atoms
    !  pos(:,1) to pos(:,6) denote I atoms
    !
    print*,''
    write(*,'(6x,2a8)')'R2 === > ', sym_atom(2)
    do iline = 1, 6
        print*, modulus(pos_car(:,iline) - pos_car(:,8),3)
    end do
   
    ! now we find the atoms between two transition metal atoms
    print*,''
    itmp = 0
    do iline = 1, 6
     if( abs( modulus(pos_car(:,iline) - pos_car(:,8),3) - modulus(pos_car(:,iline) - pos_car(:,7),3))<my_prec) then
        print*, 'index of ligand atom',iline
        itmp = itmp + 1
     end if
    end do
    iloop = 1
    if(itmp /= 2) then
       stop 'something about ligland atoms must be wrong, please check!'
    else 
       allocate(indx_nest(itmp))
       do iline = 1, 6
         if( abs( modulus(pos_car(:,iline) - pos_car(:,8),3) - modulus(pos_car(:,iline) - pos_car(:,7),3))<my_prec) then
           indx_nest(iloop) = iline
           iloop = iloop + 1
         end if
       end do
    end if
    print*, 'index of nearest of ligand atom ', indx_nest(:)


! Cr1 - I5 - I6 - I1 form right hand coordination
! I5 - Cr1 === I6 - Cr1
    print*, 'original bond '//sym_atom(1) //ligand, indx_nest(1)
    t1(:) = pos_car(:,indx_nest(1)) - pos_car(:,7)
    print*, 'original bond '//sym_atom(2) //ligand, indx_nest(2)
    t2(:) = pos_car(:,indx_nest(2)) - pos_car(:,7)
    write(*,'(6x,a15,3x, f12.6,a10)')'angle = ', angle_vec(t1,t2,3) * rad2deg,'degree'
    call cross_product(t1, t2, rtmp)
    call normalize(rtmp,kvnew,3)
    write(*,'(6x,a20,6x,3f15.8)') 'new k direction is : ',kvnew
    R2R1 = pos_car(:,8) - pos_car(:,7)
    
    print*, 'kvnew . R2R1 , which should be zero!'
    write(*,'(f12.5)') dot_product(kvnew, R2R1)

    if( abs(dot_product(kvnew,R2R1)) < my_prec)then 
        lvert = 1
    else
        lvert = 0
    end if

    if(lvert == 1) then
        write(*,1000) sym_atom(1), '==', ligand, indx_nest(1),' == ', sym_atom(2),'==',&
                &ligand, indx_nest(2),'  is in the same plane'
    else
        write(*,1000) sym_atom(1), '==', ligand, indx_nest(1),' == ', sym_atom(2),'==',&
                &ligand, indx_nest(2),' is not in the same plane!'
    end if
1000 format(6x,3a4,i2,4a4,i2,a25)

!
! now we testify the rotation matrix
!
    taxis(1) = 0.0d0
    taxis(2) = 0.0d0
    taxis(3) = 1.0d0
    alpha = 90.0d0
    print*,''
    print*,'======== test rot_mat subroutine ========'
    call rot_mat(alpha,taxis,rotmat) ! counterclockwise rotation about axis with angle alpha
    do iloop = 1,3
        write(*,'(3f12.6)') rotmat(:,iloop)
    end do
    print*,'======== test rot_mat subroutine end ========'
    print*, ''
!
! testing end ...
!   
    if( angle_vec(t1,t2,3) * rad2deg <= 90.0d0) then
        alpha = (90 - angle_vec(t1,t2,3) * rad2deg) /2.0d0 ! 
    else
        alpha = (90 - angle_vec(t1,t2,3) * rad2deg) /2.0d0 !
    end if 
    write(*,'(a15,f12.6,a10)')'alpha = ', alpha,'degree'
    !
    !  alpha is the deviation angle ....
    !  no rotation is performed until now ....
    !


    ! =====
    call rot_mat(-alpha,kvnew,rotmat)
    ! t1 = I5 - Cr1
    print*, ''
    print*, 'Rotation matrix :'
    do iloop = 1, 3
        write(*,'(3f12.8)') rotmat(:,iloop)
    end do
    print*,''
    print*, 'old t1', t1
    t1 = matmul(transpose(rotmat),t1)
    print*, 'new t1', t1
    !
        
        rescale1 = (modulus(R2R1,3)*dsqrt(2.0d0)/2.0d0) / modulus(t1,3)
        print*,'rescale1 = ', rescale1
    !
    ! ...  is it 5 ?
    ! 
        pos_car(:,5) = pos_car(:,7) + rescale1 * t1(:) ! old I5
        new_pos_car(:,5) = pos_car(:,5)! new I5
        new_pos_car(:,7) = pos_car(:,7)
        new_pos_car(:,8) = pos_car(:,8)
        new_pos_car(:,2) = 2.0d0* new_pos_car(:,7) - new_pos_car(:,5)
    
    
        call normalize(t1,ivnew,3)
    ! ======
        call rot_mat(alpha,kvnew,rotmat)
    ! t2 = I6 - Cr1
        print*, 'old t2 ', t2
        t2 = matmul(transpose(rotmat),t2)
        print*, 'new t2', t2
        rescale2 = (modulus(R2R1,3)*dsqrt(2.0d0)/2.0d0) / modulus(t2,3)
        print*,'rescale2 = ', rescale2

        pos_car(:,6) = pos_car(:,7) + rescale2 * t2(:) ! old I6
        new_pos_car(:,6) = pos_car(:,6)
        new_pos_car(:,3) = 2.0d0* new_pos_car(:,7) - new_pos_car(:,6)

        new_pos_car(:,1) = new_pos_car(:,7) + kvnew * modulus(R2R1,3)*dsqrt(2.0d0)/2.0d0
        new_pos_car(:,4) = 2.0d0 * new_pos_car(:,7) - new_pos_car(:,1)

 
        call normalize(t2,jvnew,3)
    print*, '  ===== check the new unit direction x,y,z are pendicular to each other! ======'
    print*, modulus(ivnew,3)
    print*, modulus(jvnew,3)
    print*, modulus(kvnew,3)
    print*,'<<i|j>>', angle_vec(ivnew,jvnew,3) * rad2deg,'degree'
    print*,'<<j|k>>', angle_vec(jvnew,kvnew,3) * rad2deg,'degree'
    print*, angle_vec(ivnew,kvnew,3) * rad2deg,'degree'
    print*, angle_vec(kvnew,jvnew,3) * rad2deg,'degree'
    print*, angle_vec(jvnew,ivnew,3) * rad2deg,'degree'
    print*,'ivnew', ivnew
    print*,'jvnew', jvnew
    print*,'kvnew', kvnew
    print*, '   ========  check over!  ========= '
    Rmat(:,1) = ivnew(:)
    Rmat(:,2) = jvnew(:)
    Rmat(:,3) = kvnew(:)
    do iloop = 1 ,3
       write(*,'(3f12.7)') hd(:,iloop)
    end do
    ltest = .FALSE.
    if(ltest) then
        print*,' ==== now we check the instrinsic matmul function ===== '
    
        do iloop = 1, 3
                do jloop  = 1, 3
                        matrix1(jloop,iloop) = iloop * jloop + 2 * iloop
                end do
                write(*,'(3f12.7)') matrix1(:,iloop)
        end do

        do iloop = 1, 3
                do jloop  = 1, 3
                        matrix2(jloop,iloop) = iloop * jloop + 5 * iloop
                end do
                write(*,'(3f12.7)') matrix2(:,iloop)
        end do
        matrix3 = matmul(matrix2,matrix1)
        do iloop = 1, 3   
                write(*,'(3f12.7)') matrix3(:,iloop)
        end do
        print*, ' ==== check the instinsic matmul function over ===='
    end if
    hd_new = transpose(matmul(transpose(Rmat),hd)) ! R.hd^T
    do iloop  = 1,3
        write(*,'(3f15.8)') hd_new(:,iloop)
    end do
    !
    ! now we check the rotation of axis
    !
    hd = transpose(hd_new)
    print*,' new head of POSCAR!'
    do iloop  = 1,3
       write(*,'(3f15.7)') hd(:,iloop)
    end do
    !call cart2frac(hd)
    !
    !  rotating coordination is over!
    !
    print*,''
    ! pos_car(:,7) 
    call cart2frac(hd, matmul(transpose(Rmat),pos_car(:,7)),rtmp)
    write(*,'(6x,a4,3f12.6)') sym_atom(1), rtmp
    pos(:,7) = rtmp
    call cart2frac(hd, matmul(transpose(Rmat),pos_car(:,8)),rtmp)
    write(*,'(6x,a4,3f12.6)') sym_atom(2), rtmp
    pos(:,8) = rtmp
   
    ! 
    ! now we write new poscar for slater-koster
    !

    do iloop = 1, num_atom
        call cart2frac(hd, matmul(transpose(Rmat),new_pos_car(:,iloop)),new_pos_frac(:,iloop))
        write(*,'(6x, 3f12.6)') new_pos_frac(:,iloop)
    end do
    
    !id = 10
    rtmp = 1.0
    open(id, file = 'POSCAR_SK')
    write(id,*) 'ICr'
    write(id,*) rtmp(1)
    do iloop = 1, 3
        write(id,'(4x,3f15.8)') hd(:,iloop)
    end do
    write(id, '(2a8)') ligand,sym_atom(1)
    write(id, *)  num_atom - size(sym_atom), size(sym_atom)
    write(id, *) 'Direct'
    do iloop = 1, 8
        write(id,'(4x,3f15.8)') new_pos_frac(:,iloop)
    end do
          
        
    close(id)


    !
    !  now we should find the hopping matrix 
    !
    print*, ''
    print*, ''
    print*, ' ======== now we will set the hopping matrix along three path! ==========='
    
    if( sym_atom(1) == sym_atom(2) ) then
        print*, 'there is one class transition metal atoms in unit cell !'
        allocate( et2g(2) )
        allocate( eeg(2) )
        read(iu,*) npara
        !print*, npara
        if(npara == 7) then
                read(iu,*) et2g(1) 
                read(iu,*) eeg(1)
                et2g(2) = et2g(1)
                eeg(2) = eeg(1)
                !read(id,*) deltaE
                read(iu,*) vpds
                read(iu,*) vpdp
                read(iu,*) vdds
                read(iu,*) vddp
                read(iu,*) vddd
                read(iu,*) lsimp
         else
                stop 'parameter is set incorrectly!'
         end if
    else 
        print*, 'there is two class transition metal atoms in unit cell !'
        allocate( et2g(2) )
        allocate( eeg(2) )
        read(iu,*) npara
        if( npara == 9 ) then
                read(iu,*) et2g(1) 
                read(iu,*) eeg(1)
                read(iu,*) et2g(2) 
                read(iu,*) eeg(2)
                !read(id,*) deltaE
                read(iu,*) vpds
                read(iu,*) vpdp
                read(iu,*) vdds
                read(iu,*) vddp
                read(iu,*) vddd
                read(iu,*) lsimp
        else
                stop 'parameter is set incorrectly!'
        end if
    end if 
    close(iu)
    ! read para.dat over
    print*, 'The first path',sym_atom(1), '==',ligand, indx_nest(1),'==',sym_atom(2), '==',ligand, indx_nest(2),&
             &    'new coordinations'
   

    select case (lsimp)
    case(0) ! no simplify
         print*, 'nothing to do, we count the parameters Vpdp, Vpds, Vdds, Vddp, Vddd, DeltaE, Et2g and eeg!'
         vpdp_s = vpdp**2
         vpds_s = vpds**2
    case(1) ! partial simplify
         print*, 'we will set Vddd and Vddp to zero'
         vpdp_s = vpdp**2
         vpds_s = vpds**2
         vddp = 0.0d0
         vddd = 0.0d0
    case(2) ! full simplify
         print*, 'we will set Vpdp**2 and Vpds**2 to zero except case 1'
         vpdp_s = 0.0d0
         vpds_s = 0.0d0
         vddp = 0.0d0
         vddd = 0.0d0
    case default
          print*,'there is no actions related to simplification!'
    end select
    print*
    ! Total matrix  H(10, 10)
     write(*,'(6x,a4,a2,3f12.6)') sym_atom(1),'1', pos(:,7)
     write(*,'(6x,a4,a2,3f12.6)') sym_atom(2),'2', pos(:,8)
     write(*,'(6x,a20)')'R2 - R1 = :'
     rtmp = pos(:,8) - pos(:,7)
     write(*,'(3f15.10)') rtmp
     write(*,'(3f15.10)') rtmp(1) * hd(:,1) + rtmp(2) * hd(:,2) + rtmp(3) * hd(:,3)
     diff_vec(:,1) = rtmp(1) * hd(:,1) + rtmp(2) * hd(:,2) + rtmp(3) * hd(:,3)
     write(*,'(f15.10)') modulus(rtmp(1) * hd(:,1) + rtmp(2) * hd(:,2) + rtmp(3) * hd(:,3),3)

     allocate(hmat1(10,10))
     allocate(hmat2(10,10))
     allocate(hmat3(10,10))
     !
     ! here we will deal with the general condition...
     !  et2g(2) and eeg(2)
     !  et2g(1)  === > R1
     !  et2g(2)  === > R2
     !  eeg(1)   === > R1
     !  eeg(2)   === > R2 
     !  for all magnetic states, i.e., spin up and spin down, the hopping occurs
     !  in  the same spin state, hence, we need not deal with the Hamiltonian
     !  with all spin states simulataneously.

     hmat1= 0.0d0

     hmat1(1,1) = 2.0d0 *vpdp_s/deltaE + et2g(1)*1.d0 
     hmat1(6,1) = (vddd + 3.0d0*vdds)/4.0d0  
     hmat1(10,1) = dsqrt(3.0d0) *(vddd - vdds)/4.d0 + vpdp * vpds/deltaE 
     hmat1(2,2) = vpdp_s /deltaE + et2g(1)
     hmat1(7,2) = (vddd+ vddp)/2.0d0
     hmat1(8,2) = - vpdp_s/ deltaE -(vddd - vddp)/2.0d0
     hmat1(3,3) = hmat1(2,2)
     hmat1(7,3) = - vpdp_s/ deltaE -(vddd- vddp)/2.0d0
     hmat1(8,3) = (vddd+ vddp)/2.0d0
     hmat1(4,4) = 3.0d0*vpds_s/ (2.0d0*deltaE) + eeg(1)
     hmat1(9,4) = vddp
     hmat1(5,5) = vpds_s /(2.0d0* deltaE) + eeg(1)
     hmat1(6,5) = dsqrt(3.0d0)*(vddd - vdds)/4.0d0 + vpdp * vpds/deltaE
     hmat1(10,5) = (3.0d0 * vddd + vdds) /4.0d0
     !
     !do iloop = 1, 5
     !   hmat1(iloop+5,iloop+5) = hmat1(iloop,iloop)
     !end do
     !
     hmat1(6,6) = 2.0d0 * vpdp_s/deltaE + et2g(2)
     hmat1(7,7) = vpdp_s/deltaE + et2g(2)
     hmat1(8,8) = hmat1(7,7)
     hmat1(9,9) = 3.0d0 * vpds_s/2.0d0/deltaE + eeg(2)
     hmat1(10,10) = vpds_s /2.0d0/deltaE + eeg(2)
     
     ! here is real number, do not use conjugate ........
     do iloop = 6,10
         do jloop =  1,5
             hmat1(jloop,iloop) = hmat1(iloop,jloop)   
         end do
     end do
     print*, ' === Hamiltonian along first path R2-R1 === '
     do iloop = 1, 10
        write(*,'(10f8.4)') hmat1(:,iloop)
     end do
     print*,''
    !
    !
    ! 
    print*, ' now, we find the second nearest neighbour Cr1 atom! '

    !do iloop = -1, 1
    !    do jloop = -1 ,1
    !	    rtmp(1)  = pos(1,8) + iloop 
    !	    rtmp(2)  = pos(2,8) + jloop
    !	    rtmp(3)  = pos(3,8) 
    !        write(*,'(3f15.8)') 	    
    !	end do
    !end do
    
    cr3(1) = pos(1,8) 
    cr3(2) = pos(2,8) + 1.0d0
    cr3(3) = pos(3,8)
    print*, 'Cr3', cr3
    !call frac2cart(hd,cr3,cr3_cart)
    !print*, 'cr3_cart', cr3_cart
    rtmp = cr3(:) - pos(:,7)
    print*, 'Cr3 - Cr1'
    write(*,'(3f15.10)') rtmp
    write(*,'(3f15.10)') rtmp(1) * hd(:,1) + rtmp(2) * hd(:,2) + rtmp(3) * hd(:,3)
    diff_vec(:,2) = rtmp(1) * hd(:,1) + rtmp(2) * hd(:,2) + rtmp(3) * hd(:,3)
    write(*,'(f15.10)') modulus(rtmp(1) * hd(:,1) + rtmp(2) * hd(:,2) + rtmp(3) * hd(:,3),3)
    !
    !
    !
    print*,  'now we set the hopping matrix along the second path !'
    hmat2 = 0.0d0
    hmat2(1,1) = vpdp_s/deltaE + et2g(1)
    hmat2(6,1) = (vddd + vddp)/2.0d0
    hmat2(7,1) =  vpdp_s/deltaE + (vddd-vddp)/2.0d0
    hmat2(2,2) = hmat2(1,1)
    hmat2(6,2) =  vpdp_s/deltaE +(vddd - vddp)/2.0d0
    hmat2(7,2) = (vddd + vddp)/2.0d0
    hmat2(3,3) = 2* vpdp_s /deltaE + et2g(1)    
    hmat2(8,3) = (vddd + 3* vdds)/4.0d0
    hmat2(9,3) = 3*(vddd - vdds)/8.0d0 + dsqrt(3.0d0) * vpdp* vpds/2.0d0 / deltaE
    hmat2(10,3) = dsqrt(3.0d0)*(vddd - vdds)/8.0d0 + vpdp * vpds/2.0d0 / deltaE
    hmat2(4,4) = 3.0d0 * vpds_s/ 4.0d0/deltaE + eeg(1)
    hmat2(5,4) = -dsqrt(3.0d0) * vpds_s /4.0d0/deltaE 
    hmat2(8,4) = hmat2(9,3)
    hmat2(9,4) = 9.0*vddd/16.0d0 + vddp/4.0d0 + 3* vdds/16.0d0
    hmat2(10,4) = 3*dsqrt(3.0d0) * vddd/ 16.0d0-dsqrt(3.0d0)* vddp /4.0d0 + dsqrt(3.0d0)* vdds /16.0d0
    hmat2(4,5) = hmat2(5,4)
    hmat2(5,5) = 5.0d0 * vpds_s/4.0d0/deltaE + eeg(1)
    hmat2(8,5) = hmat2(10,3)
    hmat2(9,5) = hmat2(10,4)
    hmat2(10,5) = 3.0d0 * vddd / 16.0d0 + 3.0d0 * vddp / 4.0d0 + vdds / 16.0d0
    hmat2(10,9) = -dsqrt(3.0d0) * vpds_s /4.0d0/deltaE
    hmat2(9,10) = hmat2(10,9)
    !
    !do iloop = 1, 5
    !   hmat2(iloop+5,iloop+5) = hmat2(iloop,iloop)
    !end do
    hmat2(6,6) = vpdp_s /deltaE + et2g(2)
    hmat2(7,7) = hmat2(6,6)
    hmat2(8,8) =  2.0d0 * vpdp_s / deltaE + et2g(2)
    hmat2(9,9) = 3.0d0 * vpds_s / 4.0d0 / deltaE + eeg(2)
    hmat2(10,10) = 5.0d0 * vpds_s / 4.0d0 /deltaE + eeg(2)
    do iloop = 6,10
        do jloop =  1,5
            hmat2(jloop,iloop) = hmat2(iloop,jloop)  
    ! here is real number, do not use conjugate ........
        end do
    end do
    !
    print*, ' === Hamiltonian along second path Cr4 - Cr1 === '
    do iloop = 1, 10
        write(*,'(10f8.4)') hmat2(:,iloop)
    end do
    print*, ''
    cr4(1) = pos(1,8) - 1.0
    cr4(2) = pos(2,8) 
    cr4(3) = pos(3,8)
    print*, 'Cr4 ', cr4
    rtmp = cr4(:) -  pos(:,7)
    print*, 'Cr4 -Cr1'
    write(*,'(3f15.10)') rtmp
    write(*,'(3f15.10)') rtmp(1) * hd(:,1) + rtmp(2) * hd(:,2) + rtmp(3) * hd(:,3)
    diff_vec(:,3) = rtmp(1) * hd(:,1) + rtmp(2) * hd(:,2) + rtmp(3) * hd(:,3)
    write(*,'(f15.10)') modulus(rtmp(1) * hd(:,1) + rtmp(2) * hd(:,2) + rtmp(3) * hd(:,3),3)
    print*, ''

    hmat3 = 0.0d0
    hmat3(1,1) = vpdp_s / deltaE + et2g(1)
    hmat3(6,1) = (vddd + vddp)/2.0d0 
    hmat3(8,1) = - vpdp_s / deltaE -(vddd - vddp)/2.0d0
    hmat3(2,2) = 2* vpdp_s /deltaE + et2g(1)
    hmat3(7,2) = (vddd + 3 * vdds)/4.0d0
    hmat3(9,2) = 3*(vddd - vdds)/8.0d0 + dsqrt(3.0d0) *vpdp*vpds/2.0d0 / deltaE
    hmat3(10,2) = -dsqrt(3.0d0)*(vddd - vdds)/8.0d0 - vpdp * vpds/2.0d0/ deltaE
    hmat3(3,3) = vpdp_s / deltaE + et2g(1)
    hmat3(6,3) = hmat3(8,1)
    hmat3(8,3) = hmat3(6,1)
    hmat3(4,4) = 3.0d0 * vpds_s /4.0d0/ deltaE + eeg(1)
    hmat3(5,4) = dsqrt(3.0d0)* vpds_s/4.0d0/ deltaE
    hmat3(4,5) = hmat3(5,4)
    hmat3(7,4) = hmat3(9,2)
    hmat3(9,4) = 9.0d0 * vddd /16.0d0 + vddp/4.0d0 + 3.0d0 * vdds/ 16.0d0
    hmat3(10,4) = -3.0d0 *dsqrt(3.0d0) * vddd / 16.0 + dsqrt(3.0d0)* vddp/4.0d0 -dsqrt(3.0d0) * vdds /16.0d0
    hmat3(5,5) = 5.0d0 * vpds_s /4.0d0/ deltaE + eeg(1)
    hmat3(7,5) = hmat3(10,2)
    hmat3(9,5) = hmat3(10,4)
    hmat3(10,5) = 3.0d0 * vddd/16.0d0 + 3* vddp/ 4.0d0 + vdds/16.0d0
    
    hmat3(10,9) = dsqrt(3.0d0) * vpds_s / 4.0d0 / deltaE 
    hmat3(9,10) = hmat3(10,9)

    !do iloop = 1, 5
    !   hmat3(iloop+5,iloop+5) = hmat3(iloop,iloop)
    !end do
    hmat3(6,6) = vpdp_s / deltaE + et2g(2)
    hmat3(7,7) = 2.0d0 *vpdp_s / deltaE + et2g(2)
    hmat3(8,8) = hmat3(6,6)
    hmat3(9,9) =3.0d0* vpds_s / 4.0d0 / deltaE + eeg(2)
    hmat3(10,10) =5.0d0* vpds_s / 4.0d0 / deltaE + eeg(2)

    do iloop = 6,10
        do jloop =  1,5
            hmat3(jloop,iloop) = hmat3(iloop,jloop)   ! here is real number, do not use conjugate ........
        end do
    end do
    !
    print*, ' === Hamiltonian along third path Cr4-Cr1 =='
    do iloop = 1, 10
       write(*,'(10f8.4)') hmat3(:,iloop)
    end do
    !
    ! now we will write the parameters into explicit.fdf
    !
    file_out = 'explicit.fdf'
!%block Hop_Parameters
! 1 1 .. .. .. # is1, is2, distance vecotr(in bohr) H(i_alpha, j_beta)
!
!......
!  blank space
!1 1
!%endblock Hop_Parameters
   if( sym_atom(1) == sym_atom(2) )then
        is(:) = 1
   else
        is(1) = 1
        is(2) = 2
   end if
   diff_vec(:,4) = pos(:,7)-pos(:,7)
   rtmp(:) = pos(:,8)-pos(:,7)
   diff_vec(:,5) = rtmp(1) * hd(:,1) + rtmp(2) * hd(:,2) + rtmp(3) * hd(:,3)
   diff_vec = diff_vec/autoa ! bohr 
   open(id, file = file_out )
   Hcr1site(:,:) = 0.0
   Hcr2site(:,:) = 0.0

   do iloop = 1, 5
      Hcr1site(iloop, iloop ) = hmat1(iloop,iloop) + hmat2(iloop,iloop) &
                & + hmat3(iloop,iloop) 
      Hcr2site(iloop, iloop ) = hmat1(iloop+5,iloop+5) + hmat2(iloop+5,iloop+5)&
                & + hmat3(iloop+5,iloop+5)
   end do 

   write(id,*) '%block Hop_Parameters'
   !! onsite term
   !write(id,'(i2, i2, 3f12.8)') is(1), is(2), diff_vec(:,4)
   !do iloop = 1, 5
   !    write(id,'(5f12.8)') Hcr1site(1:5, iloop)
   !end do
   !write(id, *) ''

   !write(id,'(i2, i2, 3f12.8)') is, is, diff_vec(:,5)
   !do iloop = 1, 5
   !    write(id,'(5f12.8)') Hcr2site(1:5, iloop)
   !end do
   !write(id, *) ''

   write(id,'(i2, i2, 3f12.8)') is(1), is(2), diff_vec(:,1)
   do iloop = 1, 5
       write(id,'(5f12.8)') hmat1(6:10, iloop)
   end do
   write(id, *) ''
   write(id,'(i2, i2, 3f12.8)') is(1), is(2), diff_vec(:,2)
   do iloop = 1, 5
       write(id,'(5f12.8)') hmat2(6:10, iloop)
   end do
   write(id, *) ''
   write(id,'(i2, i2, 3f12.8)') is(1), is(2), diff_vec(:,3)
   do iloop = 1, 5
       write(id,'(5f12.8)') hmat3(6:10, iloop)
   end do
   write(id,*) '%endblock Hop_Parameters'
   if( sym_atom(1) == sym_atom(2) ) then
        norb = 5
        write(id,*) ''
        write(id,*) '%block Orbitals'
        write(id,'(2i4)') is(1), norb
        !write(id,'(5f6.2)') et2g, et2g,et2g, eeg, eeg
        tzero = 0.0d0   
        write(id,'(5f6.2)') tzero, tzero,tzero, tzero, tzero
        write(id,*) '%endblock Orbitals'
   else
        norb = 5
        write(id,*) ''
        write(id,*) '%block Orbitals'
        write(id,'(2i4)') is(1), norb
        !write(id,'(5f6.2)') et2g, et2g,et2g, eeg, eeg
        tzero = 0.0d0   
        write(id,'(5f6.2)') tzero, tzero,tzero, tzero, tzero
        write(id,'(2i4)') is(2), norb
        !write(id,'(5f6.2)') et2g, et2g,et2g, eeg, eeg
        tzero = 0.0d0   
        write(id,'(5f6.2)') tzero, tzero,tzero, tzero, tzero
        write(id,*) '%endblock Orbitals'
        
   end if
   !
   !
   ! lsimp = 0,1,2,3=== > 9, 7, 7, 6
   !                === > 7, 5, 5, 4
   !   
   if(sym_atom(1) == sym_atom(2) ) then
        write(id,*) '%block RPARA'
        write(id,'(i2)') npara
        if(npara == 7 ) then  
        !write(id,'(2i4)') is, norb
                write(id,'(7f6.2)') et2g(1), eeg(1), vpds, vpdp, vdds, vddp, vddd
        else if( npara  == 5) then

                write(id,'(5f6.2)') et2g(1), eeg(1), vpds, vpdp, vdds !, vddp, vddd
        else if( npara == 4 ) then
                   
                write(id,'(4f6.2)') et2g(1), eeg(1), vpds*vpdp, vdds !, vddp, vddd
        else
                stop 'error, actual parameter is not right!, please check!'
        end if
        write(id,*) '%endblock RPARA'
   else
        write(id,*) '%block RPARA'
        write(id,'(i2)') npara
          
        !write(id,'(2i4)') is, norb
        if( npara == 9 ) then        
                write(id,'(9f6.2)') et2g, eeg, vpds, vpdp, vdds, vddp, vddd
        else if ( npara == 7 ) then
                write(id,'(7f6.2)') et2g, eeg, vpds, vpdp, vdds !, vddp, vddd
        else if ( npara == 6 ) then
                write(id,'(6f6.2)') et2g, eeg, vpds * vpdp, vdds !, vddp, vddd
        else
                stop 'error, actual parameter is not right!, please check'

        end if

        write(id,*) '%endblock RPARA'
   end if
!
!%block Orbitals
! 1 5 #is,norb
! -3 -3 -3 -2 -2 #   t2g,t2g,t2g,eg,eg (on-site energy)
!%endblock Orbitals

   close(id)

   deallocate(pos,pos_car,indx_nest,et2g,eeg)
   deallocate(hmat1, hmat2, hmat3) 
   ! now we will postprocess the hunding coupling U and onsite energy
   ! approximately.

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

