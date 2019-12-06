
!: MARCH 19 2019 TUESDAY

!: to compile 
!: pgf95 -Minfo -Mneginfo -time -fast -Mconcur=allcores -mp=allcores -Munroll -Mvect ANALYSIS.F90
!: gfortran -fopenmp ANALYSIS.F90

!: contains
!: subroutine read_tdci   : reads time-dependent CI coefficients
!: subroutine read_ci_vec : reads CI expansion coefficients (eigenvectors)
!: subroutine read_mo     : reads in MO matrix elements
!: subroutine get_indices : computes hole-particle excitation pair
!: subroutine get_psi_det : converts wavefunction in CI basis to determinantal basis
!: subroutine get_expectation
!: subroutine get_mo_den  : computes MO density
!: subroutine get_hole_den( 'ground' or 't0' )
!: subroutine get_rate    : computes rate of norm decrease based on Vabs
!: subroutine get_ratios
!: subroutine get_overlap
!: subroutine name_files

  module global_variables 

    save    
    
    implicit none
    
    integer(8),    parameter :: start_dir = 1
    integer(8),    parameter :: ndir = 37
    character(5),  parameter :: e_d = '-e1-d'
   
    character(50),  parameter :: ofile = '../OUTPUT_FILES/OUTPUT'
    character(100), parameter :: f_restart = trim(ofile)//'_RESTART.bin'
    character(100), parameter :: f_mo_restart = trim(ofile)//'_RESTART_MO.bin'

    !: outfile
    integer(8), parameter :: iout = 78
    character(20), parameter :: log = 'analysis.out'

    !: only write last density
    logical, parameter :: mo_write_all = .False.

    
    !: assigned in name_files
    character(100) :: input 
    !: outfiles
    character(100) :: f_mo_out  
    character(100) :: f_pop_out 
    character(100) :: f_rate_out 
    character(100) :: f_hole_out 
    character(100) :: f_hole_pop_out
    character(100) :: f_part_pop_out
    character(100) :: f_ratios_out
    character(100) :: f_part_out 
    character(100) :: f_overlap_out
    character(100) :: f_cumediff_out
    
    real(8), parameter :: au2fs = 2.418884326509d-2

    !: define me before moving on
    integer(8) :: ntimes, nstuse, nstates
    integer(8) :: noa, nva, nob, nvb
    integer(8) :: nrorb, norb
    logical :: unrestricted 

    real(8), allocatable :: ci_vec(:,:)
    real(8), allocatable :: time(:), norm_sq(:)
    real(8), allocatable :: efield(:), dirx(:), diry(:), dirz(:) 
    complex(8), allocatable :: psi(:,:)

    !: state indices
    integer(8), allocatable :: hole1(:), part1(:)
    integer(8), allocatable :: hole(:,:), part(:,:)
      
    !: read in matrix elements in MO
    real(8), allocatable :: dipx_a(:,:), dipy_a(:,:), dipz_a(:,:), vabs_a(:,:)
    real(8), allocatable :: dipx_b(:,:), dipy_b(:,:), dipz_b(:,:), vabs_b(:,:)
    
    
  end module global_variables
  !:--------------------------:!
  !:--------------------------:!
  program analysis
    
    use global_variables
    implicit none
    
    integer(8) :: idir


    open( unit=iout, file=trim(log) )


    call read_ci_vec  !: read in CI eigenvectors
    call read_mo      !: read in MO elements
    call get_indices  !: get indices       


    do idir=start_dir, ndir
       
       write(iout,"(' *******************')" )  
       write(iout,"(' working on direction',i0)" )  idir
       flush(iout)
       
       call name_files( idir )
       call read_tdci              !: read in time-dependent CI coefficients       
       !call get_rate              !: get rate 
       !call get_mo_den            !: get MO density and print out for orbkit analysis
       call get_hole_den('t0')     !: get hole density and print out for orbkit analysis
       !call get_ratios
       !call get_overlap
       !call get_cumulative_diff


    end do

    close(iout)


    
  end program analysis
  !:--------------------------:!
  !:--------------------------:!
  subroutine read_tdci

    
    use global_variables 
    implicit none
    
    integer(8) :: itime, i, j
    real(8), allocatable :: readtmp(:), creadtmp(:)
    

    !: read input
    open( unit=100,file=trim(input), form='unformatted' )
    read(100) ntimes, nstuse, nstates    

    if ( .not.allocated( time) ) then
       allocate( time(ntimes), norm_sq(ntimes) )
       allocate( efield(ntimes), dirx(ntimes), diry(ntimes), dirz(ntimes) )
       allocate( psi(nstuse,ntimes) )
    else
       time = 0.d0 ; norm_sq = 0.d0
       efield = 0.d0 ; dirx = 0.d0 ; diry = 0.d0 ; dirz = 0.d0
       psi = dcmplx( 0.d0,0.d0 )
    end if

    allocate( readtmp(nstuse), creadtmp(nstuse) )    

    do itime=1, ntimes
       
       read(100) time(itime), efield(itime), dirx(itime), diry(itime), dirz(itime), norm_sq(itime)  
       read(100) (  readtmp(i), i=1, nstuse )
       read(100) ( creadtmp(i), i=1, nstuse )
       psi(:,itime) = dcmplx( readtmp, creadtmp )
       
    end do
    
    deallocate( readtmp, creadtmp )
    
    write(iout,*) ' --> allocated time, norm_sq, efield, dirx, diry, dirz '
    write(iout,*) ' --> allocated psi(nstuse,ntimes) '
    write(iout,*) ' --> finished reading psi(t)'
    flush(iout)
    

  end subroutine read_tdci
  !:--------------------------:!
  !:--------------------------:!
  subroutine read_ci_vec

    
    use global_variables
    implicit none

    integer(8) :: i, j

    
    open( unit=100, file=trim(f_restart), form='unformatted' )
    read(100) nstates, nstuse
    
    allocate( ci_vec(nstates,nstates) )
    read(100) ( ( ci_vec(j,i), j=1, nstates ), i=1, nstates )

    close(100) 
    
    write(iout,*) ' --> allocated ci_vec '
    write(iout,*) ' --> finished reading ci_vec '
    flush(iout)

    
  end subroutine read_ci_vec
  !:--------------------------:!
  !:--------------------------:!
  subroutine read_mo

    use global_variables
    implicit none

    integer(8) :: i
    integer(8) :: na, nb


    open( unit=100, file=trim(f_mo_restart), form='unformatted' )
    read(100) noa, nva, nob, nvb
    if( nob.ne.0 ) then
       na = noa + nva
       nb = nob + nvb
       nrorb = na ; norb = nrorb + nb
       unrestricted = .True.
       allocate( dipx_a(na,na), dipy_a(na,na), dipz_a(na,na), vabs_a(na,na) )
       allocate( dipx_b(nb,nb), dipy_b(nb,nb), dipz_b(nb,nb), vabs_b(nb,nb) )
    else
       na = noa + nva
       nrorb = na ; norb = nrorb
       unrestricted = .False.
       allocate( dipx_a(na,na), dipy_a(na,na), dipz_a(na,na), vabs_a(na,na) ) 
    end if
    
    read(100) ( vabs_a(:,i), i=1, na )
    read(100) ( dipx_a(:,i), i=1, na )
    read(100) ( dipy_a(:,i), i=1, na )
    read(100) ( dipz_a(:,i), i=1, na )
    if ( unrestricted ) then
       read(100) ( vabs_b(:,i), i=1, nb )
       read(100) ( dipx_b(:,i), i=1, nb )
       read(100) ( dipy_b(:,i), i=1, nb )
       read(100) ( dipz_b(:,i), i=1, nb )
    end if

    write(iout,*) ' --> allocated dipx, dipy, dipz, vabs '
    write(iout,*) ' --> finished reading MO elements'
    flush(iout)

    
  end subroutine read_mo
  !:--------------------------:!
  !:--------------------------:!
  subroutine get_indices
    

    use global_variables
    implicit none

    integer(8) :: i, a, ia 


    allocate( hole1(nstates), part1(nstates) )
    hole1 = 0 ; part1 = 0
       
    !: indices for CIS for now
    ia = 1
    do i=1, noa
       do a=1, nva
          ia = ia + 1 
          hole1(ia) = -i 
          part1(ia) = -a - noa
       end do
    end do
       
    if ( unrestricted ) then
       do i=1, nob
          do a=1, nvb
             ia = ia + 1 
             hole1(ia) = i 
             part1(ia) = a + nob
          end do
       end do
    end if
    
    write(iout,*) ' --> allocated hole, part '
    write(iout,*) ' --> finished assigning hole and part indices'
    flush(iout)
    

  end subroutine get_indices
  !:--------------------------:!
  !:--------------------------:!
  subroutine get_psi_det( psi_itime, psi_det )


    use global_variables
    implicit none
    
    complex(8), intent(in) :: psi_itime(nstuse)
    complex(8), intent(inout) :: psi_det(nstates)

    integer(8) :: i
    complex(8) :: psi_i


    psi_det = dcmplx(0.d0,0.d0) 
    do i=1, nstuse
       psi_i = psi_itime(i)
       psi_det(:) = psi_det(:) + psi_i * ci_vec(:,i)
    end do
    

  end subroutine get_psi_det
  !:--------------------------:!
  !:--------------------------:!
  subroutine get_expectation( rho, obsv, expect_value ) 

    
    use global_variables
    implicit none

    real(8), intent(in) :: rho(nrorb,nrorb), obsv(nrorb,nrorb)
    real(8), intent(inout) :: expect_value
    
    integer(8) :: i, j
    real(8) :: rdum


    !: note expect_value not initialized to zero
    !: obsv(i,j) = obsv(j,i)
    do i=1, nrorb
       rdum = dot_product( rho(:,i), obsv(:,i) )
       expect_value = expect_value + rdum
    end do

    
  end subroutine get_expectation
  !:--------------------------:!
  !:--------------------------:!
  subroutine get_mo_den


    use global_variables 
    use omp_lib
    implicit none
    
    !: for writing
    real(8), parameter :: tol = 1.d-6 
    real(8), parameter :: jb_tol = 1.d-10
    
    integer(8) :: itime, dwrite, funit
    integer(8) :: ia, jb, i, a , j, b, ii, aa, jj, bb
    real(8)    :: const, i_norm2
    real(8)    :: rho(norb,norb,ntimes), rho2(norb,norb)
    complex(8) :: psi0, psi_det(nstates)

    character(20) :: corb


    const = dsqrt(2.d0)
    if ( unrestricted ) const = 1.d0 
    

    !: initialize rho
    rho = 0.d0  
    

    !$OMP PARALLEL DEFAULT(NONE), &
    !$OMP PRIVATE( itime, i, ii, a, aa, j, jj, b, bb, ia, jb, i_norm2, rho2, psi_det, psi0 ), &
    !$OMP SHARED( noa, nob, nrorb, nstates, ntimes, const, psi, norm_sq, hole1, part1, rho )
    
    !$OMP DO
    do itime=1, ntimes  

       !: get psi_det and normalize
       call get_psi_det( psi(:,itime), psi_det )

       i_norm2 = norm_sq(itime)
       psi_det = psi_det / dsqrt(i_norm2)
       rho2    = 0.d0

       !: diag occupied
       do i=1, noa
          rho2(i,i) = 1.d0
       end do
       do i=1, nob
          ii = nrorb + i
          rho2(ii,ii) = 1.d0
       end do

       !: < psi_0 | psi_ia >_(N-1)
       psi0 = dconjg( psi_det(1) )
       do ia = 2, nstates
          ii = hole1(ia) ; i = -ii
          aa = part1(ia) ; a = -aa
          if ( ii.gt.0 ) then
             i = nrorb + ii 
             a = nrorb + aa 
          end if
          rho2(a,i) = const * real( psi0*psi_det(ia) )
          rho2(i,a) = const * real( psi0*psi_det(ia) )
       end do
       
       !: the rest
       do ia = 2, nstates                    
          ii = hole1(ia) ; i = -ii
          aa = part1(ia) ; a = -aa
          if ( ii.gt.0 ) then
             i = nrorb + ii
             a = nrorb + aa 
          end if

          psi0 = dconjg( psi_det(ia) )   
          
          if ( abs(psi0).gt.jb_tol ) then
             do jb = 2, nstates 
                jj = hole1(jb) ; j = -jj
                bb = part1(jb) ; b = -bb
                if ( jj.gt.0 ) then
                   j = nrorb + jj
                   b = nrorb + bb 
                end if
                if ( a.eq.b ) rho2(j,i) = rho2(j,i) - real( psi0*psi_det(jb) )
                if ( i.eq.j ) rho2(b,a) = rho2(b,a) + real( psi0*psi_det(jb) )             
             end do
          end if

       end do

       !: record
       rho(:,:,itime) = rho2(:,:)
       
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    
    !: write for movie-making
    dwrite = ntimes-1
    if ( mo_write_all ) dwrite = 1
    
    open( unit=100, file=trim(f_mo_out)  )    
    write(100,"('ntimes ',i0)") ntimes
    write(100,"('alpha_homo ',i0)") noa
    write(100,"('beta_homo ',i0)")  nob
    write(100,"('alpha_orbitals ',i0)") noa + nva      
    do itime=1, ntimes, dwrite
       rho2(:,:) = rho(:,:,itime)
       write(100,"(f13.10)") time(itime)
       do i=1, norb
          if ( abs(rho2(i,i)).gt.tol ) write(100,"(i5,i5,1x,f13.10,',')",advance='no') i, i, rho2(i,i)
          do j=(i+1), norb
             if ( abs(rho2(j,i)).gt.tol ) write(100,"(i5,i5,1x,f13.10,',')",advance='no') i, j, 2.d0*rho2(j,i)
          end do
       end do
       write(100,'(A)') ''
    end do
    close(100)
    
    
    !: write population analysis
    funit = 100
    open( unit=funit, file=trim(f_pop_out) )
    write( funit, '(A)', advance='no' ) ' time(fs)'
    !: alpha
    do i=(noa-1), 1, -1
       write( corb, '(i0)' ) i 
       write( funit, '(A)', advance='no' ) ' aHOMO-'//trim(corb)
    end do
    write( funit,'(A)', advance='no' ) ' aHOMO aLUMO'
    do i=1, nva-1
       write( corb, '(i0)' ) i
       write( funit, '(A)', advance='no' ) ' aLUMO+'//trim(corb)
    end do
    !: beta
    do i=(nob-1), 1, -1 
       write( corb, '(i0)' ) i 
       write( funit, '(A)', advance='no' ) ' bHOMO-'//trim(corb)
    end do
    write( funit, '(A)', advance='no' ) ' bHOMO bLUMO' 
    do i=1, nvb-1
       write( corb, '(i0)' ) i
       write( funit, '(A)', advance='no' ) ' bLUMO+'//trim(corb)
    end do
    write( funit, '(A)' ) ''
    do itime=1, ntimes
       write( funit, "(2000(f13.10,1x))" ) time(itime), ( rho(i,i,itime), i=1, norb )
    end do
    close(funit)

    
    !: pretty write
    !do itime=1, ntimes
    !   call write_rhos( noa, nva, nob, nvb, norb, time(itime), rho(:,:,itime), funit  )
    !end do
    !close(funit)
    
    
    write(iout,*) ' --> finished computing mo_density'
    flush(iout)
    
    
  end subroutine get_mo_den
  !:--------------------------:!
  !:--------------------------:!
  subroutine get_hole_den( option )


    use global_variables
    !use omp_lib
    implicit none
    
    character(*), intent(in) :: option

    real(8), parameter  :: tol = 1.d-6
    
    integer(8) :: funit
    integer(8) :: itime, ia, jb, i, a , j, b, ii, aa, jj, bb
    integer(8) :: ipart, jpart, ihole, jhole
    real(8) :: rho_eh(norb,norb)
    real(8) :: const
    
    complex(8) :: psi0, psi_ia        
    complex(8) :: c_sum, ipart_ihole
    complex(8) :: rho2(norb, norb), psi_det(nstates), psi_det00(nstates), rho(norb,norb,ntimes)

    character(10) :: corb
    

    !: LEFT INDEX  = PARTICLE INDEX
    !: RIGHT INDEX = HOLE INDEX
    
    
    const = sqrt(2.d0)
    if ( unrestricted ) const = 1.d0
    

    if ( trim(option).eq.'ground' ) then
       go to 400
    else if ( trim(option).eq.'t0' ) then
       go to 500
    end if
    
    
400 continue

    !: ground state psi
    psi0 = dcmplx( 1.d0, 0.d0 )
    
    !: get hole and particle densities with respect to ground state for time 0
    !$OMP PARALLEL DEFAULT(NONE), &
    !$OMP PRIVATE( itime, i, ii, a, aa, ia, psi_ia, psi_det, rho2 ), &
    !$OMP SHARED( noa, nob, nrorb, nstates, ntimes, const, norm_sq, hole1, part1, psi, psi0, rho )
    !$OMP DO
    do itime=1, ntimes
       
       call get_psi_det( psi(:,itime), psi_det )
       psi_det = psi_det / dsqrt( norm_sq(itime) )
       
       rho2 = dcmplx( 0.d0, 0.d0 )
       do i=1, noa
          rho2(i,i) = dconjg(psi_det(1)) * psi0
       end do
       do i=1, nob
          rho2(nrorb+i, nrorb+i) = dconjg(psi_det(1)) * psi0
       end do
       do ia=2, nstates
          if ( abs(psi_det(ia)) .gt.tol ) then
             !: get indices
             ii = hole1(ia) ; i = -ii
             aa = part1(ia) ; a = -aa
             if ( ii.gt.0 ) then
                i = nrorb + ii
                a = nrorb + aa
             end if
             rho2(a,i) = const * dconjg( psi_det(ia) ) * psi0
          end if
       end do
       rho(:,:,itime) = rho2(:,:)

    end do
    !$OMP END DO
    !$OMP END PARALLEL
    
    go to 600

    
500 continue
        

    !$OMP PARALLEL DEFAULT(NONE), &
    !$OMP PRIVATE( itime, i, ii, a, aa, j, jj, b, bb, ia, jb, psi_ia, psi_det, rho2, c_sum, psi_det00 ), &
    !$OMP SHARED( noa, nob, nrorb, nstates, ntimes, const, norm_sq, hole1, part1, psi, rho )        
    !$OMP DO
    do itime=2, ntimes
       
       call get_psi_det( psi(:,1), psi_det00 ) 
       psi_det00 = psi_det00 / dsqrt( norm_sq(1) )
       
       call get_psi_det( psi(:,itime), psi_det )
       psi_det = psi_det / dsqrt( norm_sq(itime) )
       
       
       !: initialize
       rho2 = dcmplx( 0.d0, 0.d0 )       
       
       do ia=2, nstates 
          
          psi_ia = psi_det(ia)
          if( abs(psi_ia) .gt. 1.d-12 ) then          
             !: get indices
             ii = hole1(ia) ; i = -ii
             aa = part1(ia) ; a = -aa
             if ( ii.gt.0 ) then
                i = nrorb + ii
                a = nrorb + aa
             end if
             
             rho2(a,i) = const * psi_ia * psi_det00(1)
             rho2(i,a) = const * dconjg( psi_det(1) ) * psi_det00(ia)          
             
             do jb=2, nstates                
                if ( jb.ne.ia ) then
                   !: get indices
                   jj = hole1(jb) ; j = -jj
                   bb = part1(jb) ; b = -bb
                   if ( jj.gt.0 ) then
                      j = nrorb + jj
                      b = nrorb + bb
                   end if
                   if ( ii.eq.jj ) rho2(a,b) = rho2(a,b) + psi_ia * psi_det00(jb)   
                   if ( bb.eq.aa ) rho2(j,i) = rho2(j,i) - psi_ia * psi_det00(jb)
                end if
             end do

          end if
       end do
       
       rho(:,:,itime) = rho2(:,:)
       
    end do
    
    !$OMP END DO
    !$OMP END PARALLEL
    

600 continue

    !: hole
    funit = 100
    open( unit=funit, file=trim(f_hole_pop_out) )
    write( funit, '(A)', advance='no' ) ' time(fs)'
    do i=(noa-1), 1, -1
       write( corb, '(i0)' ) i 
       write( funit, '(A)', advance='no' ) ' aHOMO-'//trim(corb)
    end do
    write( funit,'(A)', advance='no' ) ' aHOMO aLUMO'
    do i=1, nva-1
       write( corb, '(i0)' ) i
       write( funit, '(A)', advance='no' ) ' aLUMO+'//trim(corb)
    end do
    do i=(nob-1), 1, -1 
       write( corb, '(i0)' ) i 
       write( funit, '(A)', advance='no' ) ' bHOMO-'//trim(corb)
    end do
    write( funit, '(A)', advance='no' ) ' bHOMO bLUMO' 
    do i=1, nvb-1
       write( corb, '(i0)' ) i
       write( funit, '(A)', advance='no' ) ' bLUMO+'//trim(corb)
    end do
    write( funit, '(A)' ) ''
    do itime=1, ntimes       
       rho2(:,:) = rho(:,:,itime)
       !$OMP PARALLEL DEFAULT(NONE),PRIVATE( ihole, jhole ),SHARED( rho_eh, rho2, norb )
       !$OMP DO
       do ihole=1, norb
          do jhole=1, norb
             rho_eh(jhole,ihole) = dot_product( rho2(:,ihole), rho2(:,jhole) )
          end do
       end do
       !$OMP END DO
       !$OMP END PARALLEL
       write( funit, "(2000(f13.10,1x))" ) time(itime), ( rho_eh(i,i), i=1, norb )
       !call write_rhos( noa, nva, nob, nvb, norb, time(itime), rho_eh, funit )
    end do
    close(funit)
    
    !: particle
    open( unit=funit, file=trim(f_part_pop_out) )
    write( funit, '(A)', advance='no' ) ' time(fs)'
    do i=(noa-1), 1, -1
       write( corb, '(i0)' ) i 
       write( funit, '(A)', advance='no' ) ' aHOMO-'//trim(corb)
    end do
    write( funit,'(A)', advance='no' ) ' aHOMO aLUMO'
    do i=1, nva-1
       write( corb, '(i0)' ) i
       write( funit, '(A)', advance='no' ) ' aLUMO+'//trim(corb)
    end do
    do i=(nob-1), 1, -1 
       write( corb, '(i0)' ) i 
       write( funit, '(A)', advance='no' ) ' bHOMO-'//trim(corb)
    end do
    write( funit, '(A)', advance='no' ) ' bHOMO bLUMO' 
    do i=1, nvb-1
       write( corb, '(i0)' ) i
       write( funit, '(A)', advance='no' ) ' bLUMO+'//trim(corb)
    end do
    write( funit, '(A)' ) ''
    do itime=1, ntimes
       rho2 = rho(:,:,itime)
       !$OMP PARALLEL DEFAULT(NONE),PRIVATE( ipart, jpart ),SHARED( rho_eh, rho2, norb )
       !$OMP DO
       do ipart=1, norb
          do jpart=1, norb
             rho_eh(jpart,ipart) = dot_product( rho2(ipart,:), rho2(jpart,:) )
          end do
       end do
       !$OMP END DO
       !$OMP END PARALLEL
       write( funit, "(2000(f13.10,1x))" ) time(itime), ( rho_eh(i,i), i=1, norb )
       !call write_rhos( noa, nva, nob, nvb, norb, time(itime), rho_eh, funit )
    end do
    close(funit)
    
    
    write(iout,*) ' --> finished computing hole_density'
    flush(iout)
    
    
  end subroutine get_hole_den
  !:--------------------------:!
  !:--------------------------:!
  subroutine get_rate_diagonal


    use global_variables
    use omp_lib
    implicit none

    
    integer(8) :: itime, ii, aa, jj, bb, i, a, j, b, ia, jb
    real(8)    :: dt, const, total_loss
    real(8)    :: rate(norb,ntimes), total_rate(ntimes), cume_loss(norb,ntimes) 
    real(8)    :: rate2(norb)
    complex(8) :: psi0, psi_det(nstates)
    
    !: write
    character(5)   :: c0
    character(1000) :: cformat, cformat2


    !: for writing
    if ( .not.unrestricted ) then
       cformat  = '     loss_HOMO'
       cformat2 = '     frac_HOMO'
       do i=1,noa-1
          write( c0, '(i0)') i
          cformat  = trim(cformat)//'    loss_HOMO-'//trim(c0)
          cformat2 = trim(cformat2)//'    frac_HOMO-'//trim(c0)
       end do
    else
       cformat  = '    loss_aHOMO'
       cformat2 = '    frac_aHOMO'
       do i=1,noa-1
          write( c0, '(i0)') i
          cformat  = trim(cformat)//'  loss_aHOMO-'//trim(c0)
          cformat2 = trim(cformat2)//'  frac_aHOMO-'//trim(c0)
       end do
       cformat  = trim(cformat)//'    loss_bHOMO'
       cformat2 = trim(cformat2)//'    frac_bHOMO'
       do i=1,nob-1
          write( c0, '(i0)') i
          cformat  = trim(cformat)//'  loss_bHOMO-'//trim(c0)
          cformat2 = trim(cformat2)//'  frac_bHOMO-'//trim(c0)
       end do
    end if
    cformat2 = '#     time'//'    total_loss'//'     pred_loss'//trim(cformat2)
    cformat  = '#     time'//'         norm2'//'    rate(1/fs)'//trim(cformat)//'    pred_norm2'

    
    !: initialize
    rate = 0.d0 ; cume_loss = 0.d0 ; total_rate = 0.d0
    dt   = time(2) - time(1)
    
    const = dsqrt(2.d0)
    if ( unrestricted ) const = 1.d0
    
    
    !$OMP PARALLEL DEFAULT(NONE), &
    !$OMP PRIVATE( itime, i, ii, a, aa, j, jj, b, bb, ia, jb, &
    !$OMP          psi0, rate2, psi_det ), &
    !$OMP SHARED( ntimes, nrorb, nstates, hole1, part1, vabs_a, vabs_b, psi, rate )

    !$OMP DO
    do itime=1, ntimes
       
       call get_psi_det( psi(:,itime), psi_det )
       
       !: initialize temporary array
       rate2 = 0.d0
       
       !: diagonal
       do ia=2, nstates
          
          ii = hole1(ia) 
          aa = part1(ia)
          psi0 = dconjg( psi_det(ia) )
          
          if ( ii.lt.0 ) then
             i = -ii
             a = -aa
             rate2(i) = rate2(i) + real(dconjg(psi0)*psi0) * vabs_a(a,a) 
          else if ( ii.gt.0 ) then
             i = ii
             a = aa
             rate2(i+nrorb) = rate2(i+nrorb) + real(dconjg(psi0)*psi0) * vabs_b(a,a)
          end if
          
          !: off diagonal
          do jb=(ia+1), nstates
             
             jj = hole1(jb)
             bb = part1(jb)             
             if ( jj.eq.ii ) then                
                if ( bb.lt.0 ) then
                   b = -bb
                   rate2(i) = rate2(i) + 2.d0 * real(psi0*psi_det(jb)) * vabs_a(b,a)
                else
                   b = bb
                   rate2(i+nrorb) = rate2(i+nrorb) + 2.d0 * real(psi0*psi_det(jb)) * vabs_b(b,a)
                end if
             end if
             
          end do
       end do

       rate(:,itime) = rate2(:)       

    end do
    !$OMP END DO
    !$OMP END PARALLEL
    

    !: get cume loss
    total_rate(1)  = -2.d0 * sum( rate(:,1) ) / au2fs !: dN/dt = -2 * SUM V(i)
    cume_loss(:,1) = 0.d0
    do itime=2, ntimes
       cume_loss(:,itime) = cume_loss(:,itime-1) - 2.d0 * dt * rate(:,itime-1) / au2fs !* norm_sq(itime)
       total_rate(itime) = -2.d0 * sum( rate(:,itime) ) / au2fs
    end do
    
    !: write
    open( unit=100,file=trim(f_rate_out) )
    write(100,"('#noa = ',i0)") noa
    write(100,"('#nob = ',i0)") nob
    write(100, '(A)') trim(cformat)
    do itime=1, ntimes
       if ( unrestricted ) then
          write(100,"(f10.5,1x,2(f13.10,1x),100(f13.10,1x))") &
               time(itime),   &
               norm_sq(itime),   &
               total_rate(itime), &
               (cume_loss(i,itime), i=noa,1,-1), &
               (cume_loss(i,itime), i=nrorb+nob,nrorb+1,-1 ), 1.d0+sum( cume_loss(:,itime) )  
       else
          write(100,"(f10.5,1x,2(f13.10,1x),100(f13.10,1x))") &
               time(itime),   &
               norm_sq(itime),   &
               total_rate(itime),&
               (cume_loss(i,itime), i=noa,1,-1), 1.d0+sum( cume_loss(:,itime) )
       end if
    end do
    
    total_loss = sum( cume_loss(:,ntimes) )

    write(100, '(A)') trim(cformat2)
    if ( unrestricted ) then
       write(100,"('#',f9.5,1x,2(f13.10,1x),100(f13.10,1x))") &
            time(ntimes),   &
            -( 1.d0 - norm_sq(ntimes) ),&
            total_loss,    &
            (cume_loss(i,ntimes)/total_loss, i=noa,1,-1), &
            (cume_loss(i,ntimes)/total_loss, i=nrorb+nob,nrorb+1,-1 )
    else
       write(100,"('#',f9.5,1x,2(f13.10,1x),100(f13.10,1x))") &
            time(ntimes),   &
            -( 1.d0 - norm_sq(ntimes) ),&
            total_loss,    &
            (cume_loss(i,ntimes)/total_loss, i=noa,1,-1)
    end if
    close(100)

    write(iout,*) ' --> finished computing rate'
    flush(iout)
    

  end subroutine get_rate_diagonal
  !:--------------------------:!
  !:--------------------------:!
  subroutine get_rate


    use global_variables
    use omp_lib
    implicit none

    
    integer(8) :: itime, ii, aa, jj, bb, i, a, j, b, ia, jb
    real(8)    :: dt, const, total_loss
    real(8)    :: rate(norb,norb,ntimes), total_rate(ntimes), cume_loss(norb,norb,ntimes) 
    real(8)    :: rate2(norb,norb)
    complex(8) :: psi0, psi_det(nstates)
    
    !: write
    character(5)   :: c0
    character(1000) :: cformat, cformat2


    !: for writing
    if ( .not.unrestricted ) then
       cformat  = 'HOMO'
       cformat2 = 'HOMO'
       do i=1,noa-1
          write( c0, '(i0)') i
          cformat  = trim(cformat)//'loss_HOMO-'//trim(c0)
          cformat2 = trim(cformat2)//'HOMO-'//trim(c0)
       end do
    else
       cformat  = '    loss_aHOMO'
       cformat2 = '    frac_aHOMO'
       do i=1,noa-1
          write( c0, '(i0)') i
          cformat  = trim(cformat)//'  loss_aHOMO-'//trim(c0)
          cformat2 = trim(cformat2)//'  frac_aHOMO-'//trim(c0)
       end do
       cformat  = trim(cformat)//'    loss_bHOMO'
       cformat2 = trim(cformat2)//'    frac_bHOMO'
       do i=1,nob-1
          write( c0, '(i0)') i
          cformat  = trim(cformat)//'  loss_bHOMO-'//trim(c0)
          cformat2 = trim(cformat2)//'  frac_bHOMO-'//trim(c0)
       end do
    end if
    cformat2 = '#     time'//'    total_loss'//'     pred_loss'//trim(cformat2)
    cformat  = '#     time'//'         norm2'//'    rate(1/fs)'//trim(cformat)//'    pred_norm2'

    
    !: initialize
    rate = 0.d0 ; cume_loss = 0.d0 ; total_rate = 0.d0
    dt   = time(2) - time(1)
    
    const = dsqrt(2.d0)
    if ( unrestricted ) const = 1.d0
    
    
    !$OMP PARALLEL DEFAULT(NONE), &
    !$OMP PRIVATE( itime, i, ii, a, aa, j, jj, b, bb, ia, jb, &
    !$OMP          psi0, rate2, psi_det ), &
    !$OMP SHARED( ntimes, nrorb, nstates, hole1, part1, vabs_a, vabs_b, psi, rate )

    !$OMP DO
    do itime=1, ntimes
       
       call get_psi_det( psi(:,itime), psi_det )
       
       !: initialize temporary array
       rate2 = 0.d0
       
       !: diagonal
       do ia=2, nstates
          
          ii = hole1(ia) 
          aa = part1(ia)
          psi0 = dconjg( psi_det(ia) )
          
          if ( ii.lt.0 ) then
             i = -ii
             a = -aa
             rate2(i) = rate2(i) + real(dconjg(psi0)*psi0) * vabs_a(a,a) 
          else if ( ii.gt.0 ) then
             i = ii
             a = aa
             rate2(i+nrorb) = rate2(i+nrorb) + real(dconjg(psi0)*psi0) * vabs_b(a,a)
          end if
          
          !: off diagonal
          do jb=(ia+1), nstates
             
             jj = hole1(jb)
             bb = part1(jb)             
             if ( jj.eq.ii ) then                
                if ( bb.lt.0 ) then
                   b = -bb
                   rate2(i) = rate2(i) + 2.d0 * real(psi0*psi_det(jb)) * vabs_a(b,a)
                else
                   b = bb
                   rate2(i+nrorb) = rate2(i+nrorb) + 2.d0 * real(psi0*psi_det(jb)) * vabs_b(b,a)
                end if
             end if
             
          end do
       end do

       rate(:,itime) = rate2(:)       

    end do
    !$OMP END DO
    !$OMP END PARALLEL
    

    !: get cume loss
    total_rate(1)  = -2.d0 * sum( rate(:,1) ) / au2fs !: dN/dt = -2 * SUM V(i)
    cume_loss(:,1) = 0.d0
    do itime=2, ntimes
       cume_loss(:,itime) = cume_loss(:,itime-1) - 2.d0 * dt * rate(:,itime-1) / au2fs !* norm_sq(itime)
       total_rate(itime) = -2.d0 * sum( rate(:,itime) ) / au2fs
    end do
    
    !: write
    open( unit=100,file=trim(f_rate_out) )
    write(100,"('#noa = ',i0)") noa
    write(100,"('#nob = ',i0)") nob
    write(100, '(A)') trim(cformat)
    do itime=1, ntimes
       if ( unrestricted ) then
          write(100,"(f10.5,1x,2(f13.10,1x),100(f13.10,1x))") &
               time(itime),   &
               norm_sq(itime),   &
               total_rate(itime), &
               (cume_loss(i,itime), i=noa,1,-1), &
               (cume_loss(i,itime), i=nrorb+nob,nrorb+1,-1 ), 1.d0+sum( cume_loss(:,itime) )  
       else
          write(100,"(f10.5,1x,2(f13.10,1x),100(f13.10,1x))") &
               time(itime),   &
               norm_sq(itime),   &
               total_rate(itime),&
               (cume_loss(i,itime), i=noa,1,-1), 1.d0+sum( cume_loss(:,itime) )
       end if
    end do
    
    total_loss = sum( cume_loss(:,ntimes) )

    write(100, '(A)') trim(cformat2)
    if ( unrestricted ) then
       write(100,"('#',f9.5,1x,2(f13.10,1x),100(f13.10,1x))") &
            time(ntimes),   &
            -( 1.d0 - norm_sq(ntimes) ),&
            total_loss,    &
            (cume_loss(i,ntimes)/total_loss, i=noa,1,-1), &
            (cume_loss(i,ntimes)/total_loss, i=nrorb+nob,nrorb+1,-1 )
    else
       write(100,"('#',f9.5,1x,2(f13.10,1x),100(f13.10,1x))") &
            time(ntimes),   &
            -( 1.d0 - norm_sq(ntimes) ),&
            total_loss,    &
            (cume_loss(i,ntimes)/total_loss, i=noa,1,-1)
    end if
    close(100)

    write(iout,*) ' --> finished computing rate'
    flush(iout)
    

  end subroutine get_rate
  !:--------------------------:!
  !:--------------------------:!
  subroutine get_ratios 

    use global_variables
    implicit none

    integer(8), parameter :: tstates = 10

    integer(8) :: itime, i
    real(8)    :: ground, norm2, ratio(tstates)
    complex(8) :: psi_itime(nstuse)


    open( unit=100,file=trim(f_ratios_out) )
    write( 100, '(a10,12(1x,a12))' ) 'time(fs)', 'ground/norm2', &
         'excited1/norm2', 'excited2/norm2', &
         'excited3/norm2', 'excited4/norm2', 'excited5/norm2', 'excited6/norm2', &
         'excited7/norm2', 'excited8/norm2', 'excited9/norm2', 'excited10/norm2', &
         'excited/norm2'
    
    do itime=1, ntimes
       
       psi_itime(:) = psi(:,itime)
       ground = dconjg( psi_itime(1) ) * psi_itime(1)
       norm2  = norm_sq(itime)
       !: get first ten excited to ground, last one total, excited to ground
       do i=1, tstates
          ratio(i) = dconjg(psi_itime(i+1)) * psi_itime(i+1) / norm2
       end do
       
       write(100,"( f10.5,12(1x,f12.10) )") time(itime), &
            ground/norm2, ( ratio(i), i=1, tstates ), 1.d0-ground/norm2
       
    end do

    close(100)

    write(iout,*) ' --> finished computing the ratios'
    flush(iout)
    
    

  end subroutine get_ratios
  !:--------------------------:!
  !:--------------------------:!
  subroutine get_overlap

    use global_variables
    implicit none

    integer(8) :: itime

    real(8) :: overlap(ntimes) 
    complex(8) :: overlap2
    complex(8) :: psi1(nstuse), psi2(nstuse)

    
    !: initialize at time 1
    psi1(:) = psi(:,1) / dsqrt(norm_sq(1))
    overlap(1) = 1.d0
    
    do itime=2, ntimes
       
       psi2(:) = psi(:,itime) / dsqrt(norm_sq(itime))
       overlap2 = dot_product( psi2, psi1 )
       overlap(itime) = dconjg(overlap2) * overlap2 
       
       psi1 = psi2 

    end do

    open( unit=100,file=trim(f_overlap_out) )
    write(100,'(A)') 'time(fs)', 'overlap'
    do itime=1, ntimes 
       write(100,"(f15.10,f20.10)") time(itime), overlap(itime)
    end do
    close(100)


  end subroutine get_overlap
  !:--------------------------:!
  !:--------------------------:!
  subroutine get_cumulative_diff

    use global_variables
    implicit none

    integer(8) :: itime, i
    real(8) :: icume(nstates), record
    real(8) :: cume_diff(nstates,ntimes)
    
    complex(8) :: psi_det1(nstates), psi_det2(nstates)


    !: initialize psi_det1
    call get_psi_det( psi(:,1), psi_det1 )
    psi_det1 = psi_det1 / dsqrt( norm_sq(1) )

    icume = 0.d0
    cume_diff(:,1) = 0.d0
    
    do itime=2, ntimes
       
       call get_psi_det( psi(:,itime), psi_det2 ) 
       psi_det2 = psi_det2 / dsqrt( norm_sq(itime) )
    
       record = 0.d0
       icume = dconjg(psi_det2)*psi_det2 - dconjg(psi_det1)*psi_det1
       do i=1, nstates
          record = record + icume(i)
          icume(i) = record
       end do
       
       cume_diff(:,itime) = icume(:)
       psi_det1 = psi_det2


    end do

    open(unit=100, file=trim(f_cumediff_out) )
    write(100,"(a15,1x,9000(i0,1x))") 'time', (i, i=1, nstates )
    do itime=1, ntimes
       write(100,"(f15.10,1x,9000(f10.5,1x))") time(itime), ( cume_diff(i,itime),i=1,nstates )
    end do
    close(100)



  end subroutine get_cumulative_diff
  !:--------------------------:!
  !:--------------------------:!
  subroutine name_files( idir )

    use global_variables
    implicit none

    integer(8), intent(in) :: idir
    character(3) :: cdir
    
    write(cdir,'(i0)') idir
    
    input = 'CI'//trim(e_d)//trim(cdir)//'.bin'
    
    f_mo_out   =  'MO'//trim(e_d)//trim(cdir)//'.out'
    f_pop_out  = 'POP'//trim(e_d)//trim(cdir)//'.out'
    f_rate_out = 'RATE'//trim(e_d)//trim(cdir)//'.out'
    f_hole_out = 'HOLE'//trim(e_d)//trim(cdir)//'.out'
    f_part_out = 'PART'//trim(e_d)//trim(cdir)//'.out'
    f_hole_pop_out = 'HOLE_POP'//trim(e_d)//trim(cdir)//'.out'
    f_part_pop_out = 'PART_POP'//trim(e_d)//trim(cdir)//'.out'
    f_ratios_out = 'RATIOS'//trim(e_d)//trim(cdir)//'.out'
    f_overlap_out = 'OVERLAP'//trim(e_d)//trim(cdir)//'.out'
    f_cumediff_out = 'CUME_DIFF'//trim(e_d)//trim(cdir)//'.out'
    

  end subroutine name_files
  !:--------------------------:!
  !:--------------------------:!
  subroutine write_rhos( noa, nva, nob, nvb, rho_size, time, rho, funit )

    implicit none
    integer(8), intent(in) :: noa, nva, nob, nvb, rho_size, funit
    real(8),   intent(in)  :: time, rho( rho_size, rho_size )

    integer(8) :: i
    integer(8) :: nrorb, norb

    
    nrorb = noa + nva
    norb = noa + nva + nob + nvb

    
90  format( '    occ_a', 14(1x,f8.5) / 100( 15(1x,f8.5) /) ) !: alpha occupied
91  format( '    vir_a', 14(1x,f8.5) / 100( 15(1x,f8.5) /) ) !: alpha virtual
92  format( '    occ_b', 14(1x,f8.5) / 100( 15(1x,f8.5) /) ) !: beta occupied
93  format( '    vir_b', 14(1x,f8.5) / 100( 15(1x,f8.5) /) ) !: beta virtual

    write(funit,"(' time ', f13.5 )") time
    write(funit,90) ( rho(i,i), i=1, noa )             !: alpha occupied
    write(funit,91) ( rho(i,i), i=noa+1, noa+nva )     !: alpha virtual
    write(funit,92) ( rho(i,i), i=nrorb+1, nrorb+nob ) !: beta occupied
    write(funit,93) ( rho(i,i), i=nrorb+nob+1, norb )  !: beta virtual
    

  end subroutine write_rhos
  !:--------------------------:!
  !:--------------------------:!
  
