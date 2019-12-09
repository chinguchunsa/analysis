module misc

  use global_variables
  implicit none


contains 
  !:--------------------------:!
  !: SUBROUTINE GET_DOS       :!
  !:--------------------------:!
  subroutine get_dos

    implicit none

    integer(8) :: i
    
    
    open(unit=100, file='OUTPUT_DOS_INFO')
    write(100,'(a5,5(1x,a20))') '#', 'energy(eV)', 'tdx_i0(au)', 'tdy_i0(au)', 'tdz_i0(au)', 'vabs(au)'
    do i=1, nstates
       write(100,'(i5,(1x,f20.10),3(1x,e20.10),(1x,f20.10))') i, ci_eig(i)*au2eV, tdx(i,1), tdy(i,1), tdz(i,1), ci_vabs(i,i)
    end do
    close(100)
    
    write(iout,*) ' --> finished writing OUTPUT_DOS_INFO'
    flush(iout)

    
  end subroutine get_dos
  !:--------------------------:!
  !:--------------------------:!
  subroutine summarize_coefficients(idir)
    
    use utils
    use omp_lib
    implicit none

    integer(8), intent(in) :: idir
    
    integer(8) :: itime, i, ii, a, aa, ia
    real(8)    :: norm, ground, singles_a(noa), singles_b(nob)
    real(8)    :: ground_save(ntimes), singles_a_save(noa,ntimes), singles_b_save(nob,ntimes), ci_ten(10,ntimes)
    complex(8) :: psi_det(nstates) 
    !: for writing
    character(5) :: c0, ci
    character(1000) :: ofile, ofile10

    
    write( c0, '(i0)' ) idir
    ofile  = 'SUMMARY_CI'//trim(e_d)//trim(c0)//'.out'
    ofile10 = 'SUMMARY_CI10'//trim(e_d)//trim(c0)//'.out'
    
    open( unit=100, file=trim(ofile) )
    open( unit=110, file=trim(ofile10) )
    

    write(100,'(A)', advance='no') " time(fs) norm2 "
    write(110,'(A)', advance='no') 'time(fs) norm ground'
    write( c0, '(i0)' ) 0
    write(100,'(A)', advance='no') trim(c0)//' '
    write(110,'(A)', advance='no') trim(c0)//' '
    do i=1, noa
       write( c0, '(i0)' ) -i
       write(100,'(A)', advance='no') trim(c0)//' '
    end do
    do i=1, nob
       write( c0, '(i0)' ) i
       write(100,'(A)', advance='no') trim(c0)//' '
    end do
    do i=1, 10
       write(110,'(A)', advance='no') trim(c0)//' '
    end do
    write(100,'(A)') ''
    write(110,'(A)') ''
    

    !$OMP PARALLEL DEFAULT(NONE), &
    !$OMP PRIVATE( itime, singles_a, singles_b, norm, i, psi_det, ground, ii, ia), &
    !$OMP SHARED(  ntimes, norm_sq, nstuse, nstates, ci_vec, psi, hole1, time,     &
    !$OMP          ground_save, singles_a_save, singles_b_save, ci_ten )
    !$OMP DO
    do itime=1, ntimes

       singles_a = 0.d0 
       singles_b = 0.d0
       norm = norm_sq(itime)
       call get_psi_det( i, nstuse, nstates, ci_vec, psi(:,itime), psi_det )
       psi_det = psi_det / dsqrt(norm)
       
       ground = dconjg(psi_det(1)) * psi_det(1)
       
       do ia=2, 10
          ii = hole1(ia)
          i  = abs(ii)
          ci_ten(ia,itime) = dconjg(psi(ia,itime))*psi(ia,itime)/norm
          if ( ii.lt.0 ) singles_a(i) = singles_a(i) + dconjg(psi_det(ia)) * psi_det(ia)
          if ( ii.gt.0 ) singles_b(i) = singles_b(i) + dconjg(psi_det(ia)) * psi_det(ia)
       end do
       do ia=11, nstates
          ii = hole1(ia)
          i  = abs(ii)
          if ( ii.lt.0 ) singles_a(i) = singles_a(i) + dconjg(psi_det(ia)) * psi_det(ia)
          if ( ii.gt.0 ) singles_b(i) = singles_b(i) + dconjg(psi_det(ia)) * psi_det(ia)
       end do

       ground_save(itime) = ground
       singles_a_save(:,itime) = singles_a(:)
       singles_b_save(:,itime) = singles_b(:)
       
    end do
    !$OMP END DO
    !$OMP END PARALLEL
       
    
    !: write 
    do itime=1, ntimes
       write(200,"(f10.7)",advance='no') time(itime), norm_sq(itime), ground_save(itime)
       write(200,"(f10.7)",advance='no' ) ( singles_a_save(i,itime) , i=1, noa )
       write(200,"(f10.7)",advance='no' ) ( singles_b_save(i,itime) , i=1, nob )
       write(200,'(A)') ''
    end do

    !: write
    do itime=1, ntimes 
       write(110,'(3(1x,f15.10))', advance='no') time(itime), norm, ground    
       write(110,'(1x,f15.10)', advance='no' ) ( ci_ten(ia,itime), ia=1, 10 )
       write(110,'(A)') ''
    end do
    
    close(110)
    close(200)
    
    write(iout,'(A)') ' --> generated file '//trim(ofile)
    write(iout,'(A)') ' --> generated file '//trim(ofile10)
    flush(iout)
    

  end subroutine summarize_coefficients
  !:--------------------------:!
  !:--------------------------:!
  subroutine summarize_coefficients_ip(idir)
    
    use utils
    implicit none

    integer(8), intent(in) :: idir

    integer(8) :: itime, i, x, a, b, ia, ii, xx, aa
    real(8) :: singles(nob), doubles_ab(noa,nob), doubles_bb(nob,nob), norm
    real(8) :: singles_save(nob,ntimes), doubles_ab_save(noa,nob,ntimes), doubles_bb_save(nob,nob,ntimes), ci10_save(10,ntimes)
    complex(8) :: psi_det(nstates)
    character(5) :: c0, cx, ci
    character(1000) :: ofile, ofile10

    
    !: only set up for single BETA excitation

    
    write( c0, '(i0)' ) idir
    ofile   = 'SUMMARY_CI'//trim(e_d)//trim(c0)//'.out'
    ofile10 = 'SUMMARY_CI10'//trim(e_d)//trim(c0)//'.out'

    open( unit=110, file=trim(ofile10) )
    open( unit=200, file=trim(ofile) )
    
    write(200,'(A)', advance='no') " time(fs) norm2 "
    write(110,'(A)', advance='no') " time(fs) norm2 "
    do x=1, nob
       write( c0, '(i0)' ) x 
       write(200,'(A)', advance='no') trim(c0)//' '
       write(110,'(A)', advance='no') trim(c0)//' '
    end do
    do x=1, nob
       write( c0, '(i0)' ) x
       do i=1, noa
          write( ci, '(i0)' ) -i
          write(200,'(A)',advance='no') '('//trim(c0)//','//trim(ci)//') '
       end do
    end do
    do x=1, nob
       write( c0, '(i0)' ) x
       do i=(x+1), nob
          write( ci, '(i0)' ) i
          write(200,'(A)',advance='no') '('//trim(c0)//','//trim(ci)//') '
       end do
    end do
    write(200,'(A)') ''
    write(110,'(A)') ''

    
    !$OMP PARALLEL DEFAULT(NONE), &
    !$OMP PRIVATE( itime, i, singles, doubles_ab, doubles_bb, psi_det, norm, xx, ii, aa ), &
    !$OMP SHARED(  ntimes, nob, psi, ci_vec, nstuse, nstates, norm_sq, hole, part, singles_save,&
    !$OMP          doubles_ab_save, doubles_bb_save, ci10_save )
    !$OMP DO    
    do itime=1, ntimes
       
       singles    = 0.d0
       doubles_ab = 0.d0
       doubles_bb = 0.d0
       
       norm = norm_sq(itime)
       call get_psi_det( i, nstuse, nstates, ci_vec, psi(:,itime), psi_det )
       psi_det = psi_det / dsqrt(norm)   
       
       do ia=1, nob          
          xx = hole(ia,1)
          ii = hole(ia,2)
          aa = part(ia,1)
          if ( aa .eq. 0 ) singles(xx) = singles(xx) + dconjg(psi_det(ia)) * psi_det(ia) 
          ci10_save(ia,itime) = dconjg( psi(ia,itime) ) * psi(ia,itime) / norm          
       end do

       do ia=nob+1, nstates          
          xx = hole(ia,1)
          ii = hole(ia,2)
          aa = part(ia,1)
          if ( aa .eq. 0 ) then
             singles(xx) = singles(xx) + dconjg(psi_det(ia)) * psi_det(ia)
          else
             if ( ii.lt.0 ) doubles_ab(-ii,xx) = doubles_ab(-ii,xx) + dconjg(psi_det(ia)) * psi_det(ia)
             if ( ii.gt.0 ) doubles_bb(ii,xx)  = doubles_bb(ii,xx) + dconjg(psi_det(ia)) * psi_det(ia)
          end if
       end do

       singles_save(:,itime) = singles(:)
       doubles_ab_save(:,:,itime) = doubles_ab(:,:)
       doubles_bb_save(:,:,itime) = doubles_bb(:,:)

    end do
    !$OMP END DO
    !$OMP END PARALLEL
    
       
    do itime=1, ntimes
       !: write singlesa and doubles
       write(200,"(f10.7)",advance='no') time(itime), norm, ( singles_save(i,itime), i=1, nob )
       do x=1, nob
          write(200,"(f10.7)",advance='no' ) ( doubles_ab_save(i,x,itime), i=1, noa )
       end do
       do x=1, nob
          write(200,"(f10.7)",advance='no' ) ( doubles_bb_save(i,x,itime), i=x+1, nob )
       end do
       write(200,'(A)') ''
       !: write ci10_save
       write(110,"(f10.7)",advance='no') time(itime), norm, ( ci10_save(i,itime), i=1, nob )
       write(110,'(A)') ''
    end do

    close(200)
    close(110)

    write(iout,'(A)') ' --> generated file '//trim(ofile)
    write(iout,'(A)') ' --> generated file '//trim(ofile10)
    flush(iout)    
    
    
  end subroutine summarize_coefficients_ip
  !:--------------------------:!
  !:--------------------------:!
end module misc
