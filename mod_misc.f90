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
    complex(8) :: psi_det(nstates), psi2
    !: for writing
    character(5) :: c0, ci
    character(1000) :: ofile, ofile10

    
    write( c0, '(i0)' ) idir
    ofile  = 'SUMMARY_CI'//trim(e_d)//trim(c0)//'.out'
    ofile10 = 'SUMMARY_CI10'//trim(e_d)//trim(c0)//'.out'
    
    open( unit=100, file=trim(ofile) )
    open( unit=110, file=trim(ofile10) )
    

    write(100,'(A)', advance='no') " time(fs) norm2 "
    write( c0, '(i0)' ) 0
    write(100,'(A)', advance='no') trim(c0)//' '
    do i=1, noa
       write( c0, '(i0)' ) -i
       write(100,'(A)', advance='no') trim(c0)//' '
    end do
    write(110,'(A)', advance='no') 'time(fs) norm ground'
    do i=1, nob
       write( c0, '(i0)' ) i
       write(100,'(A)', advance='no') trim(c0)//' '
       write(110,'(A)', advance='no') trim(c0)//' '
    end do
    write(100,'(A)') ''
    write(110,'(A)') ''
    

    !$OMP PARALLEL DEFAULT(NONE), &
    !$OMP PRIVATE( itime, singles_a, singles_b, norm, i, psi_det, ground, ii, ia,  &
    !$OMP SHARED(  ntimes, norm_sq, nstuse, nstates, ci_vec, psi, hole1, time,     &
    !$OMP          ground_save, single_a_save, single_b_save, ci_ten )
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
          ci_ten(ia,itime) = dconjg(psi_det(ia)) * psi_det(ia)
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
       

    do itime=1, ntimes
       write(200,"(f10.7)",advance='no') time(itime), norm_sq(itime), ground_save(itime)
       write(200,"(f10.7)",advance='no' ) ( singles_a_save(i,itime) , i=1, noa )
       write(200,"(f10.7)",advance='no' ) ( singles_b_save(i,itime) , i=1, nob )
       write(200,'(A)') ''
    end do

    do itime=1, ntimes 
       write(110,'(3(1x,f15.10))', advance='no') time(itime), norm, ground    
       write(110,'(1x,f15.10)', advance='no' ) ( psi_ten(ia,itimes), ia=1, 10 )
       write(110,'(A)') ''
    end do
    
    close(110)
    close(200)
    

  end subroutine summarize_coefficients
  !:--------------------------:!
  !:--------------------------:!
  subroutine summarize_coefficients_ip(idir)
    
    use utils
    implicit none

    integer(8), intent(in) :: idir

    integer(8) :: itime, i, x, j, a, b, ia, ii, xx, jj, aa, bb
    real(8)    :: singles(nob), doubles_ab(noa,nob), doubles_bb(nob,nob), norm
    complex(8) :: psi_det(nstates), psi2
    character(5) :: c0, cx, ci
    character(1000) :: ofile, pyfile


    write( c0, '(i0)' ) idir
    ofile  = 'SUMMARY_CI'//trim(e_d)//trim(c0)//'.out'
    pyfile = 'SUMMARY_CI_PY'//trim(e_d)//trim(c0)//'.out' 

    open( unit=100, file=trim(ofile) )
    open( unit=200, file=trim(pyfile) )

    write(200,'(A)', advance='no') " time(fs) norm2 "
    do x=1, nob
       write( c0, '(i0)' ) x 
       write(200,'(A)', advance='no') trim(c0)//' '
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
    
    
    do itime=1, ntimes
       
       singles    = 0.d0
       doubles_ab = 0.d0
       doubles_bb = 0.d0
       
       norm = norm_sq(itime)
       call get_psi_det( i, nstuse, nstates, ci_vec, psi(:,itime), psi_det )
       psi_det = psi_det / dsqrt(norm)   

       do ia=1, nstates

          xx = hole(ia,1)
          ii = hole(ia,2)
          aa = part(ia,1)

          psi2 = dconjg(psi_det(ia)) * psi_det(ia)

          if ( aa .eq. 0 ) then
             singles(xx) = singles(xx) + psi2
          else
             if ( ii.lt.0 ) then
                doubles_ab(-ii,xx) = doubles_ab(-ii,xx) + psi2
             else
                doubles_bb(ii,xx) = doubles_bb(ii,xx) + psi2
                doubles_bb(xx,ii) = doubles_bb(xx,ii) + psi2
             end if
          end if
          
       end do

       write(100,"(' time(fs) ',f7.3, ' norm2 ', f15.10 )") time(itime), norm
       write(100,41) ( singles(i), i=1, nob )
       do x=1, nob
          write(100,42) x, ( doubles_ab(i,x), i=1, noa )
       end do
       do x=1, nob
          write(100,43) x, ( doubles_bb(i,x), i=1, nob ) 
       end do

       write(200,"(f10.7)",advance='no') time(itime), norm, ( singles(i), i=1, nob )
       do x=1, nob
          write(200,"(f10.7)",advance='no' ) ( doubles_ab(i,x), i=1, noa )
       end do
       do x=1, nob
          write(200,"(f10.7)",advance='no' ) ( doubles_bb(i,x), i=x+1, nob )
       end do
       write(200,'(A)') ''
       
    end do
    
41  format( '   singles ', 6x,    10(1x,f8.5) )
42  format( ' doublesAB ', i5'B', 10(1x,f8.5) )
43  format( ' doublesBB ', i5'B', 10(1x,f8.5) ) 

    close(100)
    close(200)
    
    
  end subroutine summarize_coefficients_ip
  !:--------------------------:!
  !:--------------------------:!
  subroutine write_ci_mod(idir)


    use utils
    implicit none

    integer(8), intent(in) :: idir

    integer(8) :: itime, i 
    real(8) :: norm2
    character(5) :: c0
    character(1000) :: ofile

    
    write(c0,'(i0)') idir
    ofile = 'CI_mod'//trim(e_d)//trim(c0)//'.out'

    open( unit=100,file=trim(ofile) )

    do itime=1, ntimes
       
       norm2 = norm_sq(itime)

       write(100,'(f13.5)',advance='no') time(itime)
       write(100,'(f9.5)',advance='no') ( psi(i,itime)*dconjg(psi(i,itime))/norm2, i=1, nstuse )
       write(100,'(A)')

    end do
       


  end subroutine write_ci_mod
  !:--------------------------:!
  !:--------------------------:!  
  subroutine svd_coefficients(idir) 

    implicit none

    integer(8), intent(in) :: idir


  end subroutine svd_coefficients
  !:--------------------------:!
  !:--------------------------:!  
  subroutine write_coefficients(idir)

    implicit none
    
    integer(8), intent(in) :: idir

    integer(8) :: i, itime

    character(5) :: c0

    write( c0, '(i0)' ) idir
    open( unit=144,file='CI_scfyz'//trim(e_d)//trim(c0)//'.out' )
    
    do itime=1, ntimes
       do i=1, nstuse
          write(144,'(f15.10,1x,f15.10)') real(psi(i,itime)), aimag(psi(i,itime))
       end do
    end do
    
    close(144)


  end subroutine write_coefficients
  !:--------------------------:!
  !:--------------------------:!  
end module misc
