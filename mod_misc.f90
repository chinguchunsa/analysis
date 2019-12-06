module misc

  use global_variables
  implicit none


contains 
  !:--------------------------:!
  !:--------------------------:!
  subroutine get_dos

    implicit none

    real(8), parameter :: de = 0.00001 !: in eV
    
    integer(8) :: istate, ipt, npts
    real(8)    :: eig, gamma, energy
    real(8), allocatable :: dos(:)


    open(unit=100, file='OUTPUT_DOS_INFO')
    write(100,'(a5,5(1x,a20))') '#', 'energy(eV)', 'tdx_i0(au)', 'tdy_i0(au)', 'tdz_i0(au)', 'vabs(au)'
    do ipt=1, nstates
       write(100,'(i5,(1x,f20.10),3(1x,e20.10),(1x,f20.10))') ipt, ci_eig(ipt)*au2eV, tdx(ipt,1), tdy(ipt,1), tdz(ipt,1), ci_vabs(ipt,ipt)
    end do
    close(100)

    npts = int ( (ci_eig(nstuse)-ci_eig(1))*au2eV / de ) 
    

    allocate( dos(npts) )
    dos = 0.d0
    
    
    do istate=2, nstuse
       
       gamma  = 0.5d0 * ci_vabs(istate,istate) * au2eV 
       eig    = ci_eig(istate) * au2eV
       
       energy = 0.d0
       do ipt = 1 , npts
          energy = energy + de
          dos(ipt) = dos(ipt) + gamma / ( (energy - eig)**2 + gamma**2 )
       end do

    end do
          
    
    open(unit=100, file='OUTPUT_DOS')
    do ipt=1, npts
       write(100,'(f20.10,1x,f20.10)') dble(ipt)*de , dos(ipt)
    end do
    close(100) 


  end subroutine get_dos
  !:--------------------------:!
  !:--------------------------:!
  subroutine summarize_coefficients(idir)
    
    use utils
    implicit none

    integer(8), intent(in) :: idir

    integer(8) :: itime, i
    real(8)    :: ground(ntimes), excited(ntimes)
    character(5) :: c0
    character(1000) :: ofile


    ground  = 0.d0
    excited = 0.d0

    !$OMP PARALLEL DEFAULT(NONE), &
    !$OMP PRIVATE( itime ), &
    !$OMP SHARED( ntimes, norm_sq, psi, ground, excited )
    !$OMP DO
    do itime=1, ntimes
       
       ground(itime) = psi(1,itime) * dconjg( psi(1,itime) ) / norm_sq(itime)
       excited(itime) = 1.d0 - ground(itime)
       
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    write( c0, '(i0)' ) idir
    ofile = 'SUMMARY_CI'//trim(e_d)//trim(c0)//'.out'
    open( unit=100, file=trim(ofile) )
    write(100, '(a7, a9, a9)' ) 'time(fs)', 'ground_pop', 'excited_pops'
    do itime = 1, ntimes
       write(100, '(f7.3,20(f9.5))' ) time(itime), &
            ( real(psi(i,itime)*dconjg(psi(i,itime)))/norm_sq(itime), i=1,10 )
    end do
    close(100)
       
    
    
  end subroutine summarize_coefficients
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
end module misc
