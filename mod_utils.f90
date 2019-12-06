module utils 

  !: GLOBAL VARIABLES ARE BY DEFAULT SHARED 
  !: UNLESS DECLARED AS THREAD PRIVATE

  
contains
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
  subroutine get_psi_det( i, nstuse, nstates, ci_vec, psi_itime, psi_det )
    
    
    implicit none
    
    integer(8), intent(in) :: nstuse, nstates
    real(8),    intent(in) :: ci_vec(nstates,nstates)
    complex(8), intent(in) :: psi_itime(nstuse)
    complex(8), intent(inout) :: psi_det(nstates)
    integer(8), intent(inout) :: i    

    !complex(8) :: psi_i
    
    
    psi_det = dcmplx(0.d0,0.d0) 
    do i=1, nstuse
       !psi_i = psi_itime(i)
       psi_det(:) = psi_det(:) + psi_itime(i) * ci_vec(:,i)
    end do
    
    
  end subroutine get_psi_det
  !:--------------------------:!
  !:--------------------------:!
end module utils
