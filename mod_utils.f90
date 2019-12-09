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
    write(iout,*) ' --> finished assigning cis hole and part indices'
    flush(iout)
    
    
  end subroutine get_indices
  !:--------------------------:!
  !:--------------------------:!
  subroutine get_indices_ip
    
    use global_variables
    implicit none
    
    integer(8) :: x, i, a, ia 
    
 
    allocate( hole(nstates,2), part(nstates,1), hole_indices(-noa:nob,-noa:nob) )
    hole = 0 ; part = 0
    
    !: indices for CIS for now
    ia = 0

    
    !: only beta electrons will be ionized
    singles : do x=1, nob
       ia = ia + 1
       hole(ia,1) = x
       hole(ia,2) = 0
       part(ia,1) = 0
       hole_indices(x,0) = ia
       hole_indices(0,x) = ia 
    end do singles
    
    aa : do i=1, noa
       ab_doubles : do x=1, nob
          hole_indices(-i,x) = ia + 1 
          hole_indices(x,-i) = ia + 1 
          do a=1, nva
             ia = ia + 1
             hole(ia,1) = x
             hole(ia,2) = -i
             part(ia,1) = -a - noa             
          end do
       end do ab_doubles
    end do aa

    bb_doubles : do x=1, nob
       bb : do i=x+1, nob
          hole_indices(i,x) = ia + 1 
          hole_indices(x,i) = ia + 1 
          do a=1, nvb
             ia = ia + 1
             hole(ia,1) = x
             hole(ia,2) = i
             part(ia,1) = a + nob
          end do
       end do bb
    end do bb_doubles


    write(iout,*) ' --> allocated hole, part '
    write(iout,*) ' --> finished assigning ip hole and part indices'
    flush(iout)
    
    
  end subroutine get_indices_ip
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
  subroutine get_1density_ip(psi_det, nstates, hole, part, rho_a, rho_b, nrorb, noa, nob) 
  
    implicit none

    integer(8) :: nstates, nrorb, noa, nob
    integer(8), intent(in) :: hole(nstates,2), part(nstates,1)
    complex(8), intent(in) :: psi_det(nstates)
    real(8), intent(inout) :: rho_a(nrorb,nrorb), rho_b(nrorb,nrorb)


    integer(8) :: ii, jj, aa, bb, xx, yy
    integer(8) :: ia, jb, i, a, j, b, x, y
    complex(8) :: psi_ia


    !: initialize
    rho_a = 0.d0 
    rho_b = 0.d0

    do i=1, noa
       rho_a(i,i) = 1.d0
    end do

    do i=1, nob
       rho_b(i,i) = 1.d0
    end do


    do ia=1, nstates
       
       xx = hole(ia,1)
       ii = hole(ia,2)
       aa = part(ia,1)
       psi_ia = dconjg(psi_det(ia))

       do jb=1, nstates

          yy = hole(jb,1)
          jj = hole(jb,2)
          bb = part(jb,1)

          !: <X|A|Y>
          SS: if ( ii.eq.0 .and. jj.eq.0 ) then
             rho_b(yy,xx) = rho_b(yy,xx) - real(psi_ia*psi_det(jb))
             go to 78
          end if SS
          
          !: <iX->a|A|X>  or  <IX->A|A|X>
          SD1: if ( YY.eq.XX .and. jj.eq.0 ) then
             if ( ii.lt.0 ) rho_a(-aa,-ii) = rho_a(-aa,-ii) + real(psi_ia*psi_det(jb))
             if ( ii.gt.0 ) rho_b(aa,ii)   = rho_b(aa,ii)   + real(psi_ia*psi_det(jb))
             go to 78
          end if SD1
          !: <X|A|jX->b>
          SD2: if ( YY.eq.XX .and. ii.eq.0 ) then
             if ( jj.lt.0 ) rho_a(-jj,-bb) = rho_a(-jj,-bb) + real(psi_ia*psi_det(jb))
             if ( jj.gt.0 ) rho_b(jj,bb)   = rho_b(jj,bb)   + real(psi_ia*psi_det(jb))
             go to 78
          end if SD2
          !: <YX->A|A|Y>
          SD3 : if ( YY.eq.II .and. jj.eq.0 ) then
             rho_b(aa,xx) = rho_b(aa,xx) - real(psi_ia*psi_det(jb))
             go to 78
          end if SD3
          !: <X|A|XY->B>
          SD4: if ( XX.eq.JJ .and. ii.eq.0 ) then
             rho_b(yy,bb) = rho_b(yy,bb) - real(psi_ia*psi_det(jb))
             go to 78
          end if SD4
          
          !: doubles
          
          !: <ix->a|A|jx->a>
          kdelta_xy_ab : if ( xx.eq.yy .and. aa.eq.bb ) then
             if ( ii.lt.0 ) rho_a(-jj,-ii) = rho_a(-jj,-ii) - real(psi_ia*psi_det(jb))
             if ( ii.gt.0 ) rho_b(jj,ii)   = rho_b(jj,ii)   - real(psi_ia*psi_det(jb))
          end if kdelta_xy_ab
          
          !: <ix->a|A|ix->b>
          kdelta_xy_ij : if ( xx.eq.yy .and. ii.eq.jj ) then
             if ( aa.lt.0 ) rho_a(-aa,-bb) = rho_a(-aa,-bb) + real(psi_ia*psi_det(jb))
             if ( aa.gt.0 ) rho_b(aa,bb)   = rho_b(aa,bb)   + real(psi_ia*psi_det(jb))
          end if kdelta_xy_ij
          
          !: <ix->a|A|iy->a>
          kdelta_ij_ab : if ( ii.eq.jj .and. aa.eq.bb ) then
             rho_b(yy,xx) = rho_b(yy,xx) - real(psi_ia*psi_det(jb))
          end if kdelta_ij_ab
          
          !: <ix->a|A|xy->a>
          kdelta_jx_ab : if ( JJ.eq.XX .and. AA.eq.BB ) then
             rho_b(yy,ii) = rho_b(yy,ii) + real(psi_ia*psi_det(jb))
          end if kdelta_jx_ab
          
          !: <yx->a|A|jy->a>
          kdelta_yi_ab : if (YY.eq.II .and. AA.eq.BB ) then
             rho_b(jj,xx) = rho_b(jj,xx) + real(psi_ia*psi_det(jb))
          end if kdelta_yi_ab
          
          !: <yx->a|A|xy->b>
          kdelta_iy_xj : if ( II.eq.YY .and. XX.eq.JJ ) then
             rho_b(aa,bb) = rho_b(aa,bb) - real(psi_ia*psi_det(jb))
          end if kdelta_iy_xj
          
78        continue
          
       end do
    end do
    
    

  end subroutine get_1density_ip
  !:--------------------------:!
  !:--------------------------:!
end module utils
