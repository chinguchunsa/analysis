module density 


  use global_variables
  use utils


contains
  !:--------------------------:!
  !:--------------------------:!
  subroutine get_density( idir )

    
    use omp_lib
    implicit none


    integer(8), intent(in) :: idir
    real(8), parameter :: tol  = 1.d-10  !: write elements with magnitude > 10-10

    integer(8) :: itime, ia, jb, i, a, j, b, ii, jj, aa, bb
    real(8) :: norm, const, pop0
    real(8) :: rho_a(nrorb,nrorb), rho_b(nrorb,nrorb)
    real(8) :: rho_save_a(nrorb,nrorb,start_time-1:end_time+1), rho_save_b(nrorb,nrorb,start_time-1:end_time+1)
    complex(8) :: psi_det(nstates), psi_ia, psi0
    
    !: write
    character(5) :: c0
    character(1000) :: ofile

    
    const = dsqrt(2.d0)
    pop0  = 2.d0
    if ( unrestricted ) then
       const = 1.d0
       pop0  = 1.d0
    end if

    
    rho_save_a = 0.d0
    rho_save_b = 0.d0


    !$OMP PARALLEL DEFAULT(NONE), &
    !$OMP PRIVATE( itime, i, j, a, b, ia, jb, ii, jj, aa, bb,   &
    !$OMP          norm, rho_a, rho_b, psi_det, psi_ia, psi0 ), &
    !$OMP SHARED( start_time, end_time, noa, nob, nstuse, nstates, unrestricted, const, pop0, &
    !$OMP         psi, ci_vec, hole1, part1, rho_save_a, rho_save_b, norm_sq )
    !$OMP DO
    do itime=start_time, end_time

       !: save norm
       norm = norm_sq(itime)
       
       !: get psi into determinantal basis
       call get_psi_det( i, nstuse, nstates, ci_vec, psi(:,itime), psi_det )
       psi_det = psi_det / dsqrt(norm)
       

       psi0 = psi_det(1)
       
       
       !:initialize
       rho_a = 0.d0
       do i=1, noa
          rho_a(i,i) = pop0
       end do
       if( unrestricted ) then
          rho_b = 0.d0
          do i=1, nob
             rho_b(i,i) = 1.d0
          end do
       end if
       
       
       do ia=2, nstates

          
          psi_ia = dconjg( psi_det(ia) )
          

          ii = hole1(ia)  ;  i = abs(ii)
          aa = part1(ia)  ;  a = abs(aa)
          
          
          if ( ii.lt.0 ) then
             !: ground excited
             rho_a(a,i) = rho_a(a,i) + const*real( psi_ia * psi0 )
             rho_a(i,a) = rho_a(a,i)
          else
             !: ground excited
             rho_b(a,i) = rho_b(a,i) + const*real( psi_ia * psi0 )
             rho_b(i,a) = rho_b(a,i)
          end if
          
          
          do jb=2, nstates
             
             
             jj = hole1(jb)  ;  j = abs(jj)
             bb = part1(jb)  ;  b = abs(bb)
             
             
             if ( ii.eq.jj ) then
                if ( ii.lt.0 ) then
                   rho_a(b,a) = rho_a(b,a) + real( psi_ia * psi_det(jb) )
                else
                   rho_b(b,a) = rho_b(b,a) + real( psi_ia * psi_det(jb) )
                end if
             end if
             
             if ( aa.eq.bb ) then
                if ( aa.lt.0 ) then
                   rho_a(j,i) = rho_a(j,i) - real( psi_ia * psi_det(jb) )
                else
                   rho_b(j,i) = rho_b(j,i) - real( psi_ia * psi_det(jb) )
                end if
             end if
             
             
          end do
       end do

       
       rho_save_a(:,:,itime) = rho_a(:,:)
       if( unrestricted ) rho_save_b(:,:,itime) = rho_b(:,:)

       
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    
    
    write( c0, '(i0)' ) idir
    ofile = 'POP'//trim(e_d)//trim(c0)//'.out'
    open( unit=100, file=trim(ofile) )
    
    write( 100,'(A)',advance='no' ) 'time(fs) '
    write( 100,'(100(i0,1x))', advance='no' ) ( i, i=1, nrorb )
    if ( unrestricted ) write( 100, '(100(i0,1x))', advance='no' ) ( i, i=nrorb+1, norb )
    write( 100, '(A)' ) ''
    
    do itime=start_time, end_time
       write( 100,'(f15.10)', advance='no' ) time(itime), ( rho_save_a(i,i,itime), i=1, nrorb )
       if ( unrestricted ) write(100,'(f15.10)',advance='no') ( rho_save_b(i,i,itime), i=1, nrorb )
       write(100,'(A)') ''
    end do

    write(iout,*) ' --> generated file '//trim(ofile)
    flush(iout)
    close(100)

    if ( Qprint_movie ) then
       
       write( c0, '(i0)' ) idir
       ofile = 'DENSITY_MOVIE'//trim(e_d)//trim(c0)//'.out'
       open( unit=100, file=trim(ofile) )
       
       write( 100, '(i0)' ) ntimes
       write( 100, '(i0)' ) noa
       write( 100, '(i0)' ) nob
       write( 100, '(i0)' ) noa + nva


       do itime=start_time, end_time
          
          write( 100, '(f13.5)' ) time(itime)
          
          !: alpha
          rho_a(:,:) = rho_save_a(:,:,itime)          
          do i=1, nrorb
             if ( abs(rho_a(i,i)).gt.tol ) write( 100,"(2(i4,1x),f14.10,', ')",advance='no' ) i,i, rho_a(i,i)
             do j=i+1, nrorb
                if( abs(rho_a(j,i)) .gt. tol ) &
                     write( 100,"( 2(i4,1x),f14.10,', ' )",advance='no' ) i,j, 2.d0*rho_a(j,i)
             end do
          end do
          
          !: beta
          if( unrestricted ) then 
             rho_b(:,:) = rho_save_b(:,:,itime)
             do i=1, nrorb
                if ( abs(rho_b(i,i)).gt.tol ) write( 100, "(2(i4,1x),f14.10,', ')", advance='no' ) i,i, rho_b(i,i)
                do j=i+1, nrorb
                   if ( abs(rho_b(j,i)) .gt. tol ) & 
                        write( 100,"( 2(i4,1x),f14.10,', ')",advance='no' ) i+nrorb, j+nrorb, 2.d0*rho_b(j,i)
                end do
             end do
          end if
          write( 100,'(A)' ) ''
          
       end do
       
       close(100)
       write(iout,*) ' --> generated file '//trim(ofile)
       flush(iout)
       
    end if

       

  end subroutine get_density
  !:--------------------------:!
  !:--------------------------:!   
  subroutine get_density_ip(idir)
    
    use omp_lib
    implicit none

    integer(8), intent(in) :: idir
    
    integer(8) :: itime, ia, jb, i, a, j, b, x, y, ii, jj, aa, bb, xx, yy
    real(8) :: norm, const, pop0
    real(8) :: rho_a(nrorb,nrorb), rho_b(nrorb,nrorb)
    real(8) :: rho_save_a(nrorb,nrorb,start_time-1:end_time+1), rho_save_b(nrorb,nrorb,start_time-1:end_time+1)
    complex(8) :: psi_det(nstates), psi_ia, psi0
    
    !: write
    character(5) :: c0
    character(1000) :: ofile

    
    rho_save_a = 0.d0
    rho_save_b = 0.d0

    
    !$OMP PARALLEL DEFAULT(NONE), &
    !$OMP PRIVATE( itime, i, j, a, b, x, y, ia, jb, ii, jj, aa, bb, xx, yy, &
    !$OMP          norm, rho_a, rho_b, psi_det, psi_ia, psi0 ), &
    !$OMP SHARED( start_time, end_time, noa, nob, nstuse, nstates, unrestricted, const, pop0, &
    !$OMP         psi, ci_vec, hole, part, rho_save_a, rho_save_b, norm_sq )
    !$OMP DO
    do itime=start_time, end_time

       !: save norm
       norm = norm_sq(itime)
       
       !: get psi into determinantal basis
       call get_psi_det( i, nstuse, nstates, ci_vec, psi(:,itime), psi_det )
       psi_det = psi_det / dsqrt(norm)
       
       
       !:initialize
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
          psi0 = dconjg(psi_det(ia))
          
          
          do jb=1, nstates
             
             yy = hole(jb,1)
             jj = hole(jb,2)
             bb = part(jb,1)
             

             !: <X|A|Y>
             SS: if ( ii.eq.0 .and. jj.eq.0 ) then
                rho_b(yy,xx) = rho_b(yy,xx) - real(psi0*psi_det(jb))
                go to 78
             end if SS

             
             !: <iX->a|A|X>  or  <IX->A|A|X>
             SD1: if ( YY.eq.XX .and. jj.eq.0 ) then
                if ( ii.lt.0 ) rho_a(-aa,-ii) = rho_a(-aa,-ii) + real(psi0*psi_det(jb))
                if ( ii.gt.0 ) rho_b(aa,ii)   = rho_b(aa,ii)   + real(psi0*psi_det(jb))
                go to 78
             end if SD1
             !: <X|A|jX->b>  
             SD2: if ( YY.eq.XX .and. ii.eq.0 ) then
                if ( jj.lt.0 ) rho_a(-jj,-bb) = rho_a(-jj,-bb) + real(psi0*psi_det(jb))
                if ( jj.gt.0 ) rho_b(jj,bb)   = rho_b(jj,bb)   + real(psi0*psi_det(jb))
                go to 78
             end if SD2
             !: <YX->A|A|Y>
             SD3 : if ( YY.eq.II .and. jj.eq.0 ) then
                rho_b(aa,xx) = rho_b(aa,xx) - real(psi0*psi_det(jb))
                go to 78
             end if SD3
             !: <X|A|XY->B>
             SD4: if ( XX.eq.JJ .and. ii.eq.0 ) then
                rho_b(yy,bb) = rho_b(yy,bb) - real(psi0*psi_det(jb))
                go to 78
             end if SD4

             
             !: doubles 

             !: <ix->a|A|jx->a>
             kdelta_xy_ab : if ( xx.eq.yy .and. aa.eq.bb ) then
                if ( ii.lt.0 ) rho_a(-jj,-ii) = rho_a(-jj,-ii) - real(psi0*psi_det(jb))
                if ( ii.gt.0 ) rho_b(jj,ii)   = rho_b(jj,ii)   - real(psi0*psi_det(jb))
             end if kdelta_xy_ab

             !: <ix->a|A|ix->b>
             kdelta_xy_ij : if ( xx.eq.yy .and. ii.eq.jj ) then
                if ( aa.lt.0 ) rho_a(-aa,-bb) = rho_a(-aa,-bb) + real(psi0*psi_det(jb))
                if ( aa.gt.0 ) rho_b(aa,bb)   = rho_b(aa,bb)   + real(psi0*psi_det(jb))
             end if kdelta_xy_ij

             !: <ix->a|A|iy->a>
             kdelta_ij_ab : if ( ii.eq.jj .and. aa.eq.bb ) then
                rho_b(yy,xx) = rho_b(yy,xx) - real(psi0*psi_det(jb))
             end if kdelta_ij_ab

             !: <ix->a|A|xy->a>
             kdelta_jx_ab : if ( JJ.eq.XX .and. AA.eq.BB ) then
                rho_b(yy,ii) = rho_b(yy,ii) + real(psi0*psi_det(jb))
             end if kdelta_jx_ab

             !: <yx->a|A|jy->a>
             kdelta_yi_ab : if (YY.eq.II .and. AA.eq.BB ) then
                rho_b(jj,xx) = rho_b(jj,xx) + real(psi0*psi_det(jb))
             end if kdelta_yi_ab

             !: <yx->a|A|xy->b>
             kdelta_iy_xj : if ( II.eq.YY .and. XX.eq.JJ ) then
                rho_b(aa,bb) = rho_b(aa,bb) - real(psi0*psi_det(jb))
             end if kdelta_iy_xj
             
             
78           continue
             
             
          end do
       end do
       
       rho_save_a(:,:,itime) = rho_a(:,:)
       rho_save_b(:,:,itime) = rho_b(:,:)
       
       
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    

89  format( ' time(fs)', f10.5, 2x, f10.5 )
90  format( '    occ_a', 14(1x,f8.5) / 100( 15(1x,f8.5) /) ) !: alpha occupied
91  format( '    vir_a', 14(1x,f8.5) / 100( 15(1x,f8.5) /) ) !: alpha virtual
92  format( '    occ_b', 14(1x,f8.5) / 100( 15(1x,f8.5) /) ) !: beta occupied
93  format( '    vir_b', 14(1x,f8.5) / 100( 15(1x,f8.5) /) ) !: beta virtual
    
    
    write( c0, '(i0)' ) idir
    ofile = 'POP_IP'//trim(e_d)//trim(c0)//'.out'
    open( unit=100, file=trim(ofile) )
    
    do itime=start_time, end_time
       
       write( 100, 89 ) time(itime) , norm_sq(itime)
       write( 100, 90 ) ( rho_save_a(i,i,itime), i=1, noa )
       write( 100, 91 ) ( rho_save_a(i,i,itime), i=noa+1, nrorb )
       write( 100, 92 ) ( rho_save_b(i,i,itime), i=1, nob )
       write( 100, 93 ) ( rho_save_b(i,i,itime), i=nob+1, nrorb )
       
    end do
    
    close(100)


  end subroutine get_density_ip
  !:--------------------------:!
  !:--------------------------:!   
end module density
