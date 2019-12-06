module rates

  use global_variables
  use utils

contains
  !:--------------------------:!
  !:--------------------------:!
  subroutine get_rate( idir )
    

    !: WARNING:  technically, rho is complex.  
    !: WARNING:  but the imaginary component doesn't matter here.
    !: WARNING:  it eventually cancels out 
    !: WARNING:  so will declare rho real to save memory
    
    
    use omp_lib
    implicit none
    
    integer(8), intent(in) :: idir
    real(8), parameter :: tol = 1.d-10
    
    integer(8) :: itime, i, j, a, b, ia, jb, ii, jj, aa, bb
    real(8)  :: norm, dt, from_occ, from_vir, tloss, const, pop0
    real(8)  :: rhoV_a(nrorb,nrorb), rhoV_b(nrorb,nrorb), loss_save( norb, start_time-1:end_time+1 )
    real(8)  :: rho_saveA(nrorb,nrorb,start_time-1:end_time+1), rho_saveB(nrorb,nrorb,start_time-1:end_time+1)
    real(8)  :: temp(nrorb,nrorb)
    real(8)  :: trate(start_time-1:end_time+1)
    complex(8) :: psi_det(nstates), psi_ia, psi0
    
    
    !: write
    character(5)    :: c0
    character(1000) :: ofile
    
    
    loss_save = 0.d0
    trate  = 0.d0
    rho_saveA = 0.d0
    rho_saveB = 0.d0
    dt = time(2) - time(1)    

    const = dsqrt(2.d0)
    pop0  = 2.d0
    if ( unrestricted ) then
       const = 1.d0
       pop0  = 1.d0
    end if
    
    
    !$OMP PARALLEL DEFAULT(NONE), &
    !$OMP PRIVATE( itime, i, j, a, b, ia, jb, ii, jj, aa, bb,     &
    !$OMP          norm, rhoV_a, rhoV_b, temp, psi_det, psi_ia, psi0 ), &
    !$OMP SHARED(  start_time, end_time, noa, nob, nrorb, nstuse, nstates, dt, au2fs, unrestricted, &
    !$OMP          loss_save, trate, vabs_b, vabs_a, norm_sq, psi, ci_vec, hole1, part1,            &
    !$OMP          rho_saveA, rho_saveB, const, pop0 )    
    !$OMP DO
    do itime = start_time, end_time
       
       !: save norm
       norm = norm_sq(itime)
       
       !: get psi into determinantal basis
       call get_psi_det( i, nstuse, nstates, ci_vec, psi(:,itime), psi_det )
       psi_det = psi_det / dsqrt(norm)
       
       !: initialize
       rhoV_a = 0.d0
       do i=1, noa
          rhoV_a(i,i) = pop0
       end do
       if( unrestricted ) then
          rhoV_b = 0.d0
          do i=1, nob
             rhoV_b(i,i) = 1.d0
          end do
       end if
       
       psi0 = psi_det(1)
       
       !: excited states
       do ia=2, nstates
          
          psi_ia = dconjg( psi_det(ia) )
          ii = hole1(ia)  ;  i = abs(ii)
          aa = part1(ia)  ;  a = abs(aa)
          
          if ( ii.lt.0 ) then
             !: ground to excited 
             rhoV_a(a,i) = const * psi_ia * psi0
             rhoV_a(i,a) = const * psi_ia * psi0
          else
             !: ground to excited
             rhoV_b(a,i) = const * psi_ia * psi0
             rhoV_b(i,a) = const * psi_ia * psi0
          end if
          
          do jb=2, nstates
             
             jj = hole1(jb)  ;  j = abs(jj)
             bb = part1(jb)  ;  b = abs(bb)
             
             if( jj.eq.ii ) then
                if( ii.lt.0 ) then
                   rhoV_a(b,a) = rhoV_a(b,a) + real( psi_ia * psi_det(jb) )
                else
                   rhoV_b(b,a) = rhoV_b(b,a) + real( psi_ia * psi_det(jb) )
                end if
             end if
             
             if( aa .eq. bb ) then
                if( aa.lt.0 ) then
                   rhoV_a(j,i) = rhoV_a(j,i) - real( psi_ia * psi_det(jb) )
                else
                   rhoV_b(j,i) = rhoV_b(j,i) - real( psi_ia * psi_det(jb) )
                end if
             end if
             
          end do
          
       end do
       
       
       !: get matrix product, expectation value
       temp = matmul( rhoV_a, vabs_a )
       rhoV_a = temp
       rho_saveA(:,:,itime) = 2.d0 * norm * rhoV_a(:,:) / au2fs * dt
       do i=1, nrorb
          loss_save(i,itime) = -2.d0 * norm * rhoV_a(i,i) / au2fs * dt
          trate(itime) = trate(itime) - 2.d0 * norm * rhoV_a(i,i) / au2fs
       end do

       if( unrestricted ) then
          temp = matmul( rhoV_b, vabs_b )
          rhoV_b = temp
          rho_saveB(:,:,itime) = 2.d0 * norm * rhoV_b(:,:) / au2fs * dt
          do i=1, nrorb
             loss_save(nrorb+i,itime) = -2.d0 * norm * rhoV_b(i,i)  / au2fs * dt
             trate(itime) = trate(itime) - 2.d0 * norm * rhoV_b(i,i) / au2fs
          end do
       end if
       
       
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    
    !: get cumulative
    do itime=start_time, end_time
       loss_save(:,itime) = loss_save(:,itime-1) + loss_save(:,itime)
       rho_saveA(:,:,itime) = rho_saveA(:,:,itime-1) + rho_saveA(:,:,itime)
       rho_saveB(:,:,itime) = rho_saveB(:,:,itime-1) + rho_saveB(:,:,itime)
    end do


    
89  format( ' time(fs)', f10.5, 1x, 'rate(1/fs)', f13.9, 1x, 'predicted_loss', f13.9, 1x, &
         'actual_loss', f13.9, 1x, 'from_occ', f13.9, 1x, 'from_vir', f13.9 )
90  format( '    occ_a', 14(1x,f8.5) / 100( 15(1x,f8.5) /) ) !: alpha occupied
91  format( '    vir_a', 14(1x,f8.5) / 100( 15(1x,f8.5) /) ) !: alpha virtual
92  format( '    occ_b', 14(1x,f8.5) / 100( 15(1x,f8.5) /) ) !: beta occupied
93  format( '    vir_b', 14(1x,f8.5) / 100( 15(1x,f8.5) /) ) !: beta virtual
    
    
    if ( Qprint_pretty ) then
       
       write( c0, '(i0)' ) idir
       ofile = 'LOSS'//trim(e_d)//trim(c0)//'.out'
       open( unit=100, file=trim(ofile) )       
       
       do itime=start_time, end_time
          
          tloss = sum( loss_save(:,itime-1) )
          from_occ = 0.d0
          do i=1, noa 
             from_occ = from_occ + loss_save(i,itime-1)
          end do
          if ( unrestricted ) then
             do i=1, nob
                from_occ = from_occ + loss_save(nrorb+i,itime-1)
             end do
          end if
          
          write( 100, 89)  time(itime), trate(itime), tloss, norm_sq(itime)-1.d0, & 
               from_occ, tloss-from_occ
          write( 100, 90 ) ( loss_save(i,itime-1), i=1, noa )       !: alpha occupied
          write( 100, 91 ) ( loss_save(i,itime-1), i=noa+1, nrorb ) !: alpha virtual
          if( unrestricted ) then
             write( 100, 92 ) ( loss_save(i,itime-1), i=nrorb+1, nrorb+nob ) !: beta occupied
             write( 100, 93 ) ( loss_save(i,itime-1), i=nrorb+nob+1, norb )  !: beta virtual
          end if
       end do
       
       close(100)


    end if
    
100 format( f10.5, 1000(f15.10) )

    
    if ( Qprint_python ) then

       write( c0, '(i0)' ) idir
       ofile = 'LOSS_PY'//trim(e_d)//trim(c0)//'.out' 
       open( unit=100, file=trim(ofile) )
       
       !: write key
       write( 100, '(A)', advance='no' ) ' time(fs) rate(1/fs) sum_loss '
       write( 100,"(100(i0,1x))", advance='no' ) ( i, i=1, norb )
       write( 100, '(A)' ) ''
       
       do itime=start_time, end_time
          write( 100,100 ) time(itime), trate(itime), sum(loss_save(:,itime-1)), ( loss_save(i,itime-1), i=1, norb )
       end do
       
       close(100)
       
    end if

    if ( Qprint_movie ) then

       write( c0, '(i0)' ) idir
       ofile = 'LOSS_MOVIE_FILE'//trim(e_d)//trim(c0)//'.out'
       open( unit=100, file=trim(ofile) )
       write(100,'(i0)') 1
       write(100,'(i0)') noa
       write(100,'(i0)') nob
       write(100,'(i0)') noa+nva
       
       !: only for last timestep
       rhoV_a(:,:) = rho_saveA(:,:,end_time)
       do i=noa+1, nrorb
          if ( abs(rhoV_a(i,i)).gt.tol ) write( 100, "(2(i4,1x),f14.10,', ')", advance='no' ) i,i, rhoV_a(i,i)
          do j=i+1, nrorb
             if ( abs(rhoV_a(j,i)).gt. tol ) write( 100,"( 2(i4,1x),f14.10,', ')",advance='no' ) i, j, 2.d0*rhoV_a(j,i)
          end do
       end do
       if ( unrestricted ) then
          rhoV_b(:,:) = rho_saveB(:,:,end_time)
          do  i=nob+1, nrorb
             if ( abs(rhoV_b(i,i)).gt.tol ) write( 100, "(2(i4,1x),f14.10,', ')", advance='no' ) i+nrorb,i+nrorb, rhoV_b(i,i)
             do j=i+1, nrorb
                if ( abs(rhoV_b(j,i)).gt. tol ) write( 100,"( 2(i4,1x),f14.10,', ')",advance='no' ) i+nrorb, j+nrorb, 2.d0*rhoV_b(j,i)
             end do
          end do
       end if
       write(100,'(f13.10)') time(end_time)
       
       close(100)

    end if

    
  end subroutine get_rate
  !:--------------------------:!
  !:--------------------------:!
  subroutine get_rate_occ( idir )

    use omp_lib
    implicit none


    integer(8), intent(in) :: idir

    integer(8) :: itime, i, j, a, b, ia, jb, ii, jj, aa, bb
    real(8)    :: norm, dt, const, const2  
    real(8)    :: rhoV_a(noa), rhoV_b(nob), loss_save( noa+nob, start_time-1:end_time+1 )
    real(8)    :: rate_save(noa+nob,start_time-1:end_time+1), trate(start_time-1:end_time+1)
    complex(8) :: psi_det(nstates), psi_ia, psi0
    
    !: write
    character(5)    :: c0
    character(1000) :: ofile
    
    
    trate = 0.d0
    loss_save = 0.d0
    rate_save = 0.d0
    dt = time(2) - time(1) 


    const  = dsqrt(2.d0) 
    const2 = 2.d0
    if ( unrestricted ) then
       const  = 1.d0
       const2 = 1.d0
    end if
    
    

    !$OMP PARALLEL DEFAULT(NONE), &
    !$OMP PRIVATE( itime, i, j, a, b, ia, jb, ii, jj, aa, bb,     &
    !$OMP          norm, psi_det, rhoV_a, rhoV_b, psi_ia, psi0 ), &
    !$OMP SHARED(  start_time, end_time, noa, nob, nstuse, nstates, dt, au2fs, unrestricted, &
    !$OMP          loss_save, rate_save, trate, vabs_b, vabs_a, norm_sq, psi, ci_vec, hole1, part1, &
    !$OMP          const, const2 ) 
    !$OMP DO
    do itime = start_time, end_time
       
       !: save norm
       norm = norm_sq(itime)
       
       !: get psi into determinantal basis
       call get_psi_det( i, nstuse, nstates, ci_vec, psi(:,itime), psi_det )
       psi_det = psi_det / dsqrt(norm)
       
       rhoV_a = 0.d0
       rhoV_b = 0.d0
       
       !: from ground state
       do i=1, noa
          rhoV_a(i) = const2 * vabs_a(i,i)
       end do
       if ( unrestricted ) then
          do i=1, nob
             rhoV_b(i) = vabs_b(i,i)
          end do
       end if

       
       psi0 = psi_det(1)
       
       
       do ia=2, nstates
          
          psi_ia = dconjg( psi_det(ia) )
          
          !: diagonal
          ii = hole1(ia)  ;  i = abs(ii)
          aa = part1(ia)  ;  a = abs(aa)
          
          if ( ii.lt.0 ) then
             rhoV_a(i) = rhoV_a(i) +  2.d0*const*real( psi_ia * psi0 ) * vabs_a(a,i)
          else
             rhoV_b(i) = rhoV_b(i) +  2.d0*const*real( psi_ia * psi0 ) * vabs_b(a,i)
          end if
          

          do jb=2, nstates
             
             jj = hole1(jb)  ;  j = abs(jj)
             bb = part1(jb)  ;  b = abs(bb)
             
             if ( ii.eq.jj ) then
                if ( ii.lt.0 ) then
                   rhoV_a(i) = rhoV_a(i) + real(psi_ia*psi_det(jb)) * vabs_a(b,a)
                else
                   rhoV_b(i) = rhoV_b(i) + real(psi_ia*psi_det(jb)) * vabs_b(b,a)
                end if
             end if
             
             if ( aa.eq.bb ) then
                if ( ii.lt.0 ) then
                   rhoV_a(i) = rhoV_a(i) - real( psi_ia * psi_det(jb) ) * vabs_a(j,i)
                else
                   rhoV_b(i) = rhoV_b(i) - real( psi_ia * psi_det(jb) ) * vabs_b(j,i)
                end if
             end if
             
          end do
          
       end do
       
       !: save
       do i=1, noa
          loss_save(i,itime) = - 2.d0 * rhoV_a(i) * norm * dt / au2fs
          rate_save(i,itime) = - 2.d0 * rhoV_a(i) * norm / au2fs
          trate(itime) = trate(itime) - 2.d0 * rhoV_a(i) * norm / au2fs
       end do
       do i=1, nob
          loss_save(noa+i,itime) = - 2.d0 * rhoV_b(i) * norm * dt / au2fs
          rate_save(noa+i,itime) = - 2.d0 * rhoV_b(i) * norm / au2fs
          trate(itime) = trate(itime) - 2.d0 * rhoV_b(i) * norm / au2fs
       end do

    end do
    !$OMP END DO
    !$OMP END PARALLEL
    
    
    !: get cumulative loss
    do itime=start_time, end_time
       loss_save(:, itime) = loss_save(:, itime-1) + loss_save(:, itime)
    end do

        
90  format( 4(f10.5), 100(f10.5) )
91  format( 4(10x),   100(f10.5) )

    if ( Qprint_pretty .or. Qprint_python ) then
    
       write( c0, '(i0)' ) idir
       ofile = 'LOSS2'//trim(e_d)//trim(c0)//'.out'
       open( unit=100, file=trim(ofile) )
    
       write( 100,'(a10,a10,a10,a10)',advance='no' ) 'time(fs)', 'rate(fs)', 'normloss', 'predloss'
       write( 100,'(100(i10))',advance='no' ) ( i, i=1, noa )
       if( unrestricted ) then
          write( 100,'(100(i10))',advance='no' ) ( i, i=nrorb+1, nrorb+nob )
       end if
       write( 100,'(A)') ''
       
       do itime=start_time, end_time
          if ( unrestricted ) then
             write(100,90) time(itime), &
                  trate(itime), norm_sq(itime)-1.d0, sum( loss_save(:,itime-1) ), &
                  ( loss_save(i,itime-1), i=1, noa ),  &
                  ( loss_save(i,itime-1), i=noa+1, noa+nob )
          else
             write(100,90) time(itime), &
                  trate(itime), norm_sq(itime)-1.d0, sum( loss_save(:,itime-1) ), &
                  ( loss_save(i,itime-1), i=1, noa )
          end if
       end do
       
       
       write(100,'(A)') '# summary_fractional_loss_from_each_occ'
       norm = sum( loss_save(:,end_time) )
       write(100,91) ( loss_save(i,end_time)/norm, i=1, noa+nob )    

       close(100 )



93     format( 4(f15.10), 100(f15.10) )
       !:  write rates
       write( c0, '(i0)' ) idir
       ofile = 'RATE2'//trim(e_d)//trim(c0)//'.out'
       open( unit=100, file=trim(ofile) )
       write( 100,'(a15,a15)',advance='no' ) 'time(fs)', 'trate(fs)'
       write( 100,'(100(i15))',advance='no' ) ( i, i=1, noa )
       if( unrestricted ) then
          write( 100,'(100(i15))',advance='no' ) ( i, i=nrorb+1, nrorb+nob )
       end if
       write( 100,'(A)') ''

       
       do itime=start_time, end_time
          if ( unrestricted ) then
             write(100,93) time(itime), sum( rate_save(:,itime) ), &
                  ( rate_save(i,itime), i=1, noa ),  &
                  ( rate_save(i,itime), i=noa+1, noa+nob )
          else
             write(100,93) time(itime), sum( rate_save(:,itime) ), &
                  ( rate_save(i,itime), i=1, noa )
          end if
       end do

       close(100)
       

    end if

    close(100)


  end subroutine get_rate_occ
  !:--------------------------:!
  !:--------------------------:!
end module rates
