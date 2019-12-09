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
          rate_save(i,itime) = - 2.d0 * rhoV_a(i) / au2fs          
          trate(itime) = trate(itime) - 2.d0 * rhoV_a(i) / au2fs
       end do
       do i=1, nob
          loss_save(noa+i,itime) = - 2.d0 * rhoV_b(i) * norm * dt / au2fs
          rate_save(noa+i,itime) = - 2.d0 * rhoV_b(i) / au2fs
          trate(itime) = trate(itime) - 2.d0 * rhoV_b(i) / au2fs
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
  subroutine get_rate_ia( idir )

    !: always for unrestricted

    use omp_lib
    implicit none


    integer(8), intent(in) :: idir

    integer(8) :: itime, i, j, a, b, ia, jb, ii, jj, aa, bb
    real(8)    :: norm, dt, const, v00 
    real(8)    :: rhoV_0, rhoV_ia(nstates) 
    real(8)    :: rate_save(nstates,start_time-1:end_time+1), trate(start_time-1:end_time+1)
    complex(8) :: psi_det(nstates), psi_ia, psi0
    
    !: write
    character(5)    :: c0
    character(1000) :: ofile
    
    
    trate = 0.d0
    rate_save = 0.d0
    dt = time(2) - time(1) 
    
    v00 = 0.d0
    do i=1, noa
       v00 = v00 + vabs_a(i,i) 
    end do
    do i=1, nob
       v00 = v00 + vabs_b(i,i)
    end do
    

    !$OMP PARALLEL DEFAULT(NONE), &
    !$OMP PRIVATE( itime, i, j, a, b, ia, jb, ii, jj, aa, bb,     &
    !$OMP          norm, psi_det, rhoV_0, rhoV_ia, psi_ia, psi0, const ), &
    !$OMP SHARED(  start_time, end_time, noa, nob, nva, nvb, nstuse, nstates, dt, au2fs, unrestricted, &
    !$OMP          trate, rate_save, vabs_b, vabs_a, norm_sq, psi, ci_vec, v00 ) 
    !$OMP DO
    do itime = start_time, end_time
       
       !: save norm
       norm = norm_sq(itime)
       
       !: get psi into determinantal basis
       call get_psi_det( i, nstuse, nstates, ci_vec, psi(:,itime), psi_det )
       psi_det = psi_det / dsqrt(norm)
       
       rhoV_0  = 0.d0
       rhoV_ia = 0.d0
       
       psi0   = dconjg(psi_det(1))

       !: from ground state
       rhoV_0 = dconjg( psi0 )*psi0 * v00
       
       !: ground and excited state coherence
       ia = 1 
       !: alpha
       do i=1, noa
          do a=1, nva
             ia = ia + 1 
             const  = real(psi0*psi_det(ia)) * vabs_a(noa+a,i)
             rhoV_0 = rhoV_0 + const 
             rhoV_ia(ia) = const
          end do
       end do
       !: beta
       do i=1, nob
          do a=1, nvb
             ia = ia + 1 
             const  = real(psi0*psi_det(ia)) * vabs_b(nob+a,i)
             rhoV_0 = rhoV_0 + const
             rhoV_ia(ia) = const
          end do
       end do

       
       !: excited states (hole transfer states from beta to beta)
       !: alpha
       !: due to coherence with excited states
       do i=1, noa
          do a=1, nva
             ia = 1 + (i-1)*nva + a  
             jb = 1 
             const  = rhoV_ia(ia)
             psi_ia = dconjg(psi_det(ia))
             do j=1, noa
                do b=1, nva
                   jb = jb + 1 
                   if ( i.eq.j ) const = const + real(psi_ia*psi_det(jb)) * vabs_a(noa+b,noa+a)
                   if ( a.eq.b ) const = const - real(psi_ia*psi_det(jb)) * vabs_a(j,i)
                end do
             end do
             !: diagonal terms
             rhoV_ia(ia) = const + dconjg(psi_ia)*psi_ia * v00
          end do
       end do


       !: excited states (hole transfer states from beta to beta)
       !: beta
       !: due to coherence with excited states
       do i=1, nob
          do a=1, nvb
             ia = 1 + noa*nva + (i-1)*nvb + a 
             jb = 1 + noa*nva
             const  = rhoV_ia(ia)
             psi_ia = dconjg(psi_det(ia))
             do j=1, nob
                do b=1, nvb
                   jb = jb + 1 
                   if ( i.eq.j ) const = const + real(psi_ia*psi_det(jb)) * vabs_b(nob+b,nob+a)
                   if ( a.eq.b ) const = const - real(psi_ia*psi_det(jb)) * vabs_b(j,i)
                end do
             end do
             !: diagonal terms
             rhoV_ia(ia) = const + dconjg(psi_ia)*psi_ia * v00
          end do
       end do
       
       
       !: save
       trate(itime) = 2.d0*sum(rhoV_ia(:))/au2fs + 2.d0*rhoV_0/au2fs
       rate_save(:,itime) = 2.d0 * rhoV_ia(:) / au2fs
       rate_save(1,itime) = 2.d0 * rhoV_0 / au2fs
       
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    
    
90  format( f10.5, 100(f15.10) )

    if ( Qprint_pretty .or. Qprint_python ) then    
       
       write( c0, '(i0)' ) idir
       ofile = 'RATE_ia'//trim(e_d)//trim(c0)//'.out'
       open( unit=100, file=trim(ofile) )
       
       write( 100,'(a10,a10,a10,a10)',advance='no' ) 'time(fs)', 'rate(fs)' 
       write( 100,'(100(i10))',advance='no' ) ( i, i=0, nob )
       write( 100,'(A)') ''
       
       do itime=start_time, end_time
          !: offset
          const = rate_save(1,itime)
          do i=1, noa
             const = const + rate_save( 1+(i-1)*nva+1, itime ) + rate_save( 1+(i-1)*nva+2, itime )
          end do
          do i=1, nob
             const = const + rate_save( 1+noa*nva+(i-1)*nvb+1, itime ) + rate_save( 1+noa*nva+(i-1)*nvb+2, itime )
          end do
          write(100,90) time(itime), trate(itime), rate_save(1,itime), &
               ( rate_save( 1+noa*nva+(i-1)*nvb+1, itime ), i=1, nob ), const
       end do
       
       close(100)

    end if


  end subroutine get_rate_ia
  !:--------------------------:!
  !:--------------------------:!
  subroutine get_rate_ground( idir )

    !: always for unrestricted

    use omp_lib
    implicit none


    integer(8), intent(in) :: idir
    
    integer(8) :: i, j, itime

    real(8) :: ratex, ratey, ratez, tmp, field0
    real(8) :: tratex(10), tratey(10), tratez(10)
    real(8) :: ionization_rate(nstuse)
    
    complex(8) :: psi0
    complex(8) :: psi_tmp(nstuse) 
        
    !: write
    character(5)    :: c0
    character(1000) :: ofile, ofile2


    !: get excitation rate and ionization rate
    itime      = 243
    psi_tmp(:) = psi(:,itime)  / dsqrt( norm_sq(itime) )
    field0     = efield(itime)    

    tratex = 0.d0
    tratey = 0.d0
    tratez = 0.d0
    ionization_rate = 0.d0
       
    
    !: ionization rate
    do i=1, nstuse
       psi0 = dconjg( psi_tmp(i) )
       tmp = 0.d0
       do j=1, nstuse
          tmp = tmp + 2.d0 * real( psi0*psi_tmp(j) ) * ci_vabs(j,i) / au2fs
       end do
       ionization_rate(i) = tmp
    end do
    
    !: excitation rate
    do i=1, 10
       psi0 = dconjg( psi_tmp(i) )
       ratex = 0.d0
       ratey = 0.d0
       ratez = 0.d0
       do j=1, nstuse
          tmp = 2.d0 * field0 * aimag( psi0* psi_tmp(j) ) / au2fs
          ratex = ratex + tmp * tdx(j,i) 
          ratey = ratey + tmp * tdy(j,i)
          ratez = ratez + tmp * tdz(j,i)             
       end do
       tratex(i) = ratex
       tratey(i) = ratey
       tratez(i) = ratez
    end do
    
    
    !: X
    write( c0, '(i0)' ) idir 
    ofile = 'EXCTATION_RATE_X'//trim(e_d)//trim(c0)//'.out'
    open( unit=100,file=trim(ofile) )
    write( 100, '(a10,1x,a15,10(1x,a15))' ) 'state', 'ion_rate', &
         'exrate1', 'exrate2','exrate3','exrate4','exrate5',    &
         'exrate6', 'exrate7','exrate8','exrate9','exrate10'
    do i=1, nstuse
       psi0 = dconjg( psi_tmp(i) )
       write(100,'(i10,1x,f15.10)',advance='no') i, ionization_rate(i) !/sum(ionization_rate)
       do j=1, 10
          write(100, '(12(1x,f15.10))', advance='no' ) 2.d0*field0/au2fs * aimag(psi0*psi_tmp(j))*tdx(j,i) !/tratex(j)
       end do
       write(100,'(A)') ''
    end do
    write(100,"('#'9x,12(1x,f15.10))") sum(ionization_rate), ( tratex(j),j=1,10 )
    close(100)
    

    !: Y
    write( c0, '(i0)' ) idir 
    ofile = 'EXCTATION_RATE_Y'//trim(e_d)//trim(c0)//'.out'
    open( unit=100,file=trim(ofile) )
    write( 100, '(a10,1x,a15,10(1x,a15))' ) 'state', 'ion_rate', &
         'exrate1', 'exrate2','exrate3','exrate4','exrate5',    &
         'exrate6', 'exrate7','exrate8','exrate9','exrate10'
    do i=1, nstuse
       psi0 = dconjg( psi_tmp(i) )
       write(100,'(i10,1x,f15.10)',advance='no') i, ionization_rate(i) !/sum(ionization_rate)
       do j=1, 10
          write(100, '(12(1x,f15.10))', advance='no' ) 2.d0*field0/au2fs * aimag(psi0*psi_tmp(j))*tdy(j,i) !/tratey(j)
       end do
       write(100,'(A)') ''
    end do
    write(100,"('#',9x,12(1x,f15.10))") sum(ionization_rate), ( tratey(j),j=1,10 )
    close(100)
    

    !: Z
    write( c0, '(i0)' ) idir 
    ofile = 'EXCTATION_RATE_Z'//trim(e_d)//trim(c0)//'.out'
    open( unit=100,file=trim(ofile) )
    write( 100, '(a10,1x,a15,10(1x,a15))' ) 'state', 'ion_rate', &
         'exrate1', 'exrate2','exrate3','exrate4','exrate5',    &
         'exrate6', 'exrate7','exrate8','exrate9','exrate10'
    do i=1, nstuse
       psi0 = dconjg( psi_tmp(i) )
       write(100,'(i10,1x,f15.10)',advance='no') i, ionization_rate(i) !/sum(ionization_rate)
       do j=1, 10
          write(100, '(12(1x,f15.10))', advance='no' ) 2.d0*field0/au2fs * aimag(psi0*psi_tmp(j))*tdz(j,i) !/tratez(j)
       end do
       write(100,'(A)') ''
    end do
    write(100,"('#',9x,12(1x,f15.10))") sum(ionization_rate), ( tratez(j),j=1,10 )
    close(100)
    
    
    
    
    

  end subroutine get_rate_ground
  !:--------------------------:!
  !:--------------------------:!
  subroutine get_rate_ip( idir )

    use omp_lib
    use utils
    implicit none

    integer(8), intent(in) :: idir

    integer(8) :: itime, i, j, a, b, ia, jb, ii, jj, x,y, aa, bb, xx, yy
    real(8) :: norm, rdum
    real(8) :: rate_ab( noa, nob, start_time:end_time+1 ), rate_bb( nob, nob, start_time:end_time+1 ), trate(start_time:end_time+1)
    real(8) :: rho_a(nrorb,nrorb), rho_b(nrorb,nrorb)
    complex(8) :: psi_det(nstates), psi_ia, psi_jb
    
    !: write stuff
    character(5)    :: c0, ci
    character(1000) :: ofile
    

    !: initialize
    trate = 0.d0
    rate_ab = 0.d0
    rate_bb = 0.d0           

    itime : do itime=start_time, end_time
       
       !: save norm
       norm = norm_sq(itime)       

       !: get psi into determinantal basis
       call get_psi_det( i, nstuse, nstates, ci_vec, psi(:,itime), psi_det )
       psi_det = psi_det / dsqrt(norm)

       !: only take into account doubles
       do ia=nob+1, nstates
          
          psi_ia = dconjg(psi_det(ia))
          xx = hole(ia,1)
          ii = hole(ia,2)
          aa = part(ia,1)

          do jb=nob+1, nstates
             
             yy = hole(jb,1)
             jj = hole(jb,2)
             bb = part(jb,1)

             if ( xx.eq.yy .and. ii.eq.jj ) then
                if ( ii.lt.0 ) rate_ab(-ii,xx,itime) = rate_ab(-ii,xx,itime) + real(psi_ia*psi_det(jb)) * vabs_a(-bb,-aa)
                if ( ii.gt.0 ) then
                   rate_bb(ii,xx,itime) = rate_bb(ii,xx,itime) + real(psi_ia*psi_det(jb)) * vabs_b(bb,aa)
                   rate_bb(xx,ii,itime) = rate_bb(xx,ii,itime) + real(psi_ia*psi_det(jb)) * vabs_b(bb,aa)
                end if
             end if

          end do
       end do
             
       call get_1density_ip(psi_det,nstates,hole,part,rho_a,rho_b,nrorb,noa,nob)
       rho_a = matmul( rho_a, vabs_a )
       rho_b = matmul( rho_b, vabs_b )
       rdum = 0.d0
       do i=1, nrorb
          rdum = rdum + rho_a(i,i) + rho_b(i,i) 
       end do
       trate(itime) = - 2.d0 * rdum / au2fs
       
    end do itime


    rate_ab = rate_ab * 2.d0 / au2fs
    rate_bb = rate_bb * 2.d0 / au2fs


42  format( ' doublesAB ', i5'B', 10(1x,f8.5) )
43  format( ' doublesBB ', i5'B', 10(1x,f8.5) )       
    if ( Qprint_pretty ) then
       write( c0, '(i0)' ) idir
       ofile = 'TEST_RATE'//trim(e_d)//trim(c0)//'.out'
       open( unit=100, file=trim(ofile) )
       do itime=start_time, end_time
          write(100,"('time(fs) ',f10.5,' rate(fs-1) ',f15.10, 'unnormrate ', f15.10 )") time(itime), trate(itime), trate(itime)*norm_sq(itime)
          do x=1, nob
             write(100,42 ) x, ( rate_ab(i,x,itime), i=1, noa )
          end do
          do x=1, nob
             write(100,43 ) x, ( rate_bb(i,x,itime), i=1, nob )
          end do
       end do
       close(100)
    end if
    
    
    !: print python
    ofile='TEST_RATE_PY'//trim(e_d)//trim(c0)//'.out'
    open( unit=100, file=trim(ofile) )
    write(100,'(A)', advance='no') " time(fs) norm2 trate "
    do x=1, nob
       write( c0, '(i0)' ) x
       do i=1, noa
          write( ci, '(i0)' ) -i
          write(100,'(A)',advance='no') '('//trim(c0)//','//trim(ci)//') '
       end do
    end do
    do x=1, nob
       write( c0, '(i0)' ) x
       do i=1, nob
          write( ci, '(i0)' ) i
          write(100,'(A)',advance='no') '('//trim(c0)//','//trim(ci)//') '
       end do
    end do
    write(100,'(A)') ''


    do itime=start_time, end_time
       write(100,"(f10.7,1x)",advance='no') time(itime), norm_sq(itime), trate(itime)
       do x=1, nob
          write(100,"(f10.7,1x)",advance='no' ) ( rate_ab(i,x,itime), i=1, noa )
       end do
       do x=1, nob
          write(100,"(f10.7,1x)",advance='no' ) ( rate_bb(i,x,itime), i=1, nob )
       end do
       write(100,'(A)') ''
    end do
    close(100)
    


  end subroutine get_rate_ip
  !:--------------------------:!
  !:--------------------------:!
  subroutine get_rateAB_ip( idir )

    use omp_lib
    use utils
    implicit none

    integer(8), intent(in) :: idir

    integer(8) :: itime, i, j, a, b, ia, jb, ii, jj, x,y, aa, bb, xx, yy
    real(8) :: norm, rdum, vsum, dt
    real(8) :: rate_aa( noa, nob, start_time:end_time+1 ), rate_bb( nob, nob, start_time:end_time+1 ), trate(start_time:end_time+1)
    real(8) :: tmp_aa(noa,nob), tmp_bb(nob,nob)
    complex(8) :: psi_det(nstates), psi_ia, psi_jb
    
    !: write stuff
    character(5)    :: c0, ci
    character(1000) :: ofile, pyfile
    
    
    !: initialize
    trate = 0.d0
    rate_aa = 0.d0
    rate_bb = 0.d0           
    
    vsum = 0.d0
    do i=1, noa
       vsum = vsum + vabs_a(i,i)
    end do
    do i=1, nob
       vsum = vsum + vabs_b(i,i)
    end do
    
    dt = time(start_time+1) - time(start_time)

    itime : do itime=start_time, end_time
       
       tmp_aa = 0.d0
       tmp_bb = 0.d0
       
       !: save norm
       norm = norm_sq(itime)   

       !: get psi into determinantal basis
       call get_psi_det( i, nstuse, nstates, ci_vec, psi(:,itime), psi_det )
       psi_det = psi_det / dsqrt(norm)

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
                tmp_bb(yy,xx) = tmp_bb(yy,xx) - real(psi_ia*psi_det(jb))*vabs_a(yy,xx)
                go to 78
             end if SS
             
             !: <iX->a|A|X>  or  <IX->A|A|X>
             SD1: if ( YY.eq.XX .and. jj.eq.0 ) then
                if ( ii.lt.0 ) tmp_aa(-ii,xx) = tmp_aa(-ii,xx) + real(psi_ia*psi_det(jb))*vabs_a(-aa,-ii)
                if ( ii.gt.0 ) tmp_bb(ii,xx)  = tmp_bb(ii,xx)  + real(psi_ia*psi_det(jb))*vabs_b(aa,ii)
                go to 78
             end if SD1
             !: <X|A|jX->b>
             SD2: if ( YY.eq.XX .and. ii.eq.0 ) then
                if ( jj.lt.0 ) tmp_aa(-jj,yy) = tmp_aa(-jj,yy) + real(psi_ia*psi_det(jb))*vabs_a(-jj,-bb)
                if ( jj.gt.0 ) tmp_bb(jj,yy)  = tmp_bb(jj,yy)  + real(psi_ia*psi_det(jb))*vabs_b(jj,bb)
                go to 78
             end if SD2
             !: <YX->A|A|Y>
             SD3 : if ( YY.eq.II .and. jj.eq.0 ) then
                tmp_bb(ii,xx) = tmp_bb(ii,xx) - real(psi_ia*psi_det(jb))*vabs_b(aa,xx)
                go to 78
             end if SD3
             !: <X|A|XY->B>
             SD4: if ( XX.eq.JJ .and. ii.eq.0 ) then
                tmp_bb(jj,yy) = tmp_bb(jj,yy) - real(psi_ia*psi_det(jb))*vabs_b(yy,bb)
                go to 78
             end if SD4
             
             !: doubles
             
             !: <ix->a|A|jx->a>
             kdelta_xy_ab : if ( xx.eq.yy .and. aa.eq.bb ) then
                if ( ii.lt.0 ) tmp_aa(-ii,xx) = tmp_aa(-ii,xx) - real(psi_ia*psi_det(jb))*vabs_a(-jj,-ii)
                if ( ii.gt.0 ) tmp_bb(ii,xx)  = tmp_bb(ii,xx)  - real(psi_ia*psi_det(jb))*vabs_b(jj,ii)
             end if kdelta_xy_ab
             
             !: <ix->a|A|ix->b>
             kdelta_xy_ij : if ( xx.eq.yy .and. ii.eq.jj ) then
                if ( aa.lt.0 ) tmp_aa(-ii,xx) = tmp_aa(-ii,xx) + real(psi_ia*psi_det(jb))*vabs_a(-aa,-bb)
                if ( aa.gt.0 ) tmp_bb(ii,xx)  = tmp_bb(ii,xx)  + real(psi_ia*psi_det(jb))*vabs_b(aa,bb)
             end if kdelta_xy_ij
             
             !: <ix->a|A|iy->a>
             kdelta_ij_ab : if ( ii.eq.jj .and. aa.eq.bb ) then
                if ( aa.lt.0 ) tmp_aa(-ii,xx) = tmp_aa(-ii,xx) - real(psi_ia*psi_det(jb))*vabs_b(yy,xx)
                if ( aa.gt.0 ) tmp_bb(ii,xx)  = tmp_bb(ii,xx)  - real(psi_ia*psi_det(jb))*vabs_b(yy,xx)
             end if kdelta_ij_ab
             
             !: <ix->a|A|xy->a>
             kdelta_jx_ab : if ( JJ.eq.XX .and. AA.eq.BB ) then
                tmp_bb(ii,xx) = tmp_bb(ii,xx) + real(psi_ia*psi_det(jb))*vabs_b(yy,ii)
             end if kdelta_jx_ab
             
             !: <yx->a|A|jy->a>
             kdelta_yi_ab : if (YY.eq.II .and. AA.eq.BB ) then
                tmp_bb(ii,xx) = tmp_bb(ii,xx) + real(psi_ia*psi_det(jb))*vabs_b(jj,xx)
             end if kdelta_yi_ab
             
             !: <yx->a|A|xy->b>
             kdelta_iy_xj : if ( II.eq.YY .and. XX.eq.JJ ) then
                tmp_bb(ii,xx) = tmp_bb(ii,xx) - real(psi_ia*psi_det(jb))*vabs_b(aa,bb)
             end if kdelta_iy_xj
             
78           continue
          
          end do

          !: get diagonal
          if ( ii.eq.0 ) tmp_bb(ii,ii)  = tmp_bb(ii,ii)  + dconjg(psi_ia)*psi_ia * vsum
          if ( ii.lt.0 ) tmp_aa(-ii,xx) = tmp_aa(-ii,xx) + dconjg(psi_ia)*psi_ia * vsum
          if ( ii.gt.0 ) tmp_bb(ii,xx)  = tmp_bb(ii,xx)  + dconjg(psi_ia)*psi_ia * vsum                       
          
       end do

       trate(itime) = ( sum( tmp_aa ) + sum( tmp_bb ) ) * 2.d0 / au2fs
       rate_aa(:,:,itime) = tmp_aa(:,:)
       rate_bb(:,:,itime) = tmp_bb(:,:)
       
    end do itime

    rate_aa = rate_aa * 2.d0 / au2fs
    rate_bb = rate_bb * 2.d0 / au2fs    
        
42  format( ' from_alpha ', i5'B', 10(1x,f8.5) )
43  format( ' from_beta  ', i5'B', 10(1x,f8.5) )       
    if ( Qprint_pretty ) then
       write( c0, '(i0)' ) idir
       ofile = 'TEST_RATE'//trim(e_d)//trim(c0)//'.out'
       open( unit=100, file=trim(ofile) )
       do itime=start_time, end_time
          write(100,"('time(fs) ',f10.5,' rate(fs-1) ',f15.10, ' comp_rate ', f15.10 )") time(itime), trate(itime), trate(itime)*norm_sq(itime)
          do x=1, nob
             write(100,42 ) x, ( rate_aa(i,x,itime), i=1, noa )
          end do
          do x=1, nob
             write(100,43 ) x, ( rate_bb(i,x,itime), i=1, nob )
          end do
       end do
       close(100)
    end if
    
    
    !: print python
    ofile='TEST_RATE_PY'//trim(e_d)//trim(c0)//'.out'
    pyfile='TEST_LOSS2_PY'//trim(e_d)//trim(c0)//'.out'
    open( unit=100, file=trim(ofile) )
    open( unit=200, file=trim(pyfile) )
    write(100,'(A)', advance='no') " time(fs) norm2 trate "
    write(200,'(A)', advance='no') " time(fs) norm2 pred_norm2 "
    do x=1, nob
       write( c0, '(i0)' ) x
       do i=1, noa
          write( ci, '(i0)' ) -i
          write(100,'(A)',advance='no') '('//trim(c0)//','//trim(ci)//') '
          write(200,'(A)',advance='no') '('//trim(c0)//','//trim(ci)//') '
       end do
    end do
    do x=1, nob
       write( c0, '(i0)' ) x
       do i=x+1, nob
          write( ci, '(i0)' ) i
          write(100,'(A)',advance='no') '('//trim(c0)//','//trim(ci)//') '
          write(200,'(A)',advance='no') '('//trim(c0)//','//trim(ci)//') '
       end do
    end do
    write(100,'(A)') ''
    write(200,'(A)') ''
    
    do itime=start_time, end_time
       write(100,"(f10.7,1x)",advance='no') time(itime), norm_sq(itime), trate(itime)
       do x=1, nob
          write(100,"(f10.7,1x)",advance='no' ) ( rate_aa(i,x,itime), i=1, noa )
       end do
       do x=1, nob
          write(100,"(f10.7,1x)",advance='no' ) ( rate_bb(i,x,itime), i=x+1, nob )
       end do
       write(100,'(A)') ''
    end do
    close(100)
    
    !: loss
    tmp_aa = 0.d0 
    tmp_bb = 0.d0
    vsum = 1.d0
    write( c0, '(i0)' ) idir
    ofile='TEST_LOSS2'//trim(e_d)//trim(c0)//'.out'
    open( unit=100,file=trim(ofile) )
    do itime=start_time, end_time
       vsum = vsum - trate(itime) * norm_sq(itime) * dt
       write(100,"('time(fs) ',f10.5,' norm2 ',f15.10, ' pred_norm2 ', f15.10 )") time(itime), norm_sq(itime), vsum
       write(200,"(f15.10,1x)",advance='no') time(itime), norm_sq(itime), vsum
       !: compute
       do x=1, nob
          do i=1, noa
             tmp_aa(i,x) = tmp_aa(i,x) + rate_aa(i,x,itime) * dt * norm_sq(itime)
          end do
          do i=1, nob
             tmp_bb(i,x) = tmp_bb(i,x) + rate_bb(i,x,itime) * dt * norm_sq(itime)
          end do
       end do
       do x=1, nob
          write(100,42 ) x, ( tmp_aa(i,x), i=1, noa )
          write(200,"(f10.7,1x)",advance='no' ) ( tmp_aa(i,x), i=1, noa )
       end do
       do x=1, nob
          write(100,43 ) x, ( tmp_bb(i,x), i=1, nob )
          write(200,"(f10.7,1x)",advance='no' ) ( tmp_bb(i,x), i=x+1, noa )
       end do
       write(200,'(A)') ''
    end do
    close(100)
    close(200)
    

  end subroutine get_rateAB_ip
  !:--------------------------:!
  !:--------------------------:!
end module rates
