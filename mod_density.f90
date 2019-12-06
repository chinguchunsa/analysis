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

    
89  format( ' time(fs)', f10.5, 2x, f10.5 )
90  format( '    occ_a', 14(1x,f8.5) / 100( 15(1x,f8.5) /) ) !: alpha occupied
91  format( '    vir_a', 14(1x,f8.5) / 100( 15(1x,f8.5) /) ) !: alpha virtual
92  format( '    occ_b', 14(1x,f8.5) / 100( 15(1x,f8.5) /) ) !: beta occupied
93  format( '    vir_b', 14(1x,f8.5) / 100( 15(1x,f8.5) /) ) !: beta virtual
    
    
    if ( Qprint_pretty ) then       

       write( c0, '(i0)' ) idir
       ofile = 'POP'//trim(e_d)//trim(c0)//'.out'
       open( unit=100, file=trim(ofile) )
       
       do itime=start_time, end_time

          write( 100, 89 ) time(itime) , norm_sq(itime)
          write( 100, 90 ) ( rho_save_a(i,i,itime), i=1, noa )
          write( 100, 91 ) ( rho_save_a(i,i,itime), i=noa+1, nrorb )
          if( unrestricted ) then
             write( 100, 92 ) ( rho_save_b(i,i,itime), i=1, nob )
             write( 100, 93 ) ( rho_save_b(i,i,itime), i=nob+1, nrorb )
          end if

       end do

       close(100)

    end if

100 format( f10.5, 1000(f15.10) )

    if ( Qprint_python ) then
       
       write( c0, '(i0)' ) idir
       ofile = 'POP_PY'//trim(e_d)//trim(c0)//'.out'
       open( unit=100, file=trim(ofile) )
       
       write( 100,'(A)',advance='no' ) 'time(fs) '
       write( 100,'(100(i0,1x))', advance='no' ) ( i, i=1, nrorb )
       if ( unrestricted ) write( 100, '(100(i0,1x))', advance='no' ) ( i, i=nrorb+1, norb )
       write( 100, '(A)' ) ''

       do itime=start_time, end_time
          write( 100,100, advance='no' ) time(itime), ( rho_save_a(i,i,itime), i=1, nrorb )
          if ( unrestricted ) write(100,100) ( rho_save_b(i,i,itime), i=1, nrorb )
          write(100,'(A)') ''
       end do
       
       
    end if
       

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
       
    end if
       

  end subroutine get_density
  !:--------------------------:!
  !:--------------------------:!
  subroutine get_hp_density( idir )


    use omp_lib
    implicit none

    
    integer(8), intent(in) :: idir
    
    real(8), parameter :: tol2 = 1.d-10 
    
    integer(8) :: itime, ia, jb, i, a, j, b, ii, jj, aa, bb
    real(8) :: norm, const
    real(8) :: rho_a(nrorb,nrorb), rho_b(nrorb,nrorb)
    real(8) :: rho_save_a(nrorb,nrorb,start_time-1:end_time+1), rho_save_b(nrorb,nrorb,start_time-1:end_time+1)
    complex(8) :: psi_det(nstates), psi_det00(nstates), psi_ia, psi_0
    
    !: write
    character(5) :: c0
    character(1000) :: ofile


    const = dsqrt(2.d0)
    if ( unrestricted ) const = 1.d0
    

    rho_save_a = 0.d0
    rho_save_b = 0.d0


    !$OMP PARALLEL DEFAULT(SHARED), &
    !$OMP PRIVATE( itime, i, j, a, b, ia, jb, ii, jj, aa, bb,              &
    !$OMP          norm, psi_det, psi_det00, rho_a, rho_b, psi_ia, psi_0 ),&
    !$OMP SHARED( psi, rho_save_a, rho_save_b, norm_sq, psi, ci_vec )
    !$OMP DO
    do itime=start_time+1, end_time
       
       rho_a = 0.d0
       rho_b = 0.d0

       !: get psi_det at itime - 1 
       norm = norm_sq(itime-1)
       call get_psi_det ( i, nstuse, nstates, ci_vec, psi(:,itime-1), psi_det00 )
       psi_det00 = psi_det00 / dsqrt(norm)
       
       !: get psi_det at itime
       norm = norm_sq(itime)
       call get_psi_det( i, nstuse, nstates, ci_vec, psi(:,itime), psi_det )
       psi_det = psi_det / dsqrt(norm)

       
       psi_0 = psi_det00(1)

       
       do ia=2, nstates
          
          ii = hole1(ia)  ;  i = abs(ii)
          aa = part1(ia)  ;  a = abs(aa)
          
          psi_ia = psi_det(ia)
          
          if( abs(psi_ia).gt.tol2 ) then
          
             !: ground -> excited
             if ( ii.lt.0 ) then
                rho_a(a,i) = rho_a(a,i) + const * real( dconjg(psi_ia) * psi_0 )
             else
                rho_b(a,i) = rho_b(a,i) + const * real( dconjg(psi_ia) * psi_0 )
             end if
             
             !: excited -> ground
             if ( ii.lt.0 ) then
                rho_a(i,a) = rho_a(i,a) + const * real( dconjg(psi_0) * psi_ia )
             else
                rho_b(i,a) = rho_b(i,a) + const * real( dconjg(psi_0) * psi_ia )
             end if

             psi_ia = dconjg(psi_ia)
             
             do jb=2, nstates
                
                jj = hole1(jb)  ;  j = abs(jj)
                bb = part1(jb)  ;  b = abs(bb)
                
                !: a --> b excitation or dexcitation
                if ( jj.eq.ii .and. aa.ne.bb ) then
                   if ( ii.lt.0 ) then
                      rho_a(a,b) = rho_a(a,b) + real( psi_ia*psi_det00(jb) )
                   else
                      rho_b(a,b) = rho_b(a,b) + real( psi_ia*psi_det00(jb) )
                   end if
                end if

                !: i --> j hole transfer 
                if ( aa.eq.bb .and. jj.ne.ii ) then
                   if ( aa.lt.0 ) then
                      rho_a(j,i) = rho_a(j,i) + real( psi_ia*psi_det00(jb) )
                   else
                      rho_b(j,i) = rho_b(j,i) + real( psi_ia*psi_det00(jb) )
                   end if
                end if
                
             end do
          end if
       end do

       !: save 
       rho_save_a(:,:,itime) = rho_a(:,:)
       if( unrestricted ) rho_save_b(:,:,itime) = rho_b(:,:)       
       
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    
    !: HOLE
    write( c0, '(i0)' ) idir 
    ofile = 'HOLE_PY'//trim(e_d)//trim(c0)//'.out'
    open( unit=100,file=trim(ofile) )
    
    !: only prints python 
    write( 100, '(A)', advance='no' ) 'time(fs) '
    write( 100, '(100(i0,1x))', advance='no' ) ( i, i=1, nrorb )
    if ( unrestricted ) write( 100, '(100(i0,1x))', advance='no' ) ( i, i=nrorb+1, norb )
    write( 100, '(A)' ) ''
    
    do itime=start_time, end_time

       !: compute hole
       rho_a(:,:) = rho_save_a(:,:,itime)
       rho_a(:,:) = matmul( transpose(rho_a), rhO_a )
       
       if( unrestricted ) then
          rho_b = rho_save_b(:,:,itime) 
          rho_b = matmul( transpose(rho_b), rho_b )
       end if
       
       !: alpha
       write( 100, '(1000(f14.10,1x))', advance='no' ) time(itime), ( rho_a(i,i), i=1, nrorb )
       if ( unrestricted ) write( 100, '(1000(f14.10,1x))', advance='no' ) ( rho_b(i,i), i=1, nrorb )
       write( 100, '(A)' ) ''
       
    end do
    
    close(100)
    
    
    !PARTICLE
    write( c0, '(i0)' ) idir 
    ofile = 'PART_PY'//trim(e_d)//trim(c0)//'.out'
    open( unit=100, file=trim(ofile) )
    
    !: only prints python 
    write( 100, '(A)', advance='no' ) 'time(fs) '
    write( 100, '(100(i0,1x))', advance='no' ) ( i, i=1, nrorb )
    if ( unrestricted ) write( 100, '(100(i0,1x))', advance='no' ) ( i, i=nrorb+1, norb ) 
    write( 100, '(A)' ) '' 
    
    do itime=start_time, end_time

       !: compute particle 
       rho_a(:,:) = rho_save_a(:,:,itime)
       rho_a(:,:) = matmul( rho_a, transpose(rho_a) )
       
       if( unrestricted ) then
          rho_b = rho_save_b(:,:,itime)
          rho_b = matmul( rho_b, transpose(rho_b) )
       end if
       
       !: alpha
       write( 100, '(1000(f14.10,1x))', advance='no' ) time(itime), ( rho_a(i,i), i=1, nrorb ) 
       if ( unrestricted ) write( 100, '(1000(f14.10,1x))', advance='no' ) time(itime), ( rho_b(i,i), i=1, nrorb )
       write( 100, '(A)' ) ''
       
    end do
    
    close(100)

    
    
  end subroutine get_hp_density
  !:--------------------------:!
  !:--------------------------:!   
end module density
