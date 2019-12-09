module rates

  use global_variables
  use utils

contains
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
          
          if( ii.lt.0 ) rhoV_a(i) = rhoV_a(i) + 2.d0*const*real( psi_ia * psi0 ) * vabs_a(a,i)
          if( ii.gt.0 ) rhoV_b(i) = rhoV_b(i) +  2.d0*const*real( psi_ia * psi0 ) * vabs_b(a,i)
          
          do jb=2, nstates             
             jj = hole1(jb)  ;  j = abs(jj)
             bb = part1(jb)  ;  b = abs(bb)             
             if ( ii.eq.jj ) then
                if ( ii.lt.0 ) rhoV_a(i) = rhoV_a(i) + real(psi_ia*psi_det(jb)) * vabs_a(b,a)
                if ( ii.gt.0 ) rhoV_b(i) = rhoV_b(i) + real(psi_ia*psi_det(jb)) * vabs_b(b,a)
             end if
             if ( aa.eq.bb ) then
                if ( ii.lt.0 ) rhoV_a(i) = rhoV_a(i) - real( psi_ia * psi_det(jb) ) * vabs_a(j,i)
                if ( ii.gt.0 ) rhoV_b(i) = rhoV_b(i) - real( psi_ia * psi_det(jb) ) * vabs_b(j,i)
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
93  format( 4(f15.10), 100(f15.10) )
    !:  write rates
    write( c0, '(i0)' ) idir
    ofile = 'RATE'//trim(e_d)//trim(c0)//'.out'
    open( unit=100, file=trim(ofile) )
    !: write heading
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
    
    write(iout,'(A)') ' --> generated file '//trim(ofile)
    flush(iout)
    
    
  end subroutine get_rate_occ
  !:--------------------------:!
  !:--------------------------:!
  subroutine get_rateAB_ip( idir )

    use omp_lib
    use utils
    implicit none

    integer(8), intent(in) :: idir
    
    integer(8) :: itime, i, j, a, b, x, y, ia, jb
    real(8) :: norm, rdum, vsum 
    real(8) :: rate_aa( noa, nob, start_time:end_time+1 ), rate_bb( nob, nob, start_time:end_time+1 ), rate_x(nob, start_time:end_time+1)
    real(8) :: trate(start_time:end_time+1)
    real(8) :: tmp_aa(noa,nob), tmp_bb(nob,nob), vdum
    complex(8) :: psi_det(nstates), psi_ia, psi_x, psi_i, psi_y
    
    !: write stuff
    character(5)    :: c0, ci
    character(1000) :: ofile 
    
    
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
    
    !$OMP PARALLEL DEFAULT(NONE), &
    !$OMP PRIVATE( itime, tmp_aa, tmp_bb, norm, psi_det, rdum, vdum, psi_x, psi_i, psi_ia, i, x, y, j, b, ia, jb ), &
    !$OMP SHARED( start_time, end_time, vsum, norm_sq, psi, ci_vec, noa, nob, nva, nvb, nstuse, nstates, &
    !$OMP         vabs_a, vabs_b, hole_indices, trate, rate_aa, rate_bb )
    !$OMP DO
    itime : do itime=start_time, end_time
       
       tmp_aa  = 0.d0
       tmp_bb  = 0.d0
       
       !: save norm
       norm = norm_sq(itime)   
       
       !: get psi into determinantal basis
       call get_psi_det( i, nstuse, nstates, ci_vec, psi(:,itime), psi_det )
       psi_det = psi_det / dsqrt(norm)
       

       !: doubles alpha beta
       do x=1, nob
          alpha : do i=1, noa
                          
             !: diagonal
             rdum = 0.d0
             ia   = hole_indices(-i,x) - 1 
             do a=1, nva
                rdum = rdum + real(dconjg(psi_det(ia+a))*psi_det(ia+a))
             end do
             rdum = rdum*vsum 
             
             !: contributions from c_x and c_i
             rdum = rdum + dconjg(psi_det(x))*psi_det(x)*vabs_a(i,i)
             

             !: sum_a : 2c_ixa(t)c_x(t) <a|V|i>
             ia = hole_indices(-i,x) - 1 
             psi_x = dconjg( psi_det(x) )             
             do a=1, nva
                rdum = rdum + 2.d0*real(psi_x*psi_det(ia+a))*vabs_a(noa+a,i)
             end do
             
             !: sum_ja : c_ixa(t)c_jxa(t) <j|V|i>
             ia = hole_indices(-i,x) - 1 
             do j=1, noa
                jb = hole_indices(-j,x) - 1
                vdum = vabs_a(j,i)
                do a=1, nva
                   rdum = rdum - real(dconjg(psi_det(ia+a))*psi_det(jb+a))*vdum
                end do
             end do
             
             !: sum_ya : c_ixa(t)c_iya(t) <y|V|x>
             ia = hole_indices(-i,x) - 1 
             do y=1, nob
                jb   = hole_indices(-i,y) - 1 
                vdum = vabs_b(y,x)
                do a=1, nva
                   rdum = rdum - real(dconjg(psi_det(ia+a))*psi_det(jb+a))*vdum
                end do
             end do
             
             !: sum_ab : c_ixa(t)c_ixb(t) <a|V|b>
             ia = hole_indices(-i,x) - 1 
             do a=1, nva
                psi_ia = dconjg( psi_det(ia+a) )
                do b=1, nva
                   rdum = rdum + real(psi_ia*psi_det(ia+b))*vabs_a(noa+a,noa+b)
                end do
             end do
             
             !: record
             tmp_aa(i,x) = rdum
             
          end do alpha
          
          beta : do i=x+1, nob
             
             !: diagonal
             rdum = 0.d0
             ia = hole_indices(i,x) - 1 
             do a=1, nvb
                rdum = rdum + real(dconjg(psi_det(ia+a))*psi_det(ia+a))
             end do
             rdum = rdum * vsum 


             !: contributions from c_x and c_i
             rdum = rdum + dconjg(psi_det(x))*psi_det(x)*vabs_b(i,i) &
                  + dconjg(psi_det(i))*psi_det(i)*vabs_b(x,x) &
                  + real( dconjg(psi_det(x))*psi_det(i) + dconjg(psi_det(i))*psi_det(x) )*vabs_b(i,x)
             
             !: sum_a : 2cixa(t)cx(t)<a|V|i> - 2cixa(t)ci(t)<a|V|x>
             ia    = hole_indices(i,x) - 1 
             psi_x = dconjg( psi_det(x) ) 
             psi_i = dconjg( psi_det(i) )
             do a=1, nvb
                rdum = rdum + 2.d0*real(psi_x*psi_det(ia+a))*vabs_b(nob+a,i) - 2.d0*real(psi_i*psi_det(ia+a))*vabs_b(nob+a,x)
             end do

             !: sum_ja : cixa(t)cjxa(t)<j|V|i> for j < x 
             ia = hole_indices(i,x) - 1 
             do j=x+1, nob
                jb   = hole_indices(j,x) - 1 
                vdum = vabs_b(j,i)
                do a=1, nvb
                   rdum = rdum - real(dconjg(psi_det(ia+a))*psi_det(jb+a))*vdum
                end do
             end do

             !: sum_ya : c_ixa(t)*ciya(t)<y|V|x> for y>i
             ia = hole_indices(i,x) - 1 
             do y=1, (i-1)
                jb = hole_indices(i,y) - 1 
                vdum = vabs_b(y,x)
                do a=1, nvb
                   rdum = rdum - real(dconjg(psi_det(ia+a))*psi_det(jb+a))*vdum
                end do
             end do

             !: sum_ab : c_ixa(t)*cixb(t)<a|V|b>
             ia = hole_indices(i,x) - 1 
             do a=1, nvb
                do b=1, nvb
                   rdum = rdum + real(dconjg(psi_det(ia+a))*psi_det(ia+b))*vabs_b(nob+b,nob+a)
                end do
             end do

             !:sum_ya : c_ixa(t)*cxya(t)<y|V|i> for y>x
             ia = hole_indices(i,x) - 1 
             do y=1, (x-1)
                jb   = hole_indices(x,y) - 1 
                vdum = vabs_b(y,i)
                do a=1, nvb
                   rdum = rdum + real(dconjg(psi_det(ia+a))*psi_det(jb+a))*vdum
                end do
             end do

             !: sum_ja : c_ixa(t)*cjia(t)<j|V|x> for j<i
             ia = hole_indices(i,x) - 1 
             do j=(i+1), nob
                jb   = hole_indices(j,i) - 1 
                vdum = vabs_b(j,x)
                do a=1, nvb
                   rdum = rdum + real(dconjg(psi_det(ia+a))*psi_det(jb+a))*vdum
                end do
             end do
             
             tmp_bb(i,x) = rdum

          end do beta
       end do
       
       trate(itime) = ( sum( tmp_aa ) + sum( tmp_bb ) ) * 2.d0 / au2fs
       rate_aa(:,:,itime) = tmp_aa(:,:)
       rate_bb(:,:,itime) = tmp_bb(:,:)
              
    end do itime
    !$OMP END DO
    !$OMP END PARALLEL

    rate_aa = rate_aa * 2.d0 / au2fs
    rate_bb = rate_bb * 2.d0 / au2fs    

    !: print python
    ofile='RATE_IP'//trim(e_d)//trim(c0)//'.out'
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
       do i=x+1, nob
          write( ci, '(i0)' ) i
          write(100,'(A)',advance='no') '('//trim(c0)//','//trim(ci)//') '
       end do
    end do
    write(100,'(A)') ''
    
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

    write(iout,'(A)') ' --> generated file '//trim(ofile)
    flush(iout)
    
  end subroutine get_rateAB_ip
  !:--------------------------:!
  !:--------------------------:!
end module rates
