module read_files

  use global_variables
  implicit none

contains

  !:--------------------------:!
  !:--------------------------:!
  subroutine read_tdci( idir )
    

    implicit none
    
    integer(8), intent(in) :: idir
    
    integer(8) :: ntimes1, itime, i, j
    character(5)    :: cdir
    character(1000) :: input
    real(8) :: rdum !: will read efield2 in later
    real(8), allocatable :: readtmp(:), creadtmp(:)
    
    
    write(cdir, '(i0)' ) idir
    input = trim(f_cidir)//'CI'//trim(e_d)//trim(cdir)//'.bin'
    
    
    !: read input
    open( unit=100,file=trim(input), form='unformatted' )
    read(100) ntimes1, nstuse, nstates    
    
    
    !: allocate + prep
    if ( .not.allocated(time) ) then
       allocate( time(ntimes1), norm_sq(ntimes1) )
       allocate( efield(ntimes1), dirx(ntimes1), diry(ntimes1), dirz(ntimes1) )
       allocate( psi(nstuse,ntimes1) )
    else
       time    = 0.d0
       norm_sq = 0.d0
       efield  = 0.d0 
       dirx = 0.d0 
       diry = 0.d0 
       dirz = 0.d0
       psi  = dcmplx( 0.d0,0.d0 )
    end if
    
    
    allocate( readtmp(nstuse), creadtmp(nstuse) )  
    do itime=1, ntimes1
       
       read(100) time(itime), efield(itime), rdum, dirx(itime), diry(itime), dirz(itime), rdum, rdum, rdum, norm_sq(itime) 
       read(100) (  readtmp(i), i=1, nstuse )
       read(100) ( creadtmp(i), i=1, nstuse )
       psi(:,itime) = dcmplx( readtmp, creadtmp )

    end do

    close(100)
    deallocate( readtmp, creadtmp )
    
    
    write(iout,*) ' --> allocated time, norm_sq, efield, dirx, diry, dirz '
    write(iout,*) ' --> finished reading psi(t)'
    flush(iout)
    
    
  end subroutine read_tdci
  !:--------------------------:!
  !:--------------------------:!
  subroutine read_ci_vec
    
    
    implicit none
    
    integer(8) :: i, j, ij
    real(8), allocatable :: rdum(:)
    

    open( unit=100, file=trim(f_restart), form='unformatted' )
    read(100) nstates, nstuse
    
    allocate( ci_vec(nstates,nstates) )
    read(100) ( ( ci_vec(j,i), j=1, nstates ), i=1, nstates )
    
    
    if ( Qget_dos ) then       

       allocate( ci_vabs(nstates,nstates), ci_eig(nstates) )
       allocate( rdum(nstates*nstates) )
       allocate( tdx(nstates,nstates), tdy(nstates,nstates), tdz(nstates,nstates) )

       ci_vabs = 0.d0 ; ci_eig = 0.d0
       tdx = 0.d0 ; tdy = 0.d0 ; tdz = 0.d0

       read(100) ( ci_eig(i),  i=1, nstates )

       !: tdx
       read(100) ( rdum(i), i=1, nstates*nstates )
       do i=1, nstuse
          ij = (i-1)*nstuse !+ j
          tdx(:,i) = rdum(ij+1:ij+nstuse)
       end do

       !: tdy
       read(100) ( rdum(i), i=1, nstates*nstates )
       do i=1, nstuse
          ij = (i-1)*nstuse !+ j
          tdy(:,i) = rdum(ij+1:ij+nstuse)
       end do

       !: tdz
       read(100) ( rdum(i), i=1, nstates*nstates )
       do i=1, nstuse
          ij = (i-1)*nstuse !+ j
          tdz(:,i) = rdum(ij+1:ij+nstuse)
       end do
       
       
       read(100) ( ( ci_vabs(j,i), j=1, nstuse ) , i=1, nstuse )


       !: pack down


    end if
    

    close(100) 
    
    
    write(iout,*) ' --> allocated ci_vec '
    write(iout,*) ' --> finished reading ci_vec '
    flush(iout)
    
    
  end subroutine read_ci_vec
  !:--------------------------:!
  !:--------------------------:!
  subroutine read_mo
    

    implicit none
  
    integer(8) :: i
    integer(8) :: na, nb
    integer(4) :: test1, test2
    
    open( unit=100, file=trim(f_mo_restart), form='unformatted' )
    if ( unrestricted ) then
       read(100) noa, nva, nob, nvb
    else
       read(100) noa, nva
       nob = 0
       nvb = 0
       write(iout,"('noa=',i0,' nva=',i0)") noa, nva
       write(iout,"('nob=',i0,' nvb=',i0)") nob, nvb
       write(iout,*) 'unrestricted=', unrestricted
       flush(iout)
    end if

    
    if( nob.ne.0 ) then
       na = noa + nva
       nb = nob + nvb
       nrorb = na 
       norb = nrorb + nb
       !: index-ing for beta elements to be consistent
       allocate( dipx_a(na,na), dipy_a(na,na), dipz_a(na,na), vabs_a(na,na) )
       allocate( dipx_b(nb,nb), dipy_b(nb,nb), dipz_b(nb,nb), vabs_b(nb,nb) )
    else
       na = noa + nva
       nrorb = na 
       norb = nrorb
       allocate( dipx_a(na,na), dipy_a(na,na), dipz_a(na,na), vabs_a(na,na) )
    end if
    

    read(100) ( vabs_a(:,i), i=1, na )
    read(100) ( dipx_a(:,i), i=1, na )
    read(100) ( dipy_a(:,i), i=1, na )
    read(100) ( dipz_a(:,i), i=1, na )
    if ( nob.gt.0 ) then
       read(100) ( vabs_b(:,i), i=1, nb )
       read(100) ( dipx_b(:,i), i=1, nb )
       read(100) ( dipy_b(:,i), i=1, nb )
       read(100) ( dipz_b(:,i), i=1, nb )
    end if
    
    close(100)
    write(iout,*) ' --> allocated dipx, dipy, dipz, vabs '
    write(iout,*) ' --> finished reading MO elements'
    flush(iout)
    
    
  end subroutine read_mo
  !:--------------------------:!
  !:--------------------------:!
end module read_files
