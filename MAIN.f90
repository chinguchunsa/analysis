program analysis

    
  use global_variables
  use read_files
  use utils
  use rates
  use density
  use misc
  implicit none

  integer(8) :: idir
  
  
  namelist /FILES/ f_cidir, f_mo_restart, f_restart, e_d
  namelist /JOB_DETAILS/ start_dir, end_dir, start_time, end_time, unrestricted, &
       jobtype, Qprint_movie, Qget_dos, Qget_mo_den, Qget_rate_occ
  
  !: set default
  e_d = '-e1-d'
  jobtype = 'cis'
  unrestricted = .True.

  
  !: read in analysis.in
  open( unit=i_in, file=trim(input) )
  read( i_in, nml=FILES )
  rewind( i_in )
  read( i_in, nml=JOB_DETAILS )
  close( i_in )

  
  !: open log file 
  open( unit=iout, file=trim(log) )
  
  
  call read_ci_vec  !: read in CI eigenvectors
  call read_mo      !: read in MO elements
  
  
  !: initialize variables
  ntimes = end_time - start_time + 1   

  
  !call get_indices  !: get indices    
  select case( trim(jobtype) ) 
  case( 'cis' ) ; call get_indices
  case( 'ip'  ) ; call get_indices_ip
  end select
  
  
  if( Qget_dos )  call get_dos
  
  
  do idir=start_dir, end_dir
     

     write(iout,"(' *******************')" )  
     write(iout,"(' working on direction',i0)" )  idir
     flush(iout)
     

     call read_tdci(idir) !: read in time-dependent CI coefficients  
     
     
     select case( trim(jobtype) ) 
     case( 'cis' ) ; call summarize_coefficients(idir)  !: summarize ground and total excited state populations
     case( 'ip' )  ; call summarize_coefficients_ip(idir)
     end select
     
     if( Qget_rate_occ ) then
        select case( trim(jobtype) ) 
        case( 'cis' ) ; call get_rate_occ(idir) 
        case( 'ip'  ) ; call get_rateAB_ip( idir )     !: get hole in occ
        end select
     end if
     if( Qget_mo_den )  then
        select case( trim(jobtype)  )
        case( 'cis' ) ; call get_density( idir )     !: get population and density matrix
        case( 'ip' )  ; call get_density_ip( idir )  !: get population and density matrix
        end select
     end if
     
  end do
  
  close(iout)
  
  
end program analysis
