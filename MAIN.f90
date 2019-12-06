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
       Qprint_pretty, Qprint_python, Qprint_movie, &
       Qget_dos, Qget_mo_den, Qget_rate, Qget_rate_occ, Qget_trans_den
  
  
  !: set default
  e_d = '-e1-d'
  unrestricted = .True.

  open( unit=i_in, file=trim(input) )
  read( i_in, nml=FILES )
  rewind( i_in )
  read( i_in, nml=JOB_DETAILS )
  close( i_in )
  
  
  open( unit=iout, file=trim(log) )

  
  call read_ci_vec  !: read in CI eigenvectors
  call read_mo      !: read in MO elements
  
  
  !: initialize variables
  ntimes = end_time - start_time + 1   

  call get_indices  !: get indices    
  
  
  if( Qget_dos  )  call get_dos

  
  do idir=start_dir, end_dir
     
     write(iout,"(' *******************')" )  
     write(iout,"(' working on direction',i0)" )  idir
     flush(iout)
          
     call read_tdci(idir)                  !: read in time-dependent CI coefficients  
     call summarize_coefficients(idir)     !: summarize ground and total excited state populations
     if( Qget_rate )      call get_rate(idir)          !: get rate 
     if( Qget_rate_occ )  call get_rate_occ( idir )    !: get hole in occ
     if( Qget_mo_den )    call get_density( idir )     !: get population and density matrix
     if( Qget_trans_den ) call get_hp_density( idir )  !: get hole and particle density matrix
     
  end do
  
  close(iout)
  
  
end program analysis
