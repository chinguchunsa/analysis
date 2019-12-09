module global_variables 
    
    
    implicit none
    
    !: field
    !character(5),   parameter :: e_d = '-e1-d'
    
    integer(8),    parameter :: i_in  = 178
    character(12), parameter :: input = 'analysis.in'

    !: outfile    
    integer(8),    parameter :: iout = 78
    character(20), parameter :: log = 'analysis.out'
    
    
    !: namelist for files
    character(10)   :: e_d
    character(1000) :: &
         f_cidir,      & !: directory of CI coefficients
         f_restart,    & !: for CI vectors
         f_mo_restart    !: for MO matrix elements


    !: namelist for job
    integer(8) :: &
         start_dir,  & !: start analysis starting for this direction
         end_dir,    & !: end with this direction
         start_time, & !: start from this time point
         end_time      !: end with this time point


    character(8) :: jobtype !: 'ip' or 'cis'

    
    logical(8) :: &
         Qprint_movie,   & !: print density for movie making
         Qget_dos,       & !: get density of states
         Qget_rate_occ,  & !: compute total loss from occupied
         Qget_mo_den       !: get MO density and population
    
    
    real(8), parameter :: au2fs = 2.418884326509d-2
    real(8), parameter :: au2eV = 27.21138602d0
    

    !: define me before moving on
    integer(8) :: ntimes, nstuse, nstates
    integer(8) :: noa, nva, nob, nvb
    integer(8) :: nrorb, norb
    logical :: unrestricted 
    

    !: state indices
    integer(8), allocatable :: hole1(:), part1(:), hole_indices(:,:)    
    integer(8), allocatable :: hole(:,:), part(:,:)
      
    
    !: read in matrix elements in MO
    real(8), allocatable :: dipx_a(:,:), dipy_a(:,:), dipz_a(:,:), vabs_a(:,:)
    real(8), allocatable :: dipx_b(:,:), dipy_b(:,:), dipz_b(:,:), vabs_b(:,:)


    !: time-dependent data
    real(8), allocatable :: ci_vec(:,:), ci_eig(:)
    real(8), allocatable :: ci_vabs(:,:), tdx(:,:), tdy(:,:), tdz(:,:)
    real(8), allocatable :: time(:), norm_sq(:)
    real(8), allocatable :: efield(:), dirx(:), diry(:), dirz(:) 
    complex(8), allocatable :: psi(:,:)
    
    
  end module global_variables
