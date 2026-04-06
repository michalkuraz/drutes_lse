! globals.f90
module globals
  use typy
  implicit none

 
  ! Time/progress helpers
  real(kind=rkind) :: tmp, tmp1, sim_time, start_time, end_time, time, time_step
  logical          :: www
  
  

  ! ---------- Node structure ----------
  type :: nodes_str
     real(kind=rkind), allocatable :: data(:,:)     ! (n_nodes,2): x,y
     logical,       allocatable    :: watershed(:)  ! mask if needed
     integer(kind=ikind)           :: kolik = 0
     real(kind=rkind), allocatable :: altitude(:)   ! elevation z
  end type nodes_str
  
  type, public :: stream_str
    integer(kind=ikind), dimension(3) :: els = 0
    real(kind=rkind), dimension(3) :: slopes = 0
    real(kind=rkind), dimension(3) :: widths = 0
  end type stream_str

  ! ---------- Hydrological balance (per element, per step) ----------
  type, public :: hydrobal_str
     real(kind=rkind) :: deltas  = 0.0_rkind
     real(kind=rkind) :: inflow  = 0.0_rkind
     real(kind=rkind) :: outflow = 0.0_rkind
     real(kind=rkind) :: Li      = 0.0_rkind
     real(kind=rkind) :: ET      = 0.0_rkind
     real(kind=rkind) :: Qgw     = 0.0_rkind
     real(kind=rkind) :: Qsurf   = 0.0_rkind
     real(kind=rkind) :: storage = 0.0_rkind
  end type hydrobal_str

  ! ---------- Element structure ----------
  type :: elements_str
     integer(kind=ikind), allocatable :: data(:,:)      ! connectivity (nel,3)
     real(kind=rkind),    allocatable :: area(:)        ! element area
     integer(kind=ikind), allocatable :: material(:)    ! optional soil ID
     type(hydrobal_str),  allocatable :: hydrobal(:)    ! hydro state
     integer(kind=ikind)              :: kolik = 0      ! number of elements
     integer(kind=ikind), allocatable :: neighbours(:,:)! (nel,max_neigh)
     real(kind=rkind),    allocatable :: avgalt(:)      ! mean elevation
     real(kind=rkind),    allocatable :: overflow(:)    ! local generated overflow
     type(stream_str), allocatable    :: downstream(:)
     type(stream_str), allocatable    :: upstream(:)
     logical, allocatable             :: outlet(:), watershed(:)
     integer(kind=ikind), allocatable :: ndwatershed(:,:), ndoutlet(:,:)
  end type elements_str

  type(nodes_str)    :: nodes
  type(elements_str) :: elements
  
  
  
  real(kind=rkind), dimension(:), allocatable :: meteotime
  
  
  
  !type, public :: hydrodata_str
    !> meteorological data
   
  ! real(kind=rkind), dimension(:), allocatable :: precip, G, Tmax, Tmin, Tmean, uz, soilwcontent, RHmax, RHmin
    !> land cover data
   
  ! real(kind=rkind) :: conduct
  
  
  !end type hydrodata_str
  
  
  !type(hydrodata_str), dimension(:), allocatable :: hydrodata
  
  

  ! Simulation time discretisation
  integer(kind=ikind), parameter :: ntot_days = 10
  real(kind=rkind),   parameter :: dt_days   = 1.0_rkind
  integer(kind=ikind), parameter :: n_steps = int(ntot_days / dt_days)
  integer(kind=ikind) :: CN, Julian_day 
  real(kind=rkind) :: phi, as, bs, z, alpha, sigma, gsc, ccrop

  ! Element-based hydro inputs & outputs (time series)

  ! Element-based hydro inputs & outputs (time series)
  real(kind=rkind), allocatable :: precip(:,:), qinter(:,:), qout(:,:)
  real(kind=rkind), allocatable :: conduct(:), G(:), Tmax(:,:), Tmin(:,:), Tmean(:,:)
  real(kind=rkind), allocatable :: RHmax(:,:), RHmin(:,:), uz(:,:), soilcontent(:,:)

  real(kind=rkind), allocatable :: Qsurf_result(:,:), ET_flux(:,:), &
                                   L_result(:,:), Qgw_result(:,:), deltas(:,:)

  integer, parameter :: terminal = 6


  ! NEW: upstream connectivity
  integer(kind=ikind), allocatable :: upstream_count(:)
  integer(kind=ikind), allocatable :: upstream_list(:,:)


  ! ---------- Configuration structure (kept from DRUtES style) ----------
  type, public :: configuration
    character(len=1)    :: run_dual
    character(len=1)    :: damped_newton
    integer(1)          :: dimen = 2
    integer(kind=ikind) :: mesh_type
    logical             :: adapt_observe
    logical             :: run_from_backup
    integer(kind=ikind) :: fnc_method
    real(kind=rkind)    :: fnc_discr_length
    integer(kind=ikind) :: it_method
    character(len=256)  :: name
    character(len=4096) :: fullname
    logical             :: rotsym = .false.
    logical             :: check4mass = .false.
  end type configuration

  type(configuration), public :: drutes_config

  ! ---------- New routing / storage globals ----------
  real(kind=rkind), allocatable :: storage(:)    ! [mm] surface/depression storage
  real(kind=rkind), allocatable :: capacity(:)   ! [mm] max local storage capacity
  real(kind=rkind), allocatable :: outlet_Q(:)   ! [mm] outlet discharge per time step

  integer(kind=ikind), allocatable :: downstream(:) ! downstream element ID (0 = outlet)
  integer(kind=ikind), allocatable :: flow_order(:) ! processing order (upstream→downstream)

  ! Number of time steps (set in initvals)
  integer(kind=ikind) :: nsteps = 0

end module globals
