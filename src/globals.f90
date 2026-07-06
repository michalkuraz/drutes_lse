! globals.f90
module globals
  use typy
  use smartarray
  implicit none

  real(kind=rkind) :: tmp, tmp1, sim_time, start_time, end_time, time, time_step
  logical          :: www

  type :: nodes_str
     real(kind=rkind), allocatable :: data(:,:)
     logical,       allocatable    :: watershed(:)
     integer(kind=ikind)           :: kolik = 0
     real(kind=rkind), allocatable :: altitude(:)
  end type nodes_str

  type, public :: stream_str
    integer(kind=ikind), dimension(3) :: els = 0
    real(kind=rkind),    dimension(3) :: slopes = 0.0_rkind
    real(kind=rkind),    dimension(3) :: widths = 0.0_rkind
  end type stream_str

  type, public :: hydrobal_str
     real(kind=rkind) :: deltas  = 0.0_rkind
     real(kind=rkind) :: inflow  = 0.0_rkind
     real(kind=rkind) :: outflow = 0.0_rkind
     real(kind=rkind) :: ET      = 0.0_rkind
     real(kind=rkind) :: q1      = 0.0_rkind
     real(kind=rkind) :: q2      = 0.0_rkind
     real(kind=rkind) :: q3      = 0.0_rkind
     real(kind=rkind) :: pc      = 0.0_rkind
     real(kind=rkind) :: bf      = 0.0_rkind
     real(kind=rkind) :: storage = 0.0_rkind
  end type hydrobal_str

  type :: elements_str
     integer(kind=ikind), allocatable :: data(:,:)
     real(kind=rkind),    allocatable :: area(:)
     integer(kind=ikind), allocatable :: material(:)
     type(hydrobal_str),  allocatable :: hydrobal(:)
     integer(kind=ikind)              :: kolik = 0
     integer(kind=ikind), allocatable :: neighbours(:,:)
     real(kind=rkind),    allocatable :: avgalt(:)
     real(kind=rkind),    allocatable :: slope(:)
     real(kind=rkind),    allocatable :: overflow(:)
     type(stream_str),    allocatable :: downstream(:)
     type(stream_str),    allocatable :: upstream(:)
     logical,             allocatable :: outlet(:), watershed(:)
     integer(kind=ikind), allocatable :: ndwatershed(:,:), ndoutlet(:,:)
  end type elements_str

  type(nodes_str)    :: nodes
  type(elements_str) :: elements

  real(kind=rkind), dimension(:), allocatable :: meteotime

  real(kind=rkind), parameter    :: ntot_days  = 10.0_rkind
  real(kind=rkind), parameter    :: dt_hours   = 1.0_rkind
  real(kind=rkind), parameter    :: dt_days    = dt_hours / 24.0_rkind
  real(kind=rkind), parameter    :: dt_seconds = dt_hours * 3600.0_rkind
  integer(kind=ikind), parameter :: n_steps    = int(ntot_days / dt_days)

  integer(kind=ikind) :: CN, Julian_day
  real(kind=rkind)    :: phi, as, bs, z, alpha, sigma, gsc, ccrop

  real(kind=rkind), allocatable :: precip(:,:), qinter(:,:), qout(:,:)
  real(kind=rkind), allocatable :: conduct(:), G(:)
  real(kind=rkind), allocatable :: Tmax(:,:), Tmin(:,:), Tmean(:,:)
  real(kind=rkind), allocatable :: RHmax(:,:), RHmin(:,:), uz(:,:), soilcontent(:,:)

  real(kind=rkind), allocatable :: Qin_result(:,:), Qout_result(:,:)
  real(kind=rkind), allocatable :: Overflow_result(:,:), Storage_result(:,:)
  real(kind=rkind), allocatable :: outlet_Q_m3s(:)
  real(kind=rkind), allocatable :: deltas(:,:)
  real(kind=rkind), allocatable :: total_deltaS(:,:)

  integer, parameter :: terminal = 6

  type(smartarray_real) :: timestamps

  integer(kind=ikind), allocatable :: upstream_count(:)
  integer(kind=ikind), allocatable :: upstream_list(:,:)

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

  real(kind=rkind), allocatable :: storage(:)
  real(kind=rkind), allocatable :: capacity(:)
  real(kind=rkind), allocatable :: outlet_Q(:)

  integer(kind=ikind), allocatable :: downstream(:)
  integer(kind=ikind), allocatable :: flow_order(:)

  real(kind=rkind), allocatable :: Ksat_surf(:)
  real(kind=rkind), allocatable :: Ksat_sub(:)
  real(kind=rkind), allocatable :: Ksat_gw(:)

  real(kind=rkind) :: ksurf_exp, ksub_exp
  real(kind=rkind) :: cn_slope_coeff, infil_slope_coeff, min_edge_width
  real(kind=rkind) :: qgw_slope_coeff, storage_slope_coeff
  real(kind=rkind) :: theta_r, theta_s
  real(kind=rkind) :: Beta1, Beta2, Beta3, Beta4, Beta5, z1, z2

  integer(kind=ikind) :: nsteps = 0

  ! ============================================================
  ! ODE water balance variables, all in depth units [mm]
  ! ============================================================

  ! Storages [mm]
  real(kind=rkind), allocatable :: Ssurf(:), Ssub(:), Sgw(:)

  ! Fluxes [mm/timestep]
  real(kind=rkind), allocatable :: Pm(:,:), E_m(:,:), If_m(:,:)
  real(kind=rkind), allocatable :: q1(:,:), q2(:,:), q3(:,:), pc(:,:), bf(:,:)

  ! ODE increments [mm/timestep]
  real(kind=rkind), allocatable :: dSsurf(:,:), dSsub(:,:), dSgw(:,:)

  ! Storage histories [mm]
  real(kind=rkind), allocatable :: Ssurf_hist(:,:), Ssub_hist(:,:), Sgw_hist(:,:)

  type :: meteodata_str
  type(smartarray_real) :: time
  type(smartarray_real) :: rainfall
  type(smartarray_real) :: Tmax
  type(smartarray_real) :: Tmin
  type(smartarray_real) :: wind
  type(smartarray_real) :: RHmax
  type(smartarray_real) :: RHmin
  type(smartarray_real) :: soilcontent
end type meteodata_str

type(meteodata_str), dimension(:), allocatable :: meteodata

end module globals