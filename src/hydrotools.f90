module hydrotools
  use typy
  use globals
  use tools
  use readtools, only: comment
  implicit none

contains

  !==============================================================
  ! READ MESH
  !==============================================================
  subroutine read_mesh(filename)

  character(len=*), intent(in) :: filename

  integer :: fileid
  integer :: ios
  integer(kind=ikind) :: i
  real(kind=rkind)    :: node_id
  real(kind=rkind)    :: element_id

  open(newunit=fileid, file=filename, status='old', action='read', iostat=ios)

  if (ios /= 0_ikind) then
     error stop 'Error: could not open mesh file.'
  end if

  call comment(fileid)
  read(fileid, *, iostat=ios) nodes%kolik
  if (ios /= 0_ikind) error stop 'Error: failed to read number of nodes.'

  allocate(nodes%data(nodes%kolik, 2))
  allocate(nodes%altitude(nodes%kolik))

  do i = 1_ikind, nodes%kolik

     call comment(fileid)

     read(fileid, *, iostat=ios) node_id, nodes%data(i,1), &
                                 nodes%data(i,2), nodes%altitude(i)

     if (ios /= 0_ikind) then
        error stop 'Error: failed to read node data.'
     end if

  end do

  call comment(fileid)
  read(fileid, *, iostat=ios) elements%kolik
  if (ios /= 0_ikind) error stop 'Error: failed to read number of elements.'

  call mesh_allocater()

  do i = 1_ikind, elements%kolik

     call comment(fileid)

     read(fileid, *, iostat=ios) element_id, elements%data(i,1), &
                                    elements%data(i,2), &
                                    elements%data(i,3)

     if (ios /= 0_ikind) then
        error stop 'Error: failed to read element data.'
     end if

  end do

  close(fileid)

end subroutine read_mesh


  !==============================================================
  ! ALLOCATE MESH ARRAYS
  !==============================================================
  subroutine mesh_allocater()

  integer(kind=ikind) :: el

  if (.not. allocated(elements%data)) then
     allocate(elements%data(elements%kolik, 3))
  end if

  if (.not. allocated(elements%area)) then
     allocate(elements%area(elements%kolik))
  end if

  if (.not. allocated(elements%material)) then
     allocate(elements%material(elements%kolik))
  end if

  if (.not. allocated(elements%avgalt)) then
     allocate(elements%avgalt(elements%kolik))
  end if

  if (.not. allocated(elements%slope)) then
     allocate(elements%slope(elements%kolik))
  end if

  if (.not. allocated(elements%overflow)) then
     allocate(elements%overflow(elements%kolik))
  end if

  if (.not. allocated(elements%neighbours)) then
     allocate(elements%neighbours(elements%kolik, 3))
  end if

  if (.not. allocated(elements%upstream)) then
     allocate(elements%upstream(elements%kolik))
  end if

  if (.not. allocated(elements%downstream)) then
     allocate(elements%downstream(elements%kolik))
  end if

  if (.not. allocated(elements%hydrobal)) then
     allocate(elements%hydrobal(elements%kolik))
  end if

  elements%data       = 0_ikind
  elements%area       = 0.0_rkind
  elements%material   = 0_ikind
  elements%avgalt     = 0.0_rkind
  elements%slope      = 0.0_rkind
  elements%overflow   = 0.0_rkind
  elements%neighbours = 0_ikind

    do el = 1_ikind, elements%kolik

     elements%upstream(el)%els    = 0_ikind
     elements%upstream(el)%slopes = 0.0_rkind
     elements%upstream(el)%widths = 0.0_rkind

     elements%downstream(el)%els    = 0_ikind
     elements%downstream(el)%slopes = 0.0_rkind
     elements%downstream(el)%widths = 0.0_rkind

     elements%hydrobal(el)%inflow  = 0.0_rkind
     elements%hydrobal(el)%outflow = 0.0_rkind
     elements%hydrobal(el)%ET      = 0.0_rkind
     !elements%hydrobal(el)%Li      = 0.0_rkind
     elements%hydrobal(el)%Qgw     = 0.0_rkind
     elements%hydrobal(el)%deltas  = 0.0_rkind

  end do

end subroutine mesh_allocater


  !==============================================================
  ! INITIALIZE HYDROLOGY
  !==============================================================
  subroutine init_hydro()
    integer(kind=ikind) :: i, t, d

    real(kind=rkind), dimension(10) :: pday
    real(kind=rkind), dimension(10) :: uzday
    real(kind=rkind), dimension(10) :: tmaxday
    real(kind=rkind), dimension(10) :: tminday
    real(kind=rkind), dimension(10) :: tmeanday
    real(kind=rkind), dimension(10) :: rhmaxday
    real(kind=rkind), dimension(10) :: rhminday
    real(kind=rkind), dimension(10) :: soilday

    !------------------------------------------------------------
    ! Allocate forcing and result arrays
    !------------------------------------------------------------
    allocate(precip(elements%kolik, n_steps))
    allocate(qinter(elements%kolik, n_steps))
    allocate(qout(elements%kolik, n_steps))

    allocate(conduct(elements%kolik))
    allocate(Ksat_surf(elements%kolik))
    allocate(Ksat_sub(elements%kolik))
    allocate(Ksat_gw(elements%kolik))
    allocate(G(elements%kolik))

    allocate(Tmax(elements%kolik, n_steps))
    allocate(Tmin(elements%kolik, n_steps))
    allocate(Tmean(elements%kolik, n_steps))
    allocate(RHmax(elements%kolik, n_steps))
    allocate(RHmin(elements%kolik, n_steps))
    allocate(uz(elements%kolik, n_steps))
    allocate(soilcontent(elements%kolik, n_steps))

    allocate(Qsurf_result(elements%kolik, n_steps))
    allocate(ET_flux(elements%kolik, n_steps))
    !allocate(L_result(elements%kolik, n_steps))
    allocate(Qgw_result(elements%kolik, n_steps))
    allocate(Inf_result(elements%kolik, n_steps))
    allocate(pc(elements%kolik, n_steps))
    allocate(q2(elements%kolik, n_steps))
    allocate(q3(elements%kolik, n_steps))
    allocate(bf(elements%kolik, n_steps))
    allocate(deltas(elements%kolik, n_steps))

    allocate(Qin_result(elements%kolik, n_steps))
    allocate(Qout_result(elements%kolik, n_steps))
    allocate(Overflow_result(elements%kolik, n_steps))
    allocate(Storage_result(elements%kolik, n_steps))
    allocate(outlet_Q_m3s(n_steps))

    allocate(storage(elements%kolik))
    allocate(capacity(elements%kolik))
    allocate(outlet_Q(n_steps))

    !------------------------------------------------------------
    ! Allocate ODE water balance arrays
    !------------------------------------------------------------
    allocate(Pm(elements%kolik, n_steps))
    allocate(E_m(elements%kolik, n_steps))
    allocate(If_m(elements%kolik, n_steps))
    allocate(Qsurf_m(elements%kolik, n_steps))
    !allocate(Tv_m(elements%kolik, n_steps))
    !allocate(Qsub_m(elements%kolik, n_steps))
    !allocate(Rv_m(elements%kolik, n_steps))
    !allocate(Qgw_m(elements%kolik, n_steps))

    allocate(dVsurf(elements%kolik, n_steps))
    allocate(dVsub(elements%kolik, n_steps))
    allocate(dVgw(elements%kolik, n_steps))


    !allocate(Ssurf(elements%kolik))
    !allocate(Ssub(elements%kolik))
    !allocate(Sgw(elements%kolik))

    !allocate(Ssurf_hist(elements%kolik, n_steps))
    !allocate(Ssub_hist(elements%kolik, n_steps))
    !allocate(Sgw_hist(elements%kolik, n_steps))
    allocate(inf(elements%kolik, n_steps))

    !------------------------------------------------------------
    ! Initialize everything to zero
    !------------------------------------------------------------
    precip      = 0.0_rkind
    qinter      = 0.0_rkind
    qout        = 0.0_rkind
    conduct     = 0.0_rkind
    G           = 0.0_rkind

    Tmax        = 0.0_rkind
    Tmin        = 0.0_rkind
    Tmean       = 0.0_rkind
    RHmax       = 0.0_rkind
    RHmin       = 0.0_rkind
    uz          = 0.0_rkind
    soilcontent = 0.0_rkind

    Qsurf_result = 0.0_rkind
    ET_flux      = 0.0_rkind
    !L_result     = 0.0_rkind
    Qgw_result   = 0.0_rkind
    Inf_result   = 0.0_rkind
    pc           = 0.0_rkind
    q2           = 0.0_rkind
    q3           = 0.0_rkind
    bf           = 0.0_rkind
    storage      = 0.0_rkind
    deltas       = 0.0_rkind

    Qin_result      = 0.0_rkind
    Qout_result     = 0.0_rkind
    Overflow_result = 0.0_rkind
    Storage_result  = 0.0_rkind
    outlet_Q_m3s    = 0.0_rkind

    
    capacity = 8.0_rkind
    outlet_Q = 0.0_rkind

    Pm      = 0.0_rkind
    E_m     = 0.0_rkind
    If_m    = 0.0_rkind
    Qsurf_m = 0.0_rkind
    Tv_m    = 0.0_rkind
    Qsub_m  = 0.0_rkind
    Rv_m    = 0.0_rkind
    Qgw_m   = 0.0_rkind

    dVsurf = 0.0_rkind
    dVsub  = 0.0_rkind
    dVgw   = 0.0_rkind

    !Ssurf = 0.0_rkind
    !Ssub  = 0.0_rkind
    !Sgw   = 0.0_rkind

    !Ssurf_hist = 0.0_rkind
    !Ssub_hist  = 0.0_rkind
    !Sgw_hist   = 0.0_rkind
    inf = 0.0_rkind

    !------------------------------------------------------------
    ! Initialize hydrobal structure
    !------------------------------------------------------------
    if (allocated(elements%hydrobal)) then
      do i = 1, elements%kolik
        elements%hydrobal(i)%inflow  = 0.0_rkind
        elements%hydrobal(i)%outflow = 0.0_rkind
        elements%hydrobal(i)%Li      = 0.0_rkind
        elements%hydrobal(i)%ET      = 0.0_rkind
        elements%hydrobal(i)%Qgw     = 0.0_rkind
        elements%hydrobal(i)%Qsurf   = 0.0_rkind
        elements%hydrobal(i)%q2      = 0.0_rkind
        elements%hydrobal(i)%q3      = 0.0_rkind
        elements%hydrobal(i)%pc      = 0.0_rkind
        elements%hydrobal(i)%bf      = 0.0_rkind
        elements%hydrobal(i)%deltas  = 0.0_rkind
        elements%hydrobal(i)%storage = 0.0_rkind
      end do
    end if

    if (n_steps < 10) stop "Need at least 10 time steps for the hardcoded forcing."

    !------------------------------------------------------------
    ! Daily forcing values
    !------------------------------------------------------------
    pday = [0.0_rkind, 0.0_rkind, 17.0_rkind, 12.0_rkind, 9.0_rkind, &
            7.0_rkind, 40.0_rkind, 0.0_rkind, 0.0_rkind, 3.0_rkind]

    uzday = [4.38_rkind, 3.57_rkind, 4.026_rkind, 3.097_rkind, 4.14_rkind, &
             3.13_rkind, 3.92_rkind, 3.19_rkind, 3.98_rkind, 3.34_rkind]

    tmaxday = [19.1_rkind, 15.3_rkind, 12.8_rkind, 11.8_rkind, 10.5_rkind, &
               15.2_rkind, 11.6_rkind, 14.6_rkind, 17.2_rkind, 16.4_rkind]

    tminday = [5.4_rkind, 6.8_rkind, 8.8_rkind, 7.6_rkind, 8.4_rkind, &
               8.3_rkind, 8.8_rkind, 6.2_rkind, 4.8_rkind, 6.2_rkind]

    tmeanday = [12.25_rkind, 11.05_rkind, 10.8_rkind, 9.7_rkind, 9.45_rkind, &
                11.75_rkind, 10.2_rkind, 10.4_rkind, 11.0_rkind, 9.7_rkind]

    rhmaxday = [84.0_rkind, 85.0_rkind, 76.0_rkind, 87.0_rkind, 92.0_rkind, &
                94.0_rkind, 97.0_rkind, 92.0_rkind, 93.0_rkind, 97.0_rkind]

    rhminday = [56.0_rkind, 64.0_rkind, 64.0_rkind, 77.0_rkind, 77.0_rkind, &
                76.0_rkind, 74.0_rkind, 59.0_rkind, 62.0_rkind, 61.0_rkind]

    soilday = [0.32_rkind, 0.34_rkind, 0.33_rkind, 0.30_rkind, 0.27_rkind, &
           0.24_rkind, 0.22_rkind, 0.23_rkind, 0.26_rkind, 0.29_rkind]

    !------------------------------------------------------------
    ! Fill hourly time series from daily values
    !------------------------------------------------------------
    do i = 1, elements%kolik
      do t = 1, n_steps
        d = int((t - 1) / 24) + 1
        if (d > 10) d = 10

        precip(i,t)      = pday(d) / 24.0_rkind
        uz(i,t)          = uzday(d)
        Tmax(i,t)        = tmaxday(d)
        Tmin(i,t)        = tminday(d)
        Tmean(i,t)       = tmeanday(d)
        RHmax(i,t)       = rhmaxday(d)
        RHmin(i,t)       = rhminday(d)
        soilcontent(i,t) = soilday(d)
      end do
    end do

    !------------------------------------------------------------
    ! Model parameters
    !------------------------------------------------------------
    ! Base saturated conductivity
      conduct = 0.00002_rkind

     ! Different conductivities for different layers [m/s]
     Ksat_surf = 1.0e-6_rkind   ! unsaturated surface, lower effective K
     Ksat_sub  = 5.0e-7_rkind   ! subsurface
     Ksat_gw   = 1.0e-6_rkind    ! saturated groundwater

      ksurf_exp = 3.0_rkind
      ksub_exp  = 2.0_rkind
      G       = 0.0_rkind

      CN         = 98
      z          = 235.0_rkind
      Julian_day = 172
      phi        = 0.614_rkind      ! Latakia/Balloran area ≈ 35.2°N
      as         = 0.25_rkind
      bs         = 0.5_rkind
      alpha      = 0.23_rkind
      sigma      = 4.903e-5_rkind
      gsc        = 0.0820_rkind
      ccrop      = 0.95_rkind

    !------------------------------------------------------------
    ! Optional slope / infiltration parameters
    ! Only keep these if declared in globals.f90
    !------------------------------------------------------------
     cn_slope_coeff      = 0.15_rkind
     storage_slope_coeff = 4.0_rkind
     qgw_slope_coeff     = 0.5_rkind
     min_edge_width      = 1.0e-6_rkind
     theta_r = 0.08_rkind
     theta_s = 0.42_rkind
     Beta1 = 0.0_rkind
     Beta2 = 0.05_rkind
     Beta3 = 1.0_rkind
     Beta4 = 0.03_rkind
     Beta5 = 0.02_rkind
     z1 = 0.0_rkind
     z2 = 1.0_rkind
     infil_slope_coeff   = 8.0_rkind

    print *, "DEBUG precip(1,49) = ", precip(1,49)
    print *, "DEBUG precip(1,72) = ", precip(1,72)

  end subroutine init_hydro

end module hydrotools
