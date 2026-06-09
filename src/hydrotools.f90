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

    integer :: fileid, ios
    integer(kind=ikind) :: i
    real(kind=rkind) :: node_id, element_id

    open(newunit=fileid, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'Error: could not open mesh file.'

    call comment(fileid)
    read(fileid, *, iostat=ios) nodes%kolik
    if (ios /= 0) error stop 'Error: failed to read number of nodes.'

    allocate(nodes%data(nodes%kolik, 2))
    allocate(nodes%altitude(nodes%kolik))

    do i = 1, nodes%kolik
      call comment(fileid)
      read(fileid, *, iostat=ios) node_id, nodes%data(i,1), nodes%data(i,2), nodes%altitude(i)
      if (ios /= 0) error stop 'Error: failed to read node data.'
    end do

    call comment(fileid)
    read(fileid, *, iostat=ios) elements%kolik
    if (ios /= 0) error stop 'Error: failed to read number of elements.'

    call mesh_allocater()

    do i = 1, elements%kolik
      call comment(fileid)
      read(fileid, *, iostat=ios) element_id, elements%data(i,1), elements%data(i,2), elements%data(i,3)
      if (ios /= 0) error stop 'Error: failed to read element data.'
    end do

    close(fileid)
  end subroutine read_mesh


  !==============================================================
  ! ALLOCATE MESH ARRAYS
  !==============================================================
  subroutine mesh_allocater()
    integer(kind=ikind) :: el

    if (.not. allocated(elements%data))       allocate(elements%data(elements%kolik, 3))
    if (.not. allocated(elements%area))       allocate(elements%area(elements%kolik))
    if (.not. allocated(elements%material))   allocate(elements%material(elements%kolik))
    if (.not. allocated(elements%avgalt))     allocate(elements%avgalt(elements%kolik))
    if (.not. allocated(elements%slope))      allocate(elements%slope(elements%kolik))
    if (.not. allocated(elements%overflow))   allocate(elements%overflow(elements%kolik))
    if (.not. allocated(elements%neighbours)) allocate(elements%neighbours(elements%kolik, 3))
    if (.not. allocated(elements%upstream))   allocate(elements%upstream(elements%kolik))
    if (.not. allocated(elements%downstream)) allocate(elements%downstream(elements%kolik))
    if (.not. allocated(elements%hydrobal))   allocate(elements%hydrobal(elements%kolik))

    elements%data       = 0_ikind
    elements%area       = 0.0_rkind
    elements%material   = 0_ikind
    elements%avgalt     = 0.0_rkind
    elements%slope      = 0.0_rkind
    elements%overflow   = 0.0_rkind
    elements%neighbours = 0_ikind

    do el = 1, elements%kolik
      elements%upstream(el)%els    = 0_ikind
      elements%upstream(el)%slopes = 0.0_rkind
      elements%upstream(el)%widths = 0.0_rkind

      elements%downstream(el)%els    = 0_ikind
      elements%downstream(el)%slopes = 0.0_rkind
      elements%downstream(el)%widths = 0.0_rkind

      elements%hydrobal(el)%inflow  = 0.0_rkind
      elements%hydrobal(el)%outflow = 0.0_rkind
      elements%hydrobal(el)%ET      = 0.0_rkind
      elements%hydrobal(el)%q1      = 0.0_rkind
      elements%hydrobal(el)%q2      = 0.0_rkind
      elements%hydrobal(el)%q3      = 0.0_rkind
      elements%hydrobal(el)%pc      = 0.0_rkind
      elements%hydrobal(el)%bf      = 0.0_rkind
      elements%hydrobal(el)%deltas  = 0.0_rkind
      elements%hydrobal(el)%storage = 0.0_rkind
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

    allocate(Pm(elements%kolik, n_steps))
    allocate(E_m(elements%kolik, n_steps))
    allocate(If_m(elements%kolik, n_steps))
    allocate(q1(elements%kolik, n_steps))
    allocate(q2(elements%kolik, n_steps))
    allocate(q3(elements%kolik, n_steps))
    allocate(pc(elements%kolik, n_steps))
    allocate(bf(elements%kolik, n_steps))

    allocate(Ssurf(elements%kolik))
    allocate(Ssub(elements%kolik))
    allocate(Sgw(elements%kolik))

    allocate(dSsurf(elements%kolik, n_steps))
    allocate(dSsub(elements%kolik, n_steps))
    allocate(dSgw(elements%kolik, n_steps))

    allocate(Ssurf_hist(elements%kolik, n_steps))
    allocate(Ssub_hist(elements%kolik, n_steps))
    allocate(Sgw_hist(elements%kolik, n_steps))

    allocate(Qin_result(elements%kolik, n_steps))
    allocate(Qout_result(elements%kolik, n_steps))
    allocate(Overflow_result(elements%kolik, n_steps))
    allocate(Storage_result(elements%kolik, n_steps))
    allocate(deltas(elements%kolik, n_steps))
    allocate(total_deltaS(elements%kolik, n_steps))
    allocate(outlet_Q_m3s(n_steps))

    allocate(storage(elements%kolik))
    allocate(capacity(elements%kolik))
    allocate(outlet_Q(n_steps))

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

    Pm = 0.0_rkind
    E_m = 0.0_rkind
    If_m = 0.0_rkind
    q1 = 0.0_rkind
    q2 = 0.0_rkind
    q3 = 0.0_rkind
    pc = 0.0_rkind
    bf = 0.0_rkind

    Ssurf = 0.0_rkind
    Ssub  = 0.0_rkind
    Sgw   = 0.0_rkind

    dSsurf = 0.0_rkind
    dSsub  = 0.0_rkind
    dSgw   = 0.0_rkind

    Ssurf_hist = 0.0_rkind
    Ssub_hist  = 0.0_rkind
    Sgw_hist   = 0.0_rkind

    Qin_result      = 0.0_rkind
    Qout_result     = 0.0_rkind
    Overflow_result = 0.0_rkind
    Storage_result  = 0.0_rkind
    deltas       = 0.0_rkind
    total_deltaS = 0.0_rkind
    outlet_Q_m3s    = 0.0_rkind

    storage  = 0.0_rkind
    capacity = 8.0_rkind
    outlet_Q = 0.0_rkind

    if (n_steps < 10) stop "Need at least 10 time steps for the hardcoded forcing."

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

    conduct   = 0.00002_rkind
    Ksat_surf = 1.0e-6_rkind
    Ksat_sub  = 5.0e-7_rkind
    Ksat_gw   = 1.0e-6_rkind

    ksurf_exp = 3.0_rkind
    ksub_exp  = 2.0_rkind

    CN         = 98
    z          = 235.0_rkind
    Julian_day = 172
    phi        = 0.614_rkind
    as         = 0.25_rkind
    bs         = 0.5_rkind
    alpha      = 0.23_rkind
    sigma      = 4.903e-5_rkind
    gsc        = 0.0820_rkind
    ccrop      = 0.95_rkind

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

    infil_slope_coeff = 8.0_rkind

  end subroutine init_hydro

end module hydrotools