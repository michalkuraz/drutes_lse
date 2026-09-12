  
module hydrotools
  use typy
  use globals
  use tools
  use readtools, only: comment
  implicit none

contains

  !==============================================================
     ! READ MESH
     ! ALLOCATE MESH ARRAYS
     ! read_meteodata
     ! INITIALIZE HYDROLOGY (init_hydro)
       !call read_meteodata_csv("meteo.csv")
  !==============================================================








  !==============================================================
  ! READ MESH
  !==============================================================
  subroutine read_mesh(filename)
    character(len=*), intent(in) :: filename

    integer :: fileid, ios
    integer :: i
    real :: node_id, element_id

    open(unit=fileid, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) stop 'Error: could not open mesh file.'

    call comment(fileid)
    read(fileid, *, iostat=ios) nodes%kolik
    if (ios /= 0) stop 'Error: failed to read number of nodes.'

    allocate(nodes%data(nodes%kolik, 2))
    allocate(nodes%altitude(nodes%kolik))

    do i = 1, nodes%kolik
      call comment(fileid)
      read(fileid, *, iostat=ios) node_id, nodes%data(i,1), nodes%data(i,2), nodes%altitude(i)
      if (ios /= 0) stop 'Error: failed to read node data.'
    end do

    call comment(fileid)
    read(fileid, *, iostat=ios) elements%kolik
    if (ios /= 0) stop 'Error: failed to read number of elements.'

    call mesh_allocater()

    do i = 1, elements%kolik
      call comment(fileid)
      read(fileid, *, iostat=ios) element_id, elements%data(i,1), elements%data(i,2), elements%data(i,3)
      if (ios /= 0) stop 'Error: failed to read element data.'
    end do

    close(fileid)
  end subroutine read_mesh



  !==============================================================
  ! ALLOCATE MESH ARRAYS
  !==============================================================
  subroutine mesh_allocater()
    integer :: el

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

    elements%data       = 0
    elements%area       = 0.0
    elements%material   = 0
    elements%avgalt     = 0.0
    elements%slope      = 0.0
    elements%overflow   = 0.0
    elements%neighbours = 0

    do el = 1, elements%kolik
      elements%upstream(el)%els    = 0
      elements%upstream(el)%slopes = 0.0
      elements%upstream(el)%widths = 0.0

      elements%downstream(el)%els    = 0
      elements%downstream(el)%slopes = 0.0
      elements%downstream(el)%widths = 0.0

      elements%hydrobal(el)%inflow  = 0.0
      elements%hydrobal(el)%outflow = 0.0
      elements%hydrobal(el)%ET      = 0.0
      elements%hydrobal(el)%q1      = 0.0
      elements%hydrobal(el)%q2      = 0.0
      elements%hydrobal(el)%q3      = 0.0
      elements%hydrobal(el)%pc      = 0.0
      elements%hydrobal(el)%bf      = 0.0
      elements%hydrobal(el)%deltas  = 0.0
      elements%hydrobal(el)%storage = 0.0
    end do
  end subroutine mesh_allocater

  !==============================================================
  ! read_meteodata
  !==============================================================
  subroutine read_meteodata_csv(filename)
  character(len=*), intent(in) :: filename

  integer :: unit, ios
  character(len=1024) :: header
  real(kind=rkind) :: rowdata(8)

  if (.not. allocated(meteodata)) allocate(meteodata(1))

  open(unit=unit, file=filename, status='old', action='read', iostat=ios)

  if (ios /= 0) then
    print *, "ERROR: cannot open meteorological CSV file: ", trim(filename)
    stop
  end if

  ! Skip header line
  read(unit, '(A)', iostat=ios) header

  do
    read(unit, *, iostat=ios) rowdata
    if (ios /= 0) exit

    call meteodata(1)%time%fill(rowdata(1))
    call meteodata(1)%rainfall%fill(rowdata(2))
    call meteodata(1)%Tmax%fill(rowdata(3))
    call meteodata(1)%Tmin%fill(rowdata(4))
    call meteodata(1)%wind%fill(rowdata(5))
    call meteodata(1)%RHmax%fill(rowdata(6))
    call meteodata(1)%RHmin%fill(rowdata(7))
    call meteodata(1)%soilcontent%fill(rowdata(8))
  end do

  close(unit)

  if (meteodata(1)%time%pos < 1) then
    print *, "ERROR: meteorological file is empty."
    stop
  end if

  print *, "Meteo records read = ", meteodata(1)%time%pos

end subroutine read_meteodata_csv

  !==============================================================
  ! INITIALIZE HYDROLOGY
  !==============================================================
  subroutine init_hydro()
    integer :: i, t, d

    !real, dimension(10) :: pday
    !real, dimension(10) :: uzday
    !real, dimension(10) :: tmaxday
    !real, dimension(10) :: tminday
    !real, dimension(10) :: tmeanday
    !real, dimension(10) :: rhmaxday
    !real, dimension(10) :: rhminday
    !real, dimension(10) :: soilday

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

    precip      = 0.0
    qinter      = 0.0
    qout        = 0.0
    conduct     = 0.0
    G           = 0.0

    Tmax        = 0.0
    Tmin        = 0.0
    Tmean       = 0.0
    RHmax       = 0.0
    RHmin       = 0.0
    uz          = 0.0
    soilcontent = 0.0

    Pm = 0.0
    E_m = 0.0
    If_m = 0.0
    q1 = 0.0
    q2 = 0.0
    q3 = 0.0
    pc = 0.0
    bf = 0.0

    Ssurf = 0.0
    Ssub  = 0.0
    Sgw   = 0.0

    dSsurf = 0.0
    dSsub  = 0.0
    dSgw   = 0.0

    Ssurf_hist = 0.0
    Ssub_hist  = 0.0
    Sgw_hist   = 0.0

    Qin_result      = 0.0
    Qout_result     = 0.0
    Overflow_result = 0.0
    Storage_result  = 0.0
    deltas       = 0.0
    total_deltaS = 0.0
    outlet_Q_m3s    = 0.0

    storage  = 0.0
    capacity = 8.0
    outlet_Q = 0.0

    !if (n_steps < 10) stop "Need at least 10 time steps for the hardcoded forcing."

    !pday = [0.0, 0.0, 17.0, 12.0, 9.0, &
    !        7.0, 40.0, 0.0, 0.0, 3.0]

    !uzday = [4.38, 3.57, 4.026, 3.097, 4.14, &
    !         3.13, 3.92, 3.19, 3.98, 3.34]

    !tmaxday = [19.1, 15.3, 12.8, 11.8, 10.5, &
    !           15.2, 11.6, 14.6, 17.2, 16.4]

    !tminday = [5.4, 6.8, 8.8, 7.6, 8.4, &
    !           8.3, 8.8, 6.2, 4.8, 6.2]

    !tmeanday = [12.25, 11.05, 10.8, 9.7, 9.45, &
    !            11.75, 10.2, 10.4, 11.0, 9.7]

    !rhmaxday = [84.0, 85.0, 76.0, 87.0, 92.0, &
               !94.0, 97.0, 92.0, 93.0, 97.0]

    !rhminday = [56.0, 64.0, 64.0, 77.0, 77.0, &
    !            76.0, 74.0, 59.0, 62.0, 61.0]

    !soilday = [0.32, 0.34, 0.33, 0.30, 0.27, &
               !0.24, 0.22, 0.23, 0.26, 0.29]

    !do i = 1, elements%kolik
     ! do t = 1, n_steps
        !d = int((t - 1) / 24) + 1
       ! if (d > 10) d = 10

        !precip(i,t)      = pday(d) / 24.0
       ! uz(i,t)          = uzday(d)
       ! Tmax(i,t)        = tmaxday(d)
        !Tmin(i,t)        = tminday(d)
        !Tmean(i,t)       = tmeanday(d)
       ! RHmax(i,t)       = rhmaxday(d)
        !RHmin(i,t)       = rhminday(d)
        !soilcontent(i,t) = soilday(d)
      !end do
    !end do

  call read_meteodata_csv("meteo.csv")

     do i = 1, elements%kolik
       do t = 1, n_steps
          d = int((t - 1) / 24) + 1

         if (d > meteodata(1)%time%pos) then
          d = meteodata(1)%time%pos
          end if

    precip(i,t)      = meteodata(1)%rainfall%data(d) / 24.0
    uz(i,t)          = meteodata(1)%wind%data(d)
    Tmax(i,t)        = meteodata(1)%Tmax%data(d)
    Tmin(i,t)        = meteodata(1)%Tmin%data(d)
    Tmean(i,t)       = 0.5 * (Tmax(i,t) + Tmin(i,t))
    RHmax(i,t)       = meteodata(1)%RHmax%data(d)
    RHmin(i,t)       = meteodata(1)%RHmin%data(d)
    soilcontent(i,t) = meteodata(1)%soilcontent%data(d)
    end do
      end do

    conduct   = 0.00002
    Ksat_surf = 1.0e-6
    Ksat_sub  = 5.0e-7
    Ksat_gw   = 1.0e-6

    ksurf_exp = 3.0
    ksub_exp  = 2.0

    CN         = 98
    z          = 235.0
    Julian_day = 172
    phi        = 0.614
    as         = 0.25
    bs         = 0.5
    alpha      = 0.23
    sigma      = 4.903e-5
    gsc        = 0.0820
    ccrop      = 0.95

    cn_slope_coeff      = 0.15
    storage_slope_coeff = 4.0
    qgw_slope_coeff     = 0.5
    min_edge_width      = 1.0e-6

    theta_r = 0.08
    theta_s = 0.42

    Beta1 = 0.0
    Beta2 = 0.05
    Beta3 = 1.0
    Beta4 = 0.03
    Beta5 = 0.02

    z1 = 0.0
    z2 = 1.0

    infil_slope_coeff = 8.0

  end subroutine init_hydro

end module hydrotools