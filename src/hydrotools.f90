module hydrotools
 use typy
  use globals
  use tools
   use readtools, only: comment
  

contains


  !==============================================================
  !   READ MESH (nodes + elements)
  !==============================================================
   subroutine read_mesh(filename)
    character(len=*), intent(in) :: filename
    integer :: unit, ios, n_nodes, n_elem
    integer :: id, n1, n2, n3, i
    real(kind=rkind) :: x, y, zloc
    character(len=256) :: line
    logical :: file_exists

    inquire(file=filename, exist=file_exists)
    if (.not. file_exists) then
       print *, " Error: File ", trim(filename), " not found."
       stop
    end if

    open(newunit=unit, file=filename, status="old", action="read", iostat=ios)
    if (ios /= 0) stop " Could not open mesh file."

    ! --- number of nodes ---
    do
       read(unit,'(A)',iostat=ios) line
       if (ios < 0) stop "Unexpected EOF before node count"
       if (len_trim(line) > 0 .and. line(1:1) /= "#") exit
    end do
    read(line,*) n_nodes
    nodes%kolik = n_nodes

    allocate(nodes%data(n_nodes,2))
    allocate(nodes%altitude(n_nodes))

    i = 0
    do while (i < n_nodes)
       read(unit,'(A)',iostat=ios) line
       if (ios < 0) exit
       if (len_trim(line) == 0 .or. line(1:1) == "#") cycle
       ! Expect: id  x   y   z
       read(line,*) id, x, y, zloc
       i = i + 1
       nodes%data(i,1)   = x
       nodes%data(i,2)   = y
       nodes%altitude(i) = zloc
    end do
    print *, " Read", i, "nodes."

    ! --- number of elements ---
    do
       read(unit,'(A)',iostat=ios) line
       if (ios < 0) stop "Unexpected EOF before element count"
       if (len_trim(line) > 0 .and. line(1:1) /= "#") exit
    end do
    read(line,*) n_elem
    elements%kolik = n_elem

    allocate(elements%data(n_elem,3))
    allocate(elements%area(n_elem))
    allocate(elements%material(n_elem))
    allocate(elements%avgalt(n_elem))
    allocate(elements%overflow(n_elem))
    allocate(elements%neighbours(n_elem,3))

    i = 0
    do while (i < n_elem)
       read(unit,'(A)',iostat=ios) line
       if (ios < 0) exit
       if (len_trim(line) == 0 .or. line(1:1) == "#") cycle
       read(line,*) id, n1, n2, n3
       i = i + 1
       elements%data(i,:)   = [n1, n2, n3]
       elements%material(i) = 1_ikind
    end do
    print *, " Read", i, "triangular elements."
    close(unit)

    ! Allocate hydrological structure if needed
    if (.not. allocated(elements%hydrobal)) then
       allocate(elements%hydrobal(n_elem))
    end if
  end subroutine read_mesh


 ! -----------------------------------------------------
  ! initialize and allocate
  ! -----------------------------------------------------
  subroutine init_hydro()
    integer(kind=ikind) :: i

    allocate(precip(elements%kolik, n_steps))
    allocate(qinter(elements%kolik, n_steps))
    allocate(qout(elements%kolik, n_steps))

    allocate(conduct(elements%kolik))
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
    allocate(L_result(elements%kolik, n_steps))
    allocate(Qgw_result(elements%kolik, n_steps))
    allocate(deltas(elements%kolik, n_steps))

    allocate(Qin_result(elements%kolik, n_steps))
    allocate(Qout_result(elements%kolik, n_steps))
    allocate(Overflow_result(elements%kolik, n_steps))
    allocate(Storage_result(elements%kolik, n_steps))
    allocate(outlet_Q_m3s(n_steps))

      Qin_result      = 0.0_rkind
      Qout_result     = 0.0_rkind
      Overflow_result = 0.0_rkind
      Storage_result  = 0.0_rkind
      outlet_Q_m3s    = 0.0_rkind

    if (elements%kolik > 0) then
      if (.not. allocated(elements%hydrobal)) then
        allocate(elements%hydrobal(elements%kolik))
      end if
      do i = 1, elements%kolik
        elements%hydrobal(i)%inflow  = 0.0_rkind
        elements%hydrobal(i)%outflow = 0.0_rkind
        elements%hydrobal(i)%Li      = 0.0_rkind
        elements%hydrobal(i)%ET      = 0.0_rkind
        elements%hydrobal(i)%Qgw     = 0.0_rkind
        elements%hydrobal(i)%Qsurf   = 0.0_rkind
        elements%hydrobal(i)%deltas  = 0.0_rkind
        elements%hydrobal(i)%storage = 0.0_rkind
      end do
    end if

    if (.not. allocated(storage))  allocate(storage(elements%kolik))
    if (.not. allocated(capacity)) allocate(capacity(elements%kolik))
    if (.not. allocated(outlet_Q)) allocate(outlet_Q(n_steps))

    storage      = 0.0_rkind
    outlet_Q     = 0.0_rkind
    capacity     = 5.0_rkind

    precip      = 0.0_rkind
    qinter      = 0.0_rkind
    qout        = 0.0_rkind
    uz          = 0.0_rkind
    Tmax        = 0.0_rkind
    Tmin        = 0.0_rkind
    Tmean       = 0.0_rkind
    RHmax       = 0.0_rkind
    RHmin       = 0.0_rkind
    soilcontent = 0.0_rkind

    Qsurf_result = 0.0_rkind
    ET_flux      = 0.0_rkind
    L_result     = 0.0_rkind
    Qgw_result   = 0.0_rkind
    deltas       = 0.0_rkind

    
     
   

    if (n_steps < 10) stop "Need at least 10 time steps for the hardcoded forcing."

    do i = 1, elements%kolik
      precip(i,1:10) = [ 0.0_rkind, 0.0_rkind, 17.0_rkind, 12.0_rkind, 9.0_rkind, &
          7.0_rkind, 40.0_rkind, 0.0_rkind, 0.0_rkind, 3.0_rkind ]


      qinter(i,1:10) = [0.0_rkind, 0.0_rkind, 0.015_rkind, 0.018_rkind, 0.02_rkind, &
                        0.022_rkind, 0.025_rkind, 0.0235_rkind, 0.03_rkind, 0.031_rkind]

      uz(i,1:10) = [4.38_rkind, 3.57_rkind, 4.026_rkind, 3.097_rkind, 4.14_rkind, &
                    3.13_rkind, 3.92_rkind, 3.19_rkind, 3.98_rkind, 3.34_rkind]

      Tmax(i,1:10) = [19.1_rkind, 15.3_rkind, 12.8_rkind, 11.8_rkind, 10.5_rkind, &
                      15.2_rkind, 11.6_rkind, 14.6_rkind, 17.2_rkind, 16.4_rkind]

      Tmin(i,1:10) = [5.4_rkind, 6.8_rkind, 8.8_rkind, 7.6_rkind, 8.4_rkind, &
                      8.3_rkind, 8.8_rkind, 6.2_rkind, 4.8_rkind, 6.2_rkind]

      Tmean(i,1:10) = [12.25_rkind, 11.05_rkind, 10.8_rkind, 9.7_rkind, 9.45_rkind, &
                       11.75_rkind, 10.2_rkind, 10.4_rkind, 11.0_rkind, 9.7_rkind]

      RHmax(i,1:10) = [84.0_rkind, 85.0_rkind, 76.0_rkind, 87.0_rkind, 92.0_rkind, &
                       94.0_rkind, 97.0_rkind, 92.0_rkind, 93.0_rkind, 97.0_rkind]

      RHmin(i,1:10) = [56.0_rkind, 64.0_rkind, 64.0_rkind, 77.0_rkind, 77.0_rkind, &
                       76.0_rkind, 74.0_rkind, 59.0_rkind, 62.0_rkind, 61.0_rkind]

      soilcontent(i,1:10) = [0.05_rkind, 0.055_rkind, 0.062_rkind, 0.060_rkind, 0.040_rkind, &
           0.070_rkind, 0.090_rkind, 0.200_rkind, 0.250_rkind, 0.265_rkind ]

      
      
    end do

  
   conduct = 0.00002_rkind
    G      = 0.0_rkind
    CN         = 98
    z          = 3.0_rkind
    Julian_day = 172
    phi        = 0.872_rkind
    as         = 0.25_rkind
    bs         = 0.05_rkind
    alpha      = 0.23_rkind
    sigma      = 4.903e-5_rkind
    gsc        = 0.0820_rkind
    ccrop      = 0.8_rkind
  end subroutine init_hydro



end module hydrotools
