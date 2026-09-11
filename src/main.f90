program nour

  use typy
  use globals
  use hydrofnc
  use tools
  use routing
  use hydrotools
  use hydroprint
  use solver

  implicit none

  integer(kind=ikind) :: i, t
  integer :: unit, ios
  character(len=256) :: mesh_file_name

  !==============================================================
  ! FILE NAME
  !==============================================================

  mesh_file_name = "mesh.txt"


  !==============================================================
  ! READ MESH
  !==============================================================

  call read_mesh(mesh_file_name)


  !==============================================================
  ! COMPUTE MESH PROPERTIES
  !==============================================================

  call compute_areas()
  call compute_avgalt()


  !==============================================================
  ! MESH DIAGNOSTICS
  !==============================================================

  print *, "-----------------------------------------------"
  print *, "Mesh diagnostics"
  print *, "nodes%kolik    = ", nodes%kolik
  print *, "elements%kolik = ", elements%kolik
  print *, "-----------------------------------------------"


  do i = 1, nodes%kolik

    print *, "Node", i, &
             "x=", nodes%data(i,1), &
             "y=", nodes%data(i,2), &
             "z=", nodes%altitude(i)

  end do


  print *, "-----------------------------------------------"
  print *, "Element Areas and Average Altitude"


  do i = 1, elements%kolik

    print *, "Element", i, &
             "area=", elements%area(i), &
             "avg z=", elements%avgalt(i)

  end do


  !==============================================================
  ! FLOW TOPOLOGY
  !==============================================================

  drutes_config%dimen = 2

  call init_flow_topology()

  call print_graph_diagnostics()


  !==============================================================
  ! TIME-STEP INFORMATION
  !==============================================================

  print *, "-----------------------------------------------"
  print *, "Time configuration"
  print *, "n_steps  = ", n_steps
  print *, "dt_hours = ", dt_hours
  print *, "dt_days  = ", dt_days
  print *, "-----------------------------------------------"


  if (n_steps <= 0) then
    print *, "ERROR: n_steps must be greater than zero."
    stop
  end if

  if (elements%kolik <= 0) then
    print *, "ERROR: number of elements must be greater than zero."
    stop
  end if


  !==============================================================
  ! INITIALIZE HYDROLOGY
  !
  ! This also reads meteo.csv and fills:
  !
  ! precip
  ! uz
  ! Tmax
  ! Tmin
  ! Tmean
  ! RHmax
  ! RHmin
  ! soilcontent
  !
  !==============================================================

  call init_hydro()


  !==============================================================
  ! RUN HYDROLOGICAL MODEL
  !==============================================================

  print *, "-----------------------------------------------"
  print *, "Starting hydrological simulation..."
  print *, "-----------------------------------------------"

  call compute_all()


  print *, "-----------------------------------------------"
  print *, "Hydrological simulation completed."
  print *, "-----------------------------------------------"


  !==============================================================
  ! UPSTREAM FLOW DIAGNOSTICS
  !
  ! These are printed to the terminal by print_upstream_flows().
  !==============================================================

  do t = 1, n_steps

    if (mod(t, int(24, kind=ikind)) == 0) then

      call print_upstream_flows(t)

    end if

  end do


  !==============================================================
  ! WRITE DETAILED WATER BALANCE CSV
  !==============================================================

  open(newunit=unit, &
       file="water_balance_detailed.csv", &
       status="replace", &
       action="write", &
       iostat=ios)

  if (ios /= 0) then

    print *, "ERROR: cannot open water_balance_detailed.csv."
    print *, "IOSTAT = ", ios
    stop

  end if


  ! CSV HEADER

  write(unit,'(A)') &
       "step,element,P,E,If,q1,q2,q3,pc,bf," // &
       "dSsurf,dSsub,dSgw,total_dS,Ssurf,Ssub,Sgw"


  ! CSV DATA

  do t = 1, n_steps

    do i = 1, elements%kolik

      write(unit,'(2(I0,","),16(ES16.8,","),ES16.8)') &
           t, &
           i, &
           Pm(i,t), &
           E_m(i,t), &
           If_m(i,t), &
           q1(i,t), &
           q2(i,t), &
           q3(i,t), &
           pc(i,t), &
           bf(i,t), &
           dSsurf(i,t), &
           dSsub(i,t), &
           dSgw(i,t), &
           total_deltaS(i,t), &
           Ssurf_hist(i,t), &
           Ssub_hist(i,t), &
           Sgw_hist(i,t)

    end do

  end do


  close(unit)


  print *, "water_balance_detailed.csv written."
  print *, "Records written = ", n_steps * elements%kolik


  !==============================================================
  ! WRITE STORAGE BALANCE CSV
  !==============================================================

  open(newunit=unit, &
       file="storage_balance.csv", &
       status="replace", &
       action="write", &
       iostat=ios)

  if (ios /= 0) then

    print *, "ERROR: cannot open storage_balance.csv."
    print *, "IOSTAT = ", ios
    stop

  end if


  ! CSV HEADER

  write(unit,'(A)') &
       "step,element,Ssurf,Ssub,Sgw"


  ! CSV DATA

  do t = 1, n_steps

    do i = 1, elements%kolik

      write(unit,'(2(I0,","),2(ES16.8,","),ES16.8)') &
           t, &
           i, &
           Ssurf_hist(i,t), &
           Ssub_hist(i,t), &
           Sgw_hist(i,t)

    end do

  end do


  close(unit)


  print *, "storage_balance.csv written."
  print *, "Records written = ", n_steps * elements%kolik


  !==============================================================
  ! WRITE ELEMENT BALANCE
  !==============================================================

  call export_element_balance("element_balance.csv")

  print *, "element_balance.csv written."


  !==============================================================
  ! WATER BALANCE DIAGNOSTICS
  !
  ! These are printed to the terminal.
  !==============================================================

  do t = 1, n_steps

    if (mod(t, int(24, kind=ikind)) == 0) then

      call print_water_balance(t)

    end if

  end do


  !==============================================================
  ! FINAL MESSAGE
  !==============================================================

  print *, "-----------------------------------------------"
  print *, "Model run completed successfully."
  print *, "-----------------------------------------------"

  print *, "Time step [hours] = ", dt_hours
  print *, "Number of steps   = ", n_steps
  print *, "Number of elements = ", elements%kolik

  print *, "Results saved to:"
  print *, "  water_balance_detailed.csv"
  print *, "  storage_balance.csv"
  print *, "  element_balance.csv"

  print *, "-----------------------------------------------"

end program nour