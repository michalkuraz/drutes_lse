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


  !==============================================================
  ! MAIN PROGRAM
  !
  ! 1. Read mesh
  ! 2. Compute mesh properties
  ! 3. Print mesh diagnostics
  ! 4. Initialize flow topology
  ! 5. Write element properties
  ! 6. Write flow graph
  ! 7. Print time configuration
  ! 8. Initialize hydrology
  ! 9. Run hydrological model
  ! 10. Print upstream-flow diagnostics
  ! 11. Write detailed water balance
  ! 12. Write routing results
  ! 13. Write storage balance
  ! 14. Write element balance
  ! 15. Print water-balance diagnostics
  ! 16. Write outlet hydrograph
  !==============================================================


  integer(kind=ikind) :: i
  integer(kind=ikind) :: t
  integer(kind=ikind) :: side

  integer :: unit
  integer :: ios

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
  ! CHECK MESH
  !==============================================================

  if (nodes%kolik <= 0) then
    print *, "ERROR: number of nodes must be greater than zero."
    stop
  end if


  if (elements%kolik <= 0) then
    print *, "ERROR: number of elements must be greater than zero."
    stop
  end if


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


  !--------------------------------------------------------------
  ! Print node coordinates and altitude
  !--------------------------------------------------------------

  do i = 1, nodes%kolik

    print *, "Node", i, &
             " x=", nodes%data(i,1), &
             " y=", nodes%data(i,2), &
             " z=", nodes%altitude(i)

  end do


  !--------------------------------------------------------------
  ! Print element area and average altitude
  !--------------------------------------------------------------

  print *, "-----------------------------------------------"
  print *, "Element Areas and Average Altitude"
  print *, "-----------------------------------------------"


  do i = 1, elements%kolik

    print *, "Element", i, &
             " area=", elements%area(i), &
             " avg z=", elements%avgalt(i)

  end do


  !==============================================================
  ! FLOW TOPOLOGY
  !==============================================================

  drutes_config%dimen = 2


  call init_flow_topology()


  !--------------------------------------------------------------
  ! Print graph diagnostics to terminal
  !--------------------------------------------------------------

  call print_graph_diagnostics()


  !==============================================================
  ! WRITE ELEMENT PROPERTIES CSV
  !
  ! Topology is initialized before this file is written so that
  ! elements%slope is available.
  !==============================================================

  open(newunit=unit, &
       file="element_properties.csv", &
       status="replace", &
       action="write", &
       iostat=ios)


  if (ios /= 0) then

    print *, "ERROR: cannot open element_properties.csv."
    print *, "IOSTAT = ", ios

    stop

  end if


  !--------------------------------------------------------------
  ! CSV HEADER
  !--------------------------------------------------------------

  write(unit,'(A)') &
       "element,area,average_altitude,slope"


  !--------------------------------------------------------------
  ! CSV DATA
  !--------------------------------------------------------------

  do i = 1, elements%kolik

    write(unit,'(I0,",",2(ES16.8,","),ES16.8)') &
         i, &
         elements%area(i), &
         elements%avgalt(i), &
         elements%slope(i)

  end do


  close(unit)


  print *, "element_properties.csv written."
  print *, "Records written = ", elements%kolik


  !==============================================================
  ! WRITE FLOW GRAPH CSV
  !
  ! For each element and each of its three sides:
  !
  ! element
  ! side
  ! downstream element
  ! downstream slope
  ! edge width
  !==============================================================

  open(newunit=unit, &
       file="flow_graph.csv", &
       status="replace", &
       action="write", &
       iostat=ios)


  if (ios /= 0) then

    print *, "ERROR: cannot open flow_graph.csv."
    print *, "IOSTAT = ", ios

    stop

  end if


  !--------------------------------------------------------------
  ! CSV HEADER
  !--------------------------------------------------------------

  write(unit,'(A)') &
       "element,side,downstream_element,slope,width"


  !--------------------------------------------------------------
  ! CSV DATA
  !--------------------------------------------------------------

  do i = 1, elements%kolik

    do side = 1, 3

      write(unit, &
           '(I0,",",I0,",",I0,",",ES16.8,",",ES16.8)') &
           i, &
           side, &
           elements%downstream(i)%els(side), &
           elements%downstream(i)%slopes(side), &
           elements%downstream(i)%widths(side)

    end do

  end do


  close(unit)


  print *, "flow_graph.csv written."
  print *, "Records written = ", 3 * elements%kolik


  !==============================================================
  ! TIME-STEP INFORMATION
  !==============================================================

  print *, "-----------------------------------------------"
  print *, "Time configuration"
  print *, "n_steps    = ", n_steps
  print *, "dt_hours   = ", dt_hours
  print *, "dt_days    = ", dt_days
  print *, "dt_seconds = ", dt_seconds
  print *, "ntot_days  = ", ntot_days
  print *, "-----------------------------------------------"


  if (n_steps <= 0) then

    print *, "ERROR: n_steps must be greater than zero."

    stop

  end if


  !==============================================================
  ! INITIALIZE HYDROLOGY
  !
  ! init_hydro() reads meteorological information and
  ! initializes the hydrological arrays.
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
  ! Current timestep is one hour, therefore every 24 steps
  ! corresponds to one simulation day.
  !==============================================================

  do t = 1, n_steps

    if (mod(t, int(24,kind=ikind)) == 0) then

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


  !--------------------------------------------------------------
  ! CSV HEADER
  !--------------------------------------------------------------

  write(unit,'(A)') &
       "step,element,P,E,If,q1,q2,q3,pc,bf," // &
       "dSsurf,dSsub,dSgw,total_dS,Ssurf,Ssub,Sgw"


  !--------------------------------------------------------------
  ! CSV DATA
  !
  ! 2 integer columns
  ! 15 real columns
  !--------------------------------------------------------------

  do t = 1, n_steps

    do i = 1, elements%kolik

      write(unit, &
           '(2(I0,","),14(ES16.8,","),ES16.8)') &
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
  ! WRITE ROUTING RESULTS CSV
  !
  ! These arrays are filled by route_step(t).
  !==============================================================

  open(newunit=unit, &
       file="routing_results.csv", &
       status="replace", &
       action="write", &
       iostat=ios)


  if (ios /= 0) then

    print *, "ERROR: cannot open routing_results.csv."
    print *, "IOSTAT = ", ios

    stop

  end if


  !--------------------------------------------------------------
  ! CSV HEADER
  !--------------------------------------------------------------

  write(unit,'(A)') &
       "step,element,Qin,Qout,Overflow,Storage,RoutingResidual"


  !--------------------------------------------------------------
  ! CSV DATA
  !--------------------------------------------------------------

  do t = 1, n_steps

    do i = 1, elements%kolik

      write(unit, &
           '(2(I0,","),4(ES16.8,","),ES16.8)') &
           t, &
           i, &
           Qin_result(i,t), &
           Qout_result(i,t), &
           Overflow_result(i,t), &
           Storage_result(i,t), &
           deltas(i,t)

    end do

  end do


  close(unit)


  print *, "routing_results.csv written."
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


  !--------------------------------------------------------------
  ! CSV HEADER
  !--------------------------------------------------------------

  write(unit,'(A)') &
       "step,element,Ssurf,Ssub,Sgw"


  !--------------------------------------------------------------
  ! CSV DATA
  !--------------------------------------------------------------

  do t = 1, n_steps

    do i = 1, elements%kolik

      write(unit, &
           '(2(I0,","),2(ES16.8,","),ES16.8)') &
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
  ! Print once per day for the current one-hour timestep.
  !==============================================================

  do t = 1, n_steps

    if (mod(t, int(24,kind=ikind)) == 0) then

      call print_water_balance(t)

    end if

  end do


  !==============================================================
  ! WRITE OUTLET HYDROGRAPH CSV
  !==============================================================

  open(newunit=unit, &
       file="outlet_hydrograph.csv", &
       status="replace", &
       action="write", &
       iostat=ios)


  if (ios /= 0) then

    print *, "ERROR: cannot open outlet_hydrograph.csv."
    print *, "IOSTAT = ", ios

    stop

  end if


  !--------------------------------------------------------------
  ! CSV HEADER
  !--------------------------------------------------------------

  write(unit,'(A)') &
       "step,time_hours,time_days,outlet_depth_mm,outlet_Q_m3s"


  !--------------------------------------------------------------
  ! CSV DATA
  !--------------------------------------------------------------

  do t = 1, n_steps

    write(unit, &
         '(I0,",",3(ES16.8,","),ES16.8)') &
         t, &
         real(t,rkind) * dt_hours, &
         real(t,rkind) * dt_days, &
         outlet_Q(t), &
         outlet_Q_m3s(t)

  end do


  close(unit)


  print *, "outlet_hydrograph.csv written."
  print *, "Records written = ", n_steps


  !==============================================================
  ! FINAL MESSAGE
  !==============================================================

  print *, "-----------------------------------------------"
  print *, "Model run completed successfully."
  print *, "-----------------------------------------------"

  print *, "Time step [hours]  = ", dt_hours
  print *, "Number of steps    = ", n_steps
  print *, "Number of elements = ", elements%kolik

  print *, "-----------------------------------------------"
  print *, "Results saved to:"
  print *, "  element_properties.csv"
  print *, "  flow_graph.csv"
  print *, "  water_balance_detailed.csv"
  print *, "  routing_results.csv"
  print *, "  outlet_hydrograph.csv"
  print *, "  storage_balance.csv"
  print *, "  element_balance.csv"
  print *, "-----------------------------------------------"


end program nour