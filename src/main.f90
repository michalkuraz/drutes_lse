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

  mesh_file_name = "mesh.txt"

  call read_mesh(mesh_file_name)
  call compute_areas()
  call compute_avgalt()

  print *, "-----------------------------------------------"
  print *, "Mesh diagnostics"
  print *, "nodes%kolik    = ", nodes%kolik
  print *, "elements%kolik = ", elements%kolik
  print *, "-----------------------------------------------"

  do i = 1, nodes%kolik
    print *, "Node", i, "x=", nodes%data(i,1), "y=", nodes%data(i,2), &
             "z=", nodes%altitude(i)
  end do

  print *, "-----------------------------------------------"
  print *, "Element Areas and Average Altitude"

  do i = 1, elements%kolik
    print *, "Element", i, "area=", elements%area(i), &
             "avg z=", elements%avgalt(i)
  end do

  drutes_config%dimen = 2

  call init_flow_topology()
  call print_graph_diagnostics()

  print *, "n_steps  = ", n_steps
  print *, "dt_hours = ", dt_hours
  print *, "dt_days  = ", dt_days

  call init_hydro()
  call compute_all()

  do t = 1, n_steps
    if (mod(t,24) == 0) then
      call print_upstream_flows(t)
    end if
  end do

  open(newunit=unit, file="water_balance_detailed.csv", status="replace", &
       action="write", iostat=ios)

  if (ios /= 0) then
    print *, "ERROR: cannot open water_balance_detailed.csv. IOSTAT=", ios
    stop
  end if

  print *, "step,element,P,E,If,q1,q2,q3,pc,bf,dSsurf,dSsub,dSgw,total_dS,Ssurf,Ssub,Sgw"

  do t = 1, n_steps
    do i = 1, elements%kolik
      print *, &
            t, i, &
           Pm(i,t), E_m(i,t), If_m(i,t), q1(i,t), q2(i,t), q3(i,t), pc(i,t), bf(i,t), &
           dSsurf(i,t), dSsub(i,t), dSgw(i,t), total_deltaS(i,t), &
           Ssurf_hist(i,t), Ssub_hist(i,t), Sgw_hist(i,t)
    end do
  end do

  close(unit)

  open(newunit=unit, file="storage_balance.csv", status="replace", &
       action="write", iostat=ios)

  if (ios /= 0) then
    print *, "ERROR: cannot open storage_balance.csv. IOSTAT=", ios
    stop
  end if

  print *, "step,element,Ssurf,Ssub,Sgw"

  do t = 1, n_steps
    do i = 1, elements%kolik
      print *, &
        t, i, Ssurf_hist(i,t), Ssub_hist(i,t), Sgw_hist(i,t)
    end do
  end do

  close(unit)

  call export_element_balance("element_balance.csv")

  do t = 1, n_steps
    if (mod(t,24) == 0) then
      call print_water_balance(t)
    end if
  end do

  print *, "-----------------------------------------------"
  print *, "Model run completed successfully."
  print *, "Time step [hours] = ", dt_hours
  print *, "Number of steps   = ", n_steps
  print *, "Results saved to:"
  print *, "  water_balance_detailed.csv"
  print *, "  storage_balance.csv"
  print *, "  element_balance.csv"
  print *, "-----------------------------------------------"

end program nour