! main.f90
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



  ! ----------------------------------------------------
  ! 1) Read mesh and compute geometry
  ! ----------------------------------------------------
  call read_mesh(mesh_file_name)
  call compute_areas()
  call compute_avgalt()

  print *, "-----------------------------------------------"
  print *, "Mesh diagnostics"
  print *, "nodes%kolik    = ", nodes%kolik
  print *, "elements%kolik = ", elements%kolik

  if (allocated(nodes%data)) then
    print *, "nodes%data size = ", size(nodes%data,1), size(nodes%data,2)
  else
    print *, "ERROR: nodes%data is not allocated"
  end if

  if (allocated(nodes%altitude)) then
    print *, "nodes%altitude size = ", size(nodes%altitude)
  else
    print *, "ERROR: nodes%altitude is not allocated"
  end if

  if (allocated(elements%area)) then
    print *, "elements%area size = ", size(elements%area)
  else
    print *, "ERROR: elements%area is not allocated"
  end if

  print *, "-----------------------------------------------"
  print *, "Node Coordinates"

  do i = 1, nodes%kolik
    print *, "Node", i, "x=", nodes%data(i,1), "y=", nodes%data(i,2), &
             "z=", nodes%altitude(i)
  end do

  print *, "-----------------------------------------------"
  print *, "Element Areas and average altitude"

  do i = 1, elements%kolik
    print *, "Element", i, "area=", elements%area(i), &
             "avg z=", elements%avgalt(i)
  end do

  ! ----------------------------------------------------
  ! 2) Flow topology
  ! ----------------------------------------------------
  drutes_config%dimen = 2
  call init_flow_topology()
  call print_graph_diagnostics()

  print *, "n_steps  = ", n_steps
  print *, "dt_hours = ", dt_hours
  print *, "dt_days  = ", dt_days

  ! ----------------------------------------------------
  ! 3) Initialize hydrology
  ! ----------------------------------------------------
  call init_hydro()


  ! ----------------------------------------------------
  ! 4) Compute water balance and routing
  ! ----------------------------------------------------
  call compute_all()

  

  ! Print upstream flows only daily
  do t = 1, n_steps
    if (mod(t,24) == 0) then
      call print_upstream_flows(t)
    end if
  end do

  ! ----------------------------------------------------
  ! 5) Export full ODE water balance
  ! ----------------------------------------------------
  open(newunit=unit, file="water_balance_detailed.csv", status="replace", &
       action="write", iostat=ios)

  if (ios /= 0) then
    print *, "ERROR: cannot open water_balance_detailed.csv. IOSTAT=", ios
    stop
  end if

  write(unit,'(A)') &
    "step,element,Pm,If_m,E_m,Qsurf_m,Tv_m,Qsub_m,Rv_m,Qgw_m," // &
    "dVsurf_m3,dVsub_m3,dVgw_m3,Ssurf_m3,Ssub_m3,Sgw_m3"

  do t = 1, n_steps
    do i = 1, elements%kolik
      write(unit,'(I6,",",I6,14(",",ES16.8))') &
        t, i, &
        Pm(i,t), If_m(i,t), E_m(i,t), Qsurf_m(i,t), &
        Tv_m(i,t), Qsub_m(i,t), Rv_m(i,t), Qgw_m(i,t), &
        dVsurf(i,t), dVsub(i,t), dVgw(i,t), &
        Ssurf_hist(i,t), Ssub_hist(i,t), Sgw_hist(i,t)
    end do
  end do

  close(unit)
  print *, "Detailed water balance saved to: water_balance_detailed.csv"

  ! ----------------------------------------------------
  ! 6) Export storage only
  ! ----------------------------------------------------
  open(newunit=unit, file="storage_balance.csv", status="replace", &
       action="write", iostat=ios)

  if (ios /= 0) then
    print *, "ERROR: cannot open storage_balance.csv. IOSTAT=", ios
    stop
  end if

  write(unit,'(A)') "step,element,Ssurf_m3,Ssub_m3,Sgw_m3"

  do t = 1, n_steps
    do i = 1, elements%kolik
      write(unit,'(I6,",",I6,3(",",ES16.8))') &
        t, i, Ssurf_hist(i,t), Ssub_hist(i,t), Sgw_hist(i,t)
    end do
  end do

  close(unit)
  print *, "Storage balance saved to: storage_balance.csv"

  ! ----------------------------------------------------
  ! 7) Export routing balance
  ! ----------------------------------------------------
  call export_element_balance("element_balance.csv")

  ! ----------------------------------------------------
  ! 8) Print daily readable summaries
  ! ----------------------------------------------------
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