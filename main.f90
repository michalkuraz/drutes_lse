! main.f90
program nour
  use typy
  use globals
  use initvals
  use hydrofnc
  use tools
  use geom_tools
  implicit none

   integer(kind=ikind) :: i, nc
   integer :: unit
   character(len=256)  :: mesh_file_name
   logical             :: file_exists
   real(kind=rkind)    :: surplus
   real(kind=rkind)    :: qout_catch   ! catchment outflow part from this element



  ! ----------------------------------------------------
  ! 0) Get mesh file name from command line or default
  ! ----------------------------------------------------
  nc = command_argument_count()
  if (nc >= 1) then
     call get_command_argument(1, mesh_file_name)
  else
     mesh_file_name = "mesh.txt"
  end if

  print *, " Opening mesh file: ", trim(mesh_file_name)
  inquire(file=mesh_file_name, exist=file_exists)
  if (.not. file_exists) stop " Mesh file not found!"

  ! ----------------------------------------------------
  ! 1) Read mesh & compute geometry
  ! ----------------------------------------------------
  call read_mesh(mesh_file_name)
  call compute_areas()
  call compute_avgalt()

  print *, "-----------------------------------------------"
  print *, " Node Coordinates"
  do i = 1, nodes%kolik
     print '(I4,3F10.3)', i, nodes%data(i,1), nodes%data(i,2), nodes%altitude(i)
  end do

  print *, "-----------------------------------------------"
  print *, " Element Areas & avg z"
  do i = 1, elements%kolik
     print '(I4, 2F12.3)', i, elements%area(i), elements%avgalt(i)
  end do

  ! ----------------------------------------------------
  ! 2) Flow topology (neighbours, downstream, order)
  ! ----------------------------------------------------
  drutes_config%dimen = 2
  call init_flow_topology()

  print *, "-----------------------------------------------"
  print *, " Downstream Graph"
  do i = 1, elements%kolik
     if (downstream(i) > 0) then
        write(*,'(A,I3,A,I3)') " Element ", i, " drains to element ", downstream(i)
     else
        write(*,'(A,I3,A)') " Element ", i, " drains to outlet"
     end if
  end do

  ! ----------------------------------------------------
   ! 4) Initialize hydrological inputs and parameters
  !     (from your initvals module)
  call init_hydro()   ! <-- new name

    



  ! ----------------------------------------------------
  ! 5) Compute hydrological balance + routing
  ! ----------------------------------------------------
  call compute_all()

    call print_upstream_flows()

  ! ----------------------------------------------------
  ! 6) Print and export results
  ! ----------------------------------------------------
  print *, "--------------------------------------------------------------------------------------------------------------"
  print *, " Element |     P        ET      Qsurf       Li       Qgw      Qin     Surplus " // &
         " Qout_elem  Qout_catch  Overflow   deltaS"
  print *, "--------------------------------------------------------------------------------------------------------------"

  do i = 1, elements%kolik

     ! Local surplus before routing
     surplus = elements%hydrobal(i)%Qsurf + &
               elements%hydrobal(i)%Li    + &
               elements%hydrobal(i)%Qgw

     ! Part that actually leaves the catchment (only for outlet elements)
     if (downstream(i) == 0_ikind) then
        qout_catch = elements%hydrobal(i)%outflow
     else
        qout_catch = 0.0_rkind
     end if

     write(*,'(I7,11F12.6)') i, &
          precip(i,1), elements%hydrobal(i)%ET, &
          elements%hydrobal(i)%Qsurf, elements%hydrobal(i)%Li, &
          elements%hydrobal(i)%Qgw, elements%hydrobal(i)%inflow, &
          surplus, elements%hydrobal(i)%outflow, &
          qout_catch, elements%overflow(i), &
          elements%hydrobal(i)%deltas
  end do


  print *, "-----------------------------------------------"
  print *, " Element |  z_avg   Qsurf_local"
  do i = 1, elements%kolik
     write(*,'(I7,2F14.6)') i, elements%avgalt(i), Qsurf_result(i,1)
  end do

  call export_element_balance("element_balance.csv")


  print *, "-----------------------------------------------"
  print *, " Catchment outlet hydrograph (per time step):"
  do i = 1, n_days
     print '(A,I3,A,F12.6)', " step ", i, "  Q_out = ", outlet_Q(i)
  end do
    ! Export outlet hydrograph to CSV
  open(newunit=unit, file="outlet_hydrograph.csv", status="replace", action="write")
  write(unit,'(A)') "step,Q_out"

  do i = 1, n_days
     write(unit,'(I4,",",F12.6)') i, outlet_Q(i)
  end do

  close(unit)
  print *, " Outlet hydrograph saved to: outlet_hydrograph.csv"



  print *, "-----------------------------------------------"
  print *, "   Model run completed successfully."
  print *, "   Element balance results saved to element_balance.csv"

end program nour
