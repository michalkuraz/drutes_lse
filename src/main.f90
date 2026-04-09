! main.f90
program nour
  use typy
  use globals
  use hydrofnc
  use tools
  implicit none

   integer(kind=ikind) :: i, nc,t
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
  call read_mesh(mesh_file_name)   ! from tools
  call compute_areas()             ! from tools
  call compute_avgalt()            ! from tools

  print *, "-----------------------------------------------"
  print *, " Node Coordinates"
  do i = 1, nodes%kolik
     print *, i, nodes%data(i,1), nodes%data(i,2), nodes%altitude(i)
  end do

  print *, "-----------------------------------------------"
  print *, " Element Areas & avg z"
  do i = 1, elements%kolik
     print *, i, elements%area(i), elements%avgalt(i)
  end do

  ! ----------------------------------------------------
  ! 2) Flow topology (neighbours, downstream, order)
  ! ----------------------------------------------------
  drutes_config%dimen = 2
  call init_flow_topology()          !from tools
  call print_graph_diagnostics()     ! from tools
  

  print *, "-----------------------------------------------"
  print *, " Downstream Graph"
  do i = 1, elements%kolik
     if (downstream(i) > 0) then
        print*, " Element ", i, " drains to element ", downstream(i)
     else
        print*, " Element ", i, " drains to outlet"
     end if
  end do

  ! ----------------------------------------------------
   ! 3) Initialize hydrological inputs and parameters
  !     (from your initvals module)
  ! ---------------------------------------------------
  call init_hydro()    ! from tools

    
  ! ----------------------------------------------------
  ! 4) Compute hydrological balance + routing
  ! ----------------------------------------------------
  call compute_all()       ! from hydrofunction

  do t = 1, n_steps
   call print_upstream_flows(t)
end do    ! from tools

  ! ----------------------------------------------------
  ! 5) Print and export results
  ! ----------------------------------------------------
    print *, "--------------------------------------------------------------------------------"
print *, " Step Elem     P       ET      Qsurf      Li       Qgw       Qin"
print *, "           Qout   Overflow   Storage   deltaS"
print *, "--------------------------------------------------------------------------------"

do t = 1, n_steps
  print *, " "
  print *, "==================== TIME STEP ", t, " ===================="

  do i = 1, elements%kolik
     print *, t, i, &
              precip(i,t), ET_flux(i,t), Qsurf_result(i,t), &
              L_result(i,t), Qgw_result(i,t), Qin_result(i,t), &
              Qout_result(i,t), Overflow_result(i,t), &
              Storage_result(i,t), deltas(i,t)
  end do
end do

 call export_element_balance("element_balance.csv")    ! from tools

  print *, "-----------------------------------------------"
  print *, " Element |  z_avg   Qsurf_local"
  do i = 1, elements%kolik
   print *, i, elements%avgalt(i), Qsurf_result(i,n_steps)
  end do

  


  print *, "-----------------------------------------------"
print *, " Catchment outlet hydrograph"
print *, " Step    Q_out(mm/step)    Q_out(m3/s)"

    do t = 1, n_steps
      print *, t, outlet_Q(t), outlet_Q_m3s(t)
    end do
    ! Export outlet hydrograph to CSV
  open(newunit=unit, file="outlet_hydrograph.csv", status="replace", action="write")
   write(unit,'(A)') "step,Q_out_mm_per_step,Q_out_m3s"

   

  close(unit)
  print *, "Outlet hydrograph saved to: outlet_hydrograph.csv"

  




  print *, "-----------------------------------------------"
  print *, "   Model run completed successfully."
  print *, "   Element balance results saved to element_balance.csv"

end program nour
