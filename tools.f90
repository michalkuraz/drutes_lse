! tools.f90
module tools
  use typy
  use globals
  use geom_tools      ! for find_neighbours
  implicit none
contains


  function avg(r) result(avgval)
    use typy
    real(kind=rkind), dimension(:), intent(in) :: r
    real(kind=rkind) :: avgval
    
    integer(kind=ikind) :: i
    
    avgval = 0
    do i=1, ubound(r,1)
      avgval = avgval + r(i)
    end do
    
    avgval = avgval/ubound(r,1)
  end function avg
    

  ! -----------------------------------------------------
  ! Read mesh file: id, x, y, z for nodes; triangles for elements
  ! -----------------------------------------------------
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
    allocate(elements%upstream(n_elem))
    allocate(elements%downstream(n_elem))

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
  ! Compute triangle areas + total area
  ! -----------------------------------------------------
  subroutine compute_areas()
    integer(kind=ikind) :: i
    integer(kind=ikind) :: n1, n2, n3
    real(kind=rkind) :: a(2), b(2), c(2)
    real(kind=rkind) :: total_area

    total_area = 0.0_rkind

    do i = 1_ikind, elements%kolik
       n1 = elements%data(i,1)
       n2 = elements%data(i,2)
       n3 = elements%data(i,3)
       a  = nodes%data(n1,:)
       b  = nodes%data(n2,:)
       c  = nodes%data(n3,:)
       elements%area(i) = area_triangle(a,b,c)
       total_area       = total_area + elements%area(i)
    end do

    print "(A,F12.3)", " Total mesh area = ", total_area
  end subroutine compute_areas

  pure function area_triangle(a,b,c) result(area)
    real(kind=rkind), intent(in) :: a(2), b(2), c(2)
    real(kind=rkind) :: area
    area = 0.5_rkind * abs(a(1)*(b(2)-c(2)) + b(1)*(c(2)-a(2)) + c(1)*(a(2)-b(2)))
  end function area_triangle

  ! -----------------------------------------------------
  ! Compute mean elevation per element
  ! -----------------------------------------------------
  subroutine compute_avgalt()
    integer(kind=ikind) :: i, n1, n2, n3
    real(kind=rkind) :: z1, z2, z3

    do i = 1_ikind, elements%kolik
       n1 = elements%data(i,1)
       n2 = elements%data(i,2)
       n3 = elements%data(i,3)

       if (allocated(nodes%altitude)) then
          z1 = nodes%altitude(n1)
          z2 = nodes%altitude(n2)
          z3 = nodes%altitude(n3)
       else
          ! fallback: use y as pseudo-elevation
          z1 = nodes%data(n1,2)
          z2 = nodes%data(n2,2)
          z3 = nodes%data(n3,2)
       end if

       elements%avgalt(i) = (z1 + z2 + z3) / 3.0_rkind
    end do
  end subroutine compute_avgalt

  ! -----------------------------------------------------
  ! Initialise flow topology: neighbours, avgalt, downstream, order
  ! -----------------------------------------------------
  subroutine init_flow_topology()
    integer(kind=ikind) :: nel

    nel = elements%kolik
    if (nel <= 0_ikind) return

    ! Neighbour table
    call find_neighbours(elements, nodes)

    ! Average elevation already computed by compute_avgalt()

    ! Allocate downstream + flow order
    if (.not. allocated(downstream))  allocate(downstream(nel))
    if (.not. allocated(flow_order))  allocate(flow_order(nel))
    

    
    allocate(elements%outlet(nel))
    allocate(elements%watershed(nel))
    
    allocate(elements%ndwatershed(nel,2))
    allocate(elements%ndoutlet(nel,2))
    
    elements%ndwatershed = 0
    elements%ndoutlet = 0
    
    elements%outlet = .false.
    elements%watershed = .false.

    elements%overflow = 0.0_rkind
    downstream        = 0_ikind
    flow_order        = 0_ikind

    !call build_downstream_graph()
    call build_graph()
    stop
    call build_flow_order()

    ! NEW: build upstream lists
    call build_upstream_graph()
  end subroutine init_flow_topology
  

  ! -----------------------------------------------------
  ! Downstream element = neighbour with lowest avgalt
  ! -----------------------------------------------------
  subroutine build_downstream_graph()
    integer(kind=ikind) :: i, j, nb, nel
    real(kind=rkind)    :: myz, bestz
    integer(kind=ikind) :: best_nb

    nel = elements%kolik

    do i = 1, nel
       myz     = elements%avgalt(i)
       bestz   = myz
       best_nb = 0

       nodes: do j = 1, ubound(elements%neighbours,2)
          nb = elements%neighbours(i,j)
          if (nb /= 0) then
            if (elements%avgalt(nb) < bestz) then
               bestz   = elements%avgalt(nb)
               best_nb = nb
            end if
          end if
       end do nodes

       downstream(i) = best_nb  ! 0 if no lower neighbour → outlet
       
    end do
    
  end subroutine build_downstream_graph
  
  function setlines(el, nd) result(pts)
    integer(kind=ikind), intent(in) :: el, nd
    integer(kind=ikind), dimension(2) :: pts
    
    integer(kind=ikind), dimension(3) :: neighnds, elnds
    integer(kind=ikind) :: pos, i, j 
    
    elnds = elements%data(el,:)
    neighnds = elements%data(nd,:)
    
    pos = 0
    
    do i=1, 3
      do j=1,3
        if (neighnds(i) == elnds(j)) then
          pos = pos + 1
          if (pos < 3) then
            pts(pos) = neighnds(i)
          else
            print *, "exited from topotools::setlines"
            print *, "mesh is strange, you have two elements with equal nodes"
            ERROR STOP
          end if
        end if
      end do
    end do
  
  end function setlines
  
  function matchlines(l1, l2) result(true)
    use typy
    
    integer(kind=ikind), dimension(2) :: l1, l2
    logical :: true
    
    true = .false.
    
    if (l1(1) == l2(1) .and. l1(2) == l2(2) ) then
      true = .true.
      return
    end if
    
    if (l1(1) == l2(2) .and. l1(2) == l2(1) ) then
      true = .true.
      return
    end if
  
  end function matchlines
  
  function getcenter(el) result(xy)
    use typy
    use globals
    
    integer(kind=ikind), intent(in) :: el
    real(kind=rkind), dimension(2) :: xy
    
    xy(1) = avg(nodes%data(elements%data(el,:),1))
    xy(2) = avg(nodes%data(elements%data(el,:),2))
    
  end function getcenter
  
  function dist(A,B) result(l)
    use typy
    use globals
    
    real(kind=rkind), dimension(:), intent(in) :: A,B
    real(kind=rkind) :: l


    l = sqrt((A(1)-B(1))*(A(1)-B(1)) + (A(2)-B(2))*(A(2)-B(2)))


  end function dist
  
  subroutine build_graph()
    use typy
    use globals
    use debug_tools
    
    integer(kind=ikind) :: el, i
    integer(kind=ikind), dimension(3) :: ngh, nghalt
    real(kind=rkind) :: myalt
    real(kind=rkind), dimension(2) :: myel, ngel
    integer(kind=ikind), dimension(3,2) :: nghlines
    
    do el=1, elements%kolik
      myalt = elements%avgalt(el)
      ngh = elements%neighbours(el,:)
      myel = getcenter(el)
      do i=1,3
        if (ngh(i) /= 0) then
          nghlines(i,:) = setlines(ngh(i), el)
          nghalt(i) = elements%avgalt(ngh(i))
        else
          nghalt(i) = -9999
        end if
      end do
      
      do i=1,3
        if (ngh(i) /= 0) then
          ngel = getcenter(ngh(i))
          
          if (nghalt(i) <= myalt) then
            elements%downstream(el)%els(i) = ngh(i)

            elements%downstream(el)%slopes(i) = (myalt - nghalt(i))/dist(myel, ngel)
            elements%downstream(el)%widths(i) = dist(nodes%data(nghlines(i,1),:), nodes%data(nghlines(i,2),:))
            print *, "down w", elements%downstream(el)%widths(i), el, ngh(i)
          else
            elements%upstream(el)%els(i) = ngh(i)
            elements%upstream(el)%slopes(i) = (myalt - nghalt(i))/dist(myel, ngel)
            elements%upstream(el)%widths(i) = dist(nodes%data(nghlines(i,1),:), nodes%data(nghlines(i,2),:))
            print *, "up w", elements%upstream(el)%widths(i), el, ngh(i)
          end if
        end if
        call wait()
      end do
      
      
    end do
    
  
  end subroutine build_graph
  
  subroutine build_downstream_graph2
    use debug_tools
    
    integer(kind=ikind) :: el, ndmin
    integer(kind=ikind), dimension(3) :: ngh, nghalt
    integer(kind=ikind), dimension(3,2) :: nghlines
    integer(kind=ikind), dimension(1) :: ndminloc
    integer(kind=ikind) :: j, up, down
    type :: elalt_str
      integer(kind=ikind) :: nd1, nd2
      real(kind=rkind)    :: alt
    end type elalt_str
    type(elalt_str), dimension(:), allocatable :: elalt
    integer(kind=ikind), dimension(1) :: inloc, outloc
    

    allocate(elalt(elements%kolik))
    
    do el=1, elements%kolik
      ngh = elements%neighbours(el,:)
      do j=1,3
        select case(j)
          case(1)
            elalt(j)%nd1 = elements%data(el,1)
            elalt(j)%nd2 = elements%data(el,2)
          case(2)
            elalt(j)%nd1 = elements%data(el,2)
            elalt(j)%nd2 = elements%data(el,3)
          case(3)
            elalt(j)%nd1 = elements%data(el,3)
            elalt(j)%nd2 = elements%data(el,1)
        end select

        elalt(j)%alt= avg([nodes%altitude(elalt(j)%nd1), nodes%altitude(elalt(j)%nd2)])

      end do
      
      inloc = maxloc(elalt(:)%alt,1)
      outloc = minloc(elalt(:)%alt,1)
      
      print *, el, "in", elalt(inloc(1))%nd1,  elalt(inloc(1))%nd2
      print *, el, "out", elalt(outloc(1))%nd1,  elalt(outloc(1))%nd2
      print *, "-------"
      
      do j=1,3 
        if (ngh(j) > 0) then
          nghlines(j,:) = setlines(ngh(j), el)
!          print *, "el:", el
!          print *, "ngh:", ngh(j)
!          print *, "lines", nghlines(j,:)
!          print *, "-------"
!          call wait()
        end if
        
      end do
      
      
    end do
      
      
      
      
      
    stop
      
      
  
  end subroutine build_downstream_graph2

  ! -----------------------------------------------------
  ! Flow order: sort elements from highest to lowest avgalt
  ! -----------------------------------------------------
  subroutine build_flow_order()
    integer(kind=ikind) :: nel, i, j, tmp_idx
    real(kind=rkind)    :: tmp_z
    integer(kind=ikind), allocatable :: idx(:)
    real(kind=rkind),    allocatable :: z(:)

    nel = elements%kolik
    allocate(idx(nel), z(nel))

    do i = 1_ikind, nel
       idx(i) = i
       z(i)   = elements%avgalt(i)
    end do

    ! simple selection sort: descending z
    do i = 1_ikind, nel-1
       do j = i+1_ikind, nel
          if (z(j) > z(i)) then
             tmp_z   = z(i);   z(i)   = z(j);   z(j)   = tmp_z
             tmp_idx = idx(i); idx(i) = idx(j); idx(j) = tmp_idx
          end if
       end do
    end do

    flow_order = idx

    deallocate(idx, z)
  end subroutine build_flow_order


    ! -----------------------------------------------------
  ! Build upstream lists from downstream(:)
  ! upstream_count(j)  = how many elements flow into j
  ! upstream_list(j,k) = k-th upstream element of j
  ! -----------------------------------------------------
  subroutine build_upstream_graph()
    use typy
    use globals
    implicit none

    integer(kind=ikind) :: nel, i, j, max_up

    nel = elements%kolik
    if (nel <= 0_ikind) return

    !--- pass 1: count how many upstream elements each cell has
    if (allocated(upstream_count)) deallocate(upstream_count)
    allocate(upstream_count(nel))
    upstream_count = 0_ikind

    do i = 1_ikind, nel
       j = downstream(i)
       if (j > 0_ikind) then
          upstream_count(j) = upstream_count(j) + 1_ikind
       end if
    end do

    ! maximum number of direct upstream contributors
    max_up = 0_ikind
    do i = 1_ikind, nel
       if (upstream_count(i) > max_up) max_up = upstream_count(i)
    end do

    !--- pass 2: allocate and fill upstream_list(:,:)
    if (allocated(upstream_list)) deallocate(upstream_list)
    if (max_up > 0_ikind) then
       allocate(upstream_list(nel, max_up))
       upstream_list = 0_ikind

       ! reuse upstream_count as a position counter
       upstream_count = 0_ikind

       do i = 1_ikind, nel
          j = downstream(i)
          if (j > 0_ikind) then
             upstream_count(j) = upstream_count(j) + 1_ikind
             upstream_list(j, upstream_count(j)) = i
          end if
       end do
    else
       ! no upstream connections at all (e.g. single-element mesh)
       if (allocated(upstream_list)) deallocate(upstream_list)
    end if

  end subroutine build_upstream_graph


    ! -----------------------------------------------------
  ! Print inflow/outflow per element based on upstream list
  ! -----------------------------------------------------
  subroutine print_upstream_flows()
    use typy
    use globals
    implicit none

    integer(kind=ikind) :: nel, e, k
    real(kind=rkind)    :: Qin_from_up

    nel = elements%kolik

    if (.not. allocated(elements%hydrobal)) then
       print *, " No hydrological data (hydrobal not allocated)."
       return
    end if
    if (.not. allocated(upstream_count)) then
       print *, " Upstream graph not built yet (call init_flow_topology first)."
       return
    end if

    print *, "---------------------------------------------------------------"
    print *, " el | n_up | Qin_from_up    Qin(field)     Qout_elem"
    print *, "---------------------------------------------------------------"

    do e = 1_ikind, nel
       Qin_from_up = 0.0_rkind

       if (allocated(upstream_list)) then
          do k = 1_ikind, upstream_count(e)
             Qin_from_up = Qin_from_up + &
                  elements%hydrobal( upstream_list(e,k) )%outflow
          end do
       end if

       write(*,'(I3,1X,I4,3F14.6)') e, upstream_count(e), Qin_from_up, &
            elements%hydrobal(e)%inflow, elements%hydrobal(e)%outflow
    end do

  end subroutine print_upstream_flows


  ! -----------------------------------------------------
  ! Optional: print element balance (for debugging)
  ! -----------------------------------------------------
  subroutine print_element_balance()
    integer(kind=ikind) :: i
    real(kind=rkind)    :: surplus

    if (.not. allocated(elements%hydrobal)) then
       print *, " No hydrological data to display."
       return
    end if

    print *, "--------------------------------------------------------------------------------------------------------------"
    print *, " Element |     P        ET      Qsurf       Li       Qgw      Qin     Surplus   Qout_elem  Overflow    deltaS"
    print *, "--------------------------------------------------------------------------------------------------------------"

    do i = 1_ikind, elements%kolik
       surplus = elements%hydrobal(i)%Qsurf + &
                 elements%hydrobal(i)%Li    + &
                 elements%hydrobal(i)%Qgw

       write(*,'(I7,9F12.6)') i, precip(i,1), elements%hydrobal(i)%ET, &
            elements%hydrobal(i)%Qsurf, elements%hydrobal(i)%Li, &
            elements%hydrobal(i)%Qgw, elements%hydrobal(i)%inflow, &
            surplus, elements%hydrobal(i)%outflow, &
            elements%hydrobal(i)%deltas
    end do
  end subroutine print_element_balance

  ! -----------------------------------------------------
  ! Export element balance to CSV, with Qout_elem & Qout_catchment
  ! -----------------------------------------------------
     subroutine export_element_balance(filename)
    character(len=*), intent(in) :: filename
    integer(kind=ikind) :: i
    integer :: unit
    real(kind=rkind) :: surplus, qout_catch

    if (.not. allocated(elements%hydrobal)) then
       print *, " No hydrological data to export."
       return
    end if

    open(newunit=unit, file=filename, status="replace", action="write")

    ! Header now includes:
    ! Area, z_avg, Downstream, and Qsurf_local
    write(unit,'(A)') "Element,Area,z_avg,Downstream," // &
                      "P,ET,Qsurf,Li,Qgw,Qin,Surplus," // &
                      "Qout_elem,Qout_catchment,Overflow,DeltaS,Qsurf_local"

    do i = 1_ikind, elements%kolik

       ! Local surplus (production term)
       surplus = elements%hydrobal(i)%Qsurf + &
                 elements%hydrobal(i)%Li    + &
                 elements%hydrobal(i)%Qgw

       ! Only outlet elements contribute to catchment outflow
       if (downstream(i) == 0_ikind) then
          qout_catch = elements%hydrobal(i)%outflow
       else
          qout_catch = 0.0_rkind
       end if

       write(unit,'(I4, ",", F10.3, ",", F10.3, ",", I4, 12(",",F10.3))')  &
            i,                              & ! Element
            elements%area(i),               & ! Area
            elements%avgalt(i),             & ! z_avg
            downstream(i),                  & ! Downstream (0 = outlet)
            precip(i,1),                    & ! P
            elements%hydrobal(i)%ET,        & ! ET
            elements%hydrobal(i)%Qsurf,     & ! Qsurf
            elements%hydrobal(i)%Li,        & ! Li
            elements%hydrobal(i)%Qgw,       & ! Qgw
            elements%hydrobal(i)%inflow,    & ! Qin
            surplus,                        & ! Surplus
            elements%hydrobal(i)%outflow,   & ! Qout_elem
            qout_catch,                     & ! Qout_catchment
            elements%overflow(i),           & ! Overflow
            elements%hydrobal(i)%deltas,    & ! DeltaS
            Qsurf_result(i,1)                 ! Qsurf_local
    end do

    close(unit)
    print *, " Element balances saved to: ", trim(filename)
  end subroutine export_element_balance



   ! -----------------------------------------------------
  ! Routing + local mass balance for one time step
  ! -----------------------------------------------------
  subroutine route_step(tstep)
    integer(kind=ikind), intent(in) :: tstep

    integer(kind=ikind) :: el, i, dwn
    real(kind=rkind)    :: old_storage, local_input, local_losses, surplus
    real(kind=rkind), allocatable :: overflow_total(:), storage_new(:)

    ! Allocate helper arrays
    allocate(overflow_total(elements%kolik))
    allocate(storage_new(elements%kolik))

    ! Reset inflow/outflow at start of step
    elements%hydrobal(:)%inflow  = 0.0_rkind
    elements%hydrobal(:)%outflow = 0.0_rkind

    overflow_total = 0.0_rkind
    storage_new    = storage

    ! 1) Local balance before routing (element-wise)
    do el = 1_ikind, elements%kolik

       old_storage = storage(el)

       ! Inputs [mm] (precip + interflow + what was already stored)
       local_input = precip(el,tstep) + qinter(el,tstep) + old_storage

       ! Losses [mm]
       local_losses = elements%hydrobal(el)%ET + &
                      elements%hydrobal(el)%Li + &
                      elements%hydrobal(el)%Qgw

       if (local_losses >= local_input) then
          storage_new(el) = 0.0_rkind
          surplus         = 0.0_rkind
       else
          surplus = local_input - local_losses

          ! Fill storage up to capacity
          if (surplus <= capacity(el)) then
             storage_new(el) = surplus
             surplus         = 0.0_rkind
          else
             storage_new(el) = capacity(el)
             surplus         = surplus - capacity(el)
          end if
       end if

       ! Local overflow that will be routed: ponding surplus + direct surface runoff
       overflow_total(el) = surplus + Qsurf_result(el,tstep)

    end do

    ! 2) Routing along downstream graph (using flow_order: upstream -> downstream)
    do i = 1_ikind, elements%kolik
       el  = flow_order(i)
       dwn = downstream(el)

       ! Pass-through outflow at this element
       elements%hydrobal(el)%outflow = elements%hydrobal(el)%inflow + &
                                       overflow_total(el)

       if (dwn > 0_ikind) then
          ! Add all this outflow to the downstream element's inflow
          elements%hydrobal(dwn)%inflow = elements%hydrobal(dwn)%inflow + &
                                          elements%hydrobal(el)%outflow
       else
          ! This is an outlet element → contributes to catchment outlet discharge
          outlet_Q(tstep) = outlet_Q(tstep) + elements%hydrobal(el)%outflow
       end if
    end do

    ! 3) Final mass balance and state update
    do el = 1_ikind, elements%kolik
       old_storage           = storage(el)
       storage(el)           = storage_new(el)
       elements%overflow(el) = overflow_total(el)

       elements%hydrobal(el)%deltas = precip(el,tstep) + &
                                      qinter(el,tstep) + &
                                      old_storage + &
                                      elements%hydrobal(el)%inflow - &
                                      ( elements%hydrobal(el)%ET      + &
                                        elements%hydrobal(el)%Li      + &
                                        elements%hydrobal(el)%Qgw     + &
                                        storage(el)                   + &
                                        elements%hydrobal(el)%outflow )

       deltas(el,tstep) = elements%hydrobal(el)%deltas
    end do

    deallocate(overflow_total, storage_new)
  end subroutine route_step


end module tools
