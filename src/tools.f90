module tools
  use typy
  use globals
  implicit none
contains

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

    do i = 1, elements%kolik
     conduct(i) = 3.0e-6_rkind + 4.0e-8_rkind * elements%avgalt(i)
      G(i)       = 0.0_rkind
    end do

    if (n_steps < 10) stop "Need at least 10 time steps for the hardcoded forcing."

    do i = 1, elements%kolik
      precip(i,1:10) = [ &
          0.0_rkind, 0.0_rkind, 17.0_rkind, 12.0_rkind, 9.0_rkind, &
          7.0_rkind, 40.0_rkind, 0.0_rkind, 0.0_rkind, 3.0_rkind ]

      precip(i,:) = precip(i,:) * (0.90_rkind + 0.0020_rkind * (elements%avgalt(i) - 115.0_rkind))


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

      soilcontent(i,1:10) = [ &
           0.05_rkind, 0.055_rkind, 0.062_rkind, 0.060_rkind, 0.040_rkind, &
           0.070_rkind, 0.090_rkind, 0.200_rkind, 0.250_rkind, 0.265_rkind ]

      soilcontent(i,:) = soilcontent(i,:) * (1.10_rkind - 0.0015_rkind * (elements%avgalt(i) - 115.0_rkind))
      
    end do

  

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


  ! -----------------------------------------------------
  ! Find neighbours: elements that share an edge
  ! -----------------------------------------------------
  subroutine find_neighbours(elems, nodes_loc)
    type(elements_str), intent(inout) :: elems
    type(nodes_str),    intent(in)    :: nodes_loc

    integer(kind=ikind) :: nel
    integer(kind=ikind) :: i, j, k, shared
    integer(kind=ikind) :: e1(3), e2(3)

    nel = elems%kolik
    if (nel <= 0_ikind) return

    if (.not. allocated(elems%neighbours)) then
      allocate(elems%neighbours(nel,3))
    end if
    elems%neighbours = 0_ikind

    do i = 1_ikind, nel
      e1 = elems%data(i,:)
      do j = 1_ikind, nel
        if (i == j) cycle
        e2 = elems%data(j,:)

        shared = 0_ikind
        do k = 1_ikind, 3_ikind
          if (any(e2 == e1(k))) shared = shared + 1_ikind
        end do

        if (shared >= 2_ikind) then
          do k = 1_ikind, 3_ikind
            if (elems%neighbours(i,k) == 0_ikind) then
              elems%neighbours(i,k) = j
              exit
            end if
          end do
        end if
      end do
    end do
  end subroutine find_neighbours


  function avg(r) result(avgval)
    real(kind=rkind), dimension(:), intent(in) :: r
    real(kind=rkind) :: avgval
    integer(kind=ikind) :: i

    avgval = 0.0_rkind
    do i = 1, size(r)
      avgval = avgval + r(i)
    end do
    avgval = avgval / real(size(r), kind=rkind)
  end function avg


  !==============================================================
  !   READ MESH (nodes + elements)
  !==============================================================
  subroutine read_mesh(filename)
    character(len=*), intent(in) :: filename
    integer :: fileid, ios, n_nodes, n_elem
    integer :: id, n1, n2, n3, i
    real(kind=rkind) :: x, y, zloc
    character(len=256) :: line
    logical :: file_exists

    inquire(file=filename, exist=file_exists)
    if (.not. file_exists) then
      print *, " Error: File ", trim(filename), " not found."
      stop
    end if

    open(newunit=fileid, file=filename, status="old", action="read", iostat=ios)
    if (ios /= 0) stop " Could not open mesh file."

    do
      read(fileid,*,iostat=ios) line
      if (ios < 0) stop "Unexpected EOF before node count"
      if (len_trim(line) > 0 .and. line(1:1) /= "#") exit
    end do
    read(line,*) n_nodes
    nodes%kolik = n_nodes

    allocate(nodes%data(n_nodes,2))
    allocate(nodes%altitude(n_nodes))

    i = 0
    do while (i < n_nodes)
      read(fileid,'(A)',iostat=ios) line
      if (ios < 0) exit
      if (len_trim(line) == 0 .or. line(1:1) == "#") cycle
      read(line,*) id, x, y, zloc
      i = i + 1
      nodes%data(i,1)   = x
      nodes%data(i,2)   = y
      nodes%altitude(i) = zloc
    end do
    print *, " Read", i, "nodes."

    do
      read(fileid,'(A)',iostat=ios) line
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
      read(fileid,'(A)',iostat=ios) line
      if (ios < 0) exit
      if (len_trim(line) == 0 .or. line(1:1) == "#") cycle
      read(line,*) id, n1, n2, n3
      i = i + 1
      elements%data(i,:)   = [n1, n2, n3]
      elements%material(i) = 1_ikind
    end do
    print *, " Read", i, "triangular elements."
    close(fileid)

    if (.not. allocated(elements%hydrobal)) then
      allocate(elements%hydrobal(n_elem))
    end if
  end subroutine read_mesh


  !==============================================================
  !   COMPUTE ELEMENT AREAS
  !==============================================================
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

    print *, " Total mesh area = ", total_area
  end subroutine compute_areas


  pure function area_triangle(a,b,c) result(area)
    real(kind=rkind), intent(in) :: a(2), b(2), c(2)
    real(kind=rkind) :: area
    area = 0.5_rkind * abs(a(1)*(b(2)-c(2)) + b(1)*(c(2)-a(2)) + c(1)*(a(2)-b(2)))
  end function area_triangle


  !==============================================================
  !   COMPUTE average altitude per element
  !==============================================================
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
        z1 = nodes%data(n1,2)
        z2 = nodes%data(n2,2)
        z3 = nodes%data(n3,2)
      end if

      elements%avgalt(i) = (z1 + z2 + z3) / 3.0_rkind
    end do
  end subroutine compute_avgalt


  !==============================================================
  !   FLOW TOPOLOGY INITIALIZATION
  !==============================================================
  subroutine init_flow_topology()
    integer(kind=ikind) :: nel, i

    nel = elements%kolik
    if (nel <= 0_ikind) return

    call find_neighbours(elements, nodes)

    if (.not. allocated(downstream)) allocate(downstream(nel))
    if (.not. allocated(flow_order)) allocate(flow_order(nel))

    if (.not. allocated(elements%outlet))      allocate(elements%outlet(nel))
    if (.not. allocated(elements%watershed))   allocate(elements%watershed(nel))
    if (.not. allocated(elements%ndwatershed)) allocate(elements%ndwatershed(nel,2))
    if (.not. allocated(elements%ndoutlet))    allocate(elements%ndoutlet(nel,2))

    downstream           = 0_ikind
    flow_order           = 0_ikind
    elements%overflow    = 0.0_rkind
    elements%outlet      = .false.
    elements%watershed   = .false.
    elements%ndwatershed = 0_ikind
    elements%ndoutlet    = 0_ikind

    do i = 1, nel
      elements%downstream(i)%els    = 0_ikind
      elements%downstream(i)%slopes = 0.0_rkind
      elements%downstream(i)%widths = 0.0_rkind

      elements%upstream(i)%els      = 0_ikind
      elements%upstream(i)%slopes   = 0.0_rkind
      elements%upstream(i)%widths   = 0.0_rkind
    end do

    call build_downstream_graph()
    call build_graph()
    call build_flow_order()
    call build_upstream_graph()
  end subroutine init_flow_topology


  !==============================================================
  !   built downstream graph based on lowest neighbour
  !==============================================================
  subroutine build_downstream_graph()
    integer(kind=ikind) :: i, j, nb, nel
    real(kind=rkind)    :: myz, bestz
    integer(kind=ikind) :: best_nb

    nel = elements%kolik

    do i = 1, nel
      myz     = elements%avgalt(i)
      bestz   = myz
      best_nb = 0_ikind

      do j = 1, ubound(elements%neighbours,2)
        nb = elements%neighbours(i,j)
        if (nb > 0_ikind) then
          if (elements%avgalt(nb) < bestz) then
            bestz   = elements%avgalt(nb)
            best_nb = nb
          end if
        end if
      end do

      downstream(i)      = best_nb
      elements%outlet(i) = (best_nb == 0_ikind)
    end do
  end subroutine build_downstream_graph


  function setlines(el, nd) result(pts)
    integer(kind=ikind), intent(in) :: el, nd
    integer(kind=ikind), dimension(2) :: pts
    integer(kind=ikind), dimension(3) :: neighnds, elnds
    integer(kind=ikind) :: pos, i, j

    elnds    = elements%data(el,:)
    neighnds = elements%data(nd,:)

    pts = 0_ikind
    pos = 0

    do i = 1, 3
      do j = 1, 3
        if (neighnds(i) == elnds(j)) then
          pos = pos + 1
          if (pos <= 2) then
            pts(pos) = neighnds(i)
          else
            print *, "Error in setlines: two elements share all three nodes."
            error stop
          end if
        end if
      end do
    end do
  end function setlines


  function matchlines(l1, l2) result(true)
    integer(kind=ikind), dimension(2), intent(in) :: l1, l2
    logical :: true

    true = .false.
    if (l1(1) == l2(1) .and. l1(2) == l2(2)) then
      true = .true.
      return
    end if

    if (l1(1) == l2(2) .and. l1(2) == l2(1)) then
      true = .true.
      return
    end if
  end function matchlines


  function getcenter(el) result(xy)
    integer(kind=ikind), intent(in) :: el
    real(kind=rkind), dimension(2) :: xy

    xy(1) = avg(nodes%data(elements%data(el,:),1))
    xy(2) = avg(nodes%data(elements%data(el,:),2))
  end function getcenter


  function dist(A,B) result(l)
    real(kind=rkind), dimension(:), intent(in) :: A,B
    real(kind=rkind) :: l
    l = sqrt((A(1)-B(1))*(A(1)-B(1)) + (A(2)-B(2))*(A(2)-B(2)))
  end function dist


  !==============================================================
  !   build graph
  !==============================================================
  subroutine build_graph()
    integer(kind=ikind) :: el, i
    integer(kind=ikind), dimension(3) :: ngh
    real(kind=rkind),    dimension(3) :: nghalt
    real(kind=rkind) :: myalt, slopeval
    real(kind=rkind), dimension(2) :: myel, ngel
    integer(kind=ikind), dimension(3,2) :: nghlines

    do el = 1, elements%kolik
      elements%downstream(el)%els    = 0_ikind
      elements%downstream(el)%slopes = 0.0_rkind
      elements%downstream(el)%widths = 0.0_rkind

      elements%upstream(el)%els      = 0_ikind
      elements%upstream(el)%slopes   = 0.0_rkind
      elements%upstream(el)%widths   = 0.0_rkind
    end do

    do el = 1, elements%kolik
      myalt = elements%avgalt(el)
      ngh   = elements%neighbours(el,:)
      myel  = getcenter(el)

      do i = 1, 3
        if (ngh(i) /= 0) then
          nghlines(i,:) = setlines(ngh(i), el)
          nghalt(i) = elements%avgalt(ngh(i))
        else
          nghlines(i,:) = 0_ikind
          nghalt(i) = -9999.0_rkind
        end if
      end do

      do i = 1, 3
        if (ngh(i) /= 0) then
          ngel = getcenter(ngh(i))
          slopeval = abs(myalt - nghalt(i)) / max(dist(myel, ngel), 1.0e-12_rkind)

          if (nghalt(i) <= myalt) then
            elements%downstream(el)%els(i)    = ngh(i)
            elements%downstream(el)%slopes(i) = slopeval
            elements%downstream(el)%widths(i) = dist(nodes%data(nghlines(i,1),:), &
                                                     nodes%data(nghlines(i,2),:))
          else
            elements%upstream(el)%els(i)    = ngh(i)
            elements%upstream(el)%slopes(i) = slopeval
            elements%upstream(el)%widths(i) = dist(nodes%data(nghlines(i,1),:), &
                                                   nodes%data(nghlines(i,2),:))
          end if
        end if
      end do
    end do
  end subroutine build_graph


  !==============================================================
  !   Build flow order based on elevation (highest first)
  !==============================================================
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


  !==============================================================
  !   Build upstream graph
  !==============================================================
  subroutine build_upstream_graph()
    integer(kind=ikind) :: nel, i, j, max_up

    nel = elements%kolik
    if (nel <= 0_ikind) return

    if (allocated(upstream_count)) deallocate(upstream_count)
    allocate(upstream_count(nel))
    upstream_count = 0_ikind

    do i = 1_ikind, nel
      j = downstream(i)
      if (j > 0_ikind) upstream_count(j) = upstream_count(j) + 1_ikind
    end do

    max_up = 0_ikind
    do i = 1_ikind, nel
      if (upstream_count(i) > max_up) max_up = upstream_count(i)
    end do

    if (allocated(upstream_list)) deallocate(upstream_list)

    if (max_up > 0_ikind) then
      allocate(upstream_list(nel, max_up))
      upstream_list  = 0_ikind
      upstream_count = 0_ikind

      do i = 1_ikind, nel
        j = downstream(i)
        if (j > 0_ikind) then
          upstream_count(j) = upstream_count(j) + 1_ikind
          upstream_list(j, upstream_count(j)) = i
        end if
      end do
    end if
  end subroutine build_upstream_graph


  !==============================================================
  !   Print upstream flows for diagnostics
  !==============================================================
  subroutine print_upstream_flows(tstep)
  integer(kind=ikind), intent(in) :: tstep
  integer(kind=ikind) :: nel, e, k
  real(kind=rkind)    :: Qin_from_up

  nel = elements%kolik

  if (.not. allocated(Qin_result)) then
     print *, " No routed time-series data available."
     return
  end if

  if (.not. allocated(Qout_result)) then
     print *, " No routed outflow time-series data available."
     return
  end if

  if (.not. allocated(upstream_count)) then
     print *, " Upstream graph not built yet (call init_flow_topology first)."
     return
  end if

  if (tstep < 1_ikind .or. tstep > n_steps) then
     print *, " Invalid time step in print_upstream_flows: ", tstep
     return
  end if

  print *, "---------------------------------------------------------------"
  print *, " Time step = ", tstep
  print *, " el | n_up | Qin_from_up    Qin(field)     Qout_elem"
  print *, "---------------------------------------------------------------"

  do e = 1_ikind, nel
     Qin_from_up = 0.0_rkind

     if (allocated(upstream_list)) then
        do k = 1_ikind, upstream_count(e)
           Qin_from_up = Qin_from_up + Qout_result(upstream_list(e,k), tstep)
        end do
     end if

     print *, e, upstream_count(e), Qin_from_up, &
              Qin_result(e, tstep), Qout_result(e, tstep)
  end do

end subroutine print_upstream_flows

  !==============================================================
  !   Print graph diagnostics
  !==============================================================
  subroutine print_graph_diagnostics()
    integer(kind=ikind) :: el, i, nb
    real(kind=rkind)    :: width, slope, nalt

    print *, "--------------------------------------------------------------------------"
    print *, " el   nb   dir      z_el        z_nb        slope        edge_width"
    print *, "--------------------------------------------------------------------------"

    do el = 1_ikind, elements%kolik
      do i = 1_ikind, 3_ikind
        nb = elements%downstream(el)%els(i)
        if (nb > 0_ikind) then
          width = elements%downstream(el)%widths(i)
          slope = elements%downstream(el)%slopes(i)
          nalt  = elements%avgalt(nb)
          print *, el, nb, "down", elements%avgalt(el), nalt, slope, width
        end if

        nb = elements%upstream(el)%els(i)
        if (nb > 0_ikind) then
          width = elements%upstream(el)%widths(i)
          slope = elements%upstream(el)%slopes(i)
          nalt  = elements%avgalt(nb)
          print *, el, nb, "up", elements%avgalt(el), nalt, slope, width
        end if
      end do
    end do

    print *, "--------------------------------------------------------------------------"
  end subroutine print_graph_diagnostics


  !==============================================================
  !  Print element balance
  !==============================================================
  subroutine print_element_balance(tstep)
    integer(kind=ikind), intent(in) :: tstep
    integer(kind=ikind) :: i
    real(kind=rkind)    :: surplus

    if (.not. allocated(elements%hydrobal)) then
      print *, "No hydrological data to display."
      return
    end if

    print *, "----------------------------------------------------------------------------------------------------------------"
    print *, " Element |     P        ET      Qsurf       Li       Qgw      Qin     Surplus   Qout_elem  Overflow    deltaS"
    print *, "----------------------------------------------------------------------------------------------------------------"

    do i = 1_ikind, elements%kolik
      surplus = elements%hydrobal(i)%Qsurf + &
                elements%hydrobal(i)%Li    + &
                elements%hydrobal(i)%Qgw

      print *, i, precip(i,tstep), elements%hydrobal(i)%ET, &
               elements%hydrobal(i)%Qsurf, elements%hydrobal(i)%Li, &
               elements%hydrobal(i)%Qgw, elements%hydrobal(i)%inflow, &
               surplus, elements%hydrobal(i)%outflow, &
               elements%overflow(i), elements%hydrobal(i)%deltas
    end do
  end subroutine print_element_balance


  !==============================================================
  !   EXPORT ELEMENT WATER BALANCE TO CSV
  !==============================================================
  subroutine export_element_balance(filename)
  character(len=*), intent(in) :: filename
  integer(kind=ikind) :: i, t
  integer :: unit, ios

  open(newunit=unit, file=filename, status="replace", action="write", iostat=ios)
  if (ios /= 0) then
     print *, "Error opening file: ", trim(filename)
     return
  end if

  write(unit,'(A)') "step,element,area,z_avg,downstream," // &
                    "P,ET,Qsurf,Li,Qgw,Qin,Qout,Overflow,Storage,DeltaS"

  do t = 1, n_steps
     do i = 1, elements%kolik
        write(unit,'(I4, ",", I4, ",", F10.3, ",", F10.3, ",", I4, 10(",",F12.6))') &
             t, i, elements%area(i), elements%avgalt(i), downstream(i), &
             precip(i,t), ET_flux(i,t), Qsurf_result(i,t), L_result(i,t), &
             Qgw_result(i,t), Qin_result(i,t), Qout_result(i,t), &
             Overflow_result(i,t), Storage_result(i,t), deltas(i,t)
     end do
  end do

  close(unit)
  print *, "Element balances saved to: ", trim(filename)
end subroutine export_element_balance


  !==============================================================
  !   Routing
  !==============================================================
  subroutine route_step(tstep)
  integer(kind=ikind), intent(in) :: tstep

  integer(kind=ikind) :: el, i, dwn
  real(kind=rkind)    :: old_storage, local_input, local_losses, surplus
  real(kind=rkind)    :: area_total_m2, volume_m3
  real(kind=rkind), allocatable :: overflow_total(:), storage_new(:)

  allocate(overflow_total(elements%kolik))
  allocate(storage_new(elements%kolik))

  elements%hydrobal(:)%inflow  = 0.0_rkind
  elements%hydrobal(:)%outflow = 0.0_rkind

  overflow_total = 0.0_rkind
  storage_new    = storage

  do el = 1_ikind, elements%kolik
     old_storage = storage(el)

     local_input = precip(el,tstep) + qinter(el,tstep) + old_storage

     local_losses = elements%hydrobal(el)%ET + &
                    elements%hydrobal(el)%Li + &
                    elements%hydrobal(el)%Qgw

     if (local_losses >= local_input) then
        storage_new(el) = 0.0_rkind
        surplus         = 0.0_rkind
     else
        surplus = local_input - local_losses

        if (surplus <= capacity(el)) then
           storage_new(el) = surplus
           surplus         = 0.0_rkind
        else
           storage_new(el) = capacity(el)
           surplus         = surplus - capacity(el)
        end if
     end if

     overflow_total(el) = surplus
  end do

  outlet_Q(tstep) = 0.0_rkind

  do i = 1_ikind, elements%kolik
     el  = flow_order(i)
     dwn = downstream(el)

     elements%hydrobal(el)%outflow = elements%hydrobal(el)%inflow + overflow_total(el)

     if (dwn > 0_ikind) then
        elements%hydrobal(dwn)%inflow = elements%hydrobal(dwn)%inflow + &
                                        elements%hydrobal(el)%outflow
     else
        outlet_Q(tstep) = outlet_Q(tstep) + elements%hydrobal(el)%outflow
     end if
  end do

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

     deltas(el,tstep)          = elements%hydrobal(el)%deltas
     Qin_result(el,tstep)      = elements%hydrobal(el)%inflow
     Qout_result(el,tstep)     = elements%hydrobal(el)%outflow
     Overflow_result(el,tstep) = elements%overflow(el)
     Storage_result(el,tstep)  = storage(el)
  end do

  ! Convert catchment outlet from mm over catchment per step to m3/s
  area_total_m2 = sum(elements%area)
  volume_m3 = outlet_Q(tstep) / 1000.0_rkind * area_total_m2

  if (dt_days > 0.0_rkind) then
     outlet_Q_m3s(tstep) = volume_m3 / (dt_days * 86400.0_rkind)
  else
     outlet_Q_m3s(tstep) = 0.0_rkind
  end if

  deallocate(overflow_total, storage_new)
end subroutine route_step
end module tools