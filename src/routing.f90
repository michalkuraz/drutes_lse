module routing
 use typy
 use globals
 use hydrotools
 use tools
 use hydrofnc
  Implicit none

contains

  !==============================================================
  !   FLOW TOPOLOGY INITIALIZATION
  !==============================================================
  subroutine init_flow_topology()
    integer(kind=ikind) :: nel, i

    nel = elements%kolik
    if (nel <= 0_ikind) return

    call find_neighbours(elements, nodes)
    call compute_element_slopes()

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


    call build_graph()
    call finalize_routing_graph()



  end subroutine init_flow_topology
  
  
  
  
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
  !   Convert stream candidates into routing arrays
  !==============================================================
  subroutine finalize_routing_graph()
    integer(kind=ikind) :: nel, el, i, j, best_nb, tmp_el
    real(kind=rkind)    :: best_slope, tmp_z
    integer(kind=ikind), allocatable :: order(:)
    real(kind=rkind),    allocatable :: order_z(:)

    nel = elements%kolik
    if (nel <= 0_ikind) return

    downstream = 0_ikind

    do el = 1_ikind, nel
       best_nb    = 0_ikind
       best_slope = 0.0_rkind

       do i = 1_ikind, 3_ikind
          if (elements%downstream(el)%els(i) > 0_ikind .and. &
              elements%downstream(el)%slopes(i) > best_slope) then
             best_nb    = elements%downstream(el)%els(i)
             best_slope = elements%downstream(el)%slopes(i)
          end if
       end do

       downstream(el) = best_nb
    end do

    if (allocated(upstream_count)) deallocate(upstream_count)
    if (allocated(upstream_list))  deallocate(upstream_list)
    allocate(upstream_count(nel))
    allocate(upstream_list(nel, nel))

    upstream_count = 0_ikind
    upstream_list  = 0_ikind

    do el = 1_ikind, nel
       if (downstream(el) > 0_ikind) then
          upstream_count(downstream(el)) = upstream_count(downstream(el)) + 1_ikind
          upstream_list(downstream(el), upstream_count(downstream(el))) = el
       end if
    end do

    allocate(order(nel))
    allocate(order_z(nel))

    do el = 1_ikind, nel
       order(el)   = el
       order_z(el) = elements%avgalt(el)
    end do

    do i = 1_ikind, nel - 1_ikind
       do j = i + 1_ikind, nel
          if (order_z(j) > order_z(i)) then
             tmp_z      = order_z(i)
             order_z(i) = order_z(j)
             order_z(j) = tmp_z

             tmp_el   = order(i)
             order(i) = order(j)
             order(j) = tmp_el
          end if
       end do
    end do

    flow_order = order

    deallocate(order, order_z)
  end subroutine finalize_routing_graph


  
  ! ============================================================
  ! routing
  ! ============================================================

  subroutine route_step(tstep)

  integer(kind=ikind), intent(in) :: tstep

  integer(kind=ikind) :: el, i, side, dwn
  real(kind=rkind)    :: old_storage
  real(kind=rkind)    :: water_available
  real(kind=rkind)    :: losses
  real(kind=rkind)    :: available_after_losses
  real(kind=rkind)    :: storage_capacity
  real(kind=rkind)    :: routed_water
  real(kind=rkind)    :: weight_sum
  real(kind=rkind)    :: edge_weight
  real(kind=rkind)    :: routed_fraction
  real(kind=rkind)    :: catchment_area
  real(kind=rkind)    :: outlet_volume

  real(kind=rkind), allocatable :: runoff(:)
  real(kind=rkind), allocatable :: new_storage(:)
  real(kind=rkind), allocatable :: local_residual(:)
  real(kind=rkind), allocatable :: routing_residual(:)

  allocate(runoff(elements%kolik))
  allocate(new_storage(elements%kolik))
  allocate(local_residual(elements%kolik))
  allocate(routing_residual(elements%kolik))
! ============================================================
  ! 1. Initialization
  ! ============================================================

  runoff           = 0.0_rkind
  new_storage      = storage
  local_residual   = 0.0_rkind
  routing_residual = 0.0_rkind

  ! ============================================================
  !2. Reset inflow and outflow
  ! ============================================================
  elements%hydrobal(:)%inflow  = 0.0_rkind
  elements%hydrobal(:)%outflow = 0.0_rkind

  do el = 1_ikind, elements%kolik

     old_storage = storage(el)

     ! ============================================================
  ! 3. Local water balance
  ! ============================================================
     water_available = precip(el,tstep) + qinter(el,tstep) + old_storage

     losses = elements%hydrobal(el)%ET + &
              elements%hydrobal(el)%Li + &
              elements%hydrobal(el)%Qgw

     losses = min(losses, water_available)

     ! ============================================================
  ! 4. Storage capacity calculation
  ! ============================================================
     storage_capacity = capacity(el) / &
          (1.0_rkind + storage_slope_coeff * max(elements%slope(el), 0.0_rkind))

     storage_capacity = max(storage_capacity, 0.0_rkind)

     ! ============================================================
  ! 5. New storage and runoff
  ! ============================================================
     available_after_losses = water_available - losses

     new_storage(el) = min(available_after_losses, storage_capacity)

     runoff(el) = max(available_after_losses - new_storage(el), 0.0_rkind)


 ! ============================================================
  !6. Local residual check
  ! ============================================================
     local_residual(el) = water_available - losses - new_storage(el) - runoff(el)

  end do


! ============================================================
  !7. Start routing
  ! ============================================================
  outlet_Q(tstep) = 0.0_rkind

  do i = 1_ikind, elements%kolik

     el = flow_order(i)


! ============================================================
  !8. Routed water calculation
  ! ============================================================
     routed_water = elements%hydrobal(el)%inflow + runoff(el)

     elements%hydrobal(el)%outflow = routed_water

     if (routed_water > 0.0_rkind) then

        weight_sum = 0.0_rkind

        do side = 1_ikind, 3_ikind


! ============================================================
  !9. Downstream weight calculation
  ! ============================================================
           dwn = elements%downstream(el)%els(side)

           if (dwn > 0_ikind) then

              edge_weight = max(elements%downstream(el)%slopes(side), 0.0_rkind) * &
                            max(elements%downstream(el)%widths(side), min_edge_width)

              weight_sum = weight_sum + edge_weight

           end if

        end do

        if (weight_sum > 0.0_rkind) then

           do side = 1_ikind, 3_ikind

              dwn = elements%downstream(el)%els(side)

              if (dwn > 0_ikind) then

                 edge_weight = max(elements%downstream(el)%slopes(side), 0.0_rkind) * &
                               max(elements%downstream(el)%widths(side), min_edge_width)

! ============================================================
  !10. Flow distribution
  ! ============================================================
                 routed_fraction = edge_weight / weight_sum

                 elements%hydrobal(dwn)%inflow = elements%hydrobal(dwn)%inflow + &
                      routed_water * routed_fraction

              end if

           end do

        else

  ! ============================================================
  !11. Outlet flow
  ! ============================================================

           outlet_Q(tstep) = outlet_Q(tstep) + routed_water

        end if

     end if

  end do

  do el = 1_ikind, elements%kolik

  ! ============================================================
  !12. Save results
  ! ============================================================

     storage(el) = new_storage(el)

     elements%overflow(el) = runoff(el)

     routing_residual(el) = elements%hydrobal(el)%inflow + runoff(el) - &
                            elements%hydrobal(el)%outflow

     elements%hydrobal(el)%deltas = local_residual(el) + routing_residual(el)

     deltas(el,tstep)          = elements%hydrobal(el)%deltas
     Qin_result(el,tstep)      = elements%hydrobal(el)%inflow
     Qout_result(el,tstep)     = elements%hydrobal(el)%outflow
     Overflow_result(el,tstep) = elements%overflow(el)
     Storage_result(el,tstep)  = storage(el)

  end do

  ! ============================================================
  !13. Outlet discharge conversion
  ! ============================================================
  catchment_area = sum(elements%area)

  outlet_volume = outlet_Q(tstep) / 1000.0_rkind * catchment_area

  if (dt_days > 0.0_rkind) then
     outlet_Q_m3s(tstep) = outlet_volume / (dt_days * 86400.0_rkind)
  else
     outlet_Q_m3s(tstep) = 0.0_rkind
  end if

  deallocate(runoff)
  deallocate(new_storage)
  deallocate(local_residual)
  deallocate(routing_residual)

end subroutine route_step

end module routing
