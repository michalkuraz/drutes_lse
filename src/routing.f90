module routing
  use typy
  use globals
  use hydrotools
  use tools
  use hydrofnc
  implicit none

contains

!==============================================================
      ! FLOW TOPOLOGY INITIALIZATION
      ! BUILD GRAPH
      ! FINALIZE ROUTING GRAPH
      ! ROUTE SURFACE EXCESS WATER route_step(tstep)
  !==============================================================





  !==============================================================
  ! FLOW TOPOLOGY INITIALIZATION
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
  ! BUILD GRAPH
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
  ! FINALIZE ROUTING GRAPH
  !==============================================================
  subroutine finalize_routing_graph()
    integer(kind=ikind) :: nel, el, i, j, best_nb, tmp_el
    real(kind=rkind)    :: best_slope, tmp_z
    integer(kind=ikind), allocatable :: order(:)
    real(kind=rkind),    allocatable :: order_z(:)

    nel = elements%kolik
    if (nel <= 0_ikind) return

    downstream = 0_ikind

    do el = 1, nel
      best_nb    = 0_ikind
      best_slope = 0.0_rkind

      do i = 1, 3
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

    do el = 1, nel
      if (downstream(el) > 0_ikind) then
        upstream_count(downstream(el)) = upstream_count(downstream(el)) + 1_ikind
        upstream_list(downstream(el), upstream_count(downstream(el))) = el
      end if
    end do

    allocate(order(nel))
    allocate(order_z(nel))

    do el = 1, nel
      order(el)   = el
      order_z(el) = elements%avgalt(el)
    end do

    do i = 1, nel - 1
      do j = i + 1, nel
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


!==============================================================
! ROUTE SURFACE EXCESS WATER route step (tstep)
!
! Routes surface runoff between elements according to
! downstream slope/width weights.
!
! Results saved for every element and timestep:
!   Qin_result
!   Qout_result
!   Overflow_result
!   Storage_result
!   deltas
!
! Catchment outlet:
!   outlet_Q
!   outlet_Q_m3s
!==============================================================
subroutine route_step(tstep)

  implicit none

  integer(kind=ikind), intent(in) :: tstep

  integer(kind=ikind) :: el
  integer(kind=ikind) :: i
  integer(kind=ikind) :: side
  integer(kind=ikind) :: dwn

  real(kind=rkind) :: old_storage
  real(kind=rkind) :: water_available
  real(kind=rkind) :: losses
  real(kind=rkind) :: available_after_losses

  real(kind=rkind) :: storage_capacity

  real(kind=rkind) :: routed_water
  real(kind=rkind) :: routed_fraction

  real(kind=rkind) :: weight_sum
  real(kind=rkind) :: edge_weight

  real(kind=rkind) :: catchment_area
  real(kind=rkind) :: outlet_volume

  real(kind=rkind), allocatable :: runoff(:)
  real(kind=rkind), allocatable :: new_storage(:)
  real(kind=rkind), allocatable :: local_residual(:)
  real(kind=rkind), allocatable :: routing_residual(:)


  !============================================================
  ! SAFETY CHECKS
  !============================================================

  if (tstep < 1_ikind .or. tstep > n_steps) then

    write(*,*) 'ERROR in route_step'
    write(*,*) 'Invalid timestep = ', tstep
    write(*,*) 'n_steps          = ', n_steps

    error stop

  end if


  if (elements%kolik <= 0_ikind) then

    write(*,*) 'ERROR in route_step'
    write(*,*) 'No mesh elements are available.'

    error stop

  end if


  if (.not. allocated(storage)) then
    error stop 'ERROR route_step: storage not allocated.'
  end if

  if (.not. allocated(capacity)) then
    error stop 'ERROR route_step: capacity not allocated.'
  end if

  if (.not. allocated(flow_order)) then
    error stop 'ERROR route_step: flow_order not allocated.'
  end if

  if (.not. allocated(Qin_result)) then
    error stop 'ERROR route_step: Qin_result not allocated.'
  end if

  if (.not. allocated(Qout_result)) then
    error stop 'ERROR route_step: Qout_result not allocated.'
  end if

  if (.not. allocated(Overflow_result)) then
    error stop 'ERROR route_step: Overflow_result not allocated.'
  end if

  if (.not. allocated(Storage_result)) then
    error stop 'ERROR route_step: Storage_result not allocated.'
  end if

  if (.not. allocated(deltas)) then
    error stop 'ERROR route_step: deltas not allocated.'
  end if

  if (.not. allocated(outlet_Q)) then
    error stop 'ERROR route_step: outlet_Q not allocated.'
  end if

  if (.not. allocated(outlet_Q_m3s)) then
    error stop 'ERROR route_step: outlet_Q_m3s not allocated.'
  end if


  !============================================================
  ! ALLOCATE TEMPORARY ROUTING ARRAYS
  !============================================================

  allocate(runoff(elements%kolik))

  allocate(new_storage(elements%kolik))

  allocate(local_residual(elements%kolik))

  allocate(routing_residual(elements%kolik))


  runoff = 0.0_rkind

  new_storage = storage

  local_residual = 0.0_rkind

  routing_residual = 0.0_rkind


  !============================================================
  ! RESET ROUTING FLUXES FOR CURRENT TIMESTEP
  !============================================================

  elements%hydrobal(:)%inflow = 0.0_rkind

  elements%hydrobal(:)%outflow = 0.0_rkind

  elements%overflow(:) = 0.0_rkind


  !============================================================
  ! PART 1
  !
  ! CALCULATE LOCAL SURFACE STORAGE AND EXCESS RUNOFF
  !============================================================

  do el = 1, elements%kolik


    !----------------------------------------------------------
    ! Previous routing storage
    !----------------------------------------------------------

    old_storage = &
         max(storage(el), 0.0_rkind)


    !----------------------------------------------------------
    ! Surface water entering local routing storage
    !
    ! Current formulation:
    ! precipitation + previous routing storage
    !----------------------------------------------------------

    water_available = &
         max(Pm(el,tstep), 0.0_rkind) + &
         old_storage


    !----------------------------------------------------------
    ! Surface losses
    !
    ! ET and q1 are taken from the hydrological calculation.
    ! Do not allow losses to exceed available surface water.
    !----------------------------------------------------------

    losses = &
         max(elements%hydrobal(el)%ET, 0.0_rkind) + &
         max(elements%hydrobal(el)%q1, 0.0_rkind)


    losses = &
         min(losses, water_available)


    !----------------------------------------------------------
    ! Remaining water
    !----------------------------------------------------------

    available_after_losses = &
         max(water_available - losses, 0.0_rkind)


    !----------------------------------------------------------
    ! Slope-adjusted surface storage capacity
    !----------------------------------------------------------

    storage_capacity = &
         capacity(el) / &
         (1.0_rkind + &
          storage_slope_coeff * &
          max(elements%slope(el), 0.0_rkind))


    storage_capacity = &
         max(storage_capacity, 0.0_rkind)


    !----------------------------------------------------------
    ! Water retained locally
    !----------------------------------------------------------

    new_storage(el) = &
         min(available_after_losses, &
             storage_capacity)


    !----------------------------------------------------------
    ! Excess water available for routing
    !----------------------------------------------------------

    runoff(el) = &
         max(available_after_losses - &
             new_storage(el), &
             0.0_rkind)


    !----------------------------------------------------------
    ! Local numerical balance check
    !
    ! Should be approximately zero:
    !
    ! water_available
    ! - losses
    ! - storage
    ! - runoff
    !----------------------------------------------------------

    local_residual(el) = &
         water_available - &
         losses - &
         new_storage(el) - &
         runoff(el)


  end do


  !============================================================
  ! PART 2
  !
  ! ROUTE EXCESS WATER THROUGH THE ELEMENT NETWORK
  !============================================================

  outlet_Q(tstep) = 0.0_rkind


  do i = 1, elements%kolik


    !----------------------------------------------------------
    ! Process elements according to routing order
    !----------------------------------------------------------

    el = flow_order(i)


    !----------------------------------------------------------
    ! Defensive check
    !----------------------------------------------------------

    if (el < 1_ikind .or. &
        el > elements%kolik) then

      write(*,*) 'ERROR in route_step'
      write(*,*) 'Invalid flow_order value.'
      write(*,*) 'Position   = ', i
      write(*,*) 'Element    = ', el
      write(*,*) 'N elements = ', elements%kolik

      error stop

    end if


    !----------------------------------------------------------
    ! Total water routed from this element:
    !
    ! upstream inflow + local excess runoff
    !----------------------------------------------------------

    routed_water = &
         max(elements%hydrobal(el)%inflow, &
             0.0_rkind) + &
         max(runoff(el), 0.0_rkind)


    elements%hydrobal(el)%outflow = &
         routed_water


    if (routed_water > 0.0_rkind) then


      !========================================================
      ! Determine total downstream routing weight
      !========================================================

      weight_sum = 0.0_rkind


      do side = 1, 3


        dwn = &
             elements%downstream(el)%els(side)


        if (dwn > 0_ikind) then


          !----------------------------------------------------
          ! Check downstream element number
          !----------------------------------------------------

          if (dwn > elements%kolik) then

            write(*,*) 'ERROR in route_step'
            write(*,*) 'Invalid downstream element.'
            write(*,*) 'Current element    = ', el
            write(*,*) 'Side               = ', side
            write(*,*) 'Downstream element = ', dwn

            error stop

          end if


          !----------------------------------------------------
          ! Routing weight:
          !
          ! slope * effective edge width
          !----------------------------------------------------

          edge_weight = &
               max(elements%downstream(el)%slopes(side), &
                   0.0_rkind) * &
               max(elements%downstream(el)%widths(side), &
                   min_edge_width)


          weight_sum = &
               weight_sum + edge_weight


        end if


      end do


      !========================================================
      ! ROUTE TO DOWNSTREAM ELEMENTS
      !========================================================

      if (weight_sum > 0.0_rkind) then


        do side = 1, 3


          dwn = &
               elements%downstream(el)%els(side)


          if (dwn > 0_ikind) then


            edge_weight = &
                 max(elements%downstream(el)%slopes(side), &
                     0.0_rkind) * &
                 max(elements%downstream(el)%widths(side), &
                     min_edge_width)


            routed_fraction = &
                 edge_weight / weight_sum


            elements%hydrobal(dwn)%inflow = &
                 elements%hydrobal(dwn)%inflow + &
                 routed_water * routed_fraction


          end if


        end do


      else


        !======================================================
        ! NO VALID DOWNSTREAM ROUTE
        !
        ! Water leaves the modeled catchment.
        !======================================================

        outlet_Q(tstep) = &
             outlet_Q(tstep) + &
             routed_water


      end if


    end if


  end do


  !============================================================
  ! PART 3
  !
  ! SAVE ROUTING STATE AND RESULTS
  !============================================================

  do el = 1, elements%kolik


    !----------------------------------------------------------
    ! Update persistent routing storage
    !----------------------------------------------------------

    storage(el) = &
         max(new_storage(el), 0.0_rkind)


    !----------------------------------------------------------
    ! Local overflow / excess runoff
    !----------------------------------------------------------

    elements%overflow(el) = &
         max(runoff(el), 0.0_rkind)


    !----------------------------------------------------------
    ! Routing conservation check
    !
    ! Qin + local runoff - Qout
    !
    ! Should be approximately zero.
    !----------------------------------------------------------

    routing_residual(el) = &
         elements%hydrobal(el)%inflow + &
         runoff(el) - &
         elements%hydrobal(el)%outflow


    !----------------------------------------------------------
    ! Combined numerical residual
    !----------------------------------------------------------

    elements%hydrobal(el)%deltas = &
         local_residual(el) + &
         routing_residual(el)


    !----------------------------------------------------------
    ! Save routing results for output after compute_all()
    !----------------------------------------------------------

    deltas(el,tstep) = &
         elements%hydrobal(el)%deltas


    Qin_result(el,tstep) = &
         elements%hydrobal(el)%inflow


    Qout_result(el,tstep) = &
         elements%hydrobal(el)%outflow


    Overflow_result(el,tstep) = &
         elements%overflow(el)


    Storage_result(el,tstep) = &
         storage(el)


  end do


  !============================================================
  ! PART 4
  !
  ! CONVERT CATCHMENT OUTLET DEPTH TO DISCHARGE
  !============================================================

  catchment_area = &
       sum(elements%area)


  if (catchment_area < 0.0_rkind) then

    write(*,*) 'ERROR in route_step'
    write(*,*) 'Negative catchment area.'

    error stop

  end if


  ! outlet_Q is expressed as an equivalent water depth [mm].
  !
  ! mm -> m:
  !     outlet_Q / 1000
  !
  ! volume:
  !     depth [m] * catchment area [m2]
  !
  outlet_volume = &
       outlet_Q(tstep) / &
       1000.0_rkind * &
       catchment_area


  !------------------------------------------------------------
  ! Convert volume per timestep to m3/s
  !------------------------------------------------------------

  if (dt_seconds > 0.0_rkind) then

    outlet_Q_m3s(tstep) = &
         outlet_volume / dt_seconds

  else

    outlet_Q_m3s(tstep) = &
         0.0_rkind

  end if


  !============================================================
  ! OPTIONAL WATER-BALANCE WARNING
  !============================================================

  do el = 1, elements%kolik

    if (abs(deltas(el,tstep)) > &
        1.0e-8_rkind) then

      write(*,*) &
           'WARNING routing residual:', &
           ' dt=', dt_hours, &
           ' step=', tstep, &
           ' el=', el, &
           ' residual=', deltas(el,tstep)

    end if

  end do


  !============================================================
  ! CLEAN UP TEMPORARY ARRAYS
  !============================================================

  deallocate(runoff)

  deallocate(new_storage)

  deallocate(local_residual)

  deallocate(routing_residual)


end subroutine route_step

end module routing