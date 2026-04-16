module routing

contains

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


    call build_graph()


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


end module routing
