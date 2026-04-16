module hydroprint

  contains

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


end module hydroprint
