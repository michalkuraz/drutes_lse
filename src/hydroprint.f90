module hydroprint
  use typy
  use globals
  use hydrotools
  use tools
  use hydrofnc
  use solver
  use routing
  implicit none

contains

  subroutine print_upstream_flows(tstep)
    integer(kind=ikind), intent(in) :: tstep
    integer(kind=ikind) :: nel, e, k
    real(kind=rkind) :: Qin_from_up

    nel = elements%kolik

    if (.not. allocated(Qin_result)) return
    if (.not. allocated(Qout_result)) return
    if (.not. allocated(upstream_count)) return

    print *, "---------------------------------------------------------------"
    print *, " Time step = ", tstep
    print *, " el | n_up | Qin_from_up    Qin(field)     Qout_elem"
    print *, "---------------------------------------------------------------"

    do e = 1, nel
      Qin_from_up = 0.0_rkind

      do k = 1, upstream_count(e)
        Qin_from_up = Qin_from_up + Qout_result(upstream_list(e,k), tstep)
      end do

      print *, e, upstream_count(e), Qin_from_up, Qin_result(e,tstep), Qout_result(e,tstep)
    end do
  end subroutine print_upstream_flows


  subroutine print_graph_diagnostics()
    integer(kind=ikind) :: el, i, nb
    real(kind=rkind) :: width, slope, nalt

    print *, "--------------------------------------------------------------------------"
    print *, " el   nb   dir      z_el        z_nb        slope        edge_width"
    print *, "--------------------------------------------------------------------------"

    do el = 1, elements%kolik
      do i = 1, 3
        nb = elements%downstream(el)%els(i)

        if (nb > 0) then
          width = elements%downstream(el)%widths(i)
          slope = elements%downstream(el)%slopes(i)
          nalt  = elements%avgalt(nb)
          print *, el, nb, "down", elements%avgalt(el), nalt, slope, width
        end if

        nb = elements%upstream(el)%els(i)

        if (nb > 0) then
          width = elements%upstream(el)%widths(i)
          slope = elements%upstream(el)%slopes(i)
          nalt  = elements%avgalt(nb)
          print *, el, nb, "up", elements%avgalt(el), nalt, slope, width
        end if
      end do
    end do

    print *, "--------------------------------------------------------------------------"
  end subroutine print_graph_diagnostics


  subroutine print_water_balance(tstep)
    integer(kind=ikind), intent(in) :: tstep
    integer(kind=ikind) :: i

    print *
    print *, "================================================================================================="
    write(*,'(A,I6)') "WATER BALANCE AT TIME STEP = ", tstep
    print *, "================================================================================================="
    print *, "Elem        P          E          If         q1         q2         q3         pc         bf" // &
         "       dSsurf      dSsub       dSgw       total_dS    Ssurf      Ssub       Sgw"
    print *, "================================================================================================="

    do i = 1, elements%kolik
      print*,  i, &
        Pm(i,tstep), E_m(i,tstep), If_m(i,tstep), &
        q1(i,tstep), q2(i,tstep), q3(i,tstep), pc(i,tstep), bf(i,tstep), &
        dSsurf(i,tstep), dSsub(i,tstep), dSgw(i,tstep), total_deltaS(i,tstep), &
        Ssurf_hist(i,tstep), Ssub_hist(i,tstep), Sgw_hist(i,tstep)
    end do

    print *, "================================================================================================="
  end subroutine print_water_balance


  subroutine export_element_balance(filename)
    character(len=*), intent(in) :: filename
    integer(kind=ikind) :: i, t
    integer :: unit, ios

    open(newunit=unit, file=filename, status="replace", action="write", iostat=ios)

    if (ios /= 0) then
      print *, "Error opening file: ", trim(filename)
      return
    end if

    print*, "step,element,area,z_avg,downstream,P,E,If,q1,q2,q3,pc,bf," // &
                  "dSsurf,dSsub,dSgw,total_dS,Ssurf,Ssub,Sgw,Qin,Qout,Overflow,Storage,deltas"

    do t = 1, n_steps
      do i = 1, elements%kolik
        print*, &
            t, i, elements%area(i), elements%avgalt(i), downstream(i), &
           Pm(i,t), E_m(i,t), If_m(i,t), q1(i,t), q2(i,t), q3(i,t), pc(i,t), bf(i,t), &
           dSsurf(i,t), dSsub(i,t), dSgw(i,t), total_deltaS(i,t), &
           Ssurf_hist(i,t), Ssub_hist(i,t), Sgw_hist(i,t), &
           Qin_result(i,t), Qout_result(i,t), Overflow_result(i,t), &
           Storage_result(i,t), deltas(i,t)
      end do
    end do

    close(unit)
    print *, "Element balances saved to: ", trim(filename)
  end subroutine export_element_balance

end module hydroprint