module solver
  use typy
  use globals
  use hydrotools
  use tools
  use hydrofnc
  use routing, only: route_step
  implicit none

contains

  subroutine update_storage_odes(tstep)
    integer(kind=ikind), intent(in) :: tstep
    integer(kind=ikind) :: el

    do el = 1, elements%kolik

      ! ----------------------------------------------------
      ! All fluxes here are mm per timestep
      ! ----------------------------------------------------
      Pm(el,tstep)      = precip_m_step(el,tstep)
      E_m(el,tstep)     = evap_m_step(el,tstep)
      If_m(el,tstep)    = infil_m_step(el,tstep)
      Qsurf_m(el,tstep) = qsurf_m_step(el,tstep)
      q2(el,tstep)      = q2_m_step(el,tstep)
      q3(el,tstep)      = q3_m_step(el,tstep)
      pc(el,tstep)      = pc_m_step(el,tstep)
      bf(el,tstep)      = bf_m_step(el,tstep)

      ! Optional aliases used by older output/routing code
      inf(el,tstep)     = If_m(el,tstep)
      Qsub_m(el,tstep)  = q2(el,tstep) + q3(el,tstep)
      Qgw_m(el,tstep)   = bf(el,tstep)

      ! ----------------------------------------------------
      ! Storage increments as equivalent depth [mm/timestep]
      ! ----------------------------------------------------
      dVsurf(el,tstep) = Pm(el,tstep) - If_m(el,tstep) - &
                         E_m(el,tstep) - Qsurf_m(el,tstep)

      dVsub(el,tstep)  = If_m(el,tstep) - q2(el,tstep) - &
                         q3(el,tstep) - pc(el,tstep)

      dVgw(el,tstep)   = pc(el,tstep) - bf(el,tstep)

      dv_total = dVsurf(el,tstep) + dVsub(el,tstep) + dVgw(el,tstep)

      ! ----------------------------------------------------
      ! Routing variables: mm per timestep
      ! ----------------------------------------------------
      ET_flux(el,tstep)      = E_m(el,tstep)
      Inf_result(el,tstep)   = If_m(el,tstep)
      Qsurf_result(el,tstep) = Qsurf_m(el,tstep)
      Qgw_result(el,tstep)   = Qgw_m(el,tstep)

      elements%hydrobal(el)%ET    = ET_flux(el,tstep)
      elements%hydrobal(el)%Qsurf = Qsurf_result(el,tstep)
      elements%hydrobal(el)%Qgw   = Qgw_result(el,tstep)
      elements%hydrobal(el)%q2    = q2(el,tstep)
      elements%hydrobal(el)%q3    = q3(el,tstep)
      elements%hydrobal(el)%pc    = pc(el,tstep)
      elements%hydrobal(el)%bf    = bf(el,tstep)

      elements%hydrobal(el)%inflow  = 0.0_rkind
      elements%hydrobal(el)%outflow = 0.0_rkind
      elements%hydrobal(el)%deltas  = 0.0_rkind

    end do
  end subroutine update_storage_odes


  subroutine compute_all()
    integer(kind=ikind) :: t

    do t = 1, n_steps
      call update_storage_odes(t)
      call route_step(t)
    end do

    print *, "DEBUG after compute:"
    print *, "Pm(1,49)  = ", Pm(1,49)
    print *, "Pm(1,72)  = ", Pm(1,72)
    print *, "Pm(1,240) = ", Pm(1,240)

  end subroutine compute_all

end module solver
