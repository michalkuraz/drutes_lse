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
    real(kind=rkind) :: Ai, Ssub_mm, Sgw_mm

    do el = 1, elements%kolik

      Ai = elements%area(el)

      ! ----------------------------------------------------
      ! All fluxes here are mm per timestep
      ! ----------------------------------------------------
     Pm(el,tstep)      = precip_m_step(el,tstep)
     E_m(el,tstep)     = evap_m_step(el,tstep)
     If_m(el,tstep)    = min(infil_m_step(el,tstep), Pm(el,tstep))
     Qsurf_m(el,tstep) = qsurf_m_step(el,tstep)

    ! Convert current storages to mm
     Ssub_mm = Ssub(el) / Ai * 1000.0_rkind
     Sgw_mm  = Sgw(el)  / Ai * 1000.0_rkind

    ! Add this-step infiltration to available subsurface water
      Ssub_mm = Ssub_mm + If_m(el,tstep)

    ! Use small fractions per timestep, not the full potential
      Tv_m(el,tstep)   = min(0.02_rkind * Ssub_mm, Ssub_mm)
      Qsub_m(el,tstep) = min(0.05_rkind * Ssub_mm, Ssub_mm - Tv_m(el,tstep))
      Rv_m(el,tstep)   = min(0.03_rkind * Ssub_mm, Ssub_mm - Tv_m(el,tstep) - Qsub_m(el,tstep))

    ! Groundwater receives recharge, then discharges slowly
      Sgw_mm = Sgw_mm + Rv_m(el,tstep)
      Qgw_m(el,tstep) = min(0.02_rkind * Sgw_mm, Sgw_mm)

      ! ----------------------------------------------------
      ! Storage increments in m3
      ! depth mm -> m by /1000
      ! ----------------------------------------------------
      dVsurf(el,tstep) = (Pm(el,tstep) - If_m(el,tstep) - &
                          E_m(el,tstep) - Qsurf_m(el,tstep)) / &
                          1000.0_rkind * Ai

      dVsub(el,tstep) = (If_m(el,tstep) - Tv_m(el,tstep) - &
                         Qsub_m(el,tstep) - Rv_m(el,tstep)) / &
                         1000.0_rkind * Ai

      dVgw(el,tstep) = (Rv_m(el,tstep) - Qgw_m(el,tstep)) / &
                       1000.0_rkind * Ai

      Ssurf(el) = max(0.0_rkind, Ssurf(el) + dVsurf(el,tstep))
      Ssub(el)  = max(0.0_rkind, Ssub(el)  + dVsub(el,tstep))
      Sgw(el)   = max(0.0_rkind, Sgw(el)   + dVgw(el,tstep))

      Ssurf_hist(el,tstep) = Ssurf(el)
      Ssub_hist(el,tstep)  = Ssub(el)
      Sgw_hist(el,tstep)   = Sgw(el)

      ! ----------------------------------------------------
      ! Routing variables: mm per timestep
      ! ----------------------------------------------------
      ET_flux(el,tstep)      = E_m(el,tstep)
      Inf_result(el,tstep)   = If_m(el,tstep)
      Qsurf_result(el,tstep) = Qsurf_m(el,tstep)
      Qgw_result(el,tstep)   = Qgw_m(el,tstep)

       L_result(el,tstep) = If_m(el,tstep)

      elements%hydrobal(el)%ET    = ET_flux(el,tstep)
      elements%hydrobal(el)%Qsurf = Qsurf_result(el,tstep)
      elements%hydrobal(el)%Qgw   = Qgw_result(el,tstep)
      elements%hydrobal(el)%Li    = L_result(el,tstep)

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