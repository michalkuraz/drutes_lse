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

      ! ==========================================================
      ! Fluxes from the LaTeX ODE system [mm/timestep]
      ! ==========================================================
      Pm(el,tstep) = precip_m_step(el,tstep)
      E_m(el,tstep) = evap_m_step(el,tstep)
      If_m(el,tstep) = infil_m_step(el,tstep)

      q1(el,tstep) = q1_m_step(el,tstep)
      q2(el,tstep) = q2_m_step(el,tstep)
      q3(el,tstep) = q3_m_step(el,tstep)
      pc(el,tstep) = pc_m_step(el,tstep)
      bf(el,tstep) = bf_m_step(el,tstep)

      ! ==========================================================
      ! Direct ODE water balance in depth units [mm]
      !
      ! dSsurf = P - E - q1
      ! dSsub  = If - q2 - q3 - pc
      ! dSgw   = pc - bf
      ! ==========================================================

      dSsurf(el,tstep) = Pm(el,tstep) - E_m(el,tstep) - q1(el,tstep)

      dSsub(el,tstep) = If_m(el,tstep) - q2(el,tstep) - &
                        q3(el,tstep) - pc(el,tstep)

      dSgw(el,tstep) = pc(el,tstep) - bf(el,tstep)
      total_deltaS(el,tstep) = dSsurf(el,tstep) + dSsub(el,tstep) + dSgw(el,tstep)

      Ssurf(el) = max(0.0_rkind, Ssurf(el) + dSsurf(el,tstep))
      Ssub(el)  = max(0.0_rkind, Ssub(el)  + dSsub(el,tstep))
      Sgw(el)   = max(0.0_rkind, Sgw(el)   + dSgw(el,tstep))

      Ssurf_hist(el,tstep) = Ssurf(el)
      Ssub_hist(el,tstep)  = Ssub(el)
      Sgw_hist(el,tstep)   = Sgw(el)
      

      ! Values needed by routing
      elements%hydrobal(el)%ET      = E_m(el,tstep)
      elements%hydrobal(el)%q1      = q1(el,tstep)
      elements%hydrobal(el)%q2      = q2(el,tstep)
      elements%hydrobal(el)%q3      = q3(el,tstep)
      elements%hydrobal(el)%pc      = pc(el,tstep)
      elements%hydrobal(el)%bf      = bf(el,tstep)
      elements%hydrobal(el)%storage = Ssurf(el)

      elements%hydrobal(el)%inflow  = 0.0_rkind
      elements%hydrobal(el)%outflow = 0.0_rkind
      elements%hydrobal(el)%deltas  = 0.0_rkind

    end do
  end subroutine update_storage_odes


  subroutine compute_all()
    integer(kind=ikind) :: t

    time      = 0.0_rkind
    time_step = dt_seconds
    end_time  = ntot_days * 86400.0_rkind

    call timestamps%clear(.true.)

    do t = 1, n_steps
      time = time + time_step
      call timestamps%fill(time)

      call update_storage_odes(t)
      call route_step(t)

      if (time >= end_time) exit
    end do

  end subroutine compute_all

end module solver