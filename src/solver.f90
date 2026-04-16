module solver

  contains

 !==============================================================
  !  Main driver: compute all steps with routing
  !==============================================================
  subroutine compute_all()
    use tools
    integer(kind=ikind) :: el, t
    real(kind=rkind)    :: ET_pot, ET_act
    real(kind=rkind)    :: Qsurf_loc
    real(kind=rkind)    :: Qgw_pot, Qgw_act
    real(kind=rkind)    :: water_available
    real(kind=rkind)    :: water_after_ET

    do 
      t = t + 1
      time  = time + time_step

      do el = 1, elements%kolik

        ! -----------------------------------------------
        ! 1) Potential fluxes
        ! -----------------------------------------------
        ET_pot    = penman_monteith(el, t) * ccrop
        Qsurf_loc = surface_runoff(el, t)
        Qgw_pot   = ground_water(el, t)

        ! -----------------------------------------------
        ! 2) Water available after quick runoff
        ! -----------------------------------------------
        water_available = precip(el, t) + qinter(el, t) + storage(el) - Qsurf_loc

        if (water_available < 0.0_rkind) then
          water_available = 0.0_rkind
        end if

        ! -----------------------------------------------
        ! 3) Limit actual ET by available water
        ! -----------------------------------------------
        ET_act = min(ET_pot, water_available)

        ! -----------------------------------------------
        ! 4) Water left after ET
        ! -----------------------------------------------
        water_after_ET = water_available - ET_act

        if (water_after_ET < 0.0_rkind) then
          water_after_ET = 0.0_rkind
        end if

        ! -----------------------------------------------
        ! 5) Limit actual groundwater flow by remaining water
        ! -----------------------------------------------
        Qgw_act = min(Qgw_pot, water_after_ET)

        ! -----------------------------------------------
        ! 6) Store local fluxes
        ! -----------------------------------------------
        elements%hydrobal(el)%ET    = ET_act
        elements%hydrobal(el)%Qsurf = Qsurf_loc
        elements%hydrobal(el)%Qgw   = Qgw_act

        ET_flux(el, t)      = elements%hydrobal(el)%ET
        Qsurf_result(el, t) = elements%hydrobal(el)%Qsurf
        Qgw_result(el, t)   = elements%hydrobal(el)%Qgw

        ! -----------------------------------------------
        ! 7) Leakage after ET, runoff, and groundwater are known
        ! -----------------------------------------------
        elements%hydrobal(el)%Li = max(0.0_rkind, &
             qinter(el, t) + precip(el, t) - &
             ET_flux(el, t) - Qsurf_result(el, t) - Qgw_result(el, t))

        L_result(el, t) = elements%hydrobal(el)%Li

        ! Reset routed terms for this step
        elements%hydrobal(el)%inflow  = 0.0_rkind
        elements%hydrobal(el)%outflow = 0.0_rkind
        elements%hydrobal(el)%deltas  = 0.0_rkind

      end do

      call route_step(t)
      
      if (time > end_time) EXIT

    end do
  end subroutine compute_all


end module solver
