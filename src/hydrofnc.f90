module hydrofnc
  use typy
  use globals
  implicit none

  real(kind=rkind), parameter :: PI_RKIND = 4.0_rkind * atan(1.0_rkind)

contains

  pure real(kind=rkind) function esat(T)
    real(kind=rkind), intent(in) :: T
    esat = 0.6108_rkind * exp(17.27_rkind * T / (T + 237.3_rkind))
  end function esat


  pure real(kind=rkind) function slope_vp(T)
    real(kind=rkind), intent(in) :: T
    slope_vp = 4098.0_rkind * esat(T) / ((T + 237.3_rkind)**2)
  end function slope_vp


  pure real(kind=rkind) function psychro_const(z_elev)
    real(kind=rkind), intent(in) :: z_elev
    real(kind=rkind) :: P
    P = 101.3_rkind * ((293.0_rkind - 0.0065_rkind * z_elev) / 293.0_rkind)**5.26_rkind
    psychro_const = 0.000665_rkind * P
  end function psychro_const


  pure real(kind=rkind) function wind2m(uz_meas, z_meas)
    real(kind=rkind), intent(in) :: uz_meas, z_meas
    wind2m = uz_meas * 4.87_rkind / log(67.8_rkind * z_meas - 5.42_rkind)
  end function wind2m


  pure real(kind=rkind) function Ra_daily(phi_loc, jday)
    real(kind=rkind), intent(in) :: phi_loc
    integer(kind=ikind), intent(in) :: jday
    real(kind=rkind) :: dr, delta, ws

    dr    = 1.0_rkind + 0.033_rkind * cos(2.0_rkind * PI_RKIND * real(jday, rkind) / 365.0_rkind)
    delta = 0.409_rkind * sin(2.0_rkind * PI_RKIND * real(jday, rkind) / 365.0_rkind - 1.39_rkind)
    ws    = acos(-tan(phi_loc) * tan(delta))

    Ra_daily = (24.0_rkind * 60.0_rkind / PI_RKIND) * gsc * dr * &
               (ws * sin(phi_loc) * sin(delta) + cos(phi_loc) * cos(delta) * sin(ws))
  end function Ra_daily


  pure real(kind=rkind) function Rnl_daily(TmaxK, TminK, ea, Rs, Rso)
    real(kind=rkind), intent(in) :: TmaxK, TminK, ea, Rs, Rso
    real(kind=rkind) :: tmp, ratio

    tmp = 0.5_rkind * (TmaxK**4 + TminK**4)

    if (Rso > 0.0_rkind) then
      ratio = min(Rs / Rso, 1.0_rkind)
    else
      ratio = 0.0_rkind
    end if

    Rnl_daily = sigma * tmp * (0.34_rkind - 0.14_rkind * sqrt(max(ea, 0.0_rkind))) * &
                (1.35_rkind * ratio - 0.35_rkind)
  end function Rnl_daily


  !==============================================================
  !  Penman Monteith function for ET estimation
  !==============================================================
  pure real(kind=rkind) function penman_monteith(element, tstep)
    integer(kind=ikind), intent(in) :: element, tstep
    integer(kind=ikind) :: jday
    real(kind=rkind) :: es_Tmax, es_Tmin, es, ea
    real(kind=rkind) :: delta, gamma, u2, Ra, Rs, Rso
    real(kind=rkind) :: Rns, Rnl, Rn, G0
    real(kind=rkind) :: L1, L2, L3, L4, nN

    jday = Julian_day + int((real(tstep - 1, rkind)) * dt_days, kind=ikind)

    es_Tmax = esat(Tmax(element, tstep))
    es_Tmin = esat(Tmin(element, tstep))

    es = 0.5_rkind * (es_Tmax + es_Tmin)

    ea = (es_Tmin * RHmax(element, tstep) + es_Tmax * RHmin(element, tstep)) / 200.0_rkind

    delta = slope_vp(Tmean(element, tstep))
    gamma = psychro_const(z)

    u2 = wind2m(uz(element, tstep), z)

    Ra = Ra_daily(phi, jday)

    nN = 1.0_rkind
    Rs  = (as + bs * nN) * Ra
    Rso = (0.75_rkind + 2.0e-5_rkind * z) * Ra

    Rns = (1.0_rkind - alpha) * Rs

    Rnl = Rnl_daily(Tmax(element, tstep) + 273.16_rkind, &
                    Tmin(element, tstep) + 273.16_rkind, ea, Rs, Rso)

    Rn = max(Rns - Rnl, 0.0_rkind)

    G0 = G(element)

    L1 = 0.408_rkind * delta * (Rn - G0)
    L2 = gamma * (900.0_rkind / (Tmean(element, tstep) + 273.0_rkind)) * u2 * (es - ea)
    L3 = L1 + L2
    L4 = delta + gamma * (1.0_rkind + 0.34_rkind * u2)

    if (L4 > 0.0_rkind) then
      penman_monteith = max(L3 / L4, 0.0_rkind)
    else
      penman_monteith = 0.0_rkind
    end if
  end function penman_monteith


  !==============================================================
  !  Surface runoff using SCS-CN method
  !==============================================================
  pure real(kind=rkind) function surface_runoff(element, tstep)
    integer(kind=ikind), intent(in) :: element, tstep
    real(kind=rkind) :: S, I

    S = (25400.0_rkind / real(CN, rkind)) - 254.0_rkind
    I = 0.2_rkind * S

    if (precip(element, tstep) <= I) then
      surface_runoff = 0.0_rkind
    else
      surface_runoff = (precip(element, tstep) - I)**2 / &
                       (precip(element, tstep) - I + S)
    end if
  end function surface_runoff


  !==============================================================
  !  Ground water function
  !==============================================================
  pure real(kind=rkind) function ground_water(element, tstep)
    integer(kind=ikind), intent(in) :: element, tstep
    ground_water = conduct(element) * soilcontent(element, tstep)
  end function ground_water


  !==============================================================
  !  Leakage function
  !==============================================================
  pure real(kind=rkind) function leakage(element, tstep)
    integer(kind=ikind), intent(in) :: element, tstep
    leakage = max(0.0_rkind, qinter(element, tstep) + precip(element, tstep) - &
                              ET_flux(element, tstep) - Qsurf_result(element, tstep))
  end function leakage


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

  do t = 1, n_steps

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

  end do
end subroutine compute_all
end module hydrofnc