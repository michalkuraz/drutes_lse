module hydrofnc
  use typy
  use globals
  implicit none

  real(kind=rkind), parameter :: PI_RKIND = 4.0_rkind * atan(1.0_rkind)

contains

  !==============================================================
  ! Saturation vapor pressure [kPa]
  !==============================================================
  pure real(kind=rkind) function esat(T)
    real(kind=rkind), intent(in) :: T

    esat = 0.6108_rkind * exp(17.27_rkind * T / (T + 237.3_rkind))
  end function esat


  !==============================================================
  ! Slope of saturation vapor pressure curve [kPa/degC]
  !==============================================================
  pure real(kind=rkind) function slope_vp(T)
    real(kind=rkind), intent(in) :: T

    slope_vp = 4098.0_rkind * esat(T) / ((T + 237.3_rkind)**2)
  end function slope_vp


  !==============================================================
  ! Psychrometric constant [kPa/degC]
  !==============================================================
  pure real(kind=rkind) function psychro_const(z_elev)
    real(kind=rkind), intent(in) :: z_elev
    real(kind=rkind) :: P

    P = 101.3_rkind * ((293.0_rkind - 0.0065_rkind * z_elev) / 293.0_rkind)**5.26_rkind
    psychro_const = 0.000665_rkind * P
  end function psychro_const


  !==============================================================
  ! Convert measured wind speed to wind speed at 2 m [m/s]
  !==============================================================
  pure real(kind=rkind) function wind2m(uz_meas, z_meas)
    real(kind=rkind), intent(in) :: uz_meas, z_meas

    wind2m = uz_meas * 4.87_rkind / log(67.8_rkind * z_meas - 5.42_rkind)
  end function wind2m


  !==============================================================
  ! Extraterrestrial radiation Ra [MJ m-2 day-1]
  !==============================================================
  pure real(kind=rkind) function Ra_daily(phi_loc, jday)
    real(kind=rkind), intent(in) :: phi_loc
    integer(kind=ikind), intent(in) :: jday
    real(kind=rkind) :: dr, delta, ws

    dr = 1.0_rkind + 0.033_rkind * cos(2.0_rkind * PI_RKIND * real(jday, rkind) / 365.0_rkind)

    delta = 0.409_rkind * sin(2.0_rkind * PI_RKIND * real(jday, rkind) / 365.0_rkind - 1.39_rkind)

    ws = acos(-tan(phi_loc) * tan(delta))

    Ra_daily = (24.0_rkind * 60.0_rkind / PI_RKIND) * gsc * dr * &
               (ws * sin(phi_loc) * sin(delta) + &
                cos(phi_loc) * cos(delta) * sin(ws))
  end function Ra_daily


  !==============================================================
  ! Net longwave radiation Rnl [MJ m-2 day-1]
  !==============================================================
  pure real(kind=rkind) function Rnl_daily(TmaxK, TminK, ea, Rs, Rso)
    real(kind=rkind), intent(in) :: TmaxK, TminK, ea, Rs, Rso
    real(kind=rkind) :: tmp, ratio

    tmp = 0.5_rkind * (TmaxK**4 + TminK**4)

    if (Rso > 0.0_rkind) then
      ratio = min(Rs / Rso, 1.0_rkind)
    else
      ratio = 0.0_rkind
    end if

    Rnl_daily = sigma * tmp * &
                (0.34_rkind - 0.14_rkind * sqrt(max(ea, 0.0_rkind))) * &
                (1.35_rkind * ratio - 0.35_rkind)
  end function Rnl_daily


  !==============================================================
  ! Penman-Monteith reference evapotranspiration ETo [mm/day]
  !==============================================================
  real(kind=rkind) function penman_monteith(element, tstep)
    integer(kind=ikind), intent(in) :: element, tstep

    integer(kind=ikind) :: jday
    real(kind=rkind) :: es_Tmax, es_Tmin, es, ea
    real(kind=rkind) :: delta, gamma, u2, Ra, Rs, Rso
    real(kind=rkind) :: Rns, Rnl, Rn, G0
    real(kind=rkind) :: L1, L2, L3, L4, nN

    jday = Julian_day + int(real(tstep - 1, rkind) * dt_days, kind=ikind)

    es_Tmax = esat(Tmax(element, tstep))
    es_Tmin = esat(Tmin(element, tstep))
    es      = 0.5_rkind * (es_Tmax + es_Tmin)

    ea = (es_Tmin * RHmax(element, tstep) + &
          es_Tmax * RHmin(element, tstep)) / 200.0_rkind

    delta = slope_vp(Tmean(element, tstep))
    gamma = psychro_const(z)
    u2    = wind2m(uz(element, tstep), z)

    Ra = Ra_daily(phi, jday)

    nN  = 1.0_rkind
    Rs  = (as + bs * nN) * Ra
    Rso = (0.75_rkind + 2.0e-5_rkind * z) * Ra

    Rns = (1.0_rkind - alpha) * Rs

    Rnl = Rnl_daily(Tmax(element, tstep) + 273.16_rkind, &
                    Tmin(element, tstep) + 273.16_rkind, &
                    ea, Rs, Rso)

    Rn = max(Rns - Rnl, 0.0_rkind)
    G0 = G(element)

    L1 = 0.408_rkind * delta * (Rn - G0)

    L2 = gamma * (900.0_rkind / (Tmean(element, tstep) + 273.0_rkind)) * &
         u2 * (es - ea)

    L3 = L1 + L2
    L4 = delta + gamma * (1.0_rkind + 0.34_rkind * u2)

    if (L4 > 0.0_rkind) then
      penman_monteith = max(L3 / L4, 0.0_rkind)
    else
      penman_monteith = 0.0_rkind
    end if
  end function penman_monteith


  !==============================================================
  ! Effective soil saturation factor [-]
  !==============================================================
  real(kind=rkind) function theta_effective(theta)
    real(kind=rkind), intent(in) :: theta

    if (theta_s > theta_r) then
      theta_effective = (theta - theta_r) / (theta_s - theta_r)
      theta_effective = max(0.0_rkind, min(1.0_rkind, theta_effective))
    else
      theta_effective = 1.0_rkind
    end if
  end function theta_effective


  !==============================================================
  ! Unsaturated surface hydraulic conductivity [m/s]
  !==============================================================
  real(kind=rkind) function ksurf_unsat(element, tstep)
    integer(kind=ikind), intent(in) :: element, tstep

    ksurf_unsat = Ksat_surf(element) * &
                  theta_effective(soilcontent(element,tstep))**ksurf_exp
  end function ksurf_unsat


  !==============================================================
  ! Soil moisture reduction factor for infiltration [-]
  ! High moisture gives lower infiltration capacity
  !==============================================================
  real(kind=rkind) function infil_moisture_factor(theta)
    real(kind=rkind), intent(in) :: theta

    if (theta_s > theta_r) then
      infil_moisture_factor = (theta_s - theta) / (theta_s - theta_r)
      infil_moisture_factor = max(0.0_rkind, min(1.0_rkind, infil_moisture_factor))
    else
      infil_moisture_factor = 1.0_rkind
    end if
  end function infil_moisture_factor


  !==============================================================
  ! Slope reduction factor for infiltration [-]
  ! Steeper slopes reduce infiltration opportunity
  !==============================================================
  real(kind=rkind) function infil_slope_factor(element)
    integer(kind=ikind), intent(in) :: element

    infil_slope_factor = 1.0_rkind / &
         (1.0_rkind + infil_slope_coeff * max(elements%slope(element), 0.0_rkind))
  end function infil_slope_factor


  !==============================================================
  ! Slope-adjusted Curve Number [-]
  !==============================================================
  real(kind=rkind) function effective_cn(element)
    integer(kind=ikind), intent(in) :: element
    real(kind=rkind) :: slope_pct

    slope_pct = 100.0_rkind * max(elements%slope(element), 0.0_rkind)

    effective_cn = real(CN, rkind) + cn_slope_coeff * slope_pct
    effective_cn = min(98.0_rkind, max(30.0_rkind, effective_cn))
  end function effective_cn


  !==============================================================
  ! P: Precipitation [mm/timestep]
  !==============================================================
  real(kind=rkind) function precip_m_step(element, tstep)
    integer(kind=ikind), intent(in) :: element, tstep

    precip_m_step = precip(element,tstep)
  end function precip_m_step


  !==============================================================
  ! E: Actual evapotranspiration [mm/timestep]
  ! E = ETo * soil moisture factor * dt_days
  !==============================================================
  real(kind=rkind) function evap_m_step(element, tstep)
    integer(kind=ikind), intent(in) :: element, tstep

    evap_m_step = penman_monteith(element,tstep) * &
                  theta_effective(soilcontent(element,tstep)) * dt_days
  end function evap_m_step


  !==============================================================
  ! If: Infiltration [mm/timestep]
  !==============================================================
  real(kind=rkind) function infil_m_step(element, tstep)
    integer(kind=ikind), intent(in) :: element, tstep
    real(kind=rkind) :: cap

    cap = ksurf_unsat(element,tstep) * 1000.0_rkind * dt_seconds * &
          infil_moisture_factor(soilcontent(element,tstep)) * &
          infil_slope_factor(element)

    infil_m_step = min(precip_m_step(element,tstep), max(0.0_rkind, cap))
  end function infil_m_step


  !==============================================================
  ! q1: Surface runoff [mm/timestep]
  ! Slope-adjusted SCS-CN method
  !==============================================================
  real(kind=rkind) function q1_m_step(element, tstep)
    integer(kind=ikind), intent(in) :: element, tstep
    real(kind=rkind) :: CN_eff, S, Ia, Pstep_mm, denom

    CN_eff = effective_cn(element)

    if (CN_eff <= 0.0_rkind) then
      q1_m_step = 0.0_rkind
      return
    end if

    Pstep_mm = max(0.0_rkind, precip(element,tstep))

    S  = (25400.0_rkind / CN_eff) - 254.0_rkind
    Ia = 0.2_rkind * S

    if (S <= 0.0_rkind .or. Pstep_mm <= Ia) then
      q1_m_step = 0.0_rkind
    else
      denom = Pstep_mm - Ia + S

      if (denom <= 0.0_rkind) then
        q1_m_step = 0.0_rkind
      else
        q1_m_step = ((Pstep_mm - Ia)**2 / denom)
      end if
    end if
  end function q1_m_step


  !==============================================================
  ! q2: Subsurface flow component (fast interflow) [mm/timestep]
  ! Uses previous subsurface storage Ssub
  !==============================================================
  real(kind=rkind) function q2_m_step(element, tstep)
    integer(kind=ikind), intent(in) :: element, tstep
    real(kind=rkind) :: Ssub_old

    if (tstep <= 1_ikind) then
      Ssub_old = Ssub(element)
    else
      Ssub_old = Ssub_hist(element,tstep-1)
    end if

    q2_m_step = max(infil_m_step(element,tstep)**2 + Ssub_old - &
                    Beta1 * (z2 - z1), 0.0_rkind)
  end function q2_m_step


  !==============================================================
  ! q3: Subsurface drainage component(slow interflow) [mm/timestep]
  ! Uses previous subsurface storage Ssub
  !==============================================================
  real(kind=rkind) function q3_m_step(element, tstep)
    integer(kind=ikind), intent(in) :: element, tstep
    real(kind=rkind) :: Ssub_old

    if (tstep <= 1_ikind) then
      Ssub_old = Ssub(element)
    else
      Ssub_old = Ssub_hist(element,tstep-1)
    end if

    q3_m_step = Beta2 * max(Ssub_old, 0.0_rkind)**Beta3
  end function q3_m_step


  !==============================================================
  ! pc: Percolation from subsurface to groundwater [mm/timestep]
  ! Uses previous subsurface storage Ssub
  !==============================================================
  real(kind=rkind) function pc_m_step(element, tstep)
    integer(kind=ikind), intent(in) :: element, tstep
    real(kind=rkind) :: Ssub_old

    if (tstep <= 1_ikind) then
      Ssub_old = Ssub(element)
    else
      Ssub_old = Ssub_hist(element,tstep-1)
    end if

    pc_m_step = Beta4 * max(Ssub_old, 0.0_rkind)
  end function pc_m_step


  !==============================================================
  ! bf: Baseflow from groundwater [mm/timestep]
  ! Uses previous groundwater storage Sgw
  !==============================================================
  real(kind=rkind) function bf_m_step(element, tstep)
    integer(kind=ikind), intent(in) :: element, tstep
    real(kind=rkind) :: Sgw_old

    if (tstep <= 1_ikind) then
      Sgw_old = Sgw(element)
    else
      Sgw_old = Sgw_hist(element,tstep-1)
    end if

    bf_m_step = Beta5 * max(Sgw_old, 0.0_rkind)
  end function bf_m_step

end module hydrofnc