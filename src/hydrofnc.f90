! hydrofnc.f90
module hydrofnc
  use typy
  use globals
  implicit none

  real(kind=rkind), parameter :: PI_RKIND = 4*atan(1.0_rkind)

contains

  pure real(kind=rkind) function esat(T)
    real(kind=rkind), intent(in) :: T
    esat = 0.6108_rkind * exp(17.27_rkind*T / (T + 237.3_rkind))
  end function esat

  pure real(kind=rkind) function slope_vp(T)
    real(kind=rkind), intent(in) :: T
    slope_vp = 4098._rkind * esat(T) / ((T + 237.3_rkind)**2)
  end function slope_vp

  pure real(kind=rkind) function psychro_const(z_elev)
    real(kind=rkind), intent(in) :: z_elev
    real(kind=rkind) :: P
    P = 101.3_rkind * ((293._rkind - 0.0065_rkind*z_elev) / 293._rkind)**5.26_rkind
    psychro_const = 0.000665_rkind * P
  end function psychro_const

  pure real(kind=rkind) function wind2m(uz_meas, z_meas)
    real(kind=rkind), intent(in) :: uz_meas, z_meas
    wind2m = uz_meas * 4.87_rkind / log(67.8_rkind*z_meas - 5.42_rkind)
  end function wind2m

  pure real(kind=rkind) function Ra_daily(phi, Julian_day,dt_days)
    real(kind=rkind), intent(in) :: phi
    integer(kind=ikind), intent(in)          :: Julian_day,dt_days
    real(kind=rkind) :: dr, delta, ws

    dr    = 1._rkind + 0.033_rkind * cos(2._rkind*PI_RKIND*Julian_day/365._rkind)
    delta = 0.409_rkind * sin(2._rkind*PI_RKIND*Julian_day *dt_days/365._rkind - 1.39_rkind)
    ws    = acos(-tan(phi)*tan(delta))

    Ra_daily = (24._rkind*60._rkind/PI_RKIND) * gsc * dr * &
               (ws*sin(phi)*sin(delta) + cos(phi)*cos(delta)*sin(ws))
  end function Ra_daily

  pure real(kind=rkind) function Rnl_daily(TmaxK, TminK, ea, Rs, Rso)
    real(kind=rkind), intent(in) :: TmaxK, TminK, ea, Rs, Rso
    real(kind=rkind) :: tmp
    tmp = 0.5_rkind * (TmaxK**4 + TminK**4)
    Rnl_daily = sigma * tmp * (0.34_rkind - 0.14_rkind*sqrt(ea)) * &
                (1.35_rkind*min(Rs/Rso, 1._rkind) - 0.35_rkind)
  end function Rnl_daily

  !==============================================================
!  Penman Monteith function for ET estimation
!==============================================================
  pure real(kind=rkind) function penman_monteith(element, dt_days)
    integer(kind=ikind), intent(in) :: element, dt_days
    real(kind=rkind) :: es_Tmax, es_Tmin, es, ea
    real(kind=rkind) :: delta, gamma, u2, Ra, Rs, Rso
    real(kind=rkind) :: Rns, Rnl, Rn, G0
    real(kind=rkind) :: L1, L2, L3, L4, nN
    integer(kind=ikind) :: t
  

    es_Tmax = esat(Tmax(element,dt_days))
    es_Tmin = esat(Tmin(element,dt_days))

    es = 0.5_rkind*(es_Tmax + es_Tmin)

    ea = (es_Tmin*RHmax(element,dt_days) + es_Tmax*RHmin(element,dt_days)) / 200._rkind

    delta = slope_vp(Tmean(element,dt_days))
    gamma = psychro_const(z)

    u2 = wind2m(uz(element,dt_days), z)

    Ra = Ra_daily(phi , Julian_day, dt_days)

    ! actual/possible sunshine ratio (simplified)
    nN = 1._rkind

    Rs  = (as + bs*nN) * Ra
    Rso = (0.75_rkind + 2.e-5_rkind*z) * Ra

    Rns = (1._rkind - alpha) * Rs

    Rnl = Rnl_daily( Tmax(element,dt_days)+273.16_rkind, &
                     Tmin(element,dt_days)+273.16_rkind, ea, Rs, Rso )

    Rn = max(Rns - Rnl, 0.0_rkind)

    G0 = G(element)

    L1 = 0.408_rkind * delta * (Rn - G0)
    L2 = gamma * (900._rkind / (Tmean(element,dt_days) + 273._rkind)) * u2 * (es - ea)
    L3 = L1 + L2
    L4 = delta + gamma * (1._rkind + 0.34_rkind*u2)

    penman_monteith = max(L3 / L4, 0.0_rkind)
  end function penman_monteith

  !==============================================================
!  Surface runoff using SCS-CN method
!==============================================================
  pure real(kind=rkind) function surface_runoff(element, dt_days)
    integer(kind=ikind), intent(in) :: element, dt_days
    real(kind=rkind) :: S, I
    S = (25400._rkind / CN) - 254._rkind
    I = 0.2_rkind * S

    if (precip(element,dt_days) <= I) then
       surface_runoff = 0._rkind
    else
       surface_runoff = (precip(element,dt_days) - I)**2 / &
                        (precip(element,dt_days) - I + S)
    end if
  end function surface_runoff

  !==============================================================
!  Ground water function (simple linear reservoir)
!==============================================================
  pure real(kind=rkind) function ground_water(element, dt_days)
    integer(kind=ikind), intent(in) :: element, dt_days
    ground_water = conduct(element) * soilcontent(element,dt_days)
  end function ground_water

  !==============================================================
!  Leakage function (excess water that is not ET or surface runoff)
!==============================================================
  pure real(kind=rkind) function leakage(element, dt_days)
    integer(kind=ikind), intent(in) :: element, dt_days
    leakage = max(0._rkind, qinter(element,dt_days) + precip(element,dt_days) - &
                           ET_flux(element,dt_days) - Qsurf_result(element,dt_days))
  end function leakage

  ! ---- Main driver: compute all steps with routing ----
  subroutine compute_all()
    use tools
    integer(kind=ikind) :: el, t

    do t = 1, ntot_days

       ! 1) Local fluxes per element
       do el = 1, elements%kolik
          elements%hydrobal(el)%ET    = penman_monteith(el,t) * ccrop
          elements%hydrobal(el)%Qsurf = surface_runoff(el,t)
          elements%hydrobal(el)%Li    = leakage(el,t)
          elements%hydrobal(el)%Qgw   = ground_water(el,t)

          ! keep legacy arrays
          ET_flux(el,t)      = elements%hydrobal(el)%ET
          Qsurf_result(el,t) = elements%hydrobal(el)%Qsurf
          L_result(el,t)     = elements%hydrobal(el)%Li
          Qgw_result(el,t)   = elements%hydrobal(el)%Qgw

          elements%hydrobal(el)%inflow  = 0.0_rkind
          elements%hydrobal(el)%outflow = 0.0_rkind
          elements%hydrobal(el)%deltas  = 0.0_rkind
       end do

       ! 2) Routing + mass balance
       call route_step(n_steps)

    end do

  end subroutine compute_all

end module hydrofnc
