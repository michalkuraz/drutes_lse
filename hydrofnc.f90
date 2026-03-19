! hydrofnc.f90
module hydrofnc
  use typy
  use globals
  implicit none

  real(kind=rkind), parameter :: PI_RKIND = 3.14159265358979323846_rkind

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

  pure real(kind=rkind) function Ra_daily(phi, Julian_day)
    real(kind=rkind), intent(in) :: phi
    integer(kind=ikind), intent(in)          :: Julian_day
    real(kind=rkind) :: dr, delta, ws

    dr    = 1._rkind + 0.033_rkind * cos(2._rkind*PI_RKIND*Julian_day/365._rkind)
    delta = 0.409_rkind * sin(2._rkind*PI_RKIND*Julian_day/365._rkind - 1.39_rkind)
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

  pure real(kind=rkind) function penman_monteith(element, day)
    integer(kind=ikind), intent(in) :: element, day
    real(kind=rkind) :: es_Tmax, es_Tmin, es, ea
    real(kind=rkind) :: delta, gamma, u2, Ra, Rs, Rso
    real(kind=rkind) :: Rns, Rnl, Rn, G0
    real(kind=rkind) :: L1, L2, L3, L4, nN

    es_Tmax = esat(Tmax(element,day))
    es_Tmin = esat(Tmin(element,day))

    es = 0.5_rkind*(es_Tmax + es_Tmin)

    ea = (es_Tmin*RHmax(element,day) + es_Tmax*RHmin(element,day)) / 200._rkind

    delta = slope_vp(Tmean(element,day))
    gamma = psychro_const(z)

    u2 = wind2m(uz(element,day), z)

    Ra = Ra_daily(phi, Julian_day)

    ! actual/possible sunshine ratio (simplified)
    nN = 1._rkind

    Rs  = (as + bs*nN) * Ra
    Rso = (0.75_rkind + 2.e-5_rkind*z) * Ra

    Rns = (1._rkind - alpha) * Rs

    Rnl = Rnl_daily( Tmax(element,day)+273.16_rkind, &
                     Tmin(element,day)+273.16_rkind, ea, Rs, Rso )

    Rn = max(Rns - Rnl, 0.0_rkind)

    G0 = G(element,day)

    L1 = 0.408_rkind * delta * (Rn - G0)
    L2 = gamma * (900._rkind / (Tmean(element,day) + 273._rkind)) * u2 * (es - ea)
    L3 = L1 + L2
    L4 = delta + gamma * (1._rkind + 0.34_rkind*u2)

    penman_monteith = max(L3 / L4, 0.0_rkind)
  end function penman_monteith

  pure real(kind=rkind) function surface_runoff(element, day)
    integer(kind=ikind), intent(in) :: element, day
    real(kind=rkind) :: S, I
    S = (25400._rkind / CN) - 254._rkind
    I = 0.2_rkind * S

    if (precip(element,day) <= I) then
       surface_runoff = 0._rkind
    else
       surface_runoff = (precip(element,day) - I)**2 / &
                        (precip(element,day) - I + S)
    end if
  end function surface_runoff

  pure real(kind=rkind) function ground_water(element, day)
    integer(kind=ikind), intent(in) :: element, day
    ground_water = conduct(element,day) * soilcontent(element,day)
  end function ground_water

  pure real(kind=rkind) function leakage(element, day)
    integer(kind=ikind), intent(in) :: element, day
    leakage = max(0._rkind, qinter(element,day) + precip(element,day) - &
                           ET_flux(element,day) - Qsurf_result(element,day))
  end function leakage

  ! ---- Main driver: compute all steps with routing ----
  subroutine compute_all()
    use tools, only: route_step
    integer(kind=ikind) :: el

    do time_step = 1, n_days

       ! 1) Local fluxes per element
       do el = 1, elements%kolik
          elements%hydrobal(el)%ET    = penman_monteith(el,time_step) * ccrop
          elements%hydrobal(el)%Qsurf = surface_runoff(el,time_step)
          elements%hydrobal(el)%Li    = leakage(el,time_step)
          elements%hydrobal(el)%Qgw   = ground_water(el,time_step)

          ! keep legacy arrays
          ET_flux(el,time_step)      = elements%hydrobal(el)%ET
          Qsurf_result(el,time_step) = elements%hydrobal(el)%Qsurf
          L_result(el,time_step)     = elements%hydrobal(el)%Li
          Qgw_result(el,time_step)   = elements%hydrobal(el)%Qgw

          elements%hydrobal(el)%inflow  = 0.0_rkind
          elements%hydrobal(el)%outflow = 0.0_rkind
          elements%hydrobal(el)%deltas  = 0.0_rkind
       end do

       ! 2) Routing + mass balance
       call route_step(time_step)

    end do

  end subroutine compute_all

end module hydrofnc
