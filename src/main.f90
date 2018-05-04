    program main

    use altitude_maintenance_module
    use fortran_astrodynamics_toolkit, only: wp

    implicit none

    real(wp) :: et0          !! initial ephemeris time (sec)
    real(wp) :: alt0         !! initial altitude for circular orbit (km)
    real(wp) :: inc0         !! initial inclination - IAU_MOON of date (deg)
    real(wp) :: ran0         !! initial RAAN - IAU_MOON of date (deg)
    real(wp) :: deadband_alt !! altitude below initial to trigger periapsis raise (km)
    real(wp) :: dt_max       !! how long to propagate (days)
    integer  :: n_dvs        !! number of DV maneuvers performed
    real(wp) :: dv_total     !! total DV (km/s)
    real(wp),dimension(6) :: xf           !! final state - inertial frame (km, km/s)

    et0  = 0.0_wp
    alt0 = 100.0_wp
    inc0 = 100.0_wp  ! three maneuvers
    ran0 = 45.0_wp
    !inc0 = 90.0_wp   ! no maneuvers
    !ran0 = 0.0_wp
    deadband_alt = 10.0_wp
    dt_max = 10.0_wp

    write(*,*) ''
    write(*,*) 'starting...'
    write(*,*) ''

    call altitude_maintenance(et0,alt0,inc0,ran0,deadband_alt,dt_max,n_dvs,dv_total,xf)

    write(*,*) ''
    write(*,*) 'finished'
    write(*,*) ''

    end program main
