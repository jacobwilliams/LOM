!*****************************************************************************************
!>
!
    module altitude_maintenance_module

    use ddeabm_module
    use fortran_astrodynamics_toolkit

    implicit none

    private

    public :: altitude_maintenance

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Altitude maintenance for a circular lunar orbit - periapsis only control.

    subroutine altitude_maintenance(et0,sma0,inc0,ran0,deadband_alt,dt_max,n_dvs,dv_total,xf)

    implicit none

    real(wp),intent(in)               :: et0          !! initial ephemeris time (sec)
    real(wp),intent(in)               :: sma0         !! initial SMA for circular orbit (km)
    real(wp),intent(in)               :: inc0         !! initial inclination (deg)
    real(wp),intent(in)               :: ran0         !! initial RAAN (deg)
    real(wp),intent(in)               :: deadband_alt !! altitude below initial to trigger periapsis raise (km)
    real(wp),intent(in)               :: dt_max       !! how long to propagate (days)
    integer,intent(out)               :: n_dvs        !! number of DV maneuvers performed
    real(wp),intent(out)              :: dv_total     !! total DV (km/s)
    real(wp),dimension(6),intent(out) :: xf           !! final state - inertial frame (km, km/s)


    real(wp) :: etf_maxx  !! max ephemeris time (sec)
    real(wp) :: min_altitude !! minimum altitude (km)

    etf_max = et0 + dt_max*day2sec
    min_altitude = sma0-deadband_alt

    ! propagate until the altitude is less than min_altitude


    end subroutine altitude_maintenance
!*****************************************************************************************

!*****************************************************************************************
!>
!  Equations of motion for a ballistic orbit around the moon.

    subroutine ballistic_derivs(me,t,x,xdot)

    implicit none

    class(ddeabm_class),intent(inout) :: me
    real(wp),intent(in)               :: t    !! time [sec from epoch]
    real(wp),dimension(:),intent(in)  :: x    !! state [r,v] in inertial frame (moon-centered)
    real(wp),dimension(:),intent(out) :: xdot !! derivative of state (\( dx/dt \))

    real(wp),dimension(3) :: r,rb,v
    reaL(wp),dimension(6) :: rv_earth_wrt_moon,rv_sun_wrt_moon
    real(wp),dimension(3,3) :: rotmat
    real(wp),dimension(3) :: a_geopot
    real(wp),dimension(3) :: a_earth
    real(wp),dimension(3) :: a_sun
    real(wp),dimension(3) :: a_third_body
    real(wp) :: et !! ephemeris time of `t`
    logical :: status_ok

    select type (me)

    class is (segment)

        ! get state:
        r = x(1:3)
        v = x(4:6)

        ! compute ephemeris time [sec]:
        et = et_ref + t

        ! geopotential gravity:
        rotmat = icrf_to_iau_moon(et)   ! rotation matrix from inertial to body-fixed Moon frame
        rb = matmul(rotmat,r)           ! r in body-fixed frame
        call me%grav%get_acc(rb,grav_n,grav_m,a_geopot)  ! get the acc due to the geopotential
        a_geopot = matmul(transpose(rotmat),a_geopot)    ! convert acc back to inertial frame

        ! third-body state vectors (wrt the central body, which is the moon in this case):
        ! [inertial frame]
        call me%eph%get_rv(et,body_earth,body_moon,rv_earth_wrt_moon,status_ok)
        call me%eph%get_rv(et,body_sun,body_moon,rv_sun_wrt_moon,status_ok)

        ! third-body perturbation (earth & sun):
        a_third_body = 0.0_wp
        call third_body_gravity(r,rv_earth_wrt_moon(1:3),mu_earth,a_earth)
        call third_body_gravity(r,rv_sun_wrt_moon(1:3),mu_sun,a_sun)
        a_third_body = a_earth + a_sun

        !total derivative vector:
        xdot(1:3) = v
        xdot(4:6) = a_geopot + a_third_body

    class default

        error stop 'invalid class in ballistic_derivs'

    end select

    end subroutine ballistic_derivs
!*****************************************************************************************


    end module lom_module
