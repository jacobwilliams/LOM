!*****************************************************************************************
!>
!
    module altitude_maintenance_module

    use ddeabm_module, only: ddeabm_class,ddeabm_with_event_class
    use fortran_astrodynamics_toolkit

    implicit none

    private

    character(len=*),parameter :: ephemeris_file = '../eph/JPLEPH_gfortran_mac.421'
    character(len=*),parameter :: gravfile       = '../grav/gggrx_0020pm_sha.tab'

    type,extends(ddeabm_with_event_class) :: segment

        !! the integrator

        integer :: event = 0 !! event function to use

        real(wp) :: nominal_altitude = 100.0_wp  !! nominal altitude (km)
        real(wp) :: deadband         = 10.0_wp   !! altitude below nominal to trigger maneuver (km)
        real(wp) :: r_moon           = 1737.4_wp !! radius of the moon (km)

        type(geopotential_model_pines) :: grav !! central body geopotential model
        type(jpl_ephemeris)            :: eph  !! the ephemeris

        integer  :: n_eoms = 6            !! size of EOM derivative vector [x,y,z,vx,vy,vz]
        real(wp) :: integrator_tol = 1.0e-10_wp !! integrator tols
        integer  :: maxsteps = 10000        !! integrator max steps
        integer  :: grav_n = 8            !! max grav degree
        integer  :: grav_m = 8            !! max grav order
        real(wp) :: root_tol = 1.0e-3_wp  !! event tolerance for deadband (km)

        contains

        procedure,public :: initialize_seg => initialize_segment

    end type segment

    public :: altitude_maintenance

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Initialize the segment for integration.

    subroutine initialize_segment(me,alt0,deadband_alt)

    implicit none

    class(segment),intent(inout) :: me
    real(wp),intent(in) :: alt0
    real(wp),intent(in) :: deadband_alt

    logical :: status_ok

    ! set up the ephemeris:
    write(*,*) 'loading ephemeris file: '//trim(ephemeris_file)
    call me%eph%initialize(filename=ephemeris_file,status_ok=status_ok)
    if (.not. status_ok) error stop 'error initializing ephemeris'

    ! set up the force model [main body is moon]:
    call me%grav%initialize(gravfile,me%grav_n,me%grav_m,status_ok)
    if (.not. status_ok) error stop 'error initializing gravity model'

    ! set up the integrator:
    call me%initialize_event(me%n_eoms,me%maxsteps,ballistic_derivs,&
                                [me%integrator_tol],[me%integrator_tol],&
                                 event_func,me%root_tol)

    ! set class variables for event function:
    me%nominal_altitude = alt0
    me%deadband = deadband_alt

    end subroutine initialize_segment
!*****************************************************************************************

!*****************************************************************************************
!>
!  Altitude maintenance for a circular lunar orbit - periapsis only control.

    subroutine altitude_maintenance(et0,alt0,inc0,ran0,deadband_alt,dt_max,n_dvs,dv_total,xf)

    implicit none

    real(wp),intent(in)               :: et0          !! initial ephemeris time (sec)
    real(wp),intent(in)               :: alt0         !! initial altitude for circular orbit (km)
    real(wp),intent(in)               :: inc0         !! initial inclination - IAU_MOON of date (deg)
    real(wp),intent(in)               :: ran0         !! initial RAAN - IAU_MOON of date (deg)
    real(wp),intent(in)               :: deadband_alt !! altitude below initial to trigger periapsis raise (km)
    real(wp),intent(in)               :: dt_max       !! how long to propagate (days)
    integer,intent(out)               :: n_dvs        !! number of DV maneuvers performed
    real(wp),intent(out)              :: dv_total     !! total DV (km/s)
    real(wp),dimension(6),intent(out) :: xf           !! final state - inertial frame (km, km/s)

    real(wp),dimension(3) :: r,v
    real(wp) :: sma          !! circular orbit semi-major axis [km]
    real(wp) :: etf_max      !! max ephemeris time (sec)
    real(wp) :: min_altitude !! minimum altitude (km)
    real(wp),dimension(3,3) :: rotmat  !! rotation matrix from ICRF to IAU_MOON
    real(wp),dimension(3) :: dv !! periapsis raise maneuver [km/s]
    real(wp),dimension(6) :: x
    real(wp) :: t
    real(wp) :: tf
    integer :: idid
    real(wp) :: gval

    type(segment) :: seg  !! the integrator

    ! initialize the segment:
    call seg%initialize_seg(alt0,deadband_alt)

    n_dvs = 0
    dv_total = zero
    sma = seg%r_moon + alt0
    etf_max = et0 + dt_max*day2sec

    ! get initial state in J2000 - Cartesian for integration
    ! note that inc,ran are in moon-centered-of-date-frame
    call orbital_elements_to_rv(body_moon%mu,sma,zero,inc0,ran0,zero,zero,r,v)

    ! rotate from body-fixed moon of date to j2000:
    rotmat = icrf_to_iau_moon(et0)   ! rotation matrix from inertial to body-fixed Moon frame
    x(1:3) = matmul(transpose(rotmat),r)
    x(4:6) = matmul(transpose(rotmat),v)  ! because using "of date" iau_moon frame

    ! times are ephemeris time
    t  = et0
    tf = etf_max
    seg%event = 1  ! propagate until the altitude is less than min_altitude

    ! main integration loop:
    do

        call seg%integrate_to_event(t,x,tf,idid=idid,gval=gval)

        if (idid<0) then

            write(*,'(A,*(I5/))') 'idid: ',idid
            error stop 'error in integrator'

        elseif (idid==2 .or. idid==3) then
            ! if we reached the max time, then we are done, so exit

            write(*,*) 'done'
            exit

        elseif (idid == 1000) then  ! a root has been found

            select case (seg%event)
            case(1)
                ! if we hit the min altitude, then integrate to next apoapsis
                ! and raise the periapsis:
                seg%event = 2 ! propagate until apoapsis

            case (2)

                ! we have propagated to apoapsis, perform a maneuver to raise periapsis

                ! ... compute DV to raise periapsis back to initial sma

                ! ...

                !call periapsis_apoapsis(body_moon%mu,a,e,rp,ra,vp,va)

               ! dv = ...

                x(4:6) = x(4:6) + dv ! apply the maneuver
                n_dvs = n_dvs + 1
                dv_total = dv_total + norm2(dv)

                seg%event = 1     ! continue with normal mode
                call seg%first_call()  ! have to restart the integration
                                       ! since we changed the state

            case default
                error stop 'invalid event value in altitude_maintenance'
            end select

        else
            write(*,*) 'unknown exit code from integrator: idid=',idid
            error stop 'error in altitude_maintenance'
        end if

    end do

    xf = x

    end subroutine altitude_maintenance
!*****************************************************************************************

!*****************************************************************************************
!>
!  Event function: when the altitude drops below the deadband.
!
!@note Moon is assumed to be a sphere.

    subroutine event_func(me,t,x,g)

    implicit none

    class(ddeabm_with_event_class),intent(inout) :: me
    real(wp),intent(in)              :: t  !! time
    real(wp),dimension(:),intent(in) :: x  !! state -- moon centered inertial frame
    real(wp),intent(out)             :: g  !! event function

    real(wp) :: alt  !! altitude (km)
    real(wp) :: p, ecc, inc, raan, aop, tru

    select type (me)
    class is (segment)

        select case (me%event)
        case (1)
            ! offset from deadband altitude
            alt = norm2(x(1:3)) - me%r_moon
            g = alt - (me%nominal_altitude - me%deadband)
        case (2)
            ! offset from a true anomaly of 180 (apoapsis)
            call rv_to_orbital_elements(body_moon%mu,x(1:3),x(4:6),p,ecc,inc,raan,aop,tru)
            if (tru<zero) tru = tru + twopi  ! from 0 -> 360
            g = tru - pi
        case default
            error stop 'invalid event value in event_func'
        end select
    class default
        error stop 'invalid class in event_func'
    end select

    end subroutine event_func
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
        et = t   ! + et_ref

        ! geopotential gravity:
        rotmat = icrf_to_iau_moon(et)   ! rotation matrix from inertial to body-fixed Moon frame
        rb = matmul(rotmat,r)           ! r in body-fixed frame
        call me%grav%get_acc(rb,me%grav_n,me%grav_m,a_geopot)  ! get the acc due to the geopotential
        a_geopot = matmul(transpose(rotmat),a_geopot)    ! convert acc back to inertial frame

        ! third-body state vectors (wrt the central body, which is the moon in this case):
        ! [inertial frame]
        call me%eph%get_rv(et,body_earth,body_moon,rv_earth_wrt_moon,status_ok)
        call me%eph%get_rv(et,body_sun,body_moon,rv_sun_wrt_moon,status_ok)

        ! third-body perturbation (earth & sun):
        a_third_body = 0.0_wp
        call third_body_gravity(r,rv_earth_wrt_moon(1:3),body_earth%mu,a_earth)
        call third_body_gravity(r,rv_sun_wrt_moon(1:3),body_sun%mu,a_sun)
        a_third_body = a_earth + a_sun

        !total derivative vector:
        xdot(1:3) = v
        xdot(4:6) = a_geopot + a_third_body

    class default
        error stop 'invalid class in ballistic_derivs'
    end select

    end subroutine ballistic_derivs
!*****************************************************************************************


    end module altitude_maintenance_module
