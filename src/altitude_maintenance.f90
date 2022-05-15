!*****************************************************************************************
!>
!  Altitude maintenance for low lunar orbits.
!
!  Assumptions:
!
!  * Circular low-lunar orbit.
!  * Only periapsis altitude is controlled (all other elements float)
!  * The low-fidelity IAU_MOON frame is used for the elements and gravity model.

    module altitude_maintenance_module

    use fortran_astrodynamics_toolkit
    use ddeabm_module,   only: ddeabm_class,ddeabm_with_event_class
    use iso_fortran_env, only: error_unit,output_unit, wp => real64

    implicit none

    private

    type,extends(ddeabm_with_event_class),public :: segment

        !! the main class for integrating a low-lunar orbit.

        private

        integer :: event = 0 !! event function to use:
                             !!
                             !! 1 = integrate until minimum altitude
                             !! 2 = integrate to apoapsis

        real(wp) :: et_ref !! reference ephemeris time (sec)

        real(wp) :: nominal_altitude = 100.0_wp  !! nominal altitude (km)
        real(wp) :: deadband         = 10.0_wp   !! altitude below nominal to trigger maneuver (km)
        real(wp) :: r_moon           = 1737.4_wp !! radius of the moon (km)

        logical :: include_third_bodies = .false. !! to also include Earth and Sun in force model

        type(geopotential_model_pines) :: grav  !! central body geopotential model
        type(jpl_ephemeris)            :: eph   !! the ephemeris

        integer  :: n_eoms = 6                  !! size of EOM derivative vector [x,y,z,vx,vy,vz]
        real(wp) :: integrator_tol = 1.0e-12_wp !! integrator tols
        integer  :: maxsteps = 1000000          !! integrator max steps
        integer  :: grav_n = 8                  !! max grav degree
        integer  :: grav_m = 8                  !! max grav order
        real(wp) :: root_tol = 1.0e-6_wp        !! event tolerance for deadband (km)

        contains

        procedure,public :: initialize_seg => initialize_segment
        procedure,public :: altitude_maintenance

    end type segment

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Initialize the segment for integration.

    subroutine initialize_segment(me,alt0,deadband_alt,grav_n,grav_m,&
                                  ephemeris_file,gravfile)

    implicit none

    class(segment),intent(inout) :: me
    real(wp),intent(in) :: alt0
    real(wp),intent(in) :: deadband_alt
    integer,intent(in) :: grav_n
    integer,intent(in) :: grav_m
    character(len=*),intent(in) :: ephemeris_file
    character(len=*),intent(in) :: gravfile

    logical :: status_ok

    ! set up the integrator:
    call me%initialize_event(me%n_eoms,me%maxsteps,ballistic_derivs,&
                                [me%integrator_tol],[me%integrator_tol],&
                                 event_func,me%root_tol)

    if (me%include_third_bodies) then
        ! set up the ephemeris:
        write(output_unit,'(A)') 'loading ephemeris file: '//trim(ephemeris_file)
        call me%eph%initialize(filename=ephemeris_file,status_ok=status_ok)
        if (.not. status_ok) error stop 'error initializing ephemeris'
    end if

    ! set class variables for event function:
    me%nominal_altitude = alt0
    me%deadband = deadband_alt
    me%grav_n = grav_n
    me%grav_m = grav_m

    ! set up the force model [main body is moon]:
    write(output_unit,'(A)') 'loading gravity file: '//trim(gravfile)
    call me%grav%initialize(gravfile,me%grav_n,me%grav_m,status_ok)
    if (.not. status_ok) error stop 'error initializing gravity model'

    end subroutine initialize_segment
!*****************************************************************************************

!*****************************************************************************************
!>
!  Altitude maintenance for a circular lunar orbit - periapsis only control.

    subroutine altitude_maintenance(seg,et0,inc0,ran0,tru0,dt_max,n_dvs,dv_total,xf)

    implicit none

    class(segment),intent(inout)      :: seg
    real(wp),intent(in)               :: et0      !! initial ephemeris time (sec)
    real(wp),intent(in)               :: inc0     !! initial inclination - IAU_MOON of date (deg)
    real(wp),intent(in)               :: ran0     !! initial RAAN - IAU_MOON of date (deg)
    real(wp),intent(in)               :: tru0     !! initial true anomaly - IAU_MOON of date (deg)
    real(wp),intent(in)               :: dt_max   !! how long to propagate (days)
    integer,intent(out)               :: n_dvs    !! number of DV maneuvers performed
    real(wp),intent(out)              :: dv_total !! total DV (km/s)
    real(wp),dimension(6),intent(out) :: xf       !! final state - inertial frame (km, km/s)

    real(wp),dimension(3) :: r !! position vector (km)
    real(wp),dimension(3) :: v !! velocity vector (km/s)
    real(wp),dimension(3,3) :: rotmat  !! rotation matrix from ICRF to IAU_MOON
    real(wp),dimension(6) :: x  !! J2000-Moon state vector
    integer :: idid             !! integrator status flag
    real(wp) :: sma             !! circular orbit semi-major axis (km)
    real(wp) :: min_altitude    !! minimum altitude (km)
    real(wp) :: dv              !! periapsis raise maneuver magnitude (km/s)
    real(wp) :: t               !! integration time (sec from et0)
    real(wp) :: tf              !! final integration time (sec from et0)
    real(wp) :: gval            !! event function value
    real(wp) :: tru             !! true anomaly (deg)

    ! initialize
    n_dvs = 0
    dv_total = zero

    ! get initial state in J2000 - Cartesian for integration
    ! note that inc,ran are in moon-centered-of-date-frame
    ! [circular orbit so p=sma]
    sma = seg%r_moon + seg%nominal_altitude  ! initial orbit sma (circular)
    call orbital_elements_to_rv(body_moon%mu,sma,zero,&
                                inc0*deg2rad,ran0*deg2rad,zero,tru0*deg2rad,r,v)

    ! rotate from body-fixed moon of date to j2000:
    rotmat = icrf_to_iau_moon(et0)   ! rotation matrix from inertial to body-fixed Moon frame
    x(1:3) = matmul(transpose(rotmat),r)
    x(4:6) = matmul(transpose(rotmat),v)  ! because using "of date" iau_moon frame

    ! times are relative to initial epoch (sec)
    seg%et_ref = et0
    t  = zero
    tf = dt_max*day2sec

    ! propagate until the altitude is less than min_altitude
    seg%event = 1

    ! main integration loop:
    call seg%first_call()
    do

        call seg%integrate_to_event(t,x,tf,idid=idid,gval=gval)

        if (idid<0) then

            write(error_unit,'(A,*(I5/))') 'idid: ',idid
            error stop 'error in integrator'

        elseif (idid==2 .or. idid==3) then

            ! if we reached the max time, then we are done, so exit
            exit

        elseif (idid == 1000) then  ! a root has been found

            select case (seg%event)
            case(1)

                ! if we hit the min altitude, then integrate to next apoapsis
                ! and raise the periapsis:
                seg%event = 2 ! propagate until apoapsis
                call seg%first_call()  ! have to restart the integration
                                       ! since we just root solved

                ! compute the maneuver here that will be performed
                ! at next apoapsis to raise the periapsis, no matter what.
                dv = periapsis_raise_maneuver(x,sma)

            case (2)

                ! we have propagated to periapsis or apoapsis, if at apoapsis,
                ! perform a maneuver to raise periapsis
                ! if at periapsis, just continue on.

                ! compute current true anomaly:
                tru = true_anomaly(x)

                if (tru<179.0_wp .or. tru>181.0_wp ) then

                    ! we have to keep integrating, since we stopped at periapsis
                    ! note: this is inefficient since we have to
                    ! restart the integration. it would be better to prevent
                    ! to root solver from stopping here.
                    seg%event = 2  ! continue with mode 2
                    call seg%first_call()  ! have to restart integration
                    cycle

                else ! we are at apoapsis

                    ! perform the maneuver we computed when
                    ! the min altitude was triggered
                    ! apply the maneuver along the current apoapsis velocity vector:
                    x(4:6) = x(4:6) + dv * unit(x(4:6))

                    ! keep track of totals:
                    n_dvs = n_dvs + 1
                    dv_total = dv_total + dv

                    write(output_unit,'(*(A,F12.6))') 'maneuver at ', t*sec2hr, &
                                                      ' hr : TRU = ', tru, &
                                                      ' : dv = ', dv

                end if

                seg%event = 1     ! continue with normal mode
                call seg%first_call()  ! have to restart the integration
                                       ! since we root solved and/or changed
                                       ! the state

            case default
                error stop 'invalid event value in altitude_maintenance'
            end select

        else
            write(error_unit,'(A,I5)') 'unknown exit code from integrator: idid=',idid
            error stop 'error in altitude_maintenance'
        end if

    end do

    ! final state:
    xf = x

    end subroutine altitude_maintenance
!*****************************************************************************************

!*****************************************************************************************
!>
!  Compute true anomaly [0, 360] deg.

    pure function true_anomaly(rv) result(tru)

    implicit none

    real(wp),dimension(6),intent(in) :: rv !! [position,velocity] vector
    real(wp) :: tru !! true anomaly (deg)

    real(wp),dimension(3) :: r !! position vector
    real(wp),dimension(3) :: v !! velocity vector
    real(wp) :: p,ecc,inc,raan,aop !! orbital elements

    r = rv(1:3)
    v = rv(4:6)

    call rv_to_orbital_elements(body_moon%mu,r,v,p,ecc,inc,raan,aop,tru)

    tru = tru*rad2deg ! convert to deg
    if (tru<zero) tru = tru + 360.0_wp ! wrap from 0 to 369

    end function true_anomaly
!*****************************************************************************************

!*****************************************************************************************
!>
!  Compute the maneuver at apoapsis to raise periapsis to the specified value.

    pure function periapsis_raise_maneuver(rv,target_rp) result(dv)

    implicit none

    real(wp),dimension(6),intent(in) :: rv !! [position,velocity] vector
    real(wp),intent(in) :: target_rp  !! the rp value to target
    real(wp) :: dv !! the maneuver to perform at apoapsis to target `target_rp`

    real(wp),dimension(3) :: r !! position vector
    real(wp),dimension(3) :: v !! velocity vector
    real(wp) :: a,p,ecc,inc,raan,aop,tru !! orbital elements
    real(wp) :: rp1,ra1,vp1,va1,va2  !! periapsis/apoapsis pos/vel magnitudes

    r = rv(1:3)
    v = rv(4:6)

    call rv_to_orbital_elements(body_moon%mu,r,v,p,ecc,inc,raan,aop,tru)
    a = p / (one - ecc*ecc) ! compute semi-major axis
    call periapsis_apoapsis(body_moon%mu,a,ecc,rp1,ra1,vp1,va1)

    ! apoapsis velocity for periapsis radius of target_rp
    va2 = sqrt( two * body_moon%mu * target_rp/(ra1*(target_rp+ra1)) )

    ! delta-v to raise periapsis back to target rp
    dv = va2 - va1

    end function periapsis_raise_maneuver
!*****************************************************************************************

!*****************************************************************************************
!>
!  Compute radial velocity magnitude \( \dot{r} \)

    pure function rdot(rv) result(rd)

    implicit none

    real(wp),dimension(6),intent(in) :: rv !! [position,velocity] vector
    real(wp) :: rd !! \( \dot{r} \)

    real(wp),dimension(3) :: r !! position vector
    real(wp),dimension(3) :: v !! velocity vector
    real(wp) :: rmag !! position vector magnitude

    r = rv(1:3)
    v = rv(4:6)
    rmag = norm2(r)

    rd = dot_product(r,v) / rmag

    end function rdot
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

    select type (me)
    class is (segment)

        select case (me%event)
        case (1)

            ! g is offset from deadband altitude
            alt = norm2(x(1:3)) - me%r_moon
            g = alt - (me%nominal_altitude - me%deadband)

        case (2)

            ! g is rdot, which is zero and periapsis
            ! and apoapsis of an ellipse:
            g = rdot(x(1:6))

            ! TODO: figure out how to get it to ignore periapsis roots
            ! [maybe need to update DDEABM so we have more control
            ! over which roots it stops at]

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
    logical :: status_ok  !! ephemeris status flag

    select type (me)

    class is (segment)

        ! get state:
        r = x(1:3)
        v = x(4:6)

        ! compute ephemeris time [sec]:
        et = t + me%et_ref

        ! geopotential gravity:
        rotmat = icrf_to_iau_moon(et)   ! rotation matrix from inertial to body-fixed Moon frame
        rb = matmul(rotmat,r)           ! r in body-fixed frame
        call me%grav%get_acc(rb,me%grav_n,me%grav_m,a_geopot)  ! get the acc due to the geopotential
        a_geopot = matmul(transpose(rotmat),a_geopot)    ! convert acc back to inertial frame

        if (me%include_third_bodies) then
            ! third-body state vectors (wrt the central body, which is the moon in this case):
            ! [inertial frame]
            call me%eph%get_rv(et,body_earth,body_moon,rv_earth_wrt_moon,status_ok)
            call me%eph%get_rv(et,body_sun,body_moon,rv_sun_wrt_moon,status_ok)

            ! third-body perturbation (earth & sun):
            a_third_body = 0.0_wp
            call third_body_gravity(r,rv_earth_wrt_moon(1:3),body_earth%mu,a_earth)
            call third_body_gravity(r,rv_sun_wrt_moon(1:3),body_sun%mu,a_sun)
            a_third_body = a_earth + a_sun
        else
            a_third_body = zero
        end if

        !total derivative vector:
        xdot(1:3) = v
        xdot(4:6) = a_geopot + a_third_body

    class default
        error stop 'invalid class in ballistic_derivs'
    end select

    end subroutine ballistic_derivs
!*****************************************************************************************

!*****************************************************************************************
    end module altitude_maintenance_module
!*****************************************************************************************
