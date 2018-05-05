!*****************************************************************************************
!>
!  Altitude maintenance for low lunar orbits.

    module altitude_maintenance_module

    use ddeabm_module, only: ddeabm_class,ddeabm_with_event_class
    use fortran_astrodynamics_toolkit

    implicit none

    private

    character(len=*),parameter :: ephemeris_file = '../data/eph/JPLEPH_gfortran_mac.421'
    character(len=*),parameter :: gravfile       = '../data/grav/gggrx_0020pm_sha.tab'

    type,extends(ddeabm_with_event_class),public :: segment

        !! the main class for integrating a LLO

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
        integer  :: grav_n = 8              !! max grav degree
        integer  :: grav_m = 8              !! max grav order
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

    subroutine initialize_segment(me,alt0,deadband_alt)

    implicit none

    class(segment),intent(inout) :: me
    real(wp),intent(in) :: alt0
    real(wp),intent(in) :: deadband_alt

    logical :: status_ok

    ! set up the integrator:
    call me%initialize_event(me%n_eoms,me%maxsteps,ballistic_derivs,&
                                [me%integrator_tol],[me%integrator_tol],&
                                 event_func,me%root_tol)

    if (me%include_third_bodies) then
        ! set up the ephemeris:
        write(*,*) 'loading ephemeris file: '//trim(ephemeris_file)
        call me%eph%initialize(filename=ephemeris_file,status_ok=status_ok)
        if (.not. status_ok) error stop 'error initializing ephemeris'
    end if

    ! set up the force model [main body is moon]:
    write(*,*) 'loading gravity file: '//trim(gravfile)
    call me%grav%initialize(gravfile,me%grav_n,me%grav_m,status_ok)
    if (.not. status_ok) error stop 'error initializing gravity model'

    ! set class variables for event function:
    me%nominal_altitude = alt0
    me%deadband = deadband_alt

    end subroutine initialize_segment
!*****************************************************************************************

!*****************************************************************************************
!>
!  Altitude maintenance for a circular lunar orbit - periapsis only control.

    subroutine altitude_maintenance(seg,et0,inc0,ran0,dt_max,n_dvs,dv_total,xf)

    implicit none

    class(segment),intent(inout)      :: seg
    real(wp),intent(in)               :: et0          !! initial ephemeris time (sec)
    real(wp),intent(in)               :: inc0         !! initial inclination - IAU_MOON of date (deg)
    real(wp),intent(in)               :: ran0         !! initial RAAN - IAU_MOON of date (deg)
    real(wp),intent(in)               :: dt_max       !! how long to propagate (days)
    integer,intent(out)               :: n_dvs        !! number of DV maneuvers performed
    real(wp),intent(out)              :: dv_total     !! total DV (km/s)
    real(wp),dimension(6),intent(out) :: xf           !! final state - inertial frame (km, km/s)

    real(wp),dimension(3) :: r,v
    real(wp) :: sma          !! circular orbit semi-major axis [km]
    real(wp) :: min_altitude !! minimum altitude (km)
    real(wp),dimension(3,3) :: rotmat  !! rotation matrix from ICRF to IAU_MOON
    real(wp) :: dv !! periapsis raise maneuver magnitude [km/s]
    real(wp),dimension(6) :: x  !! J2000-Moon state vector
    real(wp) :: t  !! integration time (sec from et0)
    real(wp) :: tf !! final integration time (sec from et0)
    integer :: idid !! integrator status flag
    real(wp) :: gval  !! event function value
    real(wp) :: a, p, ecc, inc, raan, aop, tru, tru2
    real(wp) :: rp1,ra1,vp1,va1,rp2,va2

    n_dvs = 0
    dv_total = zero
    sma = seg%r_moon + seg%nominal_altitude  ! initial orbit sma (circular)

    ! get initial state in J2000 - Cartesian for integration
    ! note that inc,ran are in moon-centered-of-date-frame
    call orbital_elements_to_rv(body_moon%mu,sma,zero,inc0*deg2rad,ran0*deg2rad,zero,zero,r,v)

    ! rotate from body-fixed moon of date to j2000:
    rotmat = icrf_to_iau_moon(et0)   ! rotation matrix from inertial to body-fixed Moon frame
    x(1:3) = matmul(transpose(rotmat),r)
    x(4:6) = matmul(transpose(rotmat),v)  ! because using "of date" iau_moon frame

    ! times are ephemeris time
    seg%et_ref = et0
    t  = zero
    tf = dt_max*day2sec
    seg%event = 1  ! propagate until the altitude is less than min_altitude

    !write(*,*) 'starting integration loop...'
    call seg%first_call()
    ! main integration loop:
    do

        call seg%integrate_to_event(t,x,tf,idid=idid,gval=gval)

        if (idid<0) then

            write(*,'(A,*(I5/))') 'idid: ',idid
            error stop 'error in integrator'

        elseif (idid==2 .or. idid==3) then
            ! if we reached the max time, then we are done, so exit

            !write(*,*) 'done'
            exit

        elseif (idid == 1000) then  ! a root has been found

            !write(*,*) 'event found'

            select case (seg%event)
            case(1)
                write(*,*) 'min altitude at ', t*sec2hr, 'hr'
                ! if we hit the min altitude, then integrate to next apoapsis
                ! and raise the periapsis:
                seg%event = 2 ! propagate until apoapsis
                call seg%first_call()  ! have to restart the integration
                                       ! since we just root solved

            case (2)

                ! we have propagated to apoapsis, perform a maneuver to raise periapsis

                ! compute current orbit elements:
                call rv_to_orbital_elements(body_moon%mu,x(1:3),x(4:6),p,ecc,inc,raan,aop,tru)
                a = p / (one - ecc*ecc)
                tru = tru*rad2deg ! convert to deg
                if (tru<0.0_wp) tru = tru + 360.0_wp
                call periapsis_apoapsis(body_moon%mu,a,ecc,rp1,ra1,vp1,va1)
                rp2 = sma ! desired periapsis radius for new orbit

                if (rp1>(rp2-seg%deadband)) then
                    write(*,*) 'dv not necessary'
                    ! in this case, the osculating periapsis radius has not
                    ! violated the deadband altitude after all, so just
                    ! continue without doing a maneuver.

                ! elseif (tru<179.0_wp .or. tru>181.0_wp ) then  ! HACK - IGNORE PERIAPSIS ROOTS

                !     write(*,*)  'oops stopped at periapsis : TRU=',tru, ' : t=', t*sec2hr, 'gval=',gval

                !     ! we have to keep integrating, we stopped at periapsis ...
                !     seg%event = 2  ! keep this...
                !     call seg%first_call()  ! have to restart integration...
                !     cycle

                !     !... something wrong here... it's not working... getting stuck...

                ! below: if we stopped at periapsis, just reset and continue with event mode 1
                ! [not efficient, since we have to restart integration]
                !
                elseif (tru>=179.0_wp .and. tru<=181.0_wp ) then  ! HACK - IGNORE PERIAPSIS ROOTS
                    !write(*,*) 'dv to raise rp from ', rp1, ' at ', t*sec2day, 'days'

                    ! apoapsis velocity to raise periapsis radius to rp2
                    va2 = sqrt( two * body_moon%mu * ( one/ra1 - one/(rp2+ra1) ) )

                    ! delta-v to raise periapsis back to initial sma
                    dv = va2 - va1

                    ! apply the maneuver along the current apoapsis velocity vector:
                    x(4:6) = x(4:6) + dv * unit(x(4:6))

                    ! ... check results:
                    call rv_to_orbital_elements(body_moon%mu,x(1:3),x(4:6),p,ecc,inc,raan,aop,tru2)
                    a = p / (one - ecc*ecc)
                    call periapsis_apoapsis(body_moon%mu,a,ecc,rp1,ra1,vp1,va1)
                    ! write(*,*) ''
                    ! write(*,*) '-----'
                    ! write(*,*) 'after dv'
                    ! write(*,*) 'rp = ',rp1
                    ! write(*,*) 'ra = ',ra1
                    ! write(*,*) '-----'
                    ! write(*,*) ''

                    ! keep track of totals:
                    n_dvs = n_dvs + 1
                    dv_total = dv_total + dv

                    write(*,*) 'apoapsis at ', t*sec2hr, 'hr', ' : TRU = ', tru, ' : dv = ', dv

                    !write(*,*) 'DV', n_dvs, dv

                end if

                seg%event = 1     ! continue with normal mode
                call seg%first_call()  ! have to restart the integration
                                       ! since we root solved and/or changed
                                       ! the state

            case default
                error stop 'invalid event value in altitude_maintenance'
            end select

        else
            write(*,*) 'unknown exit code from integrator: idid=',idid
            error stop 'error in altitude_maintenance'
        end if

    end do

    ! final state:
    xf = x

    write(*,*) ''
    write(*,*) '=============='
    write(*,*) 't        = ', t*sec2day, 'days'
    write(*,*) 'n_dvs    = ', n_dvs
    write(*,*) 'dv_total = ', dv_total
    write(*,*) '=============='
    write(*,*) ''

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

            !... this is also catching periapsis ....  how to get only apoapsis ?????
            ! ... have to update ddeabm to allow for user-specified bracket function
            !     so we can test to see which one it's near ...

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
