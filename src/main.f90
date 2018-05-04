    program main

    use altitude_maintenance_module, only: segment
    use fortran_astrodynamics_toolkit, only: wp,km2m
    use pyplot_module, only : pyplot

    implicit none

    real(wp) :: et0              !! initial ephemeris time (sec)
    real(wp) :: inc0             !! initial inclination - IAU_MOON of date (deg)
    real(wp) :: ran0             !! initial RAAN - IAU_MOON of date (deg)
    real(wp) :: dt_max           !! how long to propagate (days)
    integer  :: n_dvs            !! number of DV maneuvers performed
    real(wp) :: dv_total         !! total DV (km/s)
    real(wp),dimension(6) :: xf  !! final state - inertial frame (km, km/s)
    integer :: i_inc             !! inclination counter
    integer :: i_raan            !! raan counter
    type(pyplot) :: plt   !! pyplot handler
    real(wp),dimension(:),allocatable :: x    !! x array for plot (raan)
    real(wp),dimension(:),allocatable :: y    !! y array for plot (inc)
    real(wp),dimension(:,:),allocatable :: z  !! z array for plot (dv)
    integer :: istat  !! pyplot status code
    character(len=10) :: istr

    !integer,parameter :: inc_start = 90
    !integer,parameter :: inc_stop  = 180
    !integer,parameter :: lan_start = -180
    !integer,parameter :: lan_stop  = 180

    ! integer,parameter :: inc_start = 90
    ! integer,parameter :: inc_stop  = 92
    ! integer,parameter :: lan_start = 0
    ! integer,parameter :: lan_stop  = 45

    !integer,parameter :: inc_start = 80
    !integer,parameter :: inc_stop  = 100
    !integer,parameter :: lan_start = -180
    !integer,parameter :: lan_stop  = 179

    integer,parameter :: inc_start = 80
    integer,parameter :: inc_stop  = 85
    integer,parameter :: lan_start = 100
    integer,parameter :: lan_stop  = 150


    ! ... test cases ...
    !inc0 = 100.0_wp  ! three maneuvers
    !ran0 = 45.0_wp
    !inc0 = 90.0_wp   ! no maneuvers
    !ran0 = 0.0_wp

    type(segment) :: seg  !! the integrator

    real(wp),parameter :: alt0 = 100.0_wp  !! initial altitude for circular orbit (km)
    real(wp),parameter :: deadband_alt = 10.0_wp !! altitude below initial to trigger periapsis raise (km)

    ! initialize the segment:
    call seg%initialize_seg(alt0,deadband_alt)

    call plt%initialize(grid=.true.,xlabel='LAN (deg)',&
                        ylabel='Inc (deg)',figsize=[10,10],&
                        title='Lunar Orbit Maintenance $\Delta v$ (m/s) : deadband = 10 km : dt = 10 days',&
                        real_fmt='(E9.3)')

    ! initialize the indep arrays:
    !x = [(real(i_raan, wp), i_raan = lan_start, lan_stop)]
    !y = [(real(i_inc, wp),  i_inc  = inc_start, inc_stop)]
    allocate(x(lan_start:lan_stop))
    allocate(y(inc_start:inc_stop))
    !allocate(z(size(x), size(y)))
    allocate(z(lan_start:lan_stop, inc_start:inc_stop))
    z = 0.0_wp

    do i_inc = inc_start, inc_stop

        inc0 = real(i_inc, wp)
        y(i_inc) = inc0

        do i_raan = lan_start, lan_stop

            ran0 = real(i_raan, wp)
            x(i_raan) = ran0

            write(*,*) ''
            write(*,*) '============'
            write(*,*) 'inc, ran', inc0, ran0
            write(*,*) '============'
            write(*,*) ''

            et0  = 0.0_wp
            dt_max = 10.0_wp

            ! write(*,*) ''
            ! write(*,*) 'starting...'
            ! write(*,*) ''

            call seg%altitude_maintenance(et0,inc0,ran0,dt_max,n_dvs,dv_total,xf)

            z(i_raan,i_inc) = dv_total * km2m ! in meters
            ! write(*,*) ''
            ! write(*,*) 'finished'
            ! write(*,*) ''

        end do

    end do

    call plt%add_contour(x, y, z, label='contour', linestyle='-', &
                         linewidth=2, filled=.true., cmap='jet', colorbar=.true.,istat=istat)
    call plt%savefig('lom.png',pyfile='lom.py',istat=istat)

    ! print some stats:
    write(*,*) ''
    write(*,*) '-- Stats --'
    do i_inc = inc_start, inc_stop
        write(*,*) ''
        write(*,*) 'max dv for inc ', y(i_inc), ' : ', maxval(z(:,i_inc)), 'm/s'
        write(*,*) 'min dv for inc ', y(i_inc), ' : ', minval(z(:,i_inc)), 'm/s'

        ! also make plots for each inc value:
        write(istr,'(I10)') int(y(i_inc))
        call plt%destroy()
        call plt%initialize(grid=.true.,xlabel='LAN (deg)',&
                            ylabel='$\Delta v$ (m/s)',figsize=[10,10],&
                            title='Lunar Orbit Maintenance : deadband = 10 km : dt = 10 days : Inc='//trim(adjustl(istr)),&
                            real_fmt='(E9.3)')
        call plt%add_plot(x, z(:,i_inc), linestyle='b-', label='dv', istat=istat)
        call plt%savefig('lom_INC='//trim(adjustl(istr))//'.png',&
                          pyfile='lom_INC='//trim(adjustl(istr))//'.py',istat=istat)
    end do

    end program main
