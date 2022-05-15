!*****************************************************************************************
!>
!  Altitude maintenance test program

    program main

    use altitude_maintenance_module,    only: segment
    use fortran_astrodynamics_toolkit,  only: km2m
    use pyplot_module,                  only: pyplot
    use iso_fortran_env,                only: error_unit,output_unit, wp => real64

    implicit none

    ! run variables, populated by the config file:
    real(wp) :: dt_max        !! how long to propagate (days)
    real(wp) :: et0           !! initial ephemeris time (sec)
                              !! (only matters if including Earth/Sun perturbations)
    real(wp) :: tru0          !! initial true anomaly (deg)
    real(wp) :: alt0          !! initial altitude for circular orbit (km)
    real(wp) :: deadband_alt  !! altitude below initial to trigger periapsis raise (km)
    real(wp) :: inc_start     !! inc/lan grid
    real(wp) :: inc_stop
    real(wp) :: inc_step
    real(wp) :: lan_start
    real(wp) :: lan_stop
    real(wp) :: lan_step
    integer  :: grav_n     !! max grav degree
    integer  :: grav_m     !! max grav order

    real(wp),dimension(:),allocatable :: x    !! x array for plot (raan)
    real(wp),dimension(:),allocatable :: y    !! y array for plot (inc)
    real(wp),dimension(:,:),allocatable :: z  !! z array for plot (dv)
    real(wp) :: inc0             !! initial inclination - IAU_MOON of date (deg)
    real(wp) :: ran0             !! initial RAAN - IAU_MOON of date (deg)
    integer  :: n_dvs            !! number of DV maneuvers performed
    real(wp) :: dv_total         !! total DV (km/s)
    real(wp),dimension(6) :: xf  !! final state - inertial frame (km, km/s)
    integer :: i_inc             !! inclination counter
    integer :: i_raan            !! raan counter
    type(pyplot) :: plt          !! pyplot handler
    integer :: istat             !! pyplot status code
    character(len=10) :: istr    !! for integer to string conversion
    character(len=10) :: dt_max_str
    character(len=10) :: deadband_alt_str
    type(segment) :: seg  !! the integrator for a ballistic moon-centered trajectory

    ! populate the run variables:
    call read_config_file()

    ! initialize the segment:
    call seg%initialize_seg(alt0,deadband_alt,grav_n,grav_m)

    write(dt_max_str,'(I10)') int(dt_max); dt_max_str = adjustl(dt_max_str)
    write(deadband_alt_str,'(I10)') int(deadband_alt); deadband_alt_str = adjustl(deadband_alt_str)

    call plt%initialize(grid=.true.,xlabel='LAN (deg)',&
                        ylabel='Inc (deg)',figsize=[20,10],&
                        title='Lunar Orbit Maintenance $\Delta v$ (m/s) : deadband = '//&
                                trim(deadband_alt_str)//' km : dt = '//trim(dt_max_str)//' days',&
                        real_fmt='(E9.3)')

    ! initialize the indep arrays:
    call size_arrays(lan_start,lan_stop,lan_step,inc_start,inc_stop,inc_step,x,y,z)

    do i_inc = 1, size(y)

        inc0 = y(i_inc) ! inclination

        do i_raan = 1, size(x)

            ran0 = x(i_raan) ! right ascension of ascending node

            write(output_unit,'(A)') ''
            write(output_unit,'(A)') '=============='
            write(output_unit,'(A,1X,F12.6,1X,F12.6)') 'inc, ran:', inc0, ran0
            write(output_unit,'(A)') ''

            call seg%altitude_maintenance(et0,inc0,ran0,tru0,dt_max,n_dvs,dv_total,xf)

            z(i_raan,i_inc) = dv_total * km2m ! delta-v in m/s

            write(output_unit,'(A)') ''
            write(output_unit,'(A,I5)') 'n_dvs    = ', n_dvs
            write(output_unit,'(A,F10.6)') 'dv_total = ', dv_total
            write(output_unit,'(A)') '=============='
            write(output_unit,'(A)') ''

        end do

    end do

    ! create a contour plot for the delta-v:
    call plt%add_contour(x, y, z, linestyle='-', &
                         linewidth=2, filled=.true., cmap='Blues',&
                         colorbar=.true.,istat=istat)
    call plt%savefig('lom.png',pyfile='lom.py',istat=istat)

    ! print some stats:
    write(*,*) ''
    write(*,*) '-- Stats --'
    do i_inc = 1, size(y)
        write(*,*) ''
        write(*,*) 'max dv for inc ', y(i_inc), ' : ', maxval(z(:,i_inc)), 'm/s'
        write(*,*) 'min dv for inc ', y(i_inc), ' : ', minval(z(:,i_inc)), 'm/s'

        ! also make plots for each inc value:
        write(istr,'(I10)') int(y(i_inc))
        call plt%destroy()
        call plt%initialize(grid=.true.,xlabel='LAN (deg)',&
                            ylabel='$\Delta v$ (m/s)',figsize=[10,10],&
                            title='Lunar Orbit Maintenance : deadband = '//&
                            trim(deadband_alt_str)//' km : dt = '//&
                            trim(dt_max_str)//' days : Inc='//trim(adjustl(istr)),&
                            real_fmt='(E9.3)')
        call plt%add_plot(x, z(:,i_inc), linestyle='b-', label='dv', istat=istat)
        call plt%savefig('lom_INC='//trim(adjustl(istr))//'.png',&
                          pyfile='lom_INC='//trim(adjustl(istr))//'.py',istat=istat)
    end do

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Read the config file and populate the global variables.

    subroutine read_config_file()

    use argv_module, only: argv
    use json_module, only: json_file

    implicit none

    character(len=:),allocatable :: filename
    type(json_file) :: json
    logical :: status_ok, found
    character(len=:),allocatable :: error_msg

    filename = argv(1) ! get the first argument
    if (filename == '') then
        error stop 'The first command line arg should be the config file name'
    else

        write(*,'(A)') 'Reading config file: '//trim(filename)

        call json%initialize() ! no optional args for now
        call json%load(filename)
        if (json%failed()) then
            call json%check_for_errors(error_msg=error_msg)
            error stop error_msg
        else
            ! populate the run variables:
            call json%get('dt_max',       dt_max,       found); if (.not. found) error stop 'dt_max not found in config file.'
            call json%get('et0',          et0,          found); if (.not. found) error stop 'et0 not found in config file.'
            call json%get('tru0',         tru0,         found); if (.not. found) error stop 'tru0 not found in config file.'
            call json%get('alt0',         alt0,         found); if (.not. found) error stop 'alt0 not found in config file.'
            call json%get('deadband_alt', deadband_alt, found); if (.not. found) error stop 'deadband_alt not found in config file.'
            call json%get('inc_start',    inc_start,    found); if (.not. found) error stop 'inc_start not found in config file.'
            call json%get('inc_stop',     inc_stop,     found); if (.not. found) error stop 'inc_stop not found in config file.'
            call json%get('inc_step',     inc_step,     found); if (.not. found) error stop 'inc_step not found in config file.'
            call json%get('lan_start',    lan_start,    found); if (.not. found) error stop 'lan_start not found in config file.'
            call json%get('lan_stop',     lan_stop,     found); if (.not. found) error stop 'lan_stop not found in config file.'
            call json%get('lan_step',     lan_step,     found); if (.not. found) error stop 'lan_step not found in config file.'
            call json%get('grav_n',       grav_n,       found); if (.not. found) error stop 'grav_n not found in config file.'
            call json%get('grav_m',       grav_m,       found); if (.not. found) error stop 'grav_m not found in config file.'
        end if

        call json%destroy()

    end if

    end subroutine read_config_file
!*****************************************************************************************

!*****************************************************************************************
!>
!  Allocate the arrays.

    subroutine size_arrays(xstart, xstop, xstep, ystart, ystop, ystep, x, y, z)

    implicit none

    real(wp),intent(in) :: xstart, xstop, xstep, ystart, ystop, ystep
    real(wp),dimension(:),allocatable,intent(out) :: x,y
    real(wp),dimension(:,:),allocatable,intent(out) :: z  !! f(x,y)

    integer :: i,j,nx,ny
    real(wp) :: tmp

    nx = 1
    tmp = xstart
    x = [tmp]
    do
        tmp = tmp + xstep
        if (tmp>xstop) exit
        nx = nx + 1
        x = [x,tmp]
    end do

    ny = 1
    tmp = ystart
    y = [tmp]
    do
        tmp = tmp + ystep
        if (tmp>ystop) exit
        ny = ny + 1
        y = [y,tmp]
    end do

    allocate(z(nx,ny))

    z = 0.0_wp

    end subroutine size_arrays
!*****************************************************************************************

!*****************************************************************************************
    end program main
!*****************************************************************************************
