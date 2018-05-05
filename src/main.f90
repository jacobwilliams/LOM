    program main

    use altitude_maintenance_module, only: segment
    use fortran_astrodynamics_toolkit, only: wp,km2m
    use pyplot_module, only : pyplot

    implicit none

    real(wp),parameter :: dt_max = 10.0_wp    !! how long to propagate (days)
    real(wp),parameter :: et0  = 0.0_wp      !! initial ephemeris time (sec)
                                             !! (only matters if including Earth/Sun perturbations)

    real(wp) :: inc0             !! initial inclination - IAU_MOON of date (deg)
    real(wp) :: ran0             !! initial RAAN - IAU_MOON of date (deg)
    integer  :: n_dvs            !! number of DV maneuvers performed
    real(wp) :: dv_total         !! total DV (km/s)
    real(wp),dimension(6) :: xf  !! final state - inertial frame (km, km/s)
    integer :: i_inc             !! inclination counter
    integer :: i_raan            !! raan counter
    type(pyplot) :: plt          !! pyplot handler
    real(wp),dimension(:),allocatable :: x    !! x array for plot (raan)
    real(wp),dimension(:),allocatable :: y    !! y array for plot (inc)
    real(wp),dimension(:,:),allocatable :: z  !! z array for plot (dv)
    integer :: istat             !! pyplot status code
    character(len=10) :: istr    !! for integer to string conversion
    character(len=10) :: dt_max_str,deadband_alt_str

    ! test case:
    ! real(wp),parameter :: inc_start = 80.0
    ! real(wp),parameter :: inc_stop  = 85.0
    ! real(wp),parameter :: inc_step  = 2.0
    ! real(wp),parameter :: lan_start = -180.0
    ! real(wp),parameter :: lan_stop  = 179.0
    ! real(wp),parameter :: lan_step  = 2.0

    ! full data set:
    ! real(wp),parameter :: inc_start = 80.0
    ! real(wp),parameter :: inc_stop  = 180.0
    ! real(wp),parameter :: inc_step  = 2.0
    ! real(wp),parameter :: lan_start = -180.0
    ! real(wp),parameter :: lan_stop  = 180.0
    ! real(wp),parameter :: lan_step  = 4.0

    real(wp),parameter :: inc_start = 80.0
    real(wp),parameter :: inc_stop  = 180.0
    real(wp),parameter :: inc_step  = 10.0
    real(wp),parameter :: lan_start = -180.0
    real(wp),parameter :: lan_stop  = 180.0
    real(wp),parameter :: lan_step  = 20.0

    type(segment) :: seg  !! the integrator

    real(wp),parameter :: alt0 = 100.0_wp  !! initial altitude for circular orbit (km)
    real(wp),parameter :: deadband_alt = 10.0_wp !! altitude below initial to trigger periapsis raise (km)

    ! initialize the segment:
    call seg%initialize_seg(alt0,deadband_alt)

    write(dt_max_str,'(I10)') int(dt_max); dt_max_str = adjustl(dt_max_str)
    write(deadband_alt_str,'(I10)') int(deadband_alt); deadband_alt_str = adjustl(deadband_alt_str)

    call plt%initialize(grid=.true.,xlabel='LAN (deg)',&
                        ylabel='Inc (deg)',figsize=[20,10],&
                        title='Lunar Orbit Maintenance $\Delta v$ (m/s) : deadband = '//&
                                trim(deadband_alt_str)//' km : dt = '//trim(dt_max_str)//' days',&
                        real_fmt='(E9.3)')

    ! initialize the indep arrays:
    call size_arrays(lan_start,lan_stop,lan_step,inc_start,inc_stop,inc_step,x,y,z)

    !allocate(x(lan_start:lan_stop))
    !allocate(y(inc_start:inc_stop))
    !allocate(z(lan_start:lan_stop, inc_start:inc_stop))

    do i_inc = 1, size(y)

        inc0 = y(i_inc)

        do i_raan = 1, size(x)

            ran0 = x(i_raan)

            write(*,*) ''
            write(*,*) '============'
            write(*,*) 'inc, ran', inc0, ran0
            write(*,*) '============'
            write(*,*) ''

            call seg%altitude_maintenance(et0,inc0,ran0,dt_max,n_dvs,dv_total,xf)

            z(i_raan,i_inc) = dv_total * km2m ! in meters

        end do

    end do

    call plt%add_contour(x, y, z, label='contour', linestyle='-', &
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

    subroutine size_arrays(xstart, xstop, xstep, ystart, ystop, ystep, x, y, z)

    implicit none

    real(wp),intent(in) :: xstart, xstop, xstep, ystart, ystop, ystep
    real(wp),dimension(:),allocatable,intent(out) :: x,y
    real(wp),dimension(:,:),allocatable,intent(out) :: z

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

    end program main
