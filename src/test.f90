!*****************************************************************************************
!> author: Jacob Williams
!  date: 4/14/2015
!  license: BSD
!
!  Unit test for [[pyplot_module]].

    program test

    use pyplot_module, only : pyplot, wp => pyplot_wp

    implicit none

    integer,parameter :: n = 100

    real(wp), dimension(:),allocatable   :: x     !! x values
    real(wp), dimension(:),allocatable   :: y     !! y values
    real(wp), dimension(:),allocatable   :: xerr  !! error values  
    real(wp), dimension(:),allocatable   :: yerr  !! error values for bar chart
    real(wp), dimension(:),allocatable   :: sx    !! sin(x) values
    real(wp), dimension(:),allocatable   :: cx    !! cos(x) values
    real(wp), dimension(:),allocatable   :: zx    !!
    real(wp), dimension(:),allocatable   :: tx    !! sin(x)*cos(x) values
    real(wp), dimension(:,:),allocatable :: z     !! z matrix for contour plot
    real(wp), dimension(:,:),allocatable :: mat   !! image values
    type(pyplot)             :: plt   !! pytplot handler
    integer                  :: i     !! counter
    integer                  :: j     !! counter
    real(wp)                 :: r2    !! temp variable
    integer                  :: istat !! status code

    real(wp),parameter :: pi = acos(-1.0_wp)
    real(wp),parameter :: deg2rad = pi/180.0_wp

    ! size arrays:
    allocate(x(n))   
    allocate(y(n))   
    allocate(yerr(n))
    allocate(sx(n))  
    allocate(cx(n))  
    allocate(zx(n))  
    allocate(tx(n))  
    allocate(z(n,n))   
    allocate(mat(n,n)) 

    !generate some data:
    x    = [(real(i,wp), i=0,size(x)-1)]/5.0_wp
    sx   = sin(x)
    cx   = cos(x)
    tx   = sx * cx
    zx   = 0.0_wp
    yerr = abs(sx*.25_wp)

    do i=1,n
        do j=1,n
            mat(i,j) = sin(real(i,wp)*real(j,wp))
        end do
    end do

    !2d line plot:
    call plt%initialize(grid=.true.,xlabel='angle (rad)',figsize=[20,10],&
                        title='plot test',legend=.true.,axis_equal=.true.,&
                        tight_layout=.true.)
    call plt%add_plot(x,sx,label='$\sin (x)$',linestyle='b-o',markersize=5,linewidth=2,istat=istat)
    call plt%add_plot(x,cx,label='$\cos (x)$',linestyle='r-o',markersize=5,linewidth=2,istat=istat)
    call plt%add_plot(x,tx,label='$\sin (x) \cos (x)$',linestyle='g-o',markersize=2,linewidth=1,istat=istat)
    call plt%savefig('plottest.png', pyfile='plottest.py',istat=istat)

    !bar chart:
    tx = 0.1_wp !for bar width
    call plt%initialize(grid=.true.,xlabel='angle (rad)',&
                        title='bar test',legend=.true.,figsize=[20,10],&
                        font_size = 20,&
                        axes_labelsize = 20,&
                        xtick_labelsize = 20,&
                        ytick_labelsize = 20,&
                        legend_fontsize = 20, raw_strings=.true. )
    call plt%add_bar(x=x,height=sx,width=tx,label='$\sin (x)$',&
                        color='r',yerr=yerr,xlim=[0.0_wp, 20.0_wp],align='center',istat=istat)
    call plt%savefig('bartest.png', pyfile='bartest.py',istat=istat)

    !contour plot:
    x = [(real(i,wp), i=0,n-1)]/100.0_wp
    y = [(real(i,wp), i=0,n-1)]/100.0_wp
    do i=1,n
        do j=1,n
            r2 = x(i)**2 + y(j)**2
            z(i,j) = sin(x(i))*cos(y(j))*sin(r2)/(1.0_wp+log(r2+1.0_wp))
        end do
    end do
    call plt%initialize(grid=.true.,xlabel='x angle (rad)',&
                        ylabel='y angle (rad)',figsize=[10,10],&
                        title='Contour plot test', real_fmt='*',&
                        axisbelow=.false.)
    call plt%add_contour(x, y, z, linestyle='-', &
                         filled=.true., cmap='bone', colorbar=.true.,&
                         istat=istat)
    call plt%savefig('contour.png',pyfile='contour.py',istat=istat)

    !image plot:
    call plt%initialize(grid=.true.,xlabel='x',ylabel='y',figsize=[20,20],&
                        title='imshow test',&
                        real_fmt='(F9.3)')
    call plt%add_imshow(mat,xlim=[0.0_wp, 100.0_wp],ylim=[0.0_wp, 100.0_wp],istat=istat)
    call plt%savefig('imshow.png', pyfile='imshow.py',istat=istat)

    !histogram chart:
    x = [0.194,0.501,-1.241,1.425,-2.217,-0.342,-0.979,0.909,0.994,0.101,       &
         -0.131,-0.627,0.463,1.404,0.036,-2.000,0.109,1.250,-1.035,-1.115,      &
          0.935,0.496,1.100,0.770,-1.307,-0.693,-0.072,-1.331,-0.701,           &
         -0.494,0.666,-0.313,-0.430,-0.654,1.638,-0.334,-0.418,0.550,-0.034,    &
          0.473,0.704,0.801,-0.157,0.055,-0.057,-1.049,-1.022,0.495,0.756,      &
          0.149,0.543,-0.813,-0.171,-0.994,-1.532,0.502,1.324,-0.593,-0.467,    &
          0.372,-0.904,1.255,0.931,-0.779,1.529,-0.036,0.783,0.292,-0.792,      &
          -0.223,-0.325,0.225,-0.492,-0.941,0.065,1.300,-1.241,-1.124,-0.499,   &
          1.233,-0.845,-0.948,-1.060,1.103,-1.154,-0.594,0.335,-1.423,0.571,    &
         -0.903,1.129,-0.372,-1.043,-1.327,0.147,1.056,1.068,-0.699,0.988,-0.630]

    call plt%initialize(grid=.true.,xlabel='x',&
                        title='hist test',&
                        legend=.true.,figsize=[20,10],&
                        font_size = 20,&
                        axes_labelsize = 20,&
                        xtick_labelsize = 20,&
                        ytick_labelsize = 20,&
                        legend_fontsize = 20 )

    call plt%add_hist(x=x, label='x',istat=istat)
    call plt%savefig('histtest1.png', pyfile='histtest1.py',istat=istat)

    call plt%initialize(grid=.true.,xlabel='x',&
                        title='cumulative hist test',&
                        legend=.true.,figsize=[20,10],&
                        font_size = 20,&
                        axes_labelsize = 20,&
                        xtick_labelsize = 20,&
                        ytick_labelsize = 20,&
                        legend_fontsize = 20 )

    call plt%add_hist(x=x, label='x', bins=8, cumulative=.true.,istat=istat)
    call plt%savefig('histtest2.png', &
                        pyfile='histtest2.py', &
                        dpi='200', &
                        transparent=.true.,istat=istat)

    !3D plot:
    call plt%initialize(grid=.true.,&
                        xlabel='x (km)',ylabel='y (km)',zlabel='z (km)',&
                        title='Orbit',&
                        axis_equal=.true.,&
                        mplot3d=.true.,&
                        figsize=[20,10] )

    x    = [(real(i,wp), i=0,size(x)-1)]/5.0_wp
    sx   = 7000.0_wp * cos(x * deg2rad)
    cx   = 7000.0_wp * sin(x * deg2rad)
    zx = 0.0_wp
    call plt%add_3d_plot(sx,cx,zx,&
                            label='orbit',linestyle='r-',&
                            linewidth=2,istat=istat)
    call plt%add_sphere(6378.0_wp,0.0_wp,0.0_wp,0.0_wp,&
                        color='Blue',n_facets=500,&
                        antialiased=.true.,istat=istat)
    call plt%savefig('orbit.png', &
                      pyfile='orbit.py', &
                      dpi='200', &
                      transparent=.true.,istat=istat)

    ! Errorbar plot:
    call plt%initialize(grid=.true.,&
                        xlabel='x',ylabel='y',&
                        title='Errorbar Plot Example',&
                        figsize=[20,10] )

    x = [(real(i,wp), i=0, 360, 10)]
    y = 10.0_wp * cos(x * deg2rad)
    xerr =  sin(x * deg2rad) * 5.0_wp
    yerr =  sin(y * deg2rad) * 10.0_wp

    call plt%add_errorbar(x, y, label='y', linestyle='.', &
                            xerr=xerr, yerr=yerr, istat=istat)
    call plt%savefig('errorbar.png', &
                     pyfile='errorbar.py', &
                     dpi='200', &
                     transparent=.true.,istat=istat)

    end program test
!*****************************************************************************************
