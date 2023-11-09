program main

    implicit none

    integer ( kind = 4 ) p
    real ( kind = 8 ) radius

    integer :: io, stat
    logical :: exists

    ! Check if the file exists
    inquire(file="sphere_cubed_grid_points.csv", exist=exists)

    ! If the file exists, open and delete it
    if (exists) then
        open(file="sphere_cubed_grid_points.csv", newunit=io, iostat=stat)
        if (stat == 0) then
            close(io, status="delete", iostat=stat)
        end if
    end if

    do p = 1, 10
      radius = (1.0D+00 * p)/10
      call sphere_cubed_grid_lines_display_test ( radius )
    end do

    stop 0
  end

  subroutine sphere_cubed_grid_lines_display_test ( radius )

    implicit none
  
    real ( kind = 8 ), allocatable :: line_data(:,:,:)
    integer ( kind = 4 ) line_num
    integer ( kind = 4 ) n
    character ( len = 255 ) prefix
    real ( kind = 8 ) radius
  ! -------------------------------------------------------------------------------------
  n = 10
  ! -------------------------------------------------------------------------------------
  
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'SPHERE_CUBED_GRID_LINES_DISPLAY_TEST'
    write ( *, '(a)' ) '  SPHERE_CUBED_GRID_LINES_DISPLAY displays the lines'
    write ( *, '(a,f6.2)' ) '  on a cubed sphere grid of radius ', radius
    write ( *, '(a,i2,a,i2,a)' ) '  Each cube face is divided into ', n, ' by ', n, ' subfaces'
  
    call sphere_cubed_grid_line_count ( n, line_num )
  
    write ( *, '(a,i4)' ) '  The number of grid lines is ', line_num
  
    allocate ( line_data(1:3,1:2,1:line_num) )
  
    call sphere_cubed_grid_lines ( n, line_num, line_data, radius )
  
    prefix = 'sphere_cubed_grid'
  
    call sphere_cubed_grid_lines_display ( line_num, line_data, prefix )
  
    deallocate ( line_data )
  
    return
  end


  subroutine get_unit ( iunit )

    implicit none
  
    integer ( kind = 4 ) i
    integer ( kind = 4 ) ios
    integer ( kind = 4 ) iunit
    logical ( kind = 4 ) lopen
  
    iunit = 0
  
    do i = 1, 99
  
      if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then
  
        inquire ( unit = i, opened = lopen, iostat = ios )
  
        if ( ios == 0 ) then
          if ( .not. lopen ) then
            iunit = i
            return
          end if
        end if
  
      end if
  
    end do
  
    return
  end

  subroutine sphere_cubed_grid_ijk_to_xyz ( n, i, j, k, xyz, radius )
  
    implicit none
  
    integer ( kind = 4 ) i
    integer ( kind = 4 ) j
    integer ( kind = 4 ) k
    integer ( kind = 4 ) n
    real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
    real ( kind = 8 ) radius
    real ( kind = 8 ) xc
    real ( kind = 8 ) xyz(3)
    real ( kind = 8 ) xyzn
    real ( kind = 8 ) yc
    real ( kind = 8 ) zc
  
    if ( i == 0 ) then
      xc = -1.0D+00
    else if ( i == n ) then
      xc = +1.0D+00
    else
      xc = tan ( real ( 2 * i - n, kind = 8 ) * 0.25D+00 * r8_pi &
        / real ( n, kind = 8 ) )
    end if
  
    if ( j == 0 ) then
      yc = -1.0D+00
    else if ( j == n ) then
      yc = +1.0D+00
    else
      yc = tan ( real ( 2 * j - n, kind = 8 ) * 0.25D+00 * r8_pi &
        / real ( n, kind = 8 ) )
    end if
  
    if ( k == 0 ) then
      zc = -1.0D+00
    else if ( k == n ) then
      zc = +1.0D+00
    else
      zc = tan ( real ( 2 * k - n, kind = 8 ) * 0.25D+00 * r8_pi &
        / real ( n, kind = 8 ) )
    end if

    xyzn = sqrt ( xc ** 2 + yc ** 2 + zc ** 2 )
  
    xyz(1) = xc / xyzn
    xyz(2) = yc / xyzn
    xyz(3) = zc / xyzn
    
    ! After computing xyz using subroutine:
    xyz(1) = xyz(1) * radius
    xyz(2) = xyz(2) * radius
    xyz(3) = xyz(3) * radius

    return
  end
  
  subroutine sphere_cubed_grid_line_count ( n, line_num )
  
    implicit none
  
    integer ( kind = 4 ) line_num
    integer ( kind = 4 ) n
  
    line_num = 0

    if ( n == 1 ) then
      line_num = 12
      return
 
    else
      line_num = line_num + 8 * 3
    end if

    if ( 2 < n ) then
      line_num = line_num + 12 * ( n - 2 )
    end if
 
    if ( 1 < n ) then
      line_num = line_num + 6 * 2 * n * ( n - 1 )
    end if
  
    return
  end

  subroutine sphere_cubed_grid_lines ( n, line_num, line_data, radius )

    implicit none
  
    integer ( kind = 4 ) line_num
  
    integer ( kind = 4 ) i
    integer ( kind = 4 ) j
    integer ( kind = 4 ) l
    real ( kind = 8 ) line_data(3,2,line_num)
    integer ( kind = 4 ) n
    real ( kind = 8 ) radius 
  
    l = 0
 
    if ( n == 1 ) then
  
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, 0, line_data(1:3,1,l), radius )
      call sphere_cubed_grid_ijk_to_xyz ( n, n, 0, 0, line_data(1:3,2,l), radius )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, n, 0, 0, line_data(1:3,1,l), radius )
      call sphere_cubed_grid_ijk_to_xyz ( n, n, n, 0, line_data(1:3,2,l), radius )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, n, n, 0, line_data(1:3,1,l), radius )
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, n, 0, line_data(1:3,2,l), radius )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, n, 0, line_data(1:3,1,l), radius )
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, 0, line_data(1:3,2,l), radius )
  
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, n, line_data(1:3,1,l), radius )
      call sphere_cubed_grid_ijk_to_xyz ( n, n, 0, n, line_data(1:3,2,l), radius )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, n, 0, n, line_data(1:3,1,l), radius )
      call sphere_cubed_grid_ijk_to_xyz ( n, n, n, n, line_data(1:3,2,l), radius )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, n, n, n, line_data(1:3,1,l), radius )
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, n, n, line_data(1:3,2,l), radius )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, n, n, line_data(1:3,1,l), radius )
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, n, line_data(1:3,2,l), radius )
  
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, 0, line_data(1:3,1,l), radius )
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, n, line_data(1:3,2,l), radius )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, n, 0, 0, line_data(1:3,1,l), radius )
      call sphere_cubed_grid_ijk_to_xyz ( n, n, 0, n, line_data(1:3,2,l), radius )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, n, n, 0, line_data(1:3,1,l), radius )
      call sphere_cubed_grid_ijk_to_xyz ( n, n, n, n, line_data(1:3,2,l), radius )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, n, 0, line_data(1:3,1,l), radius )
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, n, n, line_data(1:3,2,l), radius )
      return

    else
  
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, 0, line_data(1:3,1,l), radius )
      call sphere_cubed_grid_ijk_to_xyz ( n, 1, 0, 0, line_data(1:3,2,l), radius )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, 0, line_data(1:3,1,l), radius )
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, 1, 0, line_data(1:3,2,l), radius )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, 0, line_data(1:3,1,l), radius )
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, 1, line_data(1:3,2,l), radius )
  
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, n,   0, 0, line_data(1:3,1,l), radius )
      call sphere_cubed_grid_ijk_to_xyz ( n, n-1, 0, 0, line_data(1:3,2,l), radius )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, n, 0, 0, line_data(1:3,1,l), radius )
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, 1, 0, line_data(1:3,2,l), radius )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, n, 0, 0, line_data(1:3,1,l), radius )
      call sphere_cubed_grid_ijk_to_xyz ( n, n, 0, 1, line_data(1:3,2,l), radius )
  
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, n,   n, 0, line_data(1:3,1,l), radius )
      call sphere_cubed_grid_ijk_to_xyz ( n, n-1, n, 0, line_data(1:3,2,l), radius )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, n, n,   0, line_data(1:3,1,l), radius )
      call sphere_cubed_grid_ijk_to_xyz ( n, n, n-1, 0, line_data(1:3,2,l), radius )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, n, n, 0, line_data(1:3,1,l), radius )
      call sphere_cubed_grid_ijk_to_xyz ( n, n, n, 1, line_data(1:3,2,l), radius )
  
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, n, 0, line_data(1:3,1,l), radius )
      call sphere_cubed_grid_ijk_to_xyz ( n, 1, n, 0, line_data(1:3,2,l), radius )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, n,   0, line_data(1:3,1,l), radius )
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, n-1, 0, line_data(1:3,2,l), radius )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, n, 0, line_data(1:3,1,l), radius )
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, n, 1, line_data(1:3,2,l), radius )
  
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, n, line_data(1:3,1,l), radius )
      call sphere_cubed_grid_ijk_to_xyz ( n, 1, 0, n, line_data(1:3,2,l), radius )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, n, line_data(1:3,1,l), radius )
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, 1, n, line_data(1:3,2,l), radius )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, n,   line_data(1:3,1,l), radius )
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, n-1, line_data(1:3,2,l), radius )
  
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, n,   0, n, line_data(1:3,1,l), radius )
      call sphere_cubed_grid_ijk_to_xyz ( n, n-1, 0, n, line_data(1:3,2,l), radius )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, n, 0, n, line_data(1:3,1,l), radius )
      call sphere_cubed_grid_ijk_to_xyz ( n, n, 1, n, line_data(1:3,2,l), radius )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, n, 0, n,   line_data(1:3,1,l), radius )
      call sphere_cubed_grid_ijk_to_xyz ( n, n, 0, n-1, line_data(1:3,2,l), radius )
  
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, n,   n, n, line_data(1:3,1,l), radius )
      call sphere_cubed_grid_ijk_to_xyz ( n, n-1, n, n, line_data(1:3,2,l), radius )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, n, n,   n, line_data(1:3,1,l), radius )
      call sphere_cubed_grid_ijk_to_xyz ( n, n, n-1, n, line_data(1:3,2,l), radius )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, n, n, n,   line_data(1:3,1,l), radius )
      call sphere_cubed_grid_ijk_to_xyz ( n, n, n, n-1, line_data(1:3,2,l), radius )
  
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, n, n, line_data(1:3,1,l), radius )
      call sphere_cubed_grid_ijk_to_xyz ( n, 1, n, n, line_data(1:3,2,l), radius )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, n,   n, line_data(1:3,1,l), radius )
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, n-1, n, line_data(1:3,2,l), radius )
      l = l + 1
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, n, n,   line_data(1:3,1,l), radius )
      call sphere_cubed_grid_ijk_to_xyz ( n, 0, n, n-1, line_data(1:3,2,l), radius )
  
    end if
  
    if ( 2 < n ) then
  
      do i = 1, n - 2
        l = l + 1
        call sphere_cubed_grid_ijk_to_xyz ( n, i,   0, 0, line_data(1:3,1,l), radius )
        call sphere_cubed_grid_ijk_to_xyz ( n, i+1, 0, 0, line_data(1:3,2,l), radius )
      end do
      do i = 1, n - 2
        l = l + 1
        call sphere_cubed_grid_ijk_to_xyz ( n, n,   i, 0, line_data(1:3,1,l), radius )
        call sphere_cubed_grid_ijk_to_xyz ( n, n, i+1, 0, line_data(1:3,2,l), radius )
      end do
      do i = 1, n - 2
        l = l + 1
        call sphere_cubed_grid_ijk_to_xyz ( n, n-i,   n, 0, line_data(1:3,1,l), radius )
        call sphere_cubed_grid_ijk_to_xyz ( n, n-i-1, n, 0, line_data(1:3,2,l), radius )
      end do
      do i = 1, n - 2
        l = l + 1
        call sphere_cubed_grid_ijk_to_xyz ( n, 0, n-i,   0, line_data(1:3,1,l), radius )
        call sphere_cubed_grid_ijk_to_xyz ( n, 0, n-i-1, 0, line_data(1:3,2,l), radius )
      end do
  
      do i = 1, n - 2
        l = l + 1
        call sphere_cubed_grid_ijk_to_xyz ( n, i,   0, n, line_data(1:3,1,l), radius )
        call sphere_cubed_grid_ijk_to_xyz ( n, i+1, 0, n, line_data(1:3,2,l), radius )
      end do
      do i = 1, n - 2
        l = l + 1
        call sphere_cubed_grid_ijk_to_xyz ( n, n,   i, n, line_data(1:3,1,l), radius )
        call sphere_cubed_grid_ijk_to_xyz ( n, n, i+1, n, line_data(1:3,2,l), radius )
      end do
      do i = 1, n - 2
        l = l + 1
        call sphere_cubed_grid_ijk_to_xyz ( n, n-i,   n, n, line_data(1:3,1,l), radius )
        call sphere_cubed_grid_ijk_to_xyz ( n, n-i-1, n, n, line_data(1:3,2,l), radius )
      end do
      do i = 1, n - 2
        l = l + 1
        call sphere_cubed_grid_ijk_to_xyz ( n, 0, n-i,   n, line_data(1:3,1,l), radius )
        call sphere_cubed_grid_ijk_to_xyz ( n, 0, n-i-1, n, line_data(1:3,2,l), radius )
      end do
  
      do i = 1, n - 2
        l = l + 1
        call sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, i,   line_data(1:3,1,l), radius )
        call sphere_cubed_grid_ijk_to_xyz ( n, 0, 0, i+1, line_data(1:3,2,l), radius )
      end do
      do i = 1, n - 2
        l = l + 1
        call sphere_cubed_grid_ijk_to_xyz ( n, n, 0, i,   line_data(1:3,1,l), radius )
        call sphere_cubed_grid_ijk_to_xyz ( n, n, 0, i+1, line_data(1:3,2,l), radius )
      end do
      do i = 1, n - 2
        l = l + 1
        call sphere_cubed_grid_ijk_to_xyz ( n, n, n, i,   line_data(1:3,1,l), radius )
        call sphere_cubed_grid_ijk_to_xyz ( n, n, n, i+1, line_data(1:3,2,l), radius )
      end do
      do i = 1, n - 2
        l = l + 1
        call sphere_cubed_grid_ijk_to_xyz ( n, 0, n, i,   line_data(1:3,1,l), radius )
        call sphere_cubed_grid_ijk_to_xyz ( n, 0, n, i+1, line_data(1:3,2,l), radius )
      end do
  
    end if
 
    if ( 1 < n ) then

      do i = 1, n - 1
        do j = 0, n - 1
          l = l + 1
          call sphere_cubed_grid_ijk_to_xyz ( n, i, j,   0, line_data(1:3,1,l), radius )
          call sphere_cubed_grid_ijk_to_xyz ( n, i, j+1, 0, line_data(1:3,2,l), radius )
        end do
      end do
      do j = 1, n - 1
        do i = 0, n - 1
          l = l + 1
          call sphere_cubed_grid_ijk_to_xyz ( n, i,   j, 0, line_data(1:3,1,l), radius )
          call sphere_cubed_grid_ijk_to_xyz ( n, i+1, j, 0, line_data(1:3,2,l), radius )
        end do
      end do
  
      do i = 1, n - 1
        do j = 0, n - 1
          l = l + 1
          call sphere_cubed_grid_ijk_to_xyz ( n, i, j,   n, line_data(1:3,1,l), radius )
          call sphere_cubed_grid_ijk_to_xyz ( n, i, j+1, n, line_data(1:3,2,l), radius )
        end do
      end do
      do j = 1, n - 1
        do i = 0, n - 1
          l = l + 1
          call sphere_cubed_grid_ijk_to_xyz ( n, i,   j, n, line_data(1:3,1,l), radius )
          call sphere_cubed_grid_ijk_to_xyz ( n, i+1, j, n, line_data(1:3,2,l), radius )
        end do
      end do

      do i = 1, n - 1
        do j = 0, n - 1
          l = l + 1
          call sphere_cubed_grid_ijk_to_xyz ( n, i, 0, j,   line_data(1:3,1,l), radius   )
          call sphere_cubed_grid_ijk_to_xyz ( n, i, 0, j+1, line_data(1:3,2,l), radius )
        end do
      end do
      do j = 1, n - 1
        do i = 0, n - 1
          l = l + 1
          call sphere_cubed_grid_ijk_to_xyz ( n, i,   0, j, line_data(1:3,1,l), radius )
          call sphere_cubed_grid_ijk_to_xyz ( n, i+1, 0, j, line_data(1:3,2,l), radius )
        end do
      end do

      do i = 1, n - 1
        do j = 0, n - 1
          l = l + 1
          call sphere_cubed_grid_ijk_to_xyz ( n, i, n, j,   line_data(1:3,1,l), radius   )
          call sphere_cubed_grid_ijk_to_xyz ( n, i, n, j+1, line_data(1:3,2,l), radius )
        end do
      end do
      do j = 1, n - 1
        do i = 0, n - 1
          l = l + 1
          call sphere_cubed_grid_ijk_to_xyz ( n, i,   n, j, line_data(1:3,1,l), radius )
          call sphere_cubed_grid_ijk_to_xyz ( n, i+1, n, j, line_data(1:3,2,l), radius )
        end do
      end do
 
      do i = 1, n - 1
        do j = 0, n - 1
          l = l + 1
          call sphere_cubed_grid_ijk_to_xyz ( n, 0, i, j,   line_data(1:3,1,l), radius   )
          call sphere_cubed_grid_ijk_to_xyz ( n, 0, i, j+1, line_data(1:3,2,l), radius )
        end do
      end do
      do j = 1, n - 1
        do i = 0, n - 1
          l = l + 1
          call sphere_cubed_grid_ijk_to_xyz ( n, 0, i,   j, line_data(1:3,1,l), radius )
          call sphere_cubed_grid_ijk_to_xyz ( n, 0, i+1, j, line_data(1:3,2,l), radius )
        end do
      end do

      do i = 1, n - 1
        do j = 0, n - 1
          l = l + 1
          call sphere_cubed_grid_ijk_to_xyz ( n, n, i, j,   line_data(1:3,1,l), radius   )
          call sphere_cubed_grid_ijk_to_xyz ( n, n, i, j+1, line_data(1:3,2,l), radius )
        end do
      end do
      do j = 1, n - 1
        do i = 0, n - 1
          l = l + 1
          call sphere_cubed_grid_ijk_to_xyz ( n, n, i,   j, line_data(1:3,1,l), radius )
          call sphere_cubed_grid_ijk_to_xyz ( n, n, i+1, j, line_data(1:3,2,l), radius )
        end do
      end do
  
    end if
  
    if ( l /= line_num ) then
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'SPHERE_CUBED_GRID_LINES - Fatal error!'
      write ( *, '(a,i6)' ) '  LINE_NUM = ', line_num
      write ( *, '(a,i6)' ) '  L = ', l
      stop 1
    end if
  
    return
  end

  subroutine sphere_cubed_grid_lines_display ( line_num, line_data, prefix )

    implicit none

    integer(kind=4) :: line_num

    integer(kind=4) :: i
    integer(kind=4) :: l
    real(kind=8) :: line_data(3, 2, line_num)
    character(len=255) :: line_filename
    integer(kind=4) :: line_unit
    character(len=*) :: prefix
    logical :: file_exists

    call get_unit(line_unit)
    line_filename = trim(prefix) // '_points.csv' ! Change file extension to .csv

    inquire(file=line_filename, exist=file_exists)

    if (file_exists) then
        open(unit=line_unit, file=line_filename, status='old', position='append')
    else
        open(unit=line_unit, file=line_filename, status='new')
    end if

    do l = 1, line_num
        do i = 1, 3
            ! Check if it's the last entry in the row
            if (i == 3) then
                ! Last entry, no comma
                write(line_unit, '(f14.6)', advance='no') line_data(i, 1, l)
            else
                ! Not the last entry, add a comma
                write(line_unit, '(f14.6, ",")', advance='no') line_data(i, 1, l)
            end if
        end do
        write(line_unit, *)
    end do
    close(unit=line_unit)
    write(*, '(a)') '  Updated line file "' // trim(line_filename) // '".'

end subroutine sphere_cubed_grid_lines_display




      subroutine sphere_cubed_grid_point_count ( n, ns )
    
        implicit none
      
        integer ( kind = 4 ) n
        integer ( kind = 4 ) ns
      
        ns = ( n + 1 ) ** 3 - ( n - 1 ) ** 3
      
        return
      end
  
  
  
      subroutine sphere_cubed_grid_points ( n, ns, xyz )
     
        implicit none
      
        integer ( kind = 4 ) ns
      
        integer ( kind = 4 ) n
        integer ( kind = 4 ) ns2
        real ( kind = 8 ) xyz(3,ns)
      
        ns2 = 0
     
        call sphere_cubed_grid_points_face ( n, 0, 0, 0, n, n, 0, ns2, xyz )
     
        call sphere_cubed_grid_points_face ( n, 0, 0, 1, 0,   n-1, n-1, ns2, xyz )
        call sphere_cubed_grid_points_face ( n, 0, n, 1, n-1, n,   n-1, ns2, xyz )
        call sphere_cubed_grid_points_face ( n, n, 1, 1, n,   n,   n-1, ns2, xyz )
        call sphere_cubed_grid_points_face ( n, 1, 0, 1, n,   0,   n-1, ns2, xyz )
     
        call sphere_cubed_grid_points_face ( n, 0, 0, n, n, n, n, ns2, xyz )
        
        if ( ns2 /= ns ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SPHERE_CUBED_GRID_POINTS - Fatal error!'
          write ( *, '(a,i8,a)' ) '  Expected to generated NS = ', ns, ' points.'
          write ( *, '(a,i8,a)' ) '  Generated ', ns2, ' points.'
          stop
        end if
      
        return
      end
  
      subroutine sphere_cubed_grid_points_face ( n, i1, j1, k1, i2, j2, k2, ns, xyz )
   
        implicit none
      
        integer ( kind = 4 ) i
        integer ( kind = 4 ) i1
        integer ( kind = 4 ) i2
        integer ( kind = 4 ) j
        integer ( kind = 4 ) j1
        integer ( kind = 4 ) j2
        integer ( kind = 4 ) k
        integer ( kind = 4 ) k1
        integer ( kind = 4 ) k2
        integer ( kind = 4 ) n
        integer ( kind = 4 ) ns
        real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
        real ( kind = 8 ) xyz(3,*)
        real ( kind = 8 ) xyzn
        real ( kind = 8 ) xc
        real ( kind = 8 ) yc
        real ( kind = 8 ) zc
      
        do i = i1, i2
      
          if ( i1 < i2 ) then
            xc = tan ( real ( 2 * i - n, kind = 8 ) * 0.25D+00 * r8_pi &
              / real ( n, kind = 8 ) )
          else if ( i1 == 0 ) then
            xc = -1.0D+00
          else if ( i1 == n ) then
            xc = +1.0D+00
          else
            xc = 0.0D+00
          end if
      
          do j = j1, j2
      
            if ( j1 < j2 ) then
              yc = tan ( real ( 2 * j - n, kind = 8 ) * 0.25D+00 * r8_pi &
                / real ( n, kind = 8 ) )
            else if ( j1 == 0 ) then
              yc = -1.0D+00
            else if ( j1 == n ) then
              yc = +1.0D+00
            else
              yc = 0.0D+00
            end if
      
            do k = k1, k2
      
              if ( k1 < k2 ) then
                zc = tan ( real ( 2 * k - n, kind = 8 ) * 0.25D+00 * r8_pi &
                  / real ( n, kind = 8 ) )
              else if ( k1 == 0 ) then
                zc = -1.0D+00
              else if ( k1 == n ) then
                zc = +1.0D+00
              else
                zc = 0.0D+00
              end if
      
              xyzn = sqrt ( xc ** 2 + yc ** 2 + zc ** 2 )
      
              ns = ns + 1
              xyz(1,ns) = xc / xyzn
              xyz(2,ns) = yc / xyzn
              xyz(3,ns) = zc / xyzn
      
            end do
          end do
        end do
      
        return
      end