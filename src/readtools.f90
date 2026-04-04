! Copyright 2008 Michal Kuraz, Petr Mayer, Copyright 2016  Michal Kuraz,
! Petr Mayer, Johanna Bloecher
!
! This file is part of DRUtES.
! DRUtES is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! DRUtES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with DRUtES. If not, see <http://www.gnu.org/licenses/>.

!> \file readtools.f90
!! \brief Tools for reading input config files.
!<

module readtools
  public :: file_error
  public :: fileread
  private :: read_int, read_int_std, read_int_array, read_real,  read_real_array, read_char, read_char_array
  public :: comment, reverse_comment
!  public :: readbcvals
!  public :: set_tensor
  public :: read_sep

  interface fileread
    module procedure read_int
    module procedure read_int_std
    module procedure read_int_array
    module procedure read_real
    module procedure read_real_quad
    module procedure read_real_array
    module procedure read_real_quad_array
    module procedure read_char
    module procedure read_logical
    module procedure read_char_array
    module procedure read_logical_array
  end interface fileread

contains

  subroutine read_int(i, fileid, errmsg, ranges, minimal, maximal, noexit)
    use typy
    use globals
    integer(kind=ikind), intent(out) :: i
    integer, intent(in) :: fileid
    character(len=*), intent(in), optional :: errmsg
    integer(kind=ikind), dimension(:), intent(in),  optional :: ranges
    integer, intent(in), optional :: minimal, maximal
    logical, intent(in), optional :: noexit

    logical :: go4exit = .false.
    integer :: ierr

    call comment(fileid)
    read(unit=fileid, fmt=*, iostat=ierr) i

    if (ierr /= 0) then
      if (present(errmsg)) then
        call file_error(fileid, errmsg)
      end if
      call file_error(fileid)
    end if

    ! --- Range checks using 'ranges' OR minimal/maximal ---
    if (present(ranges)) then
      if (i < ranges(1) .or. i > ranges(2)) then
        write(unit=terminal, fmt=*) " " //achar(27)//'[43m', "incorrect ranges of input parameter(s)" //achar(27)//'[0m'
        if (.not. present(noexit)) then
          go4exit = .true.
        else
          if (.not. noexit) then
            go4exit = .true.
          end if
        end if
        if (go4exit) then
          if (present(errmsg)) then
            call file_error(fileid,errmsg)
          else
            call file_error(fileid)
          end if
        end if
      end if
    else
      ! Fallback: use minimal/maximal if supplied
      if (present(minimal)) then
        if (i < minimal) then
          write(unit=terminal, fmt=*) " " //achar(27)//'[43m', "input value is below minimal allowed" //achar(27)//'[0m'
          if (.not. present(noexit)) then
            go4exit = .true.
          else
            if (.not. noexit) go4exit = .true.
          end if
          if (go4exit) then
            if (present(errmsg)) then
              call file_error(fileid,errmsg)
            else
              call file_error(fileid)
            end if
          end if
        end if
      end if

      if (present(maximal)) then
        if (i > maximal) then
          write(unit=terminal, fmt=*) " " //achar(27)//'[43m', "input value is above maximal allowed" //achar(27)//'[0m'
          if (.not. present(noexit)) then
            go4exit = .true.
          else
            if (.not. noexit) go4exit = .true.
          end if
          if (go4exit) then
            if (present(errmsg)) then
              call file_error(fileid,errmsg)
            else
              call file_error(fileid)
            end if
          end if
        end if
      end if
    end if

  end subroutine read_int


  subroutine read_int_std(i, fileid, errmsg, ranges, noexit)
    use typy
    use globals
    integer(kind=lowikind), intent(out) :: i
    integer, intent(in) :: fileid
    character(len=*), intent(in), optional :: errmsg
    integer(kind=ikind), dimension(:), intent(in),  optional :: ranges
    logical, intent(in), optional :: noexit

    integer :: ierr

    call comment(fileid)
    read(unit=fileid, fmt=*, iostat=ierr) i

    if (ierr /= 0) then
      if (.not. present(noexit)) then
        if (present(errmsg)) then
          call file_error(fileid,errmsg)
        else
          call file_error(fileid)
        end if
      end if
    end if

    if (present(ranges)) then
      if (i < ranges(1) .or. i > ranges(2)) then
        write(unit=terminal, fmt=*) " " //achar(27)//'[43m', "incorrect ranges of input parameter(s)" //achar(27)//'[0m'
        if (.not. present(noexit)) then
          if (present(errmsg)) then
            call file_error(fileid,errmsg)
          else
            call file_error(fileid)
          end if
        end if
      end if
    end if

  end subroutine read_int_std


  subroutine read_int_array(r, fileid, ranges, errmsg, checklen, noexit)
    use typy
    use globals
    integer(kind=ikind), dimension(:), intent(out) :: r
    integer, intent(in) :: fileid
    integer(kind=ikind), dimension(:), optional :: ranges
    character(len=*), intent(in), optional :: errmsg
    logical, intent(in), optional :: checklen
    logical, intent(in), optional :: noexit

    integer :: ierr, ierr2
    integer(kind=ikind) :: i, i1, i2, arraybound, current_pos, chpos
    real(kind=rkind), dimension(:), allocatable :: tmpdata
    logical :: terminate = .false.

    call comment(fileid)

    arraybound = 1
    allocate(tmpdata(arraybound))

    if (present(checklen)) then
      if (checklen) then
        current_pos = ftell(fileid)
        do
          call comment(fileid)
          read(unit=fileid, fmt=*, iostat=ierr) tmpdata(1:arraybound-1)
          i1 = ftell(fileid)
          backspace fileid
          call comment(fileid)
          read(unit=fileid, fmt=*, iostat=ierr) tmpdata(1:arraybound)
          i2 = ftell(fileid)
          backspace fileid

          if (ierr /= 0) then
            arraybound = arraybound - 1
            do
              chpos = ftell(fileid)
              if (chpos == current_pos) then
                exit
              else
                backspace(unit=fileid, iostat=ierr2)
              end if
            end do
            deallocate(tmpdata)
            exit
          end if

          if (i1 == i2) then
            arraybound = arraybound + 1
            if (ubound(tmpdata,1) < arraybound) then
              deallocate(tmpdata)
              allocate(tmpdata(2*arraybound))
            end if
          else
            arraybound = arraybound - 1
            do
              chpos = ftell(fileid)
              if (chpos == current_pos) then
                exit
              else
                backspace(unit=fileid, iostat=ierr2)
              end if
            end do
            deallocate(tmpdata)
            exit
          end if
        end do

        if (arraybound /= ubound(r,1)) then
          write(unit=terminal, fmt=*) " " //achar(27)//'[43m', "the line in your input file has " // &
               "incorrect number of values, check your inputs" //achar(27)//'[0m'
          call file_error(fileid, errmsg)
        end if
      end if
    end if

    call comment(fileid)
    read(unit=fileid, fmt=*, iostat=ierr) r

    if (ierr /= 0) then
      if (.not. present(noexit)) then
        terminate = .true.
      else
        if (noexit) then
          terminate = .false.
        else
          terminate = .true.
        end if
      end if
    end if

    if (terminate) then
      if (present(errmsg)) then
        call file_error(fileid,errmsg)
      end if
      if (.not. present(noexit)) then
        if (present(errmsg)) then
          call file_error(fileid,errmsg)
        else
          call file_error(fileid)
        end if
      end if
    end if

    if (present(ranges)) then
      do i = 1, ubound(r,1)
        if (r(i) < ranges(1) .or. r(i) > ranges(2)) then
          write(unit=terminal, fmt=*) " " //achar(27)//'[43m', "incorrect ranges of input parameter(s)" //achar(27)//'[0m'
          call file_error(fileid, errmsg)
        end if
      end do
    end if

  end subroutine read_int_array


  subroutine read_real(r, fileid, ranges, errmsg, noexit)
    use typy
    use globals
    real(kind=rkind), intent(out) :: r
    integer, intent(in) :: fileid
    real(kind=rkind), dimension(:), optional :: ranges
    character(len=*), intent(in), optional :: errmsg
    logical, intent(in), optional :: noexit

    integer :: ierr
    logical :: terminate = .false.

    call comment(fileid)
    read(unit=fileid, fmt=*, iostat=ierr) r

    if (present(ranges)) then
      if (r < ranges(1) .or. r > ranges(2)) then
        write(unit=terminal, fmt=*) " " //achar(27)//'[43m', "incorrect ranges of input parameter(s)" //achar(27)//'[0m'
        ierr = -1
      end if
    end if

    if (ierr /= 0) then
      if (.not. present(noexit)) then
        terminate = .true.
      else
        if (noexit) then
          terminate = .false.
        else
          terminate = .true.
        end if
      end if
    end if

    if (terminate) then
      if (present(errmsg)) then
        call file_error(fileid,errmsg)
      else
        call file_error(fileid)
      end if
    end if

  end subroutine read_real


  subroutine read_real_quad(r, fileid, ranges, errmsg, noexit)
    use typy
    use globals
    real(kind=qprec), intent(out) :: r
    integer, intent(in) :: fileid
    real(kind=qprec), dimension(:), optional :: ranges
    character(len=*), intent(in), optional :: errmsg
    logical, intent(in), optional :: noexit

    integer :: ierr
    logical :: terminate = .false.

    call comment(fileid)
    read(unit=fileid, fmt=*, iostat=ierr) r

    if (present(ranges)) then
      if (r < ranges(1) .or. r > ranges(2)) then
        write(unit=terminal, fmt=*) " " //achar(27)//'[43m', "incorrect ranges of input parameter(s)" //achar(27)//'[0m'
        ierr = -1
      end if
    end if

    if (ierr /= 0) then
      if (.not. present(noexit)) then
        terminate = .true.
      else
        if (noexit) then
          terminate = .false.
        else
          terminate = .true.
        end if
      end if
    end if

    if (terminate) then
      if (present(errmsg)) then
        call file_error(fileid,errmsg)
      else
        call file_error(fileid)
      end if
    end if

  end subroutine read_real_quad


  subroutine read_real_array(r, fileid, ranges, errmsg, checklen, noexit, eof)
    use typy
    use globals
    real(kind=rkind), dimension(:), intent(out) :: r
    integer, intent(in) :: fileid
    real(kind=rkind), dimension(:), optional :: ranges
    character(len=*), intent(in), optional :: errmsg
    logical, intent(in), optional :: checklen
    logical, intent(in), optional :: noexit
    logical, intent(out), optional :: eof

    integer :: ierr, ierr2, ierr3
    integer(kind=ikind) :: i, i1, i2, arraybound, current_pos, chpos
    real(kind=rkind), dimension(:), allocatable :: tmpdata
    logical :: terminate = .false.
    logical :: noexit_local

    call comment(fileid)

    if (present(noexit)) then
      noexit_local = noexit
    else
      noexit_local = .false.
    end if

    arraybound = 1
    allocate(tmpdata(arraybound))

    if (present(checklen)) then
      if (checklen) then
        current_pos = ftell(fileid)
        do
          call comment(fileid)
          read(unit=fileid, fmt=*, iostat=ierr) tmpdata(1:arraybound-1)
          i1 = ftell(fileid)
          backspace fileid
          call comment(fileid)
          read(unit=fileid, fmt=*, iostat=ierr3) tmpdata(1:arraybound)
          i2 = ftell(fileid)
          backspace fileid

          if (ierr /= 0 .or. ierr3 /= 0) then
            arraybound = arraybound - 1
            do
              chpos = ftell(fileid)
              if (chpos == current_pos) then
                exit
              else
                backspace(unit=fileid, iostat=ierr2)
              end if
            end do
            deallocate(tmpdata)
            exit
          end if

          if (i1 == i2) then
            arraybound = arraybound + 1
            if (ubound(tmpdata,1) < arraybound) then
              deallocate(tmpdata)
              allocate(tmpdata(2*arraybound))
            end if
          else
            arraybound = arraybound - 1
            do
              chpos = ftell(fileid)
              if (chpos == current_pos) then
                exit
              else
                backspace(unit=fileid, iostat=ierr2)
              end if
            end do
            deallocate(tmpdata)
            exit
          end if
        end do

        if (arraybound /= ubound(r,1)) then
          if (.not. noexit_local) then
            write(unit=terminal, fmt=*) " " //achar(27)//'[43m', "the line in your input file has " // &
                 "incorrect number of values, check your inputs" //achar(27)//'[0m'
            call file_error(fileid, errmsg)
          else
            if (present(eof)) then
              eof = .true.
            end if
            return
          end if
        end if
      end if
    end if

    call comment(fileid)
    read(unit=fileid, fmt=*, iostat=ierr) r
    if (present(eof)) then
      eof = .false.
    end if

    if (ierr /= 0) then
      terminate = .true.
    end if

    if (terminate) then
      if (present(errmsg)) then
        if (.not. noexit_local) call file_error(fileid,errmsg)
      else
        if (.not. noexit_local) call file_error(fileid)
      end if
    end if

    if (present(ranges)) then
      do i = 1, ubound(r,1)
        if (r(i) < ranges(1) .or. r(i) > ranges(2)) then
          write(unit=terminal, fmt=*) " " //achar(27)//'[43m', "incorrect ranges of input parameter(s)" //achar(27)//'[0m'
          if (.not. noexit_local) then
            call file_error(fileid, errmsg)
          else
            if (present(eof)) eof = .true.
          end if
        end if
      end do
    end if

  end subroutine read_real_array


  subroutine read_real_quad_array(r, fileid, ranges, errmsg, checklen, noexit, eof)
    use typy
    use globals
    real(kind=qprec), dimension(:), intent(out) :: r
    integer, intent(in) :: fileid
    real(kind=qprec), dimension(:), optional :: ranges
    character(len=*), intent(in), optional :: errmsg
    logical, intent(in), optional :: checklen
    logical, intent(in), optional :: noexit
    logical, intent(out), optional :: eof

    integer :: ierr, ierr2, ierr3
    integer(kind=ikind) :: i, i1, i2, arraybound, current_pos, chpos
    real(kind=rkind), dimension(:), allocatable :: tmpdata
    logical :: terminate = .false.
    logical :: noexit_local

    call comment(fileid)

    if (present(noexit)) then
      noexit_local = noexit
    else
      noexit_local = .false.
    end if

    arraybound = 1
    allocate(tmpdata(arraybound))

    if (present(checklen)) then
      if (checklen) then
        current_pos = ftell(fileid)
        do
          call comment(fileid)
          read(unit=fileid, fmt=*, iostat=ierr) tmpdata(1:arraybound-1)
          i1 = ftell(fileid)
          backspace fileid
          call comment(fileid)
          read(unit=fileid, fmt=*, iostat=ierr3) tmpdata(1:arraybound)
          i2 = ftell(fileid)
          backspace fileid

          if (ierr /= 0) then
            arraybound = arraybound - 1
            do
              chpos = ftell(fileid)
              if (chpos == current_pos) then
                exit
              else
                backspace(unit=fileid, iostat=ierr2)
              end if
            end do
            deallocate(tmpdata)
            exit
          end if

          if (i1 == i2) then
            arraybound = arraybound + 1
            if (ubound(tmpdata,1) < arraybound) then
              deallocate(tmpdata)
              allocate(tmpdata(2*arraybound))
            end if
          else
            arraybound = arraybound - 1
            do
              chpos = ftell(fileid)
              if (chpos == current_pos) then
                exit
              else
                backspace(unit=fileid, iostat=ierr2)
              end if
            end do
            deallocate(tmpdata)
            exit
          end if
        end do

        if (arraybound /= ubound(r,1)) then
          if (.not. noexit_local) then
            write(unit=terminal, fmt=*) " " //achar(27)//'[43m', "the line in your input file has " // &
                 "incorrect number of values, check your inputs" //achar(27)//'[0m'
            call file_error(fileid, errmsg)
          else
            if (present(eof)) then
              eof = .true.
            end if
            return
          end if
        end if
      end if
    end if

    call comment(fileid)
    read(unit=fileid, fmt=*, iostat=ierr) r
    if (present(eof)) then
      eof = .false.
    end if

    if (ierr /= 0) then
      terminate = .true.
    end if

    if (terminate) then
      if (present(errmsg)) then
        if (.not. noexit_local) call file_error(fileid,errmsg)
      else
        if (.not. noexit_local) call file_error(fileid)
      end if
    end if

    if (present(ranges)) then
      do i = 1, ubound(r,1)
        if (r(i) < ranges(1) .or. r(i) > ranges(2)) then
          write(unit=terminal, fmt=*) " " //achar(27)//'[43m', "incorrect ranges of input parameter(s)" //achar(27)//'[0m'
          if (.not. noexit_local) then
            call file_error(fileid, errmsg)
          else
            if (present(eof)) eof = .true.
          end if
        end if
      end do
    end if

  end subroutine read_real_quad_array


  subroutine read_char(ch, fileid, errmsg, options, noexit)
    use debug_tools
    use core_tools

    character(len=*), intent(out) :: ch
    integer, intent(in) :: fileid
    character(len=*), intent(in), optional :: errmsg
    character(len=*), dimension(:), intent(in), optional :: options
    logical, intent(in), optional :: noexit

    integer :: ierr, i
    logical :: ok = .true.
    logical :: terminate = .false.

    call comment(fileid)
    read(unit=fileid, fmt="(a)", iostat=ierr) ch

    if (present(options)) then
      ok = .false.
      do i = 1, ubound(options, 1)
        if (cut(ch) == cut(options(i))) then
          ok = .true.
          exit
        end if
      end do
    end if

    if (ierr /= 0 .or. .not. ok) then
      if (.not. present(noexit)) then
        terminate = .true.
      else
        if (noexit) then
          terminate = .false.
        else
          terminate = .true.
        end if
      end if
    end if

    if (terminate) then
      if (.not. ok) then
        print *, "the value in your input file was: ", trim(ch)
        print *, "however the allowed options for this field are:"
        do i = 1, ubound(options,1)
          print *, "-  ", cut(options(i))
        end do
      end if
      if (present(errmsg)) then
        call file_error(fileid,errmsg)
      else
        call file_error(fileid)
      end if
    end if

  end subroutine read_char


  subroutine read_char_array(r, fileid, options, errmsg, noexit)
    use typy
    character(len=*), dimension(:), intent(out) :: r
    integer, intent(in) :: fileid
    character(len=*), dimension(:), optional :: options
    character(len=*), intent(in), optional :: errmsg
    logical, intent(in), optional :: noexit

    integer :: ierr
    integer(kind=ikind) :: i, j
    logical :: ok
    logical :: terminate = .false.

    call comment(fileid)
    read(unit=fileid, fmt=*, iostat=ierr) r

    if (ierr /= 0) then
      terminate = .true.
      if (present(noexit)) then
        if (noexit) then
          terminate = .false.
        else
          terminate = .true.
        end if
      end if
    end if

    if (terminate) then
      if (present(errmsg)) then
        call file_error(fileid,errmsg)
      else
        call file_error(fileid)
      end if
    end if

    if (present(options)) then
      do i = 1, ubound(r,1)
        ok = .false.
        optcheck: do j = 1, ubound(options,1)
          if (r(i) == options(j)) then
            ok = .true.
            exit optcheck
          end if
        end do optcheck
        if (.not. ok) then
          call file_error(fileid, errmsg)
        end if
      end do
    end if

  end subroutine read_char_array


  subroutine read_logical(l, fileid, errmsg, noexit, defaultval, success)
    use globals
    logical, intent(out) :: l
    integer, intent(in) :: fileid
    character(len=*), intent(in), optional :: errmsg
    logical, intent(in), optional :: noexit
    logical, intent(in), optional :: defaultval
    logical, intent(out), optional :: success

    integer :: ierr
    character(len=1) :: ch
    logical :: terminate = .false.

    call comment(fileid)
    read(unit=fileid, fmt=*, iostat=ierr) ch

    if (ierr == 0) then
      if (present(success)) success = .true.
      select case(ch)
        case("y")
          l = .true.
        case("n")
          l = .false.
        case default
          ierr = -1
          write(unit=terminal, fmt=*) " " //achar(27)//'[43m', "You have defined incorrect value for logical parameter" &
             //achar(27)//'[0m'
      end select
    end if

    if  (ierr /= 0 ) then
      if (.not. present(noexit)) then
        terminate = .true.
      else
        if (noexit) then
          terminate = .false.
        else
          terminate = .true.
        end if
      end if

      if (.not. present(defaultval)) then
        if (terminate) then
          if (present(errmsg)) then
            call file_error(fileid,errmsg)
          else
            call file_error(fileid)
          end if
        end if
      else
        l = defaultval
        if (present(success)) success = .false.
      end if
    end if

  end subroutine read_logical


  subroutine read_logical_array(yes, fileid,  errmsg, noexit)
    use typy
    use globals
    logical, dimension(:), intent(out) :: yes
    integer, intent(in) :: fileid
    character(len=*), intent(in), optional :: errmsg
    logical, intent(in), optional :: noexit

    integer(kind=ikind) :: i
    character(len=1), dimension(:), allocatable :: tmpdata
    logical :: noexit_loc

    if (present(noexit)) then
      noexit_loc = noexit
    else
      noexit_loc = .false.
    end if

    allocate(tmpdata(ubound(yes,1)))
    call read_char_array(tmpdata, fileid, NOEXIT=.true.)

    do i = 1, ubound(tmpdata,1)
      if (tmpdata(i) /= "y" .and. tmpdata(i) /= "n") then
        if (.not. noexit_loc) then
          if (present(errmsg)) then
            call file_error(fileid, errmsg)
          else
            call file_error(fileid, errmsg="Incorrect inputs, have you set all required [y/n] values?")
          end if
        end if
      else
        if (tmpdata(i) == "y") yes(i)=.true.
        if (tmpdata(i) == "n") yes(i)=.false.
      end if
    end do

  end subroutine read_logical_array


  subroutine read_sep(fileid)
    integer, intent(in) :: fileid
    character(len=3) :: separator

    call fileread(separator, fileid, options=(/"---"/), errmsg="Missing block separator. Check your inputs.")

  end subroutine read_sep


  subroutine file_error(iunit, errmsg)
    use core_tools
    use globals

    integer, intent(in) :: iunit
    character(len=*), intent(in), optional :: errmsg

    character(len=512) :: filename
    character(len=256) :: iaction
    integer :: line, i_err, err_read
    integer(kind=likind) :: i1, i2
    character(len=4096) :: string

    string = "last value"
    i1 = -100
    i2 = -100
    inquire(unit=iunit, name=filename, action=iaction)
    line = 1
    backspace(unit=iunit, iostat = i_err)

    read(unit=iunit, fmt=*, iostat = err_read) string

    do
      i1 = ftell(iunit)
      backspace (unit=iunit, iostat=i_err)
      line = line + 1
      i2 = ftell(iunit)
      if (i1 == i2) then
        exit
      end if
    end do

    if (line <= 2) then
      print *, "-------------------------------------------------------------------"
      print *, "WARNING: the line number was identified incorrectly, "
      print *, "likely it is a compiler error in ftell function, "
      print *, "find the input error on your own. :("
    end if

    print *, "-------------------------------------------------------------------"
    print *, "  "
    print *, " " //achar(27)//'[43m', "!!ERROR!!" //achar(27)//'[0m', &
             "  in reading from file: ", trim(filename), " near line: ", line-1, &
             " opened with status: ", trim(iaction),  " -> exiting DRUtES"
    print *, "  "

    if (err_read == 0) write(unit=terminal, fmt=*) "the value you have typed is: ", trim(string)
    write(unit=terminal, fmt=*) "--------------------------------"

    if (present(errmsg)) then
      write(unit=terminal, fmt=*) " " //achar(27)//'[43m', cut(errmsg) //achar(27)//'[0m'
    end if

    stop

  end subroutine file_error


  subroutine comment(Unit, mark)
    integer, intent(in) :: Unit
    character(len=1), intent(in), optional :: mark
    character(len=1) :: symbol
    character(len=1) :: String
    integer :: i_err

    if (present(mark)) then
      symbol = mark
    else
      symbol = "#"
    end if

    do
      read(unit=Unit,fmt = *, iostat = i_err ) String
      if (i_err /= 0) then
        return
      end if
      if (String == symbol) then
        cycle
      else
        backspace Unit
        return
      endif
    end do

  end subroutine comment


  subroutine reverse_comment(Unit,mark)
    use debug_tools
    integer, intent(in) :: Unit
    character(len=1), intent(in), optional :: mark
    character(len=1) :: symbol
    character(len=1) :: String
    integer :: i_err

    if (present(mark)) then
      symbol = mark
    else
      symbol = "#"
    end if

    do
      read(unit=Unit,fmt = *, iostat = i_err ) String
      print *, String , i_err
      print *, ftell(Unit) ; call wait()

      if (i_err /= 0) then
        return
      end if

      if (String == symbol) then
        print *, "vifu", String, symbol
        backspace Unit
        backspace Unit
      else
        backspace Unit
        return
      end if
    end do

  end subroutine reverse_comment

end module readtools
