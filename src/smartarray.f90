module smartarray
  use typy
  !> smart array, you don't need to allocate -- usefull if you don't know how many data you will write in, 
  !! it is simply smart enough, it does all allocations automatically. By default the array data allocates on dimension 2, if the data 
  !! doesn't fit into this array, then the array is reallocated for a double of its original size
  !<
  type, public :: smartarray_int
    !> array with data
    integer(kind=ikind), dimension(:), allocatable :: data
    !> put there what ever you like, sometimes it can be usefull if you have this array
    integer(kind=ikind), dimension(:), allocatable :: info
    !> current dimension of your data, e.g. your smart array is declared as smart, then don't try ubound(smart%%data) NEVER!!
    !! for the array dimension see the value smart%pos
    !<
    integer(kind=ikind) :: pos=0, infopos=0
    contains
      !> fills value into array, its dimension is automaticaly allocated
      procedure :: fill => ismartfill
      !> clear the smart array, input parameter is logical FULL, this parameter is optional
      !! by default if FULL not passed or FULL is .false., then only the value smart%pos is set to zero
      !! if FULL .true., the also the array data is deallocated
      !! it turns out that FULL .false. is faster, but less memory efficient
      !<
      procedure :: clear => ismartclear
      !> disjoint fill -- it fills value into array, only if the same value does not already exist in the array
      procedure :: nrfill => ismartfill_norepeat
      !> return value is logical, if the value already exist in the array, then the reutrn value is .true. otherwise .false.
      procedure :: exist => ismartexist
  end type smartarray_int
  
  
  !> For comments and description see smartarray_int. 
  type, public ::  smartarray_real
    real(kind=rkind), dimension(:), allocatable :: data
    integer(kind=ikind), dimension(:), allocatable :: info
    integer(kind=ikind) :: pos=0, infopos=0
    contains
      !> fills value into array, its dimension is automaticaly allocated
      procedure :: fill => rsmartfill
      procedure :: clear => rsmartclear
      !> disjoint fill -- it fills value into array, only if the same value does not exist in the array. Since we deal with real numbers in this structure, then the real number is checked for absolute value of differece between two real numbers
      !! the value is checked as 
      !!  abs(array%data(i) - input) < max(abs(array%data(i)), abs(input))*(epsilon(input)
      !<
      procedure :: nrfill => rsmartfill_norepeat
  end type smartarray_real  
  
    private :: ismartfill, ismartclear, ismartfill_norepeat, rsmartfill, rsmartclear, rsmartfill_norepeat, ismartexist
  
  contains
    !> writes into smartarray integer vectors
    subroutine ismartfill(array,input, info)
      use typy
      class(smartarray_int), intent(in out) :: array
      integer(kind=ikind), intent(in) :: input
      integer(kind=ikind), intent(in), optional :: info
      
      integer(kind=ikind) :: l
      integer(kind=ikind), dimension(:), allocatable :: itmp
      integer(kind=ikind), dimension(:), allocatable  :: logtmp
      
      ! check allocations
      if (.not. allocated(array%data)) then
        array%pos = 0
        allocate(array%data(1))
        if (present(info)) then
          array%infopos = 0
          allocate(array%info(1))
        end if
      end if
      
      array%pos = array%pos+1
      
      if (present(info)) then
        array%infopos = array%infopos + 1
        if (array%pos /= array%infopos) then
            print *, "this is a bug in the code"
            print *, "smartarray%data and smartarray%info has different amount of data written"
            print *, "called from global_objs::rsmartfill"
            print *, "contact Michal -> michalkuraz@gmail.com"
            ERROR STOP
        end if
      end if
      
      if (ubound(array%data,1) < array%pos) then
        l = ubound(array%data,1)
        allocate(itmp(l))
        itmp = array%data
        deallocate(array%data)
        allocate(array%data(2*l))
        array%data(1:l) = itmp
        deallocate(itmp)
        if (present(info)) then
          allocate (logtmp(l))
          logtmp = array%info
          deallocate(array%info)
          allocate(array%info(2*l))
          array%info(1:l) = logtmp
          deallocate(logtmp)
        end if	  
      end if
      
      array%data(array%pos) = input
      
      if (present(info)) then
        array%info = info
      end if
      
      
    end subroutine ismartfill
   
    
    
    !> writes into smartarray integer vectors, only if the value is not present in the smartvector
    subroutine ismartfill_norepeat(array, input, info)
      use typy
      class(smartarray_int), intent(in out) :: array
      integer(kind=ikind), intent(in) :: input
      integer(kind=ikind), intent(in), optional :: info
      
      integer(kind=ikind) :: i
!       logical :: exist
!       
!       exist = .false.
!       if (allocated(array%data)) then
!         do i=1, array%pos
!           if (array%data(i) == input) then
!             exist = .true.
!             EXIT
!           end if
!         end do
!        end if

      
      if (.not. ismartexist(array,input)) then
        if (present(info)) then
          call ismartfill(array, input, info)
        else
          call ismartfill(array, input)
        end if
      end if
    
    
    end subroutine ismartfill_norepeat
    
    !> clears (deallocates) smart vector
    subroutine ismartclear(array, full)
      class(smartarray_int), intent(in out) :: array
      !> if present and .true. the smart vector is completely deallocated
      logical, intent(in), optional :: full
      
      if (present(full)) then
         if (full .and. allocated(array%data)) then
           deallocate(array%data)
         end if
      end if

      array%pos = 0
      
    end subroutine ismartclear
    
    !> checks if the value is present in the vector
    function ismartexist(array, value) result(exist)
      use typy
      class(smartarray_int), intent(in) :: array
      integer(kind=ikind), intent(in) :: value
      logical :: exist
      
      integer(kind=ikind) :: i
      
      if (minval (abs(array%data(1:array%pos) - value)) == 0) then
        exist = .true.
      else
        exist = .false.
      end if
      
      
    end function ismartexist
    
    !> write real into smartarray vector
    subroutine rsmartfill(array,input, info)
      use typy
      class(smartarray_real), intent(in out) :: array
      real(kind=rkind), intent(in) :: input
      integer(kind=ikind), intent(in), optional :: info
      
      integer(kind=ikind) :: l
      real(kind=rkind), dimension(:), allocatable :: rtmp
      integer(kind=ikind), dimension(:), allocatable  :: logtmp
      
      ! check allocations
      if (.not. allocated(array%data)) then
        array%pos = 0
        allocate(array%data(1))
        if (present(info)) then
          array%infopos = 0
          allocate(array%info(1))
        end if
      end if
      
      array%pos = array%pos+1
      
      if (present(info)) then
        array%infopos = array%infopos + 1
        if (array%pos /= array%infopos) then
            print *, "this is a bug in the code"
            print *, "smartarray%data and smartarray%info has different amount of data written"
            print *, "called from global_objs::rsmartfill"
            print *, "contact Michal -> michalkuraz@gmail.com"
            ERROR STOP
        end if
      end if
      
      if (ubound(array%data,1) < array%pos) then
        l = ubound(array%data,1)
        allocate(rtmp(l))
        rtmp = array%data
        deallocate(array%data)
        allocate(array%data(2*l))
        array%data(1:l) = rtmp
        deallocate(rtmp)
        if (present(info)) then
          allocate (logtmp(l))
          logtmp = array%info
          deallocate(array%info)
          allocate(array%info(2*l))
          array%info(1:l) = logtmp
          deallocate(logtmp)
        end if	  
      end if
      
      array%data(array%pos) = input
      
      if (present(info)) then
        array%info = info
      end if
      
      
    end subroutine rsmartfill
   
    
    
    !> writes real into smartarray vector, only if value sufficiently differ from the values already present. Analogical to smartfill_norepeat for integers
    subroutine rsmartfill_norepeat(array, input, info)
      use typy
      class(smartarray_real), intent(in out) :: array
      real(kind=rkind), intent(in) :: input
      integer(kind=ikind), intent(in), optional :: info
      
      integer(kind=ikind) :: i
      logical :: exist
      
      exist = .false.
      if (allocated(array%data)) then
        do i=1, array%pos
          if (abs(array%data(i) - input) < max(abs(array%data(i)), abs(input))*(epsilon(input))) then
            exist = .true.
            EXIT
          end if
        end do
       end if
      
      if (.not. exist) then
        if (present(info)) then
          call rsmartfill(array, input, info)
        else
          call rsmartfill(array, input)
        end if
      end if
    
    
    end subroutine rsmartfill_norepeat
    
    !> clears / deallocates real smartarray vector
    subroutine rsmartclear(array, full)
      class(smartarray_real), intent(in out) :: array
      logical, intent(in), optional :: full
      
      if (present(full) ) then
        if (full .and. allocated(array%data)) then
          deallocate(array%data)
        end if
      end if

      array%pos = 0
      
    end subroutine rsmartclear


end module smartarray
