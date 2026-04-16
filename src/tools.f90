module tools
  use typy
  use globals
  implicit none
contains

 

  ! -----------------------------------------------------
  ! Find neighbours: elements that share an edge
  ! -----------------------------------------------------
  subroutine find_neighbours(elems, nodes_loc)
    type(elements_str), intent(inout) :: elems
    type(nodes_str),    intent(in)    :: nodes_loc

    integer(kind=ikind) :: nel
    integer(kind=ikind) :: i, j, k, shared
    integer(kind=ikind) :: e1(3), e2(3)

    nel = elems%kolik
    if (nel <= 0_ikind) return

    if (.not. allocated(elems%neighbours)) then
      allocate(elems%neighbours(nel,3))
    end if
    elems%neighbours = 0_ikind

    do i = 1_ikind, nel
      e1 = elems%data(i,:)
      do j = 1_ikind, nel
        if (i == j) cycle
        e2 = elems%data(j,:)

        shared = 0_ikind
        do k = 1_ikind, 3_ikind
          if (any(e2 == e1(k))) shared = shared + 1_ikind
        end do

        if (shared >= 2_ikind) then
          do k = 1_ikind, 3_ikind
            if (elems%neighbours(i,k) == 0_ikind) then
              elems%neighbours(i,k) = j
              exit
            end if
          end do
        end if
      end do
    end do
  end subroutine find_neighbours


  function avg(r) result(avgval)
    real(kind=rkind), dimension(:), intent(in) :: r
    real(kind=rkind) :: avgval
    integer(kind=ikind) :: i

    avgval = 0.0_rkind
    do i = 1, size(r)
      avgval = avgval + r(i)
    end do
    avgval = avgval / real(size(r), kind=rkind)
  end function avg





  !==============================================================
  !   COMPUTE ELEMENT AREAS
  !==============================================================
  subroutine compute_areas()
    integer(kind=ikind) :: i
    integer(kind=ikind) :: n1, n2, n3
    real(kind=rkind) :: a(2), b(2), c(2)
    real(kind=rkind) :: total_area

    total_area = 0.0_rkind

    do i = 1_ikind, elements%kolik
      n1 = elements%data(i,1)
      n2 = elements%data(i,2)
      n3 = elements%data(i,3)
      a  = nodes%data(n1,:)
      b  = nodes%data(n2,:)
      c  = nodes%data(n3,:)
      elements%area(i) = area_triangle(a,b,c)
      total_area       = total_area + elements%area(i)
    end do

    print *, " Total mesh area = ", total_area
  end subroutine compute_areas


  pure function area_triangle(a,b,c) result(area)
    real(kind=rkind), intent(in) :: a(2), b(2), c(2)
    real(kind=rkind) :: area
    area = 0.5_rkind * abs(a(1)*(b(2)-c(2)) + b(1)*(c(2)-a(2)) + c(1)*(a(2)-b(2)))
  end function area_triangle


  !==============================================================
  !   COMPUTE average altitude per element
  !==============================================================
  subroutine compute_avgalt()
    integer(kind=ikind) :: i, n1, n2, n3
    real(kind=rkind) :: z1, z2, z3

    do i = 1_ikind, elements%kolik
      n1 = elements%data(i,1)
      n2 = elements%data(i,2)
      n3 = elements%data(i,3)

      if (allocated(nodes%altitude)) then
        z1 = nodes%altitude(n1)
        z2 = nodes%altitude(n2)
        z3 = nodes%altitude(n3)
      else
        z1 = nodes%data(n1,2)
        z2 = nodes%data(n2,2)
        z3 = nodes%data(n3,2)
      end if

      elements%avgalt(i) = (z1 + z2 + z3) / 3.0_rkind
    end do
  end subroutine compute_avgalt


 




  function setlines(el, nd) result(pts)
    integer(kind=ikind), intent(in) :: el, nd
    integer(kind=ikind), dimension(2) :: pts
    integer(kind=ikind), dimension(3) :: neighnds, elnds
    integer(kind=ikind) :: pos, i, j

    elnds    = elements%data(el,:)
    neighnds = elements%data(nd,:)

    pts = 0_ikind
    pos = 0

    do i = 1, 3
      do j = 1, 3
        if (neighnds(i) == elnds(j)) then
          pos = pos + 1
          if (pos <= 2) then
            pts(pos) = neighnds(i)
          else
            print *, "Error in setlines: two elements share all three nodes."
            error stop
          end if
        end if
      end do
    end do
  end function setlines


  function matchlines(l1, l2) result(true)
    integer(kind=ikind), dimension(2), intent(in) :: l1, l2
    logical :: true

    true = .false.
    if (l1(1) == l2(1) .and. l1(2) == l2(2)) then
      true = .true.
      return
    end if

    if (l1(1) == l2(2) .and. l1(2) == l2(1)) then
      true = .true.
      return
    end if
  end function matchlines


  function getcenter(el) result(xy)
    integer(kind=ikind), intent(in) :: el
    real(kind=rkind), dimension(2) :: xy

    xy(1) = avg(nodes%data(elements%data(el,:),1))
    xy(2) = avg(nodes%data(elements%data(el,:),2))
  end function getcenter


  function dist(A,B) result(l)
    real(kind=rkind), dimension(:), intent(in) :: A,B
    real(kind=rkind) :: l
    l = sqrt((A(1)-B(1))*(A(1)-B(1)) + (A(2)-B(2))*(A(2)-B(2)))
  end function dist





 
end module tools
