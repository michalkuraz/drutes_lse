! geom_tools.f90
module geom_tools
  use typy
  use globals
  implicit none
contains

  ! -----------------------------------------------------
  ! Find neighbours: elements that share an edge
  ! -----------------------------------------------------
  subroutine find_neighbours(elems, nodes_loc)
    type(elements_str), intent(inout) :: elems
    type(nodes_str),    intent(in)    :: nodes_loc   ! not used now, but fine

    integer(kind=ikind) :: nel
    integer(kind=ikind) :: i, j, k, shared
    integer(kind=ikind) :: e1(3), e2(3)

    ! Number of elements in the mesh
    nel = elems%kolik
    if (nel <= 0_ikind) return

    ! Allocate neighbour table if needed: each triangle can have up to 3 neighbours
    if (.not. allocated(elems%neighbours)) then
       allocate(elems%neighbours(nel,3))
    end if
    elems%neighbours = 0_ikind   ! 0 = no neighbour

    ! Double loop over all element pairs (i, j)
    do i = 1_ikind, nel
       e1 = elems%data(i,:)          ! nodes of element i

       do j = 1_ikind, nel
          if (i == j) cycle          ! skip self
          e2 = elems%data(j,:)       ! nodes of element j

          ! Count how many nodes are shared between element i and j
          shared = 0_ikind
          do k = 1_ikind, 3_ikind
             if (any(e2 == e1(k))) shared = shared + 1_ikind
          end do

          ! If they share at least 2 nodes, they share an edge => neighbours
          if (shared >= 2_ikind) then
             ! Put j into the next free neighbour slot of element i
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

end module geom_tools
