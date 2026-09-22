! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   MiLaDy (Machine Learning Dynamics) 
! HND X   
! HND X
! HND X   Milady was written and designed by Mihai-Cosmin Marinica, Alexandra M. Goryaeva 
! HND X   Copyright 2015-2025.
! HND X
! HND X   Portions of MiLady were written by Wesley UnnToc, Clovis Lapointe, Anruo Zhong, 
! HND X   Jacopo Baima, Anida Khizar, Christian van Wambeke 
! HND X
! HND X   MiLaDy is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   MiLaDy is distributed at https://github.com/ai-atoms/milady in the hope that it 
! HND X   will be useful for non-commercial academic research, but WITHOUT ANY WARRANTY;  
! HND X   without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
! HND X   PURPOSE. See the ASL for more details.
! HND X
! HND X   You should have received a copy of the ASL along with this program
! HND X   (e.g. in a LICENSE.md file); if not, you can write to the original licensors,
! HND X   Mihai-Cosmin Marinica and Alexandra M. Goryaeva. The ASL is also published at
! HND X   http://github.com/ai-atoms/milady/ASL
! HND X
! HND X   When using this software, please cite the following reference:
! HND X
! HND X   A.M. Goryaeva et al. Comput. Mater. Sci. 166, 200-209 (2019)
! HND X   A.M. Goryaeva et al. Phys. Rev. Mater. 5: 103803 (2021)
! HND X
! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "../../MLD_MACROS.INC"



!>------begin_module_graph_cg------------------------------------------------------------
module module_graph_cg
    use mld_logger
    implicit none

    type :: NodeType
        integer :: id
        integer, dimension(:), allocatable :: l_values
        integer, dimension(:,:), allocatable :: L_mat 
        integer :: lini 
        integer :: parent1, parent2
        integer :: child  ! Added to track the child node
        integer :: level  ! Added to track the level of the node (how many agregations) 
        character(len=80) :: name ! Added to track the name of the node
        integer :: Lorder ! Added to track the order of the node in the level ordered list
    end type NodeType

    type :: MatLevelType
        integer, dimension(:,:), allocatable :: mat
    end type MatLevelType

    type :: GraphType
        type(NodeType), allocatable :: nodes(:)
        integer :: nu ! Added to track the order of the graph
        integer,  dimension(:), allocatable :: level_ordered_nodes 
        integer, dimension(:), allocatable :: level0_nodes
        integer, dimension(:,:), allocatable :: LmatLevel1

        type(MatLevelType), dimension(:), allocatable :: mlevel 
        !TODOace to remove 
        integer, allocatable :: edges(:,:)
    contains
        procedure :: add_node
        procedure :: initialize
        procedure :: print_graph
        procedure :: build => build_ralf
        !procedure :: build => build_genevieve
        procedure :: find_children
        procedure :: aggregation
        procedure :: deallocate_graph_cg
    end type GraphType

contains

    subroutine initialize(this, nu)
      class(GraphType), intent(inout) :: this
      integer, intent(in) :: nu
      integer :: num_nodes, num_edges
      ! Compute the number of nodes and edges
   
      
      num_nodes = calculate_num_nodes(nu)
      num_edges = calculate_num_edges(nu)  
      ! Allocate arrays for nodes and edges
      if (allocated(this%nodes)) deallocate(this%nodes) ; allocate(this%nodes(num_nodes))
      !$! allocate(this%edges(2, num_edges))  ! Two entries per edge (from, to)  
      this%nu = nu

    end subroutine initialize

    subroutine print_graph(this)
      use mld_mpi, only : mld_rank
      class(GraphType), intent(in) :: this
      integer :: i
      if (mld_rank == 0 ) then 
        write(*,'("------------- Graph of kinetic momenta coupling: info ----------------------")')
        write(*,'("  index", "  ", "Node id", "  ", "Node Order", "   ", "Par 1", "  ", "Par2", "  ", "Level", "   ", "Name")')        
        do i = 1, size(this%nodes)
          write(*,'(i8, 5i8, "     ", A)') i, this%nodes(i)%id,  this%nodes(i)%Lorder, this%nodes(i)%parent1, &
          this%nodes(i)%parent2,  this%nodes(i)%level,  trim(this%nodes(i)%name)
        end do
        write(*,'("-----------------------------------------------------------------------------")')
      end if 
    end subroutine print_graph

    function calculate_num_edges(nu) result(num_edges)
      integer, intent(in) :: nu
      integer :: num_edges

      !Numerical computation
      if (nu >= 2) then
        num_edges = 2  ! Two edges for nu = 2
        if (nu > 2) then
            num_edges = num_edges + 2 * (nu - 2)  ! Add two edges for each additional nu
        end if
      else
        num_edges = 0
      end if
      !$! ! Output the result from closed formula
      !$! print *, 'Closed formula edge computation: ', num_edges
    end function calculate_num_edges

    function calculate_num_nodes(nu) result(num_nodes)
      integer, intent(in) :: nu
      integer ::  num_nodes 

      num_nodes = 2*nu-1 

      !$! print *, 'Closed formula computation versus num: ', num_nodes
    end function calculate_num_nodes

    subroutine build_genevieve(this, nu)
      ! The build_genevieve subroutine in the module_graph_cg module follows a straightforward 
      ! iterative approach for constructing a graph representing quantum angular momentum coupling 
      ! up to a given order nu. 
      ! This method does not utilize a recursive strategy but builds the graph by sequentially 
      ! adding nodes for each order.
      !
      ! Key Steps:
      ! - Individual momentum nodes (l1, l2, ...) are added for each value of nu.
      ! - Combined momentum nodes (L12, L34, ...) are created for each pair of momenta.
      ! - The subroutine iteratively adds these nodes, ensuring that each combined 
      !    node (e.g., L1234) is linked to its constituent momenta or the previous combined nodes.
      !
      ! The subroutine does not employ recursion but rather a loop that iterates over 
      ! the range of angular momenta, systematically adding and connecting nodes
      class(GraphType), intent(inout) :: this
      integer, intent(in) :: nu
      integer :: i, node_index, edge_index
      character(len=100) :: temp_str
      integer :: initdim, idn, idcurrent, idlevel
      ! Initialize node_index and edge_index
      node_index = nu + 1  ! Start index for combined nodes
      edge_index = 1
      
      this%nodes(:)%level = 0
      ! Initialize for nu = 2
      call add_node(this, 0, 'l1', 1, 0, 0)
      call add_node(this, 0, 'l2', 2, 0, 0)
      call add_node(this, 1, 'L12', 3, 1, 2)
      if (nu == 2) then 
        call find_children(this)
        call build_level_ordered_nodes(this)
        call print_graph(this)
        return
      end if 
    
      initdim = 3 
      idn = 0 
      temp_str = '12'
      do i = 3, nu
        idn = idn + 1
        idcurrent = idn + initdim
        call add_node(this, 0, 'l' // trim(write_str(i)), idcurrent, 0 , 0)
        idn = idn + 1
        idcurrent = idn + initdim
        idlevel = max(this%nodes(idcurrent-1)%level, this%nodes(idcurrent-2)%level) + 1
        temp_str = trim(temp_str) // trim(write_str(i))
        call add_node(this, idlevel, 'L' // temp_str, idcurrent, idcurrent-1, idcurrent-2)
      end do 
      call find_children(this)
      call build_level_ordered_nodes(this)
      call print_graph(this)

    end subroutine build_genevieve

    subroutine build_ralf(this, nu)
        class(GraphType), intent(inout) :: this
        integer, intent(in) :: nu
        integer :: i, idn, idcurrent, initdim, idlevel 
    
        this%nodes(:)%level = 0
        ! Initialize for nu = 2
        call add_node(this, 0, 'l1', 1, 0, 0)
        call add_node(this, 0, 'l2', 2, 0, 0)
        call add_node(this, 1, 'L12', 3, 1, 2)
    
        if (nu == 2) then 
          call find_children(this)
          call build_level_ordered_nodes(this)
          call print_graph(this)
          return
        end if   

        ! Initialize for nu = 3
        call add_node(this, 0, 'l3', 4, 0 , 0)
        call add_node(this, 2, 'L123', 5, 3 , 4)
            
        if (nu == 3) then 
          call find_children(this)
          call build_level_ordered_nodes(this)
          call print_graph(this)
          return
        end if  
    
        !initialize for nu = 4
        call add_node(this, 0, 'l4', 5, 0 , 0)
        call add_node(this, 1, 'L34', 6, 4, 5)
        call add_node(this, 2, 'L1234', 7, 3, 6)
    
        if (nu == 4) then 
          call find_children(this)
          call build_level_ordered_nodes(this)
          call print_graph(this)
          return
        end if
    
        ! From now we are ready for recursion  
    
        initdim = 7 
        ! Handle cases for nu > 2
        idn = 0 
        do i = 5, nu
            if (mod(i, 2) == 0) then
                ! For even nu, add l_i, L(i-1)i, and L123...i
                !no update we crash the previous node 
                !idn = idn + 1
                idcurrent = idn + initdim
                call add_node(this, 0, 'l' // trim(write_str(i)), idcurrent, 0 , 0)
                idn = idn + 1
                idcurrent = idn + initdim
                idlevel = max(this%nodes(idcurrent-1)%level, this%nodes(idcurrent-2)%level) + 1
                call add_node(this, idlevel, 'L' // trim(write_str(i - 1)) // trim(write_str(i)), idcurrent, &
                              idcurrent-1, idcurrent-2)
                idn = idn + 1
                idcurrent = idn + initdim 
                idlevel = max(this%nodes(idcurrent-1)%level, this%nodes(idcurrent-4)%level) + 1          
                call add_node(this, idlevel, 'L' // trim(write_str_range(1, i)), idcurrent, idcurrent-1, idcurrent-4)
            else
                ! For odd nu, add l_i and L123...i
                idn = idn + 1
                idcurrent = idn + initdim
                call add_node(this, 0, 'l' // trim(write_str(i)), idcurrent, 0 , 0)
                idn = idn + 1 
                idcurrent = idn + initdim
                idlevel = max(this%nodes(idcurrent-1)%level, this%nodes(idcurrent-2)%level) + 1
                call add_node(this, idlevel, 'L' // trim(write_str_range(1, i)), idcurrent, idcurrent-1, idcurrent-2)
            end if
        end do

        
        call find_children(this)
        call build_level_ordered_nodes(this)
        call print_graph(this)

    end subroutine build_ralf 

    subroutine add_node(this, level, name, id, parent1, parent2)
        class(GraphType), intent(inout) :: this
        integer, intent(in) :: id, parent1, parent2, level 
        character(len=*), intent(in) :: name
    
        this%nodes(id)%id = id
        this%nodes(id)%name = name
        this%nodes(id)%parent1 = parent1
        this%nodes(id)%parent2 = parent2
        this%nodes(id)%level = level

    end subroutine add_node

    function write_str(num) result(str)
        integer, intent(in) :: num
        character(len=10) :: str
        write(str, '(I0)') num
        str = trim(adjustl(str))
    end function write_str

    function write_str_range(start, end) result(str)
        integer, intent(in) :: start, end
        character(len=100) :: str
        integer :: i
        str = ''
        do i = start, end
            str = trim(str) // trim(write_str(i))
        end do
    end function write_str_range

    subroutine find_children(this)
      !  num_nodes ^ 2 = (2*nu -1)^2 algo 
      class(GraphType), intent(inout) :: this
      integer :: num_nodes, i, parent_id
  
      num_nodes = size(this%nodes)

      do parent_id = 1, num_nodes
          do i = 1, num_nodes
              if (this%nodes(i)%parent1 == parent_id .or. this%nodes(i)%parent2 == parent_id) then
                  this%nodes(parent_id)%child = i
              end if
          end do
      end do
    end subroutine find_children

    subroutine build_adjacency_matrix(this, adj_matrix)
      class(GraphType), intent(in) :: this
      integer, allocatable, intent(out) :: adj_matrix(:,:)
      integer :: num_nodes, i
  
      num_nodes = size(this%nodes)
      allocate(adj_matrix(num_nodes, num_nodes))
      adj_matrix = 0  ! Initialize the matrix with zeros
  
      ! Fill in the adjacency matrix
      do i = 1, num_nodes
          if (this%nodes(i)%parent1 > 0) then
              adj_matrix(this%nodes(i)%parent1, i) = 1
          end if
          if (this%nodes(i)%parent2 > 0) then
              adj_matrix(this%nodes(i)%parent2, i) = 1
          end if
      end do
    end subroutine build_adjacency_matrix

    subroutine build_level_ordered_nodes(this)
    class(GraphType), intent(inout) :: this
    integer :: i, j, vector_index, node_count, level_node_count 


    ! Count the number of nodes with level = 0 and made a list 
    node_count = 0
    do i = 1, size(this%nodes)
        if (this%nodes(i)%level == 0) then
            node_count = node_count + 1
        end if
    end do
    if (allocated(this%level0_nodes)) deallocate(this%level0_nodes) ; allocate(this%level0_nodes(node_count))
    node_count = 0
    do i = 1, size(this%nodes)
        if (this%nodes(i)%level == 0) then
            node_count = node_count + 1
            this%level0_nodes(node_count) = this%nodes(i)%id
            !debug! write(*,*) 'Level 0 node ', i, this%nodes(i)%id, node_count
        end if
    end do
    !end count level 0 .........................................

    ! Count the number of nodes with level >= 1
    node_count = 0
    do i = 1, size(this%nodes)
        if (this%nodes(i)%level >= 1) then
            node_count = node_count + 1
        end if
    end do

    ! Allocate the vector with the appropriate size
    allocate(this%level_ordered_nodes(node_count))
    ! put all Lorder to zero ... 
    this%nodes(:)%Lorder = 0
    ! Populate the vector with node IDs, ordered by level
    vector_index = 1
    do i = 1, maxval(this%nodes(:)%level)
        level_node_count = 0
        do j = 1, size(this%nodes)
            if (this%nodes(j)%level == i) then
                this%level_ordered_nodes(vector_index) = this%nodes(j)%id
                this%nodes(j)%Lorder = vector_index
                vector_index = vector_index + 1
                level_node_count = level_node_count + 1
            end if
        end do
        ! Check for multiple nodes at levels greater than 1 and issue a warning
        if (i > 1 .and. level_node_count > 1) then
            call log_critical("Fatal: More than one node found at level "//vtoa(i)) 
            call log_critical("Fatal: This is not supported by the current implementation, only G. Dusson and R. Drautz topologies")
            stop 'Weird node topology' 
        end if
    
    end do

    !$! do i = 1, size(this%nodes)
    !$!     write(*,*), 'Level ', i, ' has ',  this%nodes(i)%Lorder
    !$! end do 
  end subroutine build_level_ordered_nodes

subroutine aggregation(this, lini, Lval, Lmat)
  use tuple_module, only: generate_combinations, vector_tuple
  class(GraphType), intent(inout) :: this
  integer, intent(in) :: lini(:), Lval 
  integer, dimension(:,:), allocatable, intent(out) :: Lmat
  integer :: maxLevel, iv, j, noLevel1, noLevel0, icount 
  integer :: cnode, cp1, cp2, l1, l2 , ll, ii, jj, il, tmp_size  
  integer, dimension(:,:), allocatable :: tmp_vec 
  type(vector_tuple), dimension(:), allocatable :: tmpnodes
  !integer, dimension(:), allocatable :: pvec
  integer :: pM_mat, pN_mat, lev ,Mline, Ncols  
  integer :: status 
  
  maxLevel = maxval(this%nodes(:)%level)
  !level 1 is done separatelly because in Ralf way at the level 1 can be more than one node
  noLevel1 = 0 
  noLevel0 = 0
  do j = 1, size(this%nodes)
    if (this%nodes(j)%level == 1) then  
      noLevel1 = noLevel1 + 1
    end if

    if (this%nodes(j)%level == 0) then  
      noLevel0 = noLevel0 + 1
      this%nodes(j)%lini = lini(noLevel0)
    end if
  end do

  if (noLevel0 /= this%nu ) then
     call log_critical("Fatal: Nodes at level 0 found different from nu order"//vtoa(noLevel0)//" "//vtoa(this%nu))
     stop 'Level0: Weird node topology in agregation'  
  end if 

  if (noLevel0 /= size(lini) ) then
     call log_critical("Fatal: Nodes at level 0 found different from lini size"//vtoa(noLevel0)//"  "//vtoa(size(lini))) 
     stop 'Level0: lini not right dimension'  
  end if 

  if (noLevel1 == 0) then
    call log_critical("Fatal: Zero nodes found at level  1") 
    call log_critical("Fatal: This is not supported by the current implementation, only G. Dusson and R. Drautz topologies")
    stop 'Level1: Weird node topology in agregation' 
  end if

  noLevel0 = 0
  do j = 1, size(this%nodes)
      if (this%nodes(j)%level == 0) then
        noLevel0 = noLevel0 + 1  
        !$! print '("id lini values :", I12, " ", *(I5, " " ))', j, lini(noLevel0)
      end if 
  end do


  !Level1 -------------begin--------------------------
  ! Level 1 is done separatelly because in Ralf way at the level 1 can be more than one node. Then will be 
  ! iterative process. So need of level 1 for a general case that works for Genevieve and Ralf way.  

  if (noLevel1 == 0 ) then
     call log_critical("This can be Fatal: Level1 Zero nodes found at level 1") 
  else 
     if (allocated(tmpnodes)) deallocate(tmpnodes) ;  allocate(tmpnodes(noLevel1), stat=status)
     if (status /= 0) then 
        call log_critical("This is Fatal: Allocation of tmpnodes failed") 
        stop 'Allocation failed'
     end if
  end if 

  do  iv = 1, noLevel1
    cnode = this%level_ordered_nodes(iv)
    cp1 = this%nodes(cnode)%parent1
    cp2 = this%nodes(cnode)%parent2
    l1=this%nodes(cp1)%lini
    l2=this%nodes(cp2)%lini
    allocate(tmpnodes(iv)%values(1 + l1 + l2 - abs(l1-l2) ))
    icount = 0
    do ll = abs(l1-l2),  l1 + l2
      icount = icount + 1  
      tmpnodes(iv)%values(icount) = ll
    end do
  end do

  call generate_combinations(tmpnodes, tmp_vec)

  !$! !Print all combinations
  !$! write(*, '("      combinations        ")')
  !$! allocate(pvec(size(tmp_vec,2)))
  !$! do ii = 1, size(tmp_vec, 1)
  !$!      pvec(:) = tmp_vec(ii, :)
  !$!      write(*,*) "Combinations: ", pvec(:)
  !$! end do
  if (allocated(this%LmatLevel1)) deallocate(this%LmatLevel1)
  allocate(this%LmatLevel1(size(tmp_vec,1), 2*this%nu-1))
  do iv = 1, noLevel1
    cnode = this%level_ordered_nodes(iv) 
    this%LmatLevel1(:,cnode) = tmp_vec(:,iv)
  end do

  noLevel0 = 0
  do j = 1, size(this%nodes)
      if (this%nodes(j)%level == 0) then
        noLevel0 = noLevel0 + 1  
        this%LmatLevel1(:,j) = lini(noLevel0) 
      end if 
  end do
  !$! !Print all combinations
  !$! deallocate(pvec) 
  !$! allocate(pvec(2*this%nu-1))
  !$! do ii = 1, size(tmp_vec, 1)
  !$!     pvec(:) = this%LmatLevel1(ii, :)
  !$!     print '("LmatLevel1 a ", I12, ": ", *(I5, ", "))', ii, pvec(:) 
  !$! end do
  !Level1 -------------end-----------------------------
  
  !$!write(*,*) 'maxLevel a', maxLevel, size(this%LmatLevel1,1), size(this%LmatLevel1,2)
  if (allocated(this%mlevel)) deallocate(this%mlevel) ;  allocate(this%mlevel(maxLevel))
  allocate(this%mlevel(1)%mat(size(this%LmatLevel1,1), 2*this%nu-1))
  this%mlevel(1)%mat(:,:) = this%LmatLevel1(:,:)

  !$! if (maxLevel == 1 )  then 
  !$!    deallocate(tmpnodes)
  !$!    deallocate(this%LmatLevel1)
  !$!    return
  !$! end if
  ! Starting recurence for iv = 2   
  if (maxLevel > 1 ) then 
  do iv = noLevel1 + 1,  size(this%level_ordered_nodes)
    cnode = this%level_ordered_nodes(iv)
    lev = this%nodes(cnode)%level 
    pM_mat = size(this%mlevel(lev-1)%mat, 1)
    pN_mat = size(this%mlevel(lev-1)%mat, 2)
 
    cp1 = this%nodes(cnode)%parent1
    cp2 = this%nodes(cnode)%parent2
    if (allocated(tmpnodes)) deallocate(tmpnodes)  ; allocate(tmpnodes(pM_mat))

    tmp_size = 0 
    do ii = 1, pM_mat !  size(mlevel(iv-1)%mat, 1)
      l1=this%mlevel(lev-1)%mat(ii,cp1)
      l2=this%mlevel(lev-1)%mat(ii,cp2)
      allocate(tmpnodes(ii)%values(1 + l1 + l2 - abs(l1-l2)))
      tmp_size = tmp_size + l1 + l2 - abs(l1-l2)
      icount = 0
      do ll = abs(l1-l2),  l1 + l2
        icount = icount + 1  
        tmpnodes(ii)%values(icount) = ll
      end do
    end do 

    if (allocated(this%mlevel(lev)%mat)) deallocate(this%mlevel(lev)%mat) ;
    allocate(this%mlevel(lev)%mat(tmp_size + pM_mat, pN_mat))
    il = 0 
    do ii = 1,  pM_mat !  size(mlevel(iv-1)%mat, 1)
      il = il + 1
      this%mlevel(lev)%mat(il, :) = this%mlevel(lev-1)%mat(ii, :)
      this%mlevel(lev)%mat(il, cnode) = tmpnodes(ii)%values(1)
      if (size(tmpnodes(ii)%values)>1) then 
        do jj = 2, size(tmpnodes(ii)%values)
          il = il + 1
          this%mlevel(lev)%mat(il, :) = this%mlevel(lev-1)%mat(ii, :)
          this%mlevel(lev)%mat(il, cnode) = tmpnodes(ii)%values(jj)
        end do
      end if
      !write(*,*) 'il', il, ii, cnode 
    end do

    ! some cleaning ... 
    deallocate(tmpnodes)
    deallocate(this%mlevel(lev-1)%mat)

  end do  
  else  
    ! here  maxLevel ==  1 
    deallocate(tmpnodes)
    deallocate(this%LmatLevel1) 
  end if  ! if maxLevel > 1 

  !$! do ii = 1, size(this%mlevel(maxLevel)%mat, 1)
  !$!     pvec(:) = this%mlevel(maxLevel)%mat(ii, :)
  !$!     print '("mlevel ", I12, ": ", *(I5, ", "))', ii, pvec(:) 
  !$! end do


  ! filter the final matrix ... 
  Mline = size(this%mlevel(maxLevel)%mat, 1)
  Ncols = size(this%mlevel(maxLevel)%mat, 2)

  icount = 0
  do ii = 1, Mline   
     if (this%mlevel(maxLevel)%mat(ii, Ncols) == Lval) then
        icount = icount + 1
     end if
  end do

  
  if (allocated(Lmat)) deallocate(Lmat) ; allocate(Lmat(icount, Ncols))   

  icount = 0
  do ii = 1, Mline   

     if (this%mlevel(maxLevel)%mat(ii, Ncols) == Lval) then
        icount = icount + 1
        Lmat(icount, :) = this%mlevel(maxLevel)%mat(ii, :)     
     end if
  end do

  !Final clean ... 
  deallocate(this%mlevel(maxLevel)%mat)
  deallocate(this%mlevel)

end subroutine aggregation

subroutine deallocate_graph_cg(this, nu)
    class(GraphType), intent(inout) :: this
    integer, intent(in) :: nu
    integer :: ii 

    
    if  ( this%nu /= nu ) then 
        call log_critical("Fatal: deallocate the L-graph for wrong nu order"//vtoa(this%nu)//" "//vtoa(nu))
        stop 'graph: deallocate Weird node topology in deallocate'  
    end if

    !recap! type :: NodeType
    !recap!     integer :: id
    !recap!     integer, dimension(:), allocatable :: l_values
    !recap!     integer, dimension(:,:), allocatable :: L_mat 
    !recap!     integer :: lini 
    !recap!     integer :: parent1, parent2
    !recap!     integer :: child  ! Added to track the child node
    !recap!     integer :: level  ! Added to track the level of the node (how many agregations) 
    !recap!     character(len=80) :: name ! Added to track the name of the node
    !recap!     integer :: Lorder ! Added to track the order of the node in the level ordered list
    !recap! end type NodeType
    
    do ii = 1, size(this%nodes)
      if (allocated(this%nodes(ii)%l_values)) deallocate(this%nodes(ii)%l_values)
      if (allocated(this%nodes(ii)%L_mat))    deallocate(this%nodes(ii)%L_mat)
    end do
    deallocate(this%nodes)

    !recap! type :: MatLevelType
    !recap!     integer, dimension(:,:), allocatable :: mat
    !recap! end type MatLevelType

    if (allocated(this%mlevel)) then 
      do ii = 1, size(this%mlevel)
        if(allocated(this%mlevel(ii)%mat)) deallocate(this%mlevel(ii)%mat)
      end do
      deallocate(this%mlevel) 
    end if 
    !recap! type :: GraphType
    !recap!     type(NodeType), allocatable :: nodes(:)
    !recap!     integer :: nu ! Added to track the order of the graph
    !recap!     integer,  dimension(:), allocatable :: level_ordered_nodes 
    !recap!     integer, dimension(:,:), allocatable :: LmatLevel1
    !recap!     type(MatLevelType), dimension(:), allocatable :: mlevel 
    !recap!     !TODOace to remove 
    !recap!     integer, allocatable :: edges(:,:)

    if (allocated(this%LmatLevel1)) deallocate(this%LmatLevel1)
    if (allocated(this%level_ordered_nodes))  deallocate(this%level_ordered_nodes)
    if (allocated(this%edges))  deallocate(this%edges)

  end subroutine deallocate_graph_cg       

end module module_graph_cg 
!>------end_module_graph_cg--------------------------------------------------------------


!>------begin_module_l_space_L0-----------------------------------------------------------
module module_l_space_L0
  use mld_logger
  use module_graph_cg, only: GraphType
  implicit none

  type type_base_l0 
    integer, dimension(:,:), allocatable :: Lmat
    integer :: dimr_Lmat, dimc_Lmat 
    integer, dimension(:), allocatable :: lbold
    integer, dimension(:), allocatable :: list_permutations
  end type type_base_l0

  type(type_base_l0), dimension(:), allocatable :: base_L0, tmpbase_L0
  type(GraphType)   :: graph

  contains 

    subroutine dump_classify_l(fname, tuples, map_into_classes, map_into_tuples, unique_classes)
      implicit none
      character(len=*), intent(in) :: fname
      integer, intent(in) :: tuples(:,:)                 ! (num_tuples, nu) only used for nu and num_tuples
      integer, intent(in) :: map_into_classes(:)         ! (num_tuples) 1-based
      integer, intent(in) :: map_into_tuples(:)          ! (num_classes) 1-based
      integer, intent(in) :: unique_classes(:,:)         ! (nu, num_classes)
      integer :: u, nu, num_tuples, num_classes, j
      nu = size(tuples, 2)
      num_tuples = size(tuples, 1)
      num_classes = size(unique_classes, 2)
      open(newunit=u, file=fname, status='replace', action='write', form='formatted')
      write(u,'(I0,1X,I0,1X,I0)') nu, num_tuples, num_classes
      ! map_into_classes
      write(u,"(*(I0,1X))") map_into_classes
      ! map_into_tuples
      write(u,"(*(I0,1X))") map_into_tuples
      ! unique_classes per class (columns)
      do j = 1, num_classes
        write(u,"(*(I0,1X))") unique_classes(:, j)
      end do
      close(u)
    end subroutine dump_classify_l

    subroutine generate_l_space_L0(nu, Lval, l_space_ini, l_space_L0)

        integer, dimension(:,:), intent(in) :: l_space_ini
        integer, dimension(:,:), allocatable, intent(out) :: l_space_L0
        integer, intent(in) :: Lval, nu 
        
        !local variables ... 
        integer, dimension(:,:), allocatable :: tmpLmat
        integer, dimension(:), allocatable :: ltmp 
        integer :: i, dim_l, ii, ic, dim_L0_last
        !character(len=80) :: CHFMT
        !integer, dimension(:), allocatable :: pvec 
        integer, dimension(:), allocatable :: vprint
        integer, dimension(:,:), allocatable :: tmp_l_space_L0

        dim_l = size(l_space_ini,1)
        if (dim_l == 0) then
            call log_critical("No initial l space found") 
            stop 'No initial l space found in generate_l_space_L0'
        end if
        if (allocated(tmp_l_space_L0)) deallocate(tmp_l_space_L0) ; allocate(tmp_l_space_L0(dim_L, nu))
        tmp_l_space_L0 = 0

        call graph%initialize(nu)
        call graph%build(nu)

        if (allocated(ltmp)) deallocate(ltmp) ; allocate(ltmp(graph%nu))
        if (allocated(vprint)) deallocate(vprint) ; allocate(vprint(graph%nu))
        if (allocated(tmpbase_L0)) deallocate(tmpbase_L0) ; allocate(tmpbase_L0(dim_L))
        !$! write (CHFMT, *) '("The l_bold and size_L (", ', graph%nu, '(i4), "  ) ", i4 )'

        ic = 0 
        do i = 1, dim_l
            ltmp(:) = l_space_ini(i, :)
            call graph%aggregation(ltmp, Lval,  tmpLmat)
            !$! write(*,CHFMT) ltmp(:)
            !$! if (allocated(pvec)) deallocate(pvec);  allocate(pvec(size(tmpLmat,2)))
            !$! do ii = 1, size(tmpLmat, 1)
            !$!   pvec(:) = tmpLmat(ii, :)
            !$!   print '("Lmat ", I12, ": ", *(I5, ", "))', ii, pvec(:) 
            !$! end do
            if (size(tmpLmat,1) > 0) then 
              ic = ic + 1 
              tmpbase_L0(ic)%dimc_Lmat = size(tmpLmat, 2)
              tmpbase_L0(ic)%dimr_Lmat = size(tmpLmat, 1)
              if (allocated(tmpbase_L0(ic)%Lmat)) deallocate(tmpbase_L0(ic)%Lmat)
              allocate(tmpbase_L0(ic)%Lmat(size(tmpLmat, 1), size(tmpLmat, 2)))
              tmpbase_L0(ic)%Lmat = tmpLmat

              if (allocated(tmpbase_L0(ic)%lbold)) deallocate(tmpbase_L0(ic)%lbold)
              allocate(tmpbase_L0(ic)%lbold(graph%nu))
              tmpbase_L0(ic)%lbold = ltmp
              tmp_l_space_L0(ic,:) = ltmp(:) 
            end if 
        end do

        dim_L0_last = ic 
        if (allocated(l_space_L0)) deallocate(l_space_L0) ; allocate(l_space_L0(dim_L0_last, nu))
        l_space_L0 = tmp_l_space_L0(1:dim_L0_last, :)
        deallocate(tmp_l_space_L0)

        if (allocated(base_L0)) deallocate(base_L0) ; allocate(base_L0(dim_L0_last))
        do ii = 1, dim_L0_last
             base_L0(ii)%dimc_Lmat = tmpbase_L0(ii)%dimc_Lmat
             base_L0(ii)%dimr_Lmat = tmpbase_L0(ii)%dimr_Lmat
             if (allocated(base_L0(ii)%Lmat)) deallocate(base_L0(ii)%Lmat)
             allocate(base_L0(ii)%Lmat(size(tmpbase_L0(ii)%Lmat, 1), size(tmpbase_L0(ii)%Lmat, 2)))
              base_L0(ii)%Lmat = tmpbase_L0(ii)%Lmat
              if (allocated(base_L0(ii)%lbold)) deallocate(base_L0(ii)%lbold)
              allocate(base_L0(ii)%lbold(graph%nu))
              base_L0(ii)%lbold = tmpbase_L0(ii)%lbold

              deallocate(tmpbase_L0(ii)%Lmat)
              deallocate(tmpbase_L0(ii)%lbold)
        end do 
        deallocate(tmpbase_L0)

    end subroutine generate_l_space_L0

    subroutine find_permutations_for_l(nu, l_space_L0)
        use permutation_module, only: permutation_object
        integer, intent(in) :: nu
        integer, dimension(:,:), intent(in) :: l_space_L0
        integer :: ii 
        !$! character(len=80) :: CHFMT
        integer, dimension(:), allocatable :: list_permutations, another_vector
        type(permutation_object) :: perm_obj

        call perm_obj%initialize(nu)
        allocate(another_vector(nu))

        do ii = 1, size(l_space_L0, 1)
          another_vector = l_space_L0(ii, :)
          call perm_obj%extract_permutations(another_vector, list_permutations)
          if (allocated(base_L0(ii)%list_permutations))  deallocate(base_L0(ii)%list_permutations)
          allocate(base_L0(ii)%list_permutations(size(list_permutations)))
          base_L0(ii)%list_permutations = list_permutations
          !debug  print '("Permutations for lbold ", I12, ": ", *(I5, ", "))', ii, list_permutations(:)
        end do

    end subroutine find_permutations_for_l

    subroutine classify_tuples_for_l(tuples, num_classes, map_into_classes, map_into_tuples, unique_classes)
      !-----------------------------------------------------------------------------------------------
      ! subroutine that classify the lbold up to a permuation contributions into the same class.
      ! E.g. the tuples of l (1,1,0), (0,1,1) and (1,0,1) are in the same class.
      ! tuples is the list of lbold to be classified (num_perm, nu), nu is the dimesnion of tuples. 
      ! num_classes is the number of unique classes found
      ! unique_classes is the unique classes found
      ! class_map is the mapping of each permutation to a class  
      ! Please note the trick with prime numbers in order to compute the hash of permutaions. Priceless!
      !----------------------------------------------------------------------------------------------- 
      implicit none
      integer, intent(in) :: tuples(:, :)
      integer, intent(out) :: num_classes
      !integer, intent(out) :: class_map(num_perm)
      integer, dimension(:,:), allocatable, intent(inout) :: unique_classes  ! unique_classes(dim, num_perm)
      integer, dimension(:), allocatable, intent(inout) :: map_into_classes, map_into_tuples !  hash_values(num_perm)
      integer :: nu, i, j, hash_value, num_tuples 
      integer, dimension(:), allocatable :: sorted_perm !  sorted_perm(dim)
      integer, dimension(:), allocatable :: hash_values  !  hash_values(num_perm)
      integer, dimension(:,:), allocatable :: itmpm 
      integer, dimension(:), allocatable :: itmpv
      !, dimension(:), allocatable :: map_into_perm
      logical :: is_unique
      num_tuples = size(tuples, 1)
      nu = size(tuples, 2)
      if (allocated(unique_classes)) deallocate(unique_classes) ; allocate(unique_classes(nu, num_tuples))
      if (allocated(map_into_classes)) deallocate(map_into_classes) ; allocate(map_into_classes(num_tuples))
      if (allocated(map_into_tuples)) deallocate(map_into_tuples) ; allocate(map_into_tuples(num_tuples))

      allocate(hash_values(num_tuples))
      allocate(sorted_perm(nu))
      num_classes = 0
      ! Iterate over each permutation
      do i = 1, num_tuples
         ! Copy and sort the permutation
         sorted_perm = tuples(i, :)
         call sort_vector(sorted_perm, nu)
         hash_value = calculate_hash(sorted_perm, nu)
         hash_values(i) = hash_value
         !debug write(*,*) 'vect',  i, hash_value
         ! Check if the hash value is already in unique_classes
         is_unique = .true.
         do j = 1, num_classes
            if (hash_value == hash_values(j)) then
               if (all(sorted_perm == unique_classes(:, j))) then
                  is_unique = .false.
                  map_into_classes(i) = j
                  !map_into_perm(j) = i
                  exit
               end if
            end if
         end do
         ! If unique, add to unique_classes
         if (is_unique) then
            num_classes = num_classes + 1
            unique_classes(:, num_classes) = sorted_perm
            hash_values(num_classes) = hash_value
            map_into_classes(i) = num_classes
            map_into_tuples(num_classes) = i 
         end if
      end do

      !debug! do ii = 1, num_perm
      !debug!     write(*,*) 'Permutation ', ii, ' : ', permutations(ii, :)
      !debug!     write(*,*) 'Class ', ii, ' : ', map_into_classes(ii)
      !debug!  end do  
      !debug!  write(*,*) 'Number of unique classes: ', num_classes
      !debug!  do ii = 1, num_classes
      !debug!     write(*,*) 'Class ', ii, ' : ', unique_classes(:, ii)
      !debug!  end do
      !debug!  do ii = 1, num_classes
      !debug!     write(*,*) 'Class ', ii, ' : ', map_into_perm(ii)
      !debug! end do

      ! Resize unique_classes to the actual number of classes
      allocate (itmpm(nu, num_classes))
      itmpm = unique_classes(:, 1:num_classes)
      deallocate(unique_classes)
      allocate(unique_classes(nu, num_classes))
      unique_classes = itmpm
      deallocate(itmpm)

      allocate(itmpv(num_classes))
      itmpv = map_into_tuples(1:num_classes)
      deallocate(map_into_tuples)
      allocate(map_into_tuples(num_classes))
      map_into_tuples = itmpv
      deallocate(itmpv)

      deallocate(sorted_perm)
      deallocate(hash_values)
    end subroutine classify_tuples_for_l

    subroutine classify_tuples_for_nl(lbold, tuples, num_classes, map_into_classes, map_into_tuples, unique_classes)
      !-----------------------------------------------------------------------------------------------
      ! subroutine that classify a set of tuples with respect some lbold permutation.
      ! "with respect" means than from one nbold =(n1 n2 n3)  and an lbold=(l1 l2 l3) we generate
      ! pair (n1 l1, n2 l2, n3 l3) and we classify those combinations into the same class.
      ! 
      !----------------------------------------------------------------------------------------------- 
      implicit none
      integer, dimension(:), intent(in) :: lbold
      integer, dimension(:,:), intent(in) :: tuples   ! num_perm x nu 
      integer, intent(out) :: num_classes
      !integer, intent(out) :: class_map(num_perm)
      integer, dimension(:,:), allocatable, intent(inout) :: unique_classes  ! unique_classes(dim, num_perm)
      integer, dimension(:), allocatable, intent(inout) :: map_into_classes, map_into_tuples !  hash_values(num_perm)
      integer :: nu, i,  num_tuples
      integer, dimension(:,:), allocatable :: nlbold_tuples
      !, dimension(:), allocatable :: map_into_perm
      integer, dimension(:), allocatable :: tmp_nlbold
      integer :: inu 

      num_tuples = size(tuples, 1)
      nu = size(tuples, 2)

      allocate(nlbold_tuples(num_tuples, nu))
      num_classes = 0
      !debug do ii = 1, num_perm
      !debug    write(*,*) 'Permutation ', ii, ' : ', tuples(ii, :) 
      !debug end do 
      !generate nlbold_permutations
      if (nu > 30) then
         call log_critical('ML error: vector nu size is too large for hash nbold calculation'//vtoa(nu))
         stop 
      end if

      allocate(tmp_nlbold(nu))
      !debug! write(*,*) '---------------------------------------------------'
      !debug! write(*,*) 'lbold ', lbold(:) 
      !debug! do i =1, num_perm
      !debug!     write(*,*) 'tuples ', i, ' : ', tuples(i, :)
      !debug! end do 
      do i = 1, num_tuples
         !tmp_nperm = permutations(i, :)
         do inu = 1, nu
            tmp_nlbold(inu) = 30*lbold(inu) + tuples(i, inu)
         end do
         nlbold_tuples(i, :) = tmp_nlbold(:) 
         !debug! write(*,*) 'nlbold_tuples ', i, ' : ', nlbold_tuples(i, :)
      end do

      call classify_tuples_for_l(nlbold_tuples, num_classes, map_into_classes, map_into_tuples, unique_classes)

      !debug! do i = 1, num_classes
      !debug!    write(*,*) 'Class ', i, ' : ', unique_classes(:, i)
      !debug! end do

      !debug! do i = 1, size(map_into_classes)
      !debug!    write(*,*) 'map_into_class ', i, ' : ', map_into_classes(i)
      !debug! end do

      !debug do i = 1, size(map_into_tuples)
      !debug    write(*,*) 'map_into_tuples ', i, ' : ', map_into_tuples(i)
      !debug end do

    end subroutine classify_tuples_for_nl

    integer function calculate_hash(vec, size)
      ! Function to calculate a hash value for a vector using prime multiplication
      implicit none
      integer, intent(in) :: size
      integer, intent(in) :: vec(size)
      integer :: i
      ! all first twenty prime numbers
      !integer, dimension(0:19) :: primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71]
      if (size > 20) then
         call log_critical('ML Error: vector size is too large for hash calculation')
         stop 
      end if
      !old calculate_hash = 1
      !old do i = 1, size
      !old    calculate_hash = calculate_hash * primes(vec(i))
      !old end do
      calculate_hash = 0
      do i = 1, size
         !calculate_hash = calculate_hash  +  vec(i)**3 
         calculate_hash = calculate_hash  +  vec(i)**(i-1) 
      end do

    
    end function calculate_hash

    ! Subroutine to sort a vector using bubble sort
    subroutine sort_vector(vec, size)
      implicit none
      integer, intent(in) :: size
      integer, intent(inout) :: vec(size)
      integer :: i, j, temp
      do i = 1, size - 1
         do j = 1, size - i
            if (vec(j) > vec(j + 1)) then
               temp = vec(j)
               vec(j) = vec(j + 1)
               vec(j + 1) = temp
            end if
         end do
      end do
    end subroutine sort_vector

end module module_l_space_L0
!>------end_module_l_space_L0------------------------------------------------------------


!>------begin_module_base_cnlm-----------------------------------------------------------
module module_base_cnlm
  use iso_fortran_env,  dp => real64
  use mld_logger
  implicit none 
  private 
  public :: type_base_cnlm
  public :: dump_dico_lb0
  public :: type_base_dico_ll
  
  ! ACE chemical radial descriptor types
  integer, parameter :: ACE_CHEM_RADIAL_RALF = 1
  integer, parameter :: ACE_CHEM_RADIAL_BLOCK_HSVD = 2
  integer, parameter :: ACE_CHEM_RADIAL_HSVD = 3
  integer, parameter :: ACE_CHEM_RADIAL_RANDPROJ = 5
  integer, parameter :: ACE_CHEM_RADIAL_CHEMMAP_HSVD = 4
  
  type type_base_dico_mu
    ! will be mubold vector of dimension nu
    integer, dimension(:), allocatable :: mubold
    ! stamp for a given mubold. if mustamp(i) == 0 then i (1 from mumax)  is not in anyposition in mubold(:)  
    integer, dimension(:), allocatable :: mustamp 
    integer, dimension(:), allocatable :: list_permutations 
  end type type_base_dico_mu 
  
  type type_base_dico_nn
    integer, dimension(:), allocatable :: nbold 
    integer, dimension(:), allocatable :: list_permutations
  end type type_base_dico_nn  


  type type_base_dico_ll
    integer, dimension(:), allocatable :: lbold 
    integer :: dimr_Lmat, dimc_Lmat
    ! will be   #Lbold0 vectors x (2*nu-1) 
    integer, dimension(:,:), allocatable :: Lmat 
    integer, dimension(:), allocatable :: list_permutations
    integer :: unique_nbold 
    ! will be #(mbold=0 vectors) x nu
    integer, dimension(:,:), allocatable :: M0mat
    ! will be #(mbold=0 vectors) x #Lbold0 vectors
    real(dp), dimension(:,:), allocatable :: cg
    ! for a given lbold in this matrix are all  nbold and lbold tuples in order to have a unique class 
    !        nbold = storage(1:nu, idx)    
    !        mubold = storage(nu+1:2*nu,idx)  
    !        idx = 1, 2, 3, ... #(nbold x mubold) 
    integer, dimension(:,:), allocatable :: storage_nmu 
    ! this is  #(nbold x mubold) the number of unique classes max(idx) for a given lbold size( storage_nmu, 2)
    integer :: max_storage_nmu 
  end type type_base_dico_ll

  type type_dico_cg
  ! general type of generalized CG  coefficients. 
    ! will be #(mbold=0 vectors) x #Lbold0 vectors
    real(dp), dimension(:,:), allocatable :: cg  
    integer :: dimr_cg, dimc_cg
  end type type_dico_cg 
  

  type type_base_dico_nnll
     integer :: idx_nb
     integer :: idx_lb0
  end type type_base_dico_nnll

  type type_base_dico_mu_nnll
    integer :: idx_mub
    integer :: idx_nb
    integer :: idx_lb0
    ! will be #(mbold=0 vectors) x #SVD(Lbold0) vectors
     real(dp), dimension(:,:), allocatable :: rpicg
    integer :: idx_cg  ! this map into dicoCG(idx_cg)%cg(:,:) 
  end type type_base_dico_mu_nnll  

  type type_base_dico_mu_nnll_LLi
    integer :: idx_nb       ! this map into dico_nb 
    integer :: idx_mub      ! this map into dico_mu
    integer :: idx_lb0      ! this map into dico_lb0
    integer :: idx_munblb0  ! this map into dico_mu_nnll 
    integer :: idx_Li       ! this map into dico_mub_nb_lb0%Lmat(i,:) or dico_lb0%cg(:,i)  or dico_mu_nnll%rpicg(:,i)
  end type type_base_dico_mu_nnll_LLi


  type type_base_cnlm 
    integer ::  nu 
    type(type_base_dico_mu),   dimension(:),  allocatable:: dico_mu
    type(type_base_dico_nn),   dimension(:),  allocatable:: dico_nb  
    type(type_base_dico_ll),   dimension(:),  allocatable:: dico_lb0 
    type(type_base_dico_mu_nnll), dimension(:),  allocatable :: dico_mub_nb_lb0
    type(type_base_dico_mu_nnll_LLi), dimension(:), allocatable :: dicoB ! from dico for B full  
    type(type_dico_cg), dimension(:), allocatable  :: dicoCG 
    ! number of nbold unique classes after lbold was imposed. 
    integer :: nbold_unique 
    integer, dimension(:,:), allocatable :: uniqueA_tuples
    integer, dimension(:),  allocatable  :: map_large_to_uniqueA 

    integer :: llmax 
    real(dp), dimension(:,:,:,:,:,:), allocatable :: precg 
    contains 
    procedure :: init => base_cnlm_init
    procedure :: build_Bfull => base_cnlm_build_mubnblbLbi
    procedure :: build_precg
    procedure :: base_cnlm_compute_ri_gencg
    procedure :: base_cnlm_compute_rpi_gencg
    procedure :: check_uniqueA_tuples => check_uniqueA_tuples
    procedure :: dump_dico_cg
    procedure :: dump_dico_lb0_cg
    procedure :: dump_dicoB
    procedure :: dump_uniqueA
    procedure :: compute_cg 
  end type type_base_cnlm

  type(type_base_cnlm), dimension(:), allocatable :: base_cnlm
  public :: base_cnlm

  contains 

    subroutine build_precg(this, lmax) 
      use ml_in_ndm_module, only: cg_vector
      class(type_base_cnlm), intent(inout) :: this
      integer, intent(in) :: lmax
      integer :: j0_max, j1_max, j2_max
     
     
      call compute_cg_vector(dble(lmax), 1)
      !debug call dump_cg_vector_file('cg_vector_dump.dat', 1, cg_vector)
      !debug stop 'stop after cg_vector dump'
      if (allocated(this%precg)) deallocate(this%precg)
      j0_max = int(lmax)
      j1_max = int(lmax)
      j2_max = int(lmax)
      allocate (this%precg(0:j1_max, -j1_max:j1_max, 0:j2_max, -j2_max:j2_max, 0:j0_max, -j0_max:j0_max))
      this%precg = cg_vector
      deallocate(cg_vector)

    end subroutine build_precg



    subroutine compute_cg(this, ace_gencg, nu)
      ! extra subroutine to compute the generalized CG in draft-RI version or Dusson-Ortner-RPI
      ! the outcome is agnostic dicoCG which erase the history RI or Dusson-Ortner-RPI  
      use module_ace_desc, only: GENCG_DRAFT, GENCG_DUSORT
      implicit none 
      class(type_base_cnlm), intent(inout) :: this
      integer, intent(in) :: nu, ace_gencg 
      integer :: ii 

      call this%base_cnlm_compute_ri_gencg(nu)
      
      ! Debug dump of RI CGs for Body 4
      !if (nu == 4) then
         !call dump_dico_lb0_cg(this, 'dico_lb0_cg_fortran_4.txt')
      !end if

      ! c. init the vectors of RPI generalizedCG coefficients 
      if (ace_gencg == GENCG_DUSORT) call this%base_cnlm_compute_rpi_gencg(nu)

      if (allocated(this%dicoCG)) deallocate(this%dicoCG)
      if (ace_gencg == GENCG_DRAFT)  allocate(this%dicoCG(size(this%dico_lb0)))
      if (ace_gencg == GENCG_DUSORT) allocate(this%dicoCG(size(this%dico_mub_nb_lb0))) 

      do ii = 1, size(this%dicoCG)
        if (ace_gencg == GENCG_DRAFT)  then
          this%dicoCG(ii)%dimr_cg = size(this%dico_lb0(ii)%cg, 1)
          this%dicoCG(ii)%dimc_cg = size(this%dico_lb0(ii)%cg, 2)
          allocate(this%dicoCG(ii)%cg(size(this%dico_lb0(ii)%cg,1), size(this%dico_lb0(ii)%cg,2)))
          this%dicoCG(ii)%cg = this%dico_lb0(ii)%cg
        end if
        if (ace_gencg == GENCG_DUSORT) then
          this%dicoCG(ii)%dimr_cg = size(this%dico_mub_nb_lb0(ii)%rpicg, 1)
          this%dicoCG(ii)%dimc_cg = size(this%dico_mub_nb_lb0(ii)%rpicg, 2)
          allocate(this%dicoCG(ii)%cg(size(this%dico_mub_nb_lb0(ii)%rpicg,1), size(this%dico_mub_nb_lb0(ii)%rpicg,2)))
          this%dicoCG(ii)%cg = this%dico_mub_nb_lb0(ii)%rpicg
        end if
      end do
      !debug call dump_dico_cg(this, 'dico_cg_fortran.txt')
      !debug stop 'stop after final cg dump'
    end subroutine compute_cg 


    subroutine base_cnlm_compute_ri_gencg(this, nu) 
        use module_l_space_L0, only : graph 
        use matrix_printer, only: int_printMatrix
        use module_ace_desc, only: GENCG_DRAFT,  ace_gencg
        class(type_base_cnlm), intent(inout) :: this
        integer, intent(in) :: nu
        integer, dimension(nu) :: lbold, mbold 

        integer :: ii, Msize, Lsize, Lcols, LL, mm, lev, cnode, p1, p2, id, no_nodes
        integer :: L1, L2, L12, M1, M2, M12, im 
        integer, dimension(:), allocatable :: LLbold, Mn 
        real(dp) :: tt

        if  (nu /= this%nu) then 
            call log_critical("Fatal: compute the L-graph for wrong nu order"//vtoa(this%nu)//" "//vtoa(nu))
            stop 'graph: wrird topology base_cnlm_compute_ri_gencg this nu error '  
        end if

        if  (graph%nu /= this%nu) then 
            call log_critical("Fatal: compute the L-graph for wrong nu order"//vtoa(this%nu)//" "//vtoa(nu))
            stop 'graph: weirg topology base_cnlm_compute_ri_gencg this graph error '  
        end if

        no_nodes = size(graph%nodes)
        if (allocated(Mn)) deallocate(Mn) ; allocate(Mn(no_nodes))

        Lcols = size(this%dico_lb0(1)%Lmat, 2) 
        if (allocated(LLbold)) deallocate(LLbold) ; allocate(LLbold(Lcols)) 
        if (Lcols /= (2*nu-1)) then
           call log_critical("Fatal: weird topology in the L-graph, Lmat has not 2*nu-1 size. Maybe wrong nu order"//vtoa(this%nu)//" "//vtoa(nu)) 
          stop 'stop with message: weird topology in the L-graph'
        end if 
      
        do ii = 1, size(this%dico_lb0)
            lbold = this%dico_lb0(ii)%lbold
            !call int_printMatrix(this%dico_lb0(ii)%Lmat, 'Lmat') 
            Lsize = size(this%dico_lb0(ii)%Lmat, 1)
            Msize = size(this%dico_lb0(ii)%M0mat, 1)
            ! write(*,*) 'Lsize, Msize', Lsize, Msize
            if (allocated(this%dico_lb0(ii)%cg)) deallocate(this%dico_lb0(ii)%cg) ; allocate(this%dico_lb0(ii)%cg(Msize, Lsize))
            do LL = 1, Lsize
                LLbold(:) = this%dico_lb0(ii)%Lmat(LL, :)
                do mm = 1, Msize 
                    ! put mbold componenet in the right places in Mn(1:2*nu-1)
                    Mn(:)= 0
                    mbold(:) = this%dico_lb0(ii)%M0mat(mm, :)
                    do im = 1, size(graph%level0_nodes)
                       Mn(graph%level0_nodes(im)) = mbold(im)
                    end do
                    tt = 1.0_dp 
                    do lev = 1, maxval(graph%nodes(:)%level)
                      do cnode = 1, size(graph%nodes)
                        if (graph%nodes(cnode)%level == lev) then
                          id = graph%nodes(cnode)%id
                          p1 = graph%nodes(cnode)%parent1
                          p2 = graph%nodes(cnode)%parent2
                          !L for parents and child (node) as it is given by Lmat element  
                          L1  = LLbold(p1)
                          L2  = LLbold(p2)
                          L12 = LLbold(id)
                          !M for parents and child (node) as it is given by compising selection rule 
                          !in Clebsh-Gordon coefficients m12 = m1 + m2. Otherwise is zero.  
                          Mn(id) = Mn(p1) + Mn(p2)
                          M1  = Mn(p1)
                          M2  = Mn(p2)
                          M12 = Mn(id)
                          tt = tt * this%precg(L1, M1, L2, M2, L12, M12)
                        end if
                      end do
                    end do !lev 
                    this%dico_lb0(ii)%cg(mm, LL) = tt
                    !debug! write(*,*) 'cgttt', tt, '(', mbold(:), ')', '(', LLbold(:), ')'
                end do  !mm for mbold 
            end do      !LL for LLbold (Lmat)
        end do          !ii for lbold 

        if (ace_gencg == GENCG_DRAFT)  call graph%deallocate_graph_cg(nu)              
    end subroutine base_cnlm_compute_ri_gencg
    

    !this can be MPI parallelized over ii size(this%dico_mub_nb_lb0) over all procs.
    subroutine base_cnlm_compute_rpi_gencg(this, nu)
      use module_l_space_L0, only : graph 
      use permutation_module, only: permutation_object
      use module_svd_small_matrix, only: eigen_svd
      use time_check_general, only: MY_MPI_WTIME
      class(type_base_cnlm), intent(inout) :: this
      integer, intent(in) :: nu
      !local 
      type(permutation_object) :: perm_obj 
      type(eigen_svd) :: svd_obj
      integer :: im, icg, il, nn, ll, mm, ll1, ll2, ii, pp, mu, size_Lbold0, size_M0
      integer, dimension(nu) :: nbold, lbold,  mbold 
      integer, dimension(:), allocatable :: n_perm, l_perm, mu_perm 
      
      integer, dimension(:), allocatable :: commun_perm01, commun_perm02, commun_perm, sigma_mbold, vect_Lmat
      integer, dimension(:,:), allocatable  :: M0mat 
      real(dp), dimension(:,:), allocatable :: tmpGram, mat_cg
      real(dp), dimension(:,:,:), allocatable ::  mat_cg2 
      real(dp) :: term, cg1, cg2, tmp, tmptmp 
      real(dp) :: time1, time2, time_perm, time_build, time_svd 


      if (nu > 3 ) then 
        call perm_obj%initialize(nu)

        if (allocated(vect_Lmat)) deallocate(vect_Lmat) ; allocate(vect_Lmat(2*nu-1))
        time_perm = 0.0_dp
        time_build = 0.0_dp
        time_svd = 0.0_dp
        do ii = 1, size(this%dico_mub_nb_lb0)
          nn = this%dico_mub_nb_lb0(ii)%idx_nb 
          ll = this%dico_mub_nb_lb0(ii)%idx_lb0
          mu = this%dico_mub_nb_lb0(ii)%idx_mub
          !get the nbold and lbold and permutations list 
          nbold(:) = this%dico_nb(nn)%nbold
          if (allocated(n_perm)) deallocate(n_perm) ; allocate(n_perm(size(this%dico_nb(nn)%list_permutations)))
          n_perm(:) = this%dico_nb(nn)%list_permutations
          lbold(:) = this%dico_lb0(ll)%lbold
          if (allocated(l_perm)) deallocate(l_perm) ; allocate(l_perm(size(this%dico_lb0(ll)%list_permutations)))
          l_perm(:) = this%dico_lb0(ll)%list_permutations
          if (allocated(mu_perm)) deallocate(mu_perm) ; allocate(mu_perm(size(this%dico_mu(mu)%list_permutations)))
          mu_perm(:) = this%dico_mu(mu)%list_permutations

          time1 =  MY_MPI_WTIME()
          !TODOace ... maybe is way to optimize that. 
          !match and get the common permutations ... 
          call perm_obj%extract_common_permutations(n_perm, l_perm, commun_perm01) 
          call perm_obj%extract_common_permutations(mu_perm, l_perm, commun_perm02) 
          call perm_obj%extract_common_permutations(commun_perm01, commun_perm02, commun_perm)

          time2 =  MY_MPI_WTIME()
          time_perm = time_perm + time2 - time1
          size_Lbold0 = size(this%dico_lb0(ll)%Lmat,1) 
          size_M0 = size(this%dico_lb0(ll)%M0mat,1)
          !typical for nu=9 permutations......2      131072                 40320        40320           91        110
          !write(*,*) "permutations.........", ii, size(this%dico_mub_nb_lb0), size(n_perm), size(l_perm), size_Lbold0, size_M0
        
          if (allocated(tmpGram)) deallocate(tmpGram) ; allocate(tmpGram(size_Lbold0, size_Lbold0))
          !use pointer ? 
          !write(*,*) 'TTTEST ', size_M0, size_Lbold0, size(this%dico_lb0(ll)%Lmat(:, :),1), size(this%dico_lb0(ll)%Lmat(:, :),2)
          if (allocated(mat_cg)) deallocate(mat_cg) ; allocate(mat_cg(size_M0, size_Lbold0))
          mat_cg = this%dico_lb0(ll)%cg(:, :)

          if (allocated(M0mat)) deallocate(M0mat) ; allocate(M0mat(size_M0, 2*nu-1))
          M0mat = this%dico_lb0(ll)%M0mat

          ! Here I precompute the cg coefficients for the M0mat_permuted and Lmat
          if (allocated(mat_cg2)) deallocate(mat_cg2) ; allocate(mat_cg2(size_M0,  size(commun_perm), size_Lbold0))
          do ll2 = 1, size_Lbold0
            vect_Lmat = this%dico_lb0(ll)%Lmat(ll2, :)
            do pp = 1, size(commun_perm)
              do mm = 1, size_M0
                mbold(:) = M0mat(mm, :)
                if (allocated(sigma_mbold)) deallocate(sigma_mbold) ; allocate(sigma_mbold(size(mbold)))
                call perm_obj%apply_permutation(mbold, sigma_mbold, commun_perm(pp))
                call from_mbold_lLbold_get_cg(this, nu, sigma_mbold, vect_Lmat, cg2)
                mat_cg2(mm, pp, ll2) = cg2
              end do   
            end do
          end do 

          call svd_obj%init(size_Lbold0)
          tmpGram = 0.0_dp
          do ll1 = 1, size_Lbold0
            do ll2 = 1, size_Lbold0
              !write(*,*) 'sizes', size_Lbold0, size(vect_Lmat), size(this%dico_lb0(ll)%Lmat(:, :),1), size(this%dico_lb0(ll)%Lmat(:, :),2)
              !$! vect_Lmat = this%dico_lb0(ll)%Lmat(ll2, :)
              term = 0.0_dp
              do pp = 1, size(commun_perm)
                !$! do mm = 1, size(this%dico_lb0(ll)%M0mat, 1)
                do mm = 1, size_M0 
                  !$! mbold(:) = this%dico_lb0(ll)%M0mat(mm, :)
                  !$! mbold(:) = M0mat(mm, :)
                  !$! if (allocated(sigma_mbold)) deallocate(sigma_mbold) ; allocate(sigma_mbold(size(mbold)))
                  !$! call perm_obj%apply_permutation(mbold, sigma_mbold, commun_perm(pp))
                  !$! compute the CG coefficients
                  !$! call from_mbold_lLbold_get_cg(this, nu, sigma_mbold, vect_Lmat, cg2)
                  
                  cg2 = mat_cg2(mm, pp, ll2)
                  !cg1 = this%dico_lb0(ll)%cg(mm, ll1)
                  cg1 = mat_cg(mm, ll1)
                  term = term + cg1 * cg2
                end do !mm 
              end do   !pp 
              tmpGram(ll1, ll2) = term
            end do  !ll1
          end do    !ll2 
          time1 =  MY_MPI_WTIME()
          time_build = time_build + time1 - time2

          !compute SVD of the Gram matrix and  filter-it ... 
          call svd_obj%evaluate(tmpGram)
          !d svd_obj%rank = size_Lbold0
          if (allocated(this%dico_mub_nb_lb0(ii)%rpicg))  deallocate(this%dico_mub_nb_lb0(ii)%rpicg)
          allocate(this%dico_mub_nb_lb0(ii)%rpicg(size_M0, svd_obj%rank))
          
          do icg = 1, svd_obj%rank         
            do im = 1, size_M0
              tmp = 0.0_dp
              do il = 1, size_Lbold0
                !TODOace ... it is not svd_obj%eigvecU(il,icg) ? 
                !tmptmp = svd_obj%eigvecU(icg,il) * this%dico_lb0(ll)%cg(im, il)
                tmptmp = svd_obj%eigvecU(il,icg) * this%dico_lb0(ll)%cg(im, il)
                tmp = tmp + tmptmp
              end do
              this%dico_mub_nb_lb0(ii)%rpicg(im, icg) = tmp / dsqrt(svd_obj%eigval(icg))
            end do  ! im 
          end do    !icg 

          !destroy svd ... 
          call svd_obj%destroy()
          !TODOace ... SVD 
          !debugACE! write(*,*) ii,  'SVD rank ', svd_obj%rank, ' size_Lbold0 ', size_Lbold0
          time2 = MY_MPI_WTIME()
          time_svd = time_svd + time2 - time1
          
          !debugcall repport_time(2, 0.d0, time_perm,  "ML: Permutations")
          !debugcall repport_time(2, 0.d0, time_build, "ML:   Build Gram")
          !debugcall repport_time(2, 0.d0, time_svd,   "ML:  DO FULL SVD")

        end do ! ii


      end if 

      !here is the case ...  nu < 4 
      if (nu <= 3 ) then 

        do ii = 1, size(this%dico_mub_nb_lb0)
          ll = this%dico_mub_nb_lb0(ii)%idx_lb0
          size_Lbold0 = size(this%dico_lb0(ll)%M0mat,2)
          size_M0 = size(this%dico_lb0(ll)%M0mat,1)  
          if (allocated(this%dico_mub_nb_lb0(ii)%rpicg)) deallocate(this%dico_mub_nb_lb0(ii)%rpicg)
          allocate(this%dico_mub_nb_lb0(ii)%rpicg(size_M0, size_Lbold0))
          this%dico_mub_nb_lb0(ii)%rpicg =  this%dico_lb0(ll)%cg
        end do 
 
      end if 

      ! deallocate all the graph structure. 
      call graph%deallocate_graph_cg(nu) 
    end subroutine base_cnlm_compute_rpi_gencg


    subroutine from_mbold_lLbold_get_cg(this, nu, mbold, LLbold, cgri)
        use module_l_space_L0, only : graph 
        class(type_base_cnlm), intent(inout) :: this
        integer, dimension(:), intent(in) :: mbold, LLbold
        integer, intent(in) :: nu 
        real(dp), intent(out) :: cgri
        integer :: lev, cnode, p1, p2, id
        integer :: L1, L2, L12, M1, M2, M12, im
        real(dp) :: tt
        integer, dimension(2*nu-1) ::  Mn
        
        Mn(:)= 0
        !mbold(:) = this%dico_lb0(ii)%M0mat(mm, :)
        do im = 1, size(graph%level0_nodes)
           Mn(graph%level0_nodes(im)) = mbold(im)
        end do
        tt = 1.0_dp 
        do lev = 1, maxval(graph%nodes(:)%level)
          do cnode = 1, size(graph%nodes)
            if (graph%nodes(cnode)%level == lev) then
              id = graph%nodes(cnode)%id
              p1 = graph%nodes(cnode)%parent1
              p2 = graph%nodes(cnode)%parent2
              !L for parents and child (node) as it is given by Lmat element  
              L1  = LLbold(p1)
              L2  = LLbold(p2)
              L12 = LLbold(id)
              !M for parents and child (node) as it is given by compising selection rule 
              !in Clebsh-Gordon coefficients m12 = m1 + m2. Otherwise is zero.  
              Mn(id) = Mn(p1) + Mn(p2)
              M1  = Mn(p1)
              M2  = Mn(p2)
              M12 = Mn(id)
              tt = tt * this%precg(L1, M1, L2, M2, L12, M12)
            end if
          end do
        end do !lev 
        cgri  = tt
    end subroutine from_mbold_lLbold_get_cg


    subroutine base_cnlm_init(this, nu, nmax, lmax, mumax)
      use module_l_space_first, only: generate_l_space_ini
      use module_l_space_L0, only: generate_l_space_L0, find_permutations_for_l, base_L0, classify_tuples_for_l, &
                                   dump_classify_l
      use module_n_space_first, only: generate_n_space_ini,  generate_n_space_ini_acek, find_permutations_for_n,  base_n0
      use module_m_space_first, only: generate_M0_space_for_l
      use module_mu_space_first, only: generate_mu_space_ini, generate_mu_space_ini_acek, generate_mu_space_ini_acek_block,  &
                  find_permutations_for_mu, base_mu, stamp_permutations_for_mu
      use module_tuples_for_given_lbold, only: generate_tuples_for_given_lbold
      use module_ace_desc, only: ace_radial_chem
      class(type_base_cnlm), intent(inout) :: this
      integer, intent(in) :: nu, lmax, nmax
      integer, intent(inout) :: mumax 
      !local variables ...
      integer :: dim_l, dim_n, dim_mu 
      integer, dimension(:,:), allocatable :: l_space_ini, n_space_ini, mu_space_ini, stamp_mu_space 
      integer :: Lval
      integer, dimension(:,:), allocatable :: l_space_L0
      !debugACE! character(len=80) :: CHFMT
      integer :: ii, iicc, isp,  tmpmax, llmax
      ! l-permutation 
      integer  :: num_classes_l
      integer, dimension(:,:), allocatable :: unique_classes_l  
      integer, dimension(:), allocatable :: map_into_classes_l, map_into_tuples_l 



      this%nu = nu
      
      !---------------generate full mu_space---------------------
      !maximum number of species. 
    
      !step mu1: generate mu_bold space with all possible combinations 
      !TODhea!
      !debug mumax=4 
        if (ace_radial_chem == ACE_CHEM_RADIAL_HSVD .or. ace_radial_chem == ACE_CHEM_RADIAL_BLOCK_HSVD &
          .or. ace_radial_chem == ACE_CHEM_RADIAL_RANDPROJ .or. ace_radial_chem == ACE_CHEM_RADIAL_CHEMMAP_HSVD) then
         ! HSVD, BLOCK_HSVD, and CHEMMAP_HSVD: use diagonal mu-space (same species per position)
        if (ace_radial_chem == ACE_CHEM_RADIAL_HSVD .or. ace_radial_chem == ACE_CHEM_RADIAL_RANDPROJ &
          .or. ace_radial_chem == ACE_CHEM_RADIAL_CHEMMAP_HSVD) then
          dim_mu = 1
          call generate_mu_space_ini_acek(nu, mumax, mu_space_ini, dim_mu)
        else
          call generate_mu_space_ini_acek_block(nu, mumax, mu_space_ini, dim_mu)
        end if
      else
         ! Standard ACE (RALF): generate full tensor product of species
         call generate_mu_space_ini(nu, mumax, mu_space_ini, dim_mu)
      end if
      call stamp_permutations_for_mu(nu, mumax, mu_space_ini, stamp_mu_space)
      !TODhea! 

      if (allocated(this%dico_mu)) deallocate(this%dico_mu) ; allocate(this%dico_mu(dim_mu))
      do ii = 1, dim_mu
        if (allocated(this%dico_mu(ii)%mubold)) deallocate(this%dico_mu(ii)%mubold) ; allocate(this%dico_mu(ii)%mubold(nu))
        this%dico_mu(ii)%mubold(:) = mu_space_ini(ii, :)
        if (allocated(this%dico_mu(ii)%mustamp)) deallocate(this%dico_mu(ii)%mustamp) ; allocate(this%dico_mu(ii)%mustamp(mumax))
        this%dico_mu(ii)%mustamp(:) = stamp_mu_space(ii, :)
        !debug write(*,* ) 'mubold ', this%dico_mu(ii)%mubold(:)
      end do
      deallocate(stamp_mu_space)
      !step mu2: find for each mubold (a column in mu_space_ini(ii,:) ) all the permutations 
      call find_permutations_for_mu(nu, mu_space_ini)


      !step mu3: pack into this%dico_mu
      do ii = 1, dim_mu 
        if (allocated(this%dico_mu(ii)%list_permutations)) deallocate(this%dico_mu(ii)%list_permutations) ; allocate(this%dico_mu(ii)%list_permutations(nu))
        this%dico_mu(ii)%list_permutations = base_mu(ii)%list_permutations
        deallocate(base_mu(ii)%list_permutations)
      end do
      deallocate(mu_space_ini)
      deallocate(base_mu)
      !debug_test call dump_dico_mu_raw('dico_mu_fortran.txt', nu, mumax, this%dico_mu)
      !debug_test stop ' stop after dump_dico_mu'
      !---------------generate l_space---------------------
      ! step l1: generate the initial l space with all possible combinations
      ! All possible (l_1, l_2, .. l_nu) such that
      ! \sum l_i is even
      ! 0 <= l_i <= lmax  

      !debugACE! write(*,*)  'nu, lmax, nmax = ', nu, lmax     

      call generate_l_space_ini(nu, lmax, l_space_ini, dim_l)
      !dico_lb! if (allocated(this%dico_lb)) deallocate(this%dico_lb) ; allocate(this%dico_lb(dim_l))  
      !dico_lb! do ii = 1, dim_l
      !dico_lb!     if (allocated(this%dico_lb(ii)%lbold)) deallocate(this%dico_lb(ii)%lbold) ; allocate(this%dico_lb(ii)%lbold(nu))
      !dico_lb!     this%dico_lb(ii)%lbold(:) = l_space_ini(ii, :)
      !dico_lb! end do
      call log_info("ML: ACE body "//vtoa(this%nu)//" the number of naive lbold is        :" //vtoa(dim_l))   
      !step l2: refine  l_space_ini into l_space_L0 for 
      !        which the composed L_123...\nu is 0, i.e. invariant under rotation
      ! at the same price we can generate the base_L0 too that will be used in order to compute the generalized CG coefficients. 
      Lval = 0 
      call generate_l_space_L0(nu, Lval, l_space_ini, l_space_L0)
      !step 2.5: find all the permutation of each l in the  l_space_L0
      call log_info("ML: ACE body "//vtoa(this%nu)//" the number of lbold with L=0 is     :" //vtoa(size(l_space_L0, 1)))   

      call find_permutations_for_l(nu, l_space_L0)

      !debug!  do ii = 1, size(l_space_L0, 1)
      !debug!      write(*,*) 'l_space_L0 ', ii, ' -> ', l_space_L0(ii, :)
      !debug!  end do  

      ! compress to take only the unique classes in of l_space_L0 with respect to the permutations
      call classify_tuples_for_l(l_space_L0, num_classes_l, map_into_classes_l,  map_into_tuples_l, unique_classes_l) 
      !debug_test call dump_classify_l('classify_l_fortran.txt', l_space_L0, map_into_classes_l, map_into_tuples_l, unique_classes_l)
      !debug_test stop 'stop after dump_classify_l'
      call log_info("ML: ACE body "//vtoa(this%nu)//" the final number of unique lbold is :"// vtoa(num_classes_l))    

      !debug! write(*,*) 'unique-classes', size(unique_classes_l, 1), size(unique_classes_l, 2)
      !debug! do ii = 1, num_classes_l
      !debug!     write(*,*) 'num_classes_l ', ii, ' -> ', unique_classes_l(: , ii)
      !debug! end do  
      

      !step l3: pack into this%dico_lb0
      !l_full_not_erase   if (allocated(this%dico_lb0)) deallocate(this%dico_lb0) ; allocate(this%dico_lb0(size(base_L0)))
      if (allocated(this%dico_lb0)) deallocate(this%dico_lb0) ; allocate(this%dico_lb0(num_classes_l))

      do iicc = 1, size(this%dico_lb0)
          ii = map_into_tuples_l(iicc)
          !l_full_not_erase
          !ii = iicc
          if (allocated(this%dico_lb0(iicc)%lbold)) deallocate(this%dico_lb0(iicc)%lbold) ; allocate(this%dico_lb0(iicc)%lbold(nu))
          this%dico_lb0(iicc)%lbold = l_space_L0(ii, :)
          if (allocated(this%dico_lb0(iicc)%list_permutations)) deallocate(this%dico_lb0(iicc)%list_permutations) ; 
          isp = size(base_L0(ii)%list_permutations)
          allocate(this%dico_lb0(iicc)%list_permutations(isp))
          this%dico_lb0(iicc)%list_permutations = base_L0(ii)%list_permutations
          this%dico_lb0(iicc)%dimr_Lmat = base_L0(ii)%dimr_Lmat
          this%dico_lb0(iicc)%dimc_Lmat = base_L0(ii)%dimc_Lmat
          if (allocated(this%dico_lb0(iicc)%Lmat)) deallocate(this%dico_lb0(iicc)%Lmat) ; allocate(this%dico_lb0(iicc)%Lmat(base_L0(ii)%dimr_Lmat, base_L0(ii)%dimc_Lmat))
          this%dico_lb0(iicc)%Lmat = base_L0(ii)%Lmat
      end do

      
      llmax = 0 
      do ii = 1, size(base_L0) 
        deallocate(base_L0(ii)%list_permutations)
        tmpmax = maxval(base_L0(ii)%Lmat)
        if (tmpmax > llmax) llmax = tmpmax
        deallocate(base_L0(ii)%Lmat)
      end do   
      this%llmax = llmax



      deallocate(l_space_ini)
      deallocate(l_space_L0)
      deallocate(base_L0)
      !debugACE! write (CHFMT, *) '(" l_bold and size_L (", ', nu, '(i4), "  ) ", i4, i5  )'
      do ii = 1, size(this%dico_lb0)
          !call generateM0_for_l(this%dico_lb0(ii)%lbold, this%dico_lb0(ii)%M0mat)
          call generate_M0_space_for_l(this%dico_lb0(ii)%lbold, this%dico_lb0(ii)%M0mat)
          !debugACE! write(*, CHFMT )  this%dico_lb0(ii)%lbold(:), size(this%dico_lb0(ii)%M0mat,1), ii
      end do

      !debug_test call dump_dico_lb0("dico_lb0_fortran.txt",this%nu, lmax, this%dico_lb0)
      !debug_test stop "after dump dico_lb0" 
      !----------------------------------------------------
      

      !---------------generate full n_space---------------------
      !step n1: generate n_bold space with all possible combinations
      !TODhea!
        if (ace_radial_chem == ACE_CHEM_RADIAL_HSVD .or. ace_radial_chem == ACE_CHEM_RADIAL_BLOCK_HSVD &
          .or. ace_radial_chem == ACE_CHEM_RADIAL_RANDPROJ .or. ace_radial_chem == ACE_CHEM_RADIAL_CHEMMAP_HSVD) then
        dim_n = 1 
        call generate_n_space_ini_acek(nu, nmax, n_space_ini, dim_n)
      else 
        call generate_n_space_ini(nu, nmax, n_space_ini, dim_n)
      end if
      !TODhea!
      !debug write(*,*), 'dim_n ini = ', dim_n, '  ', size(n_space_ini, 1), '  ', size(n_space_ini, 2)

      
      !step n2: find for each column nbold in n_space_ini  (dim_n columns) all the permutations 
      !         in S_\nu that leaves invariant nbold    
      call find_permutations_for_n(nu, n_space_ini)


      !step n3: pack into this%dico_nb0
      if (allocated(this%dico_nb)) deallocate(this%dico_nb) ; allocate(this%dico_nb(dim_n))
      do ii = 1, dim_n
          if (allocated(this%dico_nb(ii)%nbold)) deallocate(this%dico_nb(ii)%nbold) ; allocate(this%dico_nb(ii)%nbold(nu))
          this%dico_nb(ii)%nbold = n_space_ini(ii, :)
      end do

      do ii = 1, dim_n 
          if (allocated(this%dico_nb(ii)%list_permutations)) deallocate(this%dico_nb(ii)%list_permutations) ; allocate(this%dico_nb(ii)%list_permutations(nu))
          this%dico_nb(ii)%list_permutations = base_n0(ii)%list_permutations
          deallocate(base_n0(ii)%list_permutations)
      end do

      !still needded deallocate(n_space_ini)
      deallocate(base_n0)
      !----------------------------------------------------
      !debug_test call dump_dico_nb('dump_nb_fortran.txt',nu, nmax, this%dico_nb) 
      !debug_test stop 'dump_nb_fortran'

      deallocate(n_space_ini)

      ! new way 
      !setup the dico this%dico_mub_nb_lb0
      ! NOTE:
      ! - HSVD and CHEMMAP_HSVD collapse chemical tuples (dim_mu=1 in basis), so the acek shortcut is valid.
      ! - BLOCK_HSVD keeps only diagonal mubold tuples (mubold=(mu,...,mu)) but preserves
      !   the mu-channel index mu=1..mumax (radial is still collapsed to k).
        if (ace_radial_chem == ACE_CHEM_RADIAL_HSVD .or. ace_radial_chem == ACE_CHEM_RADIAL_RANDPROJ &
          .or. ace_radial_chem == ACE_CHEM_RADIAL_CHEMMAP_HSVD) then
        call base_cnlm_build_mubnblb_acek(this, nu, nmax, mumax)
      else if (ace_radial_chem == ACE_CHEM_RADIAL_BLOCK_HSVD) then
        call base_cnlm_build_mubnblb_block_acek(this, nu, nmax, mumax)
        !call base_cnlm_build_mubnblb_acek(this, nu, nmax, mumax)
      else
        call base_cnlm_build_mubnblb(this, nu, nmax, mumax)
      end if

      call log_info("ML: ACE body "//vtoa(this%nu)//" sizes: dico_mu="//vtoa(size(this%dico_mu))// &
                   " dico_mub_nb_lb0="//vtoa(size(this%dico_mub_nb_lb0)))
      !debug call dump_dico_mub_nb_lb0_file(this, 'dico_mub_nb_lb0_fortran.txt')
      !stop 'stop after dump dico_mub_nb_lb0'

      !debugACE! write (CHFMT, *) '(" n_bold  (", ', nu, '(i4), "  ) ", i4 )'
      !debugACE! do ii = 1, dim_n
      !debugACE!     write(*,CHFMT) this%dico_nb(ii)%nbold(:), ii
      !debugACE! end do


    end subroutine base_cnlm_init


    subroutine base_cnlm_build_mubnblb_block_acek(this, nu, nmax, mumax)
      ! BLOCK_HSVD: radial tuples are collapsed (handled elsewhere via k), but keep mu channels.
      ! Enumerate only diagonal mubold tuples: mubold=(mu,mu,...,mu) for mu=1..mumax.
      use module_ace_desc, only: ace_gencg, GENCG_DRAFT, GENCG_DUSORT
      implicit none
      class(type_base_cnlm), intent(inout) :: this
      integer, intent(in) :: nu, nmax, mumax
      integer :: ii, icount, mu
      integer :: icg

      icount = size(this%dico_lb0) * mumax
      if (allocated(this%dico_mub_nb_lb0)) deallocate(this%dico_mub_nb_lb0)
      allocate(this%dico_mub_nb_lb0(icount))

      icount = 0
      do ii = 1, size(this%dico_lb0)
        this%dico_lb0(ii)%unique_nbold = 1
        this%dico_lb0(ii)%max_storage_nmu = mumax
        do mu = 1, mumax
          icount = icount + 1
          this%dico_mub_nb_lb0(icount)%idx_lb0 = ii
          this%dico_mub_nb_lb0(icount)%idx_nb  = 1
          this%dico_mub_nb_lb0(icount)%idx_mub = mu
          if (ace_gencg == GENCG_DRAFT) icg = this%dico_mub_nb_lb0(icount)%idx_lb0
          if (ace_gencg == GENCG_DUSORT) icg = icount
          this%dico_mub_nb_lb0(icount)%idx_cg = icg
        end do
      end do

    end subroutine base_cnlm_build_mubnblb_block_acek


    subroutine base_cnlm_build_mubnblb_acek(this, nu, nmax, mumax)
      use module_ace_desc, only: ace_gencg, GENCG_DRAFT, GENCG_DUSORT
      implicit none 
      class(type_base_cnlm), intent(inout) :: this
      integer, intent(in) :: nu, nmax, mumax
      integer :: ii, icount
      integer, dimension(nu) :: tmp_lbold, tmp_nbold, tmp_mubold
      integer :: nmu
      integer :: idx, icg 
     
      ! First sorting .............. build the relevant nbold and mubold for each lbold. Avoid the repetitions by permutations.
      ! First step: build the storage_nmu for each lbold.
      icount = 0
      do ii = 1, size(this%dico_lb0)
        tmp_lbold = this%dico_lb0(ii)%lbold 
        this%dico_lb0(ii)%max_storage_nmu = 1
        do nmu = 1, this%dico_lb0(ii)%max_storage_nmu
          icount = icount + 1
        enddo   
      end do 
      
      if (allocated(this%dico_mub_nb_lb0)) deallocate(this%dico_mub_nb_lb0) ; allocate(this%dico_mub_nb_lb0(icount))

      icount = 0 
      do ii = 1, size(this%dico_lb0)
        tmp_lbold = this%dico_lb0(ii)%lbold 
        this%dico_lb0(ii)%unique_nbold = 1
        ! max_storage_nmu is 1 for acek shortcut
        do nmu = 1, this%dico_lb0(ii)%max_storage_nmu
          icount = icount + 1
          this%dico_mub_nb_lb0(icount)%idx_lb0 = ii 
          this%dico_mub_nb_lb0(icount)%idx_nb  = nmu 
          this%dico_mub_nb_lb0(icount)%idx_mub = nmu 
          if (ace_gencg == GENCG_DRAFT) icg = this%dico_mub_nb_lb0(icount)%idx_lb0
          if (ace_gencg == GENCG_DUSORT) icg = icount 
          this%dico_mub_nb_lb0(icount)%idx_cg = icg 
        end do
      end do 

    end subroutine base_cnlm_build_mubnblb_acek



    subroutine base_cnlm_build_mubnblb(this, nu, nmax, mumax)
      use module_tuples_for_given_lbold, only: generate_tuples_for_given_lbold
      use module_perm_basis_remove , only: sort_basis
      use module_ace_desc, only: ace_gencg, GENCG_DRAFT, GENCG_DUSORT, ace_radial_chem
      implicit none 
      class(type_base_cnlm), intent(inout) :: this
      integer, intent(in) :: nu, nmax, mumax
      integer :: ii, icount
      integer, dimension(nu) :: tmp_lbold, tmp_nbold, tmp_mubold
      integer :: nmu, itmp, ss
      type type_buffer 
        integer :: idx_lb0
        integer :: idx_nb
        integer :: idx_mub
      end type type_buffer
      type(type_buffer), dimension(:), allocatable :: dico_buffer 
      integer, dimension(:,:,:), allocatable :: allBases 
      integer, dimension(:), allocatable :: index_of_sorting 
      integer :: idx, icg, k, idx_nb_calc, idx_mub_calc 
     
      ! First sorting .............. build the relevant nbold and mubold for each lbold. Avoid the repetitions by permutations.
      ! First step: build the storage_nmu for each lbold.
      icount = 0
      do ii = 1, size(this%dico_lb0)
        tmp_lbold = this%dico_lb0(ii)%lbold 
        call generate_tuples_for_given_lbold(nu, nmax, mumax, tmp_lbold, this%dico_lb0(ii)%unique_nbold, this%dico_lb0(ii)%storage_nmu)
        this%dico_lb0(ii)%max_storage_nmu = size(this%dico_lb0(ii)%storage_nmu, 2) 
        do nmu = 1, this%dico_lb0(ii)%max_storage_nmu
          icount = icount + 1
        enddo   
      end do 
      
      if (allocated(this%dico_mub_nb_lb0)) deallocate(this%dico_mub_nb_lb0) ; allocate(this%dico_mub_nb_lb0(icount))

      icount = 0
      do ii = 1, size(this%dico_lb0)
        tmp_lbold = this%dico_lb0(ii)%lbold 
        this%dico_lb0(ii)%max_storage_nmu = size(this%dico_lb0(ii)%storage_nmu, 2) 
        do nmu = 1, this%dico_lb0(ii)%max_storage_nmu
          icount = icount + 1
          this%dico_mub_nb_lb0(icount)%idx_lb0 = ii
          
          ! Calculate idx_nb and idx_mub
          tmp_nbold = this%dico_lb0(ii)%storage_nmu(1:nu, nmu)
          tmp_mubold = this%dico_lb0(ii)%storage_nmu(nu+1:2*nu, nmu)

          ! 1. Calculate idx_nb (Radial Index)
          ! RALF: Full tensor product for both n and mu
          idx_nb_calc = 1
          do k = 1, nu
             idx_nb_calc = idx_nb_calc + (tmp_nbold(k) - 1) * (nmax**(nu-k))
          end do
          this%dico_mub_nb_lb0(icount)%idx_nb = idx_nb_calc

          ! 2. Calculate idx_mub (Chemical Index)
          ! RALF: Full tensor product for mu
          idx_mub_calc = 1
          do k = 1, nu
             idx_mub_calc = idx_mub_calc + (tmp_mubold(k) - 1) * (mumax**(nu-k))
          end do
          this%dico_mub_nb_lb0(icount)%idx_mub = idx_mub_calc

          if (ace_gencg == GENCG_DRAFT) icg = this%dico_mub_nb_lb0(icount)%idx_lb0
          if (ace_gencg == GENCG_DUSORT) icg = icount 
          this%dico_mub_nb_lb0(icount)%idx_cg = icg 
        end do 
      end do 
      ! ... dealloacte storage_nmu
      do ii = 1, size(this%dico_lb0) 
         if (allocated(this%dico_lb0(ii)%storage_nmu)) deallocate(this%dico_lb0(ii)%storage_nmu)
      end do 

      ! Plug the second permutations symmetries ... 
      if (allocated(allBases)) deallocate(allBases) ; allocate(allBases(3, nu, size(this%dico_mub_nb_lb0)))
      if (allocated(dico_buffer)) deallocate(dico_buffer) ; allocate(dico_buffer(size(this%dico_mub_nb_lb0)))
      do ii = 1, size(this%dico_mub_nb_lb0)
        dico_buffer(ii)%idx_lb0 = this%dico_mub_nb_lb0(ii)%idx_lb0
        dico_buffer(ii)%idx_nb = this%dico_mub_nb_lb0(ii)%idx_nb
        dico_buffer(ii)%idx_mub = this%dico_mub_nb_lb0(ii)%idx_mub
        allBases(3, :, ii) = this%dico_lb0(this%dico_mub_nb_lb0(ii)%idx_lb0)%lbold
        allBases(2, :, ii) = this%dico_nb(this%dico_mub_nb_lb0(ii)%idx_nb)%nbold
        allBases(1, :, ii) = this%dico_mu(this%dico_mub_nb_lb0(ii)%idx_mub)%mubold
      end do

      !call base_cnlm_build_sorted(this, is_unique)
      call sort_basis(this%nu, size(this%dico_mub_nb_lb0), allBases, index_of_sorting) 
      
       if (allocated(this%dico_mub_nb_lb0)) deallocate(this%dico_mub_nb_lb0) ; allocate(this%dico_mub_nb_lb0(size(index_of_sorting)))
      icount = 0
      do ii = 1, size(index_of_sorting)
          icount = icount + 1
          idx = index_of_sorting(ii)
          this%dico_mub_nb_lb0(icount)%idx_lb0 = dico_buffer(idx)%idx_lb0
          this%dico_mub_nb_lb0(icount)%idx_nb  = dico_buffer(idx)%idx_nb
          this%dico_mub_nb_lb0(icount)%idx_mub = dico_buffer(idx)%idx_mub
          if (ace_gencg == GENCG_DRAFT) icg = this%dico_mub_nb_lb0(icount)%idx_lb0
          if (ace_gencg == GENCG_DUSORT) icg = icount 
          this%dico_mub_nb_lb0(icount)%idx_cg = icg 
      end do


      if (allocated(dico_buffer)) deallocate(dico_buffer)
      if (allocated(allBases))    deallocate(allBases)
      if (allocated(index_of_sorting)) deallocate(index_of_sorting) 

    end subroutine base_cnlm_build_mubnblb


    subroutine base_cnlm_build_mubnblbLbi(this) ! the call name is build_Bfull 
    class(type_base_cnlm), intent(inout) :: this
    integer :: inl, iLLb 
    integer :: icount, icount_nunl, icg 
    !TODOace ... 
    icount = 0
    icount_nunl = 0      
    
      do inl  = 1, size(this%dico_mub_nb_lb0)
        icg  = this%dico_mub_nb_lb0(inl)%idx_cg 
        icount_nunl = icount_nunl + 1
        icount = icount + size(this%dicoCG(icg)%cg, 2)          
      end do !inl 
    
    if (allocated(this%dicoB)) deallocate(this%dicoB) ; allocate(this%dicoB(icount))
    !if (allocated(this%dico_mub_nb_lb0)) deallocate(this%dico_mub_nb_lb0) ; allocate(this%dico_mub_nb_lb0(icount_nunl))

    icount = 0      
    icount_nunl = 0
    !do imu = 1, size(this%dico_mu) 
      do inl  = 1, size(this%dico_mub_nb_lb0)
        icount_nunl = icount_nunl + 1
        icg  = this%dico_mub_nb_lb0(inl)%idx_cg 
     
        do iLLb = 1, size(this%dicoCG(icg)%cg, 2)
          icount = icount + 1
          this%dicoB(icount)%idx_nb = this%dico_mub_nb_lb0(inl)%idx_nb
          this%dicoB(icount)%idx_mub = this%dico_mub_nb_lb0(inl)%idx_mub
          this%dicoB(icount)%idx_lb0 = this%dico_mub_nb_lb0(inl)%idx_lb0
          this%dicoB(icount)%idx_munblb0 = icount_nunl
          this%dicoB(icount)%idx_Li = iLLb
        end do          
      end do !inl 
    !end do   !imu 

    !$! call log_warning("ML: ACE body "//vtoa(this%nu)//" the total number of B tuples is    "// vtoa(icount))
    !$! call log_warning("ML: ACE body "//vtoa(this%nu)//" the total number of munl tuples is "// vtoa(icount_nunl))


    end subroutine  base_cnlm_build_mubnblbLbi


    subroutine check_uniqueA_tuples(this)
      use mld_mpi, only: mld_mpi_abort
      ! use module_ace_radial, only: radialace, ACE_CHEM_RADIAL_RALF, ACE_CHEM_RADIAL_HOME, ACE_CHEM_RADIAL_HSVD
      use module_ace_desc, only: ace_radial_chem, ace_kmax
      class(type_base_cnlm), intent(inout) :: this
      integer :: ib, im, kk, start
      integer, dimension(this%nu) :: nbold, mbold, lbold, mubold
      integer :: idx_nb, idx_mub, idx_lb0, dimr_M0mat, dim_baseB, dim_munl
      integer :: icount 
      integer, dimension(:,:), allocatable :: mat_mubold
      integer, dimension(:,:), allocatable :: mat_nbold
      integer, dimension(:,:), allocatable :: mat_lbold
      integer, dimension(:,:), allocatable :: mat_mbold, tmptuples 
      integer :: acekmax_nu
      integer :: M, i, j, k, num_unique
      logical :: is_unique
      _NAMECURRENT_('check_uniqueA_tuples')
      _MLD_BEGIN_
      dim_baseB = size(this%dicoB)
      dim_munl = size(this%dico_mub_nb_lb0)
      acekmax_nu = ace_kmax(this%nu)
      !debug write(*,*)  'check', this%nu, dim_baseB, dim_munl, acekmax_nu
      
      icount = 0 
      do ib = 1, dim_munl
        !idx_nb  = this%dicoB(ib)%idx_nb
        !nbold   = this%dico_nb(idx_nb)%nbold
        !idx_mub = this%dicoB(ib)%idx_mub 
        !mubold  = this%dico_mu(idx_mub)%mubold 
        idx_lb0 = this%dico_mub_nb_lb0(ib)%idx_lb0 
        !lbold   = this%dico_lb0(idx_lb0)%lbold 
        !idx_Li  = this%dicoB(ib)%idx_Li
        dimr_M0mat = size(this%dico_lb0(idx_lb0)%M0mat, 1) 
        do im = 1, dimr_M0mat 
          !mbold = this%dico_lb0(idx_lb0)%M0mat(im, :)  
          icount = icount + 1
        end do   
      end do 

      allocate(mat_mubold(icount, this%nu))
      allocate(mat_nbold(icount, this%nu))
      allocate(mat_lbold(icount, this%nu))
      allocate(mat_mbold(icount, this%nu))

      M = icount 
      icount = 0 
      do ib = 1, dim_munl
        
        idx_nb  = this%dico_mub_nb_lb0(ib)%idx_nb
        nbold   = this%dico_nb(idx_nb)%nbold
        idx_mub = this%dico_mub_nb_lb0(ib)%idx_mub 
        mubold  = this%dico_mu(idx_mub)%mubold 
        idx_lb0 = this%dico_mub_nb_lb0(ib)%idx_lb0 
        lbold   = this%dico_lb0(idx_lb0)%lbold 
        !debug write(*,*) 'ib, idx_nb, idx_mub, idx_lb0 ', ib, nbold, mubold, lbold, idx_nb, idx_mub, idx_lb0
        !idx_Li  = this%dicoB(ib)%idx_Li
        dimr_M0mat = size(this%dico_lb0(idx_lb0)%M0mat, 1) 
        do im = 1, dimr_M0mat 
          mbold = this%dico_lb0(idx_lb0)%M0mat(im, :)  
          icount = icount + 1   
          mat_mubold(icount, :) = mubold(:)
          mat_nbold(icount, :) = nbold(:)
          mat_lbold(icount, :) = lbold(:)
          mat_mbold(icount, :) = mbold(:)
        end do 
      end do
      !debugACE write(*,*) '!!!!!!!!!all combination', icount

    if (ace_radial_chem == ACE_CHEM_RADIAL_HSVD .or. ace_radial_chem == ACE_CHEM_RADIAL_BLOCK_HSVD &
      .or. ace_radial_chem == ACE_CHEM_RADIAL_RANDPROJ .or. ace_radial_chem == ACE_CHEM_RADIAL_CHEMMAP_HSVD) then
     call log_info("ML: ACE body "//vtoa(this%nu)//" the total number of A tuples is    "// vtoa(icount*acekmax_nu))  
    end if  
    if (ace_radial_chem == ACE_CHEM_RADIAL_RALF ) then
     call log_info("ML: ACE body "//vtoa(this%nu)//" the total number of A tuples is    "// vtoa(icount))  
    end if

    if (allocated(this%uniqueA_tuples)) deallocate(this%uniqueA_tuples) ; allocate(this%uniqueA_tuples(4, M * this%nu))
    if (allocated(this%map_large_to_uniqueA)) deallocate(this%map_large_to_uniqueA) ;   allocate(this%map_large_to_uniqueA(M*this%nu))

    num_unique = 0

    if (ace_radial_chem == ACE_CHEM_RADIAL_RALF .or. ace_radial_chem == ACE_CHEM_RADIAL_HSVD &
      .or. ace_radial_chem == ACE_CHEM_RADIAL_BLOCK_HSVD .or. ace_radial_chem == ACE_CHEM_RADIAL_RANDPROJ &
      .or. ace_radial_chem == ACE_CHEM_RADIAL_CHEMMAP_HSVD) then
    ! Loop over all collections and vectors
    icount = 0
    do i = 1, M
        do j = 1, this%nu
            icount = icount + 1
            is_unique = .true.
            ! Check if the current tuple is unique
            ! For HSVD and CHEMMAP_HSVD: skip nbold check (always 1), as kk goes in slot 2 later
            do k = 1, num_unique
                if (ace_radial_chem == ACE_CHEM_RADIAL_HSVD .or. ace_radial_chem == ACE_CHEM_RADIAL_RANDPROJ &
                  .or. ace_radial_chem == ACE_CHEM_RADIAL_CHEMMAP_HSVD) then
                  ! HSVD/CHEMMAP: check only mubold, lbold, mbold (nbold always 1)
                  if (this%uniqueA_tuples(1, k) == mat_mubold(i, j) .and. &
                      this%uniqueA_tuples(3, k) == mat_lbold(i, j) .and. &
                      this%uniqueA_tuples(4, k) == mat_mbold(i, j)) then 
                    is_unique = .false.
                    this%map_large_to_uniqueA(icount) = k
                    exit
                  end if
                else
                  ! RALF/BLOCK_HSVD: check all 4 indices including nbold
                  if (this%uniqueA_tuples(1, k) == mat_mubold(i, j) .and. &
                      this%uniqueA_tuples(2, k) == mat_nbold(i, j) .and. &
                      this%uniqueA_tuples(3, k) == mat_lbold(i, j) .and. &
                      this%uniqueA_tuples(4, k) == mat_mbold(i, j)) then 
                    is_unique = .false.
                    this%map_large_to_uniqueA(icount) = k
                    exit
                  end if
                end if
            end do

            ! If unique, add to the unique_tuples array
            if (is_unique) then
                num_unique = num_unique + 1
                this%uniqueA_tuples(1, num_unique) = mat_mubold(i, j)
                if (ace_radial_chem == ACE_CHEM_RADIAL_HSVD .or. ace_radial_chem == ACE_CHEM_RADIAL_RANDPROJ &
                  .or. ace_radial_chem == ACE_CHEM_RADIAL_CHEMMAP_HSVD) then
                  ! For HSVD/CHEMMAP: don't store nbold yet (kk will go here later)
                  this%uniqueA_tuples(3, num_unique) = mat_lbold(i, j)
                  this%uniqueA_tuples(4, num_unique) = mat_mbold(i, j)
                else
                  ! For RALF/BLOCK_HSVD: store all 4 fields
                  this%uniqueA_tuples(2, num_unique) = mat_nbold(i, j)
                  this%uniqueA_tuples(3, num_unique) = mat_lbold(i, j)
                  this%uniqueA_tuples(4, num_unique) = mat_mbold(i, j)
                end if
                this%map_large_to_uniqueA(icount) = num_unique
            end if
            
        end do
    end do
    end if !ACE_CHEM_RADIAL_RALF or BLOCK_HSVD or CHEMMAP_HSVD


    if (ace_radial_chem == ACE_CHEM_RADIAL_HSVD .and. .false.) then
      ! NOTE: This branch is now disabled - HSVD follows same path as RALF/BLOCK_HSVD above
      ! HSVD A-basis design: For HSVD, dim_mu=1 in the basis, so mat_mubold(i,j)=1 always.
      ! The uniqueA_tuples(1,k) will be 1 for all k. 
      ! Loop over all collections and vectors
      ! Use same logic as BLOCK_HSVD: check all 4 indices for uniqueness
      icount = 0
      do i = 1, M
        do j = 1, this%nu
          icount = icount + 1
          is_unique = .true.
          ! Check if the current tuple is unique (including nbold)
          do k = 1, num_unique
            if (this%uniqueA_tuples(1, k) == mat_mubold(i, j) .and. &
              this%uniqueA_tuples(2, k) == mat_nbold(i, j) .and. &
              this%uniqueA_tuples(3, k) == mat_lbold(i, j) .and. &
              this%uniqueA_tuples(4, k) == mat_mbold(i, j)) then 
              is_unique = .false.
              this%map_large_to_uniqueA(icount) = k
              exit
            end if
          end do
          ! If unique, add to the unique_tuples array
          if (is_unique) then
              num_unique = num_unique + 1
              this%uniqueA_tuples(1, num_unique) = mat_mubold(i, j)
              this%uniqueA_tuples(2, num_unique) = mat_nbold(i, j)
              this%uniqueA_tuples(3, num_unique) = mat_lbold(i, j)
              this%uniqueA_tuples(4, num_unique) = mat_mbold(i, j)
              this%map_large_to_uniqueA(icount) = num_unique
          end if
            
        end do
      end do
    end if !ACE_CHEM_RADIAL_HSVD

    if (allocated(tmptuples)) deallocate(tmptuples) ; allocate(tmptuples(4, num_unique))  
    tmptuples(1:4,1:num_unique) = this%uniqueA_tuples(1:4,1:num_unique)
    ! DEBUG: Check tmptuples(1,:) values for HSVD
    ! if (ace_radial_chem == ACE_CHEM_RADIAL_HSVD .or. ace_radial_chem == ACE_CHEM_RADIAL_BLOCK_HSVD) then
    !  call log_info("ML: ACE body "//vtoa(this%nu)//" tmptuples(1,1:5)="//vtoa(tmptuples(1,1))// &
    !               " "//vtoa(tmptuples(1,min(2,num_unique)))//" "//vtoa(tmptuples(1,min(3,num_unique)))// &
    !               " "//vtoa(tmptuples(1,min(4,num_unique)))//" "//vtoa(tmptuples(1,min(5,num_unique))))
    !end if
    deallocate(this%uniqueA_tuples)

    if (ace_radial_chem == ACE_CHEM_RADIAL_HSVD .or. ace_radial_chem == ACE_CHEM_RADIAL_BLOCK_HSVD .or. &
      ace_radial_chem == ACE_CHEM_RADIAL_RANDPROJ .or. ace_radial_chem == ACE_CHEM_RADIAL_CHEMMAP_HSVD) then
      allocate(this%uniqueA_tuples(4, num_unique * acekmax_nu))
      do kk = 1, acekmax_nu
        start = num_unique*(kk-1)
        this%uniqueA_tuples(1,start+1:start+num_unique) = tmptuples(1,1:num_unique)
        this%uniqueA_tuples(2,start+1:start+num_unique) = kk
        this%uniqueA_tuples(3:4,start+1:start+num_unique) = tmptuples(3:4,1:num_unique)
      end do  
      num_unique = num_unique * acekmax_nu
    else if (ace_radial_chem == ACE_CHEM_RADIAL_RALF ) then
      allocate(this%uniqueA_tuples(4, num_unique))
      this%uniqueA_tuples(1:4,1:num_unique) = tmptuples(1:4,1:num_unique)
    else
      call log_warning("ML: in ACE the type_chem_radial is not yet implemented.")
      stop 'stop here in check_uniqueA_tuples'
    end if
    deallocate(tmptuples)

    ! call log_info("ML: ACE body "//vtoa(this%nu)//" the total number of unique A tuples is "// vtoa(num_unique))
    call log_info("ML: ACE body "//vtoa(this%nu)//" the total number of unique A tuples is "// vtoa(size(this%uniqueA_tuples, 2)))  
    call log_info("ML: ACE body "//vtoa(this%nu)//" min k "//vtoa(minval(this%map_large_to_uniqueA))//" max k "//vtoa(maxval(this%map_large_to_uniqueA)))
    if ((minval(this%map_large_to_uniqueA) /= 1) .or. &
      (ace_radial_chem == ACE_CHEM_RADIAL_RALF .and. maxval(this%map_large_to_uniqueA) /= num_unique) .or. &
      ((ace_radial_chem == ACE_CHEM_RADIAL_HSVD .or. ace_radial_chem == ACE_CHEM_RADIAL_BLOCK_HSVD .or. &
        ace_radial_chem == ACE_CHEM_RADIAL_RANDPROJ .or. ace_radial_chem == ACE_CHEM_RADIAL_CHEMMAP_HSVD) .and. &
        maxval(this%map_large_to_uniqueA)*acekmax_nu /= num_unique)) then
        call log_critical("ML: ACE body "//vtoa(this%nu)//" the map_large_to_uniqueA is not correct")
        call log_critical("This is a critical issue. Please contact the developers.")
        call log_critical("The subroutine "//NAMECURRENT//" should be re-writnen") 
        call log_critical("---- program will now stop ----")
        call mld_mpi_abort("stop: coffee and re-write")  
    end if
    !debug! do i = 1, M 
    !debug!    do j = 1, this%nu
    !debug!       ix = (i-1)*this%nu + j
    !debug!       write(*,*) ix, mat_mubold(i, j), mat_nbold(i, j), mat_lbold(i, j), mat_mbold(i, j), this%map_large_to_uniqueA(ix)
    !debug!       if (this%map_large_to_uniqueA(ix) == 0) then 
    !debug!           stop '<<<<<<<<<<<<<<BUDDY TI-AI LUAT-O IN FREZA! >>>>>>>>>>>>>>>'
    !debug!       end if     
    !debug!     end do   
    !debug! end do 
    !debug! 
    !debug! write(*,*) '-------------------------------------'
    !debug! do ii = 1, num_unique
    !debug!     write(*,*)  ii, this%uniqueA_tuples(1, ii), this%uniqueA_tuples(2, ii), this%uniqueA_tuples(3, ii), this%uniqueA_tuples(4, ii)
    !debug! end do 
    !debug! stop '<<<<<<<<<<<<<<>>>>>>>>>>>>>>>'
    _MLD_END_
    end subroutine check_uniqueA_tuples


    subroutine dump_dico_lb0(fname, nu, llmax, dico_lb0)
      use permutation_module, only: permutation_object
      implicit none
      character(len=*), intent(in) :: fname
      integer,          intent(in) :: nu, llmax
      type(type_base_dico_ll), intent(in) :: dico_lb0(:)
      integer :: u, ii, j, k, dimr, dimc
      integer :: dimr_m0, dimc_m0
      type(permutation_object) :: perm_obj
      integer, allocatable :: idv(:), perm(:), perms(:,:)

      open(newunit=u, file=fname, status='replace', action='write', form='formatted')
      write(u,'(I0,1X,I0,1X,I0)') nu, llmax, size(dico_lb0)

      call perm_obj%initialize(nu)
      allocate(idv(nu))
      do j = 1, nu
        idv(j) = j
      end do

      do ii = 1, size(dico_lb0)
        write(u,'(A,1X,I0)') 'idx', ii

        write(u,'(A,1X)', advance='no') 'lbold'
        call f_write_ints_line(u, dico_lb0(ii)%lbold, nu)

        if (allocated(dico_lb0(ii)%Lmat)) then
          dimr = size(dico_lb0(ii)%Lmat, 1)
          dimc = size(dico_lb0(ii)%Lmat, 2)
        else
          dimr = 0
          dimc = 0
        end if
        write(u,'(A,1X,I0,1X,I0)') 'dimLmat', dimr, dimc
        if (dimr > 0) then
          write(u,'(A)') 'Lmat'
          do j = 1, dimr
            call f_write_ints_line(u, dico_lb0(ii)%Lmat(j, :), dimc)
          end do
        end if

        ! Dump M0mat in the same format as Python
        if (allocated(dico_lb0(ii)%M0mat)) then
          dimr_m0 = size(dico_lb0(ii)%M0mat, 1)
          dimc_m0 = size(dico_lb0(ii)%M0mat, 2)
        else
          dimr_m0 = 0
          dimc_m0 = 0
        end if
        write(u,'(A,1X,I0,1X,I0)') 'dimM0mat', dimr_m0, dimc_m0
        if (dimr_m0 > 0) then
          write(u,'(A)') 'M0mat'
          do j = 1, dimr_m0
            call f_write_ints_line(u, dico_lb0(ii)%M0mat(j, :), dimc_m0)
          end do
        end if

        if (allocated(dico_lb0(ii)%list_permutations)) then
          k = size(dico_lb0(ii)%list_permutations)
          write(u,'(A,1X,I0)') 'raw_perm_codes', k
          if (k > 0) then
            call f_write_ints_line(u, dico_lb0(ii)%list_permutations, k)
            allocate(perms(k, nu))
            allocate(perm(nu))
            do j = 1, k
              call perm_obj%apply_permutation(idv, perm, dico_lb0(ii)%list_permutations(j))
              perms(j, :) = perm(:)
            end do
            write(u,'(A)') 'perms'
            do j = 1, k
              call f_write_ints_line(u, perms(j, :), nu)
            end do
            deallocate(perms, perm)
          else
            write(u,*) ''
            write(u,'(A)') 'perms'
          end if
        else
          write(u,'(A,1X,I0)') 'raw_perm_codes', 0
          write(u,*) ''
          write(u,'(A)') 'perms'
        end if
      end do

      deallocate(idv)
      close(u)
    end subroutine dump_dico_lb0


    subroutine dump_dico_nb(fname, nu, nmax, dico_nb)
      use permutation_module, only: permutation_object
      implicit none
      character(len=*), intent(in) :: fname
      integer, intent(in)          :: nu, nmax
      type(type_base_dico_nn), intent(in) :: dico_nb(:)
      integer :: u, ii, j, k
      type(permutation_object) :: perm_obj
      integer, allocatable :: idv(:), perm(:), perms(:,:)
      open(newunit=u, file=fname, status='replace', action='write', form='formatted')
      write(u,'(I0,1X,I0,1X,I0)') nu, nmax, size(dico_nb)
      call perm_obj%initialize(nu)
      allocate(idv(nu))
      do j = 1, nu
        idv(j) = j
      end do
      do ii = 1, size(dico_nb)
        write(u,'(A,1X,I0)') 'idx', ii
        write(u,'(A,1X)', advance='no') 'nbold'
        call f_write_ints_line(u, dico_nb(ii)%nbold, nu)
        if (allocated(dico_nb(ii)%list_permutations)) then
          k = size(dico_nb(ii)%list_permutations)
          write(u,'(A,1X,I0)') 'raw_perm_codes', k
          if (k > 0) then
            call f_write_ints_line(u, dico_nb(ii)%list_permutations, k)
            allocate(perms(k, nu))
            allocate(perm(nu))
            do j = 1, k
              call perm_obj%apply_permutation(idv, perm, dico_nb(ii)%list_permutations(j))
              perms(j, :) = perm(:)
            end do
            write(u,'(A)') 'perms'
            do j = 1, k
              call f_write_ints_line(u, perms(j, :), nu)
            end do
            deallocate(perms, perm)
          else
            write(u,*) ''
            write(u,'(A)') 'perms'
          end if
        else
          write(u,'(A,1X,I0)') 'raw_perm_codes', 0
          write(u,*) ''
          write(u,'(A)') 'perms'
        end if
      end do
      deallocate(idv)
      close(u)
    end subroutine dump_dico_nb


    subroutine dump_dico_mu_raw(fname, nu, mumax, dico_mu)
      use permutation_module, only: permutation_object
      implicit none
      character(len=*), intent(in) :: fname
      integer, intent(in) :: nu, mumax
      type(type_base_dico_mu), intent(in) :: dico_mu(:)
      integer :: u, ii, j, k
      type(permutation_object) :: perm_obj
      integer, allocatable :: idv(:), perm(:), perms(:,:)
    
      open(newunit=u, file=fname, status='replace', action='write', form='formatted')
      write(u,'(I0,1X,I0,1X,I0)') nu, mumax, size(dico_mu)
    
      call perm_obj%initialize(nu)
      allocate(idv(nu))
      do j = 1, nu
        idv(j) = j
      end do
    
      do ii = 1, size(dico_mu)
        write(u,'(A,1X,I0)') 'idx', ii
        write(u,'(A,1X)', advance='no') 'mubold'
        call f_write_ints_line(u, dico_mu(ii)%mubold, nu)
        write(u,'(A,1X)', advance='no') 'mustamp'
        call f_write_ints_line(u, dico_mu(ii)%mustamp, mumax)
    
        if (allocated(dico_mu(ii)%list_permutations)) then
          k = size(dico_mu(ii)%list_permutations)
          write(u,'(A,1X,I0)') 'raw_perm_codes', k
          if (k > 0) then
            call f_write_ints_line(u, dico_mu(ii)%list_permutations, k)
            allocate(perms(k, nu))
            allocate(perm(nu))
            do j = 1, k
              call perm_obj%apply_permutation(idv, perm, dico_mu(ii)%list_permutations(j))
              perms(j, :) = perm(:)
            end do
            write(u,'(A)') 'perms'
            do j = 1, k
              call f_write_ints_line(u, perms(j, :), nu)
            end do
            deallocate(perms, perm)
          else
            write(u,*) ''
            write(u,'(A)') 'perms'
          end if
        else
          write(u,'(A,1X,I0)') 'raw_perm_codes', 0
          write(u,*) ''
          write(u,'(A)') 'perms'
        end if
      end do
    
      deallocate(idv)
      close(u)
    end subroutine dump_dico_mu_raw

    subroutine f_write_ints_line(u, arr, n)
      implicit none
      integer, intent(in) :: u, n
      integer, intent(in) :: arr(n)
      integer :: j
      do j = 1, n
        write(u,'(I0)', advance='no') arr(j)
        if (j < n) write(u,'(A)', advance='no') ' '
      end do
      write(u,*)
    end subroutine f_write_ints_line

    subroutine sort_int_asc(arr, n)
      implicit none
      integer, intent(in) :: n
      integer, intent(inout) :: arr(n)
      integer :: i, j, key
      do i = 2, n
        key = arr(i); j = i - 1
        do while (j >= 1 .and. arr(j) > key)
          arr(j+1) = arr(j); j = j - 1
        end do
        arr(j+1) = key
      end do
    end subroutine sort_int_asc


    subroutine dump_cg_vector_file(fname, ntype, cg_vector)
      implicit none
      character(len=*), intent(in) :: fname
      integer,          intent(in) :: ntype
      real(8),          intent(in) :: cg_vector(:,:,:,:,:,:)
    
      integer :: u
      integer :: lb1,ub1,lb2,ub2,lb3,ub3,lb4,ub4,lb5,ub5,lb6,ub6
      integer :: e1,e2,e3,e4,e5,e6
      integer :: j1_max,j2_max,j0_max
      integer :: j1,j2,j,m1,m2,m
    
      lb1=lbound(cg_vector,1); ub1=ubound(cg_vector,1); e1=ub1-lb1+1
      lb2=lbound(cg_vector,2); ub2=ubound(cg_vector,2); e2=ub2-lb2+1
      lb3=lbound(cg_vector,3); ub3=ubound(cg_vector,3); e3=ub3-lb3+1
      lb4=lbound(cg_vector,4); ub4=ubound(cg_vector,4); e4=ub4-lb4+1
      lb5=lbound(cg_vector,5); ub5=ubound(cg_vector,5); e5=ub5-lb5+1
      lb6=lbound(cg_vector,6); ub6=ubound(cg_vector,6); e6=ub6-lb6+1
    
      ! Robust: derive j*_max from m-dim extents (2*j+1)
      if (mod(e2,2)==0 .or. mod(e4,2)==0 .or. mod(e6,2)==0) then
         write(*,*) 'dump_cg_vector_file: even extent in m-dim(s):', e2,e4,e6
      end if
      j1_max = (e2 - 1)/2
      j2_max = (e4 - 1)/2
      j0_max = (e6 - 1)/2
      !debug write(*,*) 'dump_cg_vector_file: derived j_max:', e1, e2, e3, e4, e5, e6
      !debug write(*,*) 'dump_cg_vector_file: derived j_max:', size(cg_vector,1)-1, size(cg_vector,3)-1, size(cg_vector,5)-1
      !debug write(*,*) 'dump_cg_vector_file: derived j_max:', j1_max, j2_max, j0_max
      
    
      ! Optional sanity checks for j-dims (should be j_max+1)
      if (e1 /= j1_max+1) write(*,*) 'warn: dim1 extent=',e1,' expected=',j1_max+1
      if (e3 /= j2_max+1) write(*,*) 'warn: dim3 extent=',e3,' expected=',j2_max+1
      if (e5 /= j0_max+1) write(*,*) 'warn: dim5 extent=',e5,' expected=',j0_max+1
    
      open(newunit=u, file=fname, status='replace', action='write')
      write(u,'(I0,1X,I0,1X,I0,1X,I0)') j1_max, j2_max, j0_max, ntype
    
      do j1 = 0, j1_max
        do j2 = 0, j2_max
          do j  = 0, j0_max
            do m1 = -j1, j1, ntype
              do m2 = -j2, j2, ntype
                do m  = -j , j , ntype
                  write(u,'(6(I0,1X),ES24.16E3)') j1, m1, j2, m2, j, m, &
         &          cg_vector( j1 + lb1, m1 + (j1_max + lb2),  &
         &                     j2 + lb3, m2 + (j2_max + lb4),  &
         &                     j  + lb5, m  + (j0_max + lb6) )
                end do
              end do
            end do
          end do
        end do
      end do
      close(u)
    end subroutine dump_cg_vector_file    

  subroutine dump_dico_mub_nb_lb0_file(this, fname)
    class(type_base_cnlm), intent(in) :: this
    character(len=*), intent(in) :: fname
    integer :: u, i, n
    if (.not. allocated(this%dico_mub_nb_lb0)) then
      write(*,*) 'dump_dico_mub_nb_lb0_file: not allocated'
      return
    end if
    n = size(this%dico_mub_nb_lb0)
    open(newunit=u, file=fname, status='replace', action='write', form='formatted')
    write(u,'(I0)') n
    do i = 1, n
      write(u,'(I0,1X,I0,1X,I0,1X,I0)') &
        this%dico_mub_nb_lb0(i)%idx_mub, &
        this%dico_mub_nb_lb0(i)%idx_nb,  &
        this%dico_mub_nb_lb0(i)%idx_lb0, &
        this%dico_mub_nb_lb0(i)%idx_cg
    end do
    close(u)
  end subroutine dump_dico_mub_nb_lb0_file

    subroutine dump_dico_cg(this, filename)
      ! Dump the dicoCG dictionary to an ASCII file matching the Python layout.
      class(type_base_cnlm), intent(in) :: this
      character(len=*), intent(in) :: filename
      integer :: unit, ios
      integer :: i, j

      if (.not. allocated(this%dicoCG)) return

      open(newunit=unit, file=trim(filename), status='replace', action='write', &
           form='formatted', iostat=ios)
      if (ios /= 0) then
        write(*,'(A,1X,A,1X,I0)') 'dump_dico_cg: unable to open', trim(filename), ios
        return
      end if

      write(unit,'(A,1X,I0)') 'dim_baseCG', size(this%dicoCG)

      do i = 1, size(this%dicoCG)
        write(unit,'(A,1X,I0)') 'entry', i
        write(unit,'(A,1X,I0)') 'dimr', this%dicoCG(i)%dimr_cg
        write(unit,'(A,1X,I0)') 'dimc', this%dicoCG(i)%dimc_cg
        write(unit,'(A)') 'cg'
        if (allocated(this%dicoCG(i)%cg)) then
          if (this%dicoCG(i)%dimr_cg > 0 .and. this%dicoCG(i)%dimc_cg > 0) then
            do j = 1, this%dicoCG(i)%dimr_cg
              write(unit,'(*(F24.16,1X))') this%dicoCG(i)%cg(j,1:this%dicoCG(i)%dimc_cg)
            end do
          end if
        else
           write(unit,'(A)') 'cg_not_allocated'
        end if
        write(unit,'(A)') 'end_entry'
      end do

      close(unit)
    end subroutine  dump_dico_cg

    subroutine dump_dico_lb0_cg(this, filename)
      class(type_base_cnlm), intent(in) :: this
      character(len=*), intent(in) :: filename
      integer :: unit, ios
      integer :: i, j

      if (.not. allocated(this%dico_lb0)) return

      open(newunit=unit, file=trim(filename), status='replace', action='write', &
           form='formatted', iostat=ios)
      if (ios /= 0) then
        write(*,'(A,1X,A,1X,I0)') 'dump_dico_lb0_cg: unable to open', trim(filename), ios
        return
      end if

      write(unit,'(A,1X,I0)') 'dim_dico_lb0', size(this%dico_lb0)

      do i = 1, size(this%dico_lb0)
        write(unit,'(A,1X,I0)') 'entry', i
        if (allocated(this%dico_lb0(i)%cg)) then
           write(unit,'(A,1X,I0)') 'dimr', size(this%dico_lb0(i)%cg, 1)
           write(unit,'(A,1X,I0)') 'dimc', size(this%dico_lb0(i)%cg, 2)
           write(unit,'(A)') 'cg'
           do j = 1, size(this%dico_lb0(i)%cg, 1)
             write(unit,'(*(F24.16,1X))') this%dico_lb0(i)%cg(j, :)
           end do
        else
           write(unit,'(A)') 'cg_not_allocated'
        end if
        write(unit,'(A)') 'end_entry'
      end do

      close(unit)
    end subroutine dump_dico_lb0_cg

    subroutine dump_dicoB(this, filename)
      class(type_base_cnlm), intent(in) :: this
      character(len=*), intent(in) :: filename
      integer :: unit, ios, i

      if (.not. allocated(this%dicoB)) return

      open(newunit=unit, file=trim(filename), status='replace', action='write', &
           form='formatted', iostat=ios)
      if (ios /= 0) then
        write(*,'(A,1X,A,1X,I0)') 'dump_dicoB: unable to open', trim(filename), ios
        return
      end if

      write(unit,'(A,1X,I0)') 'dim_baseB', size(this%dicoB)
      do i = 1, size(this%dicoB)
        write(unit,'(5(I0,1X))') this%dicoB(i)%idx_nb, this%dicoB(i)%idx_mub, this%dicoB(i)%idx_lb0, &
                                  this%dicoB(i)%idx_munblb0, this%dicoB(i)%idx_Li
      end do

      close(unit)
    end subroutine dump_dicoB


    subroutine dump_uniqueA(this, filename)
      class(type_base_cnlm), intent(in) :: this
      character(len=*), intent(in) :: filename
      integer :: unit, ios, dim_unique, dim_map, j

      open(newunit=unit, file=trim(filename), status='replace', action='write', &
           form='formatted', iostat=ios)
      if (ios /= 0) then
        write(*,'(A,1X,A,1X,I0)') 'dump_uniqueA: unable to open', trim(filename), ios
        return
      end if

      if (allocated(this%uniqueA_tuples)) then
        dim_unique = size(this%uniqueA_tuples, 2)
      else
        dim_unique = 0
      end if
      write(unit,'(A,1X,I0)') 'dim_uniqueA', dim_unique
      if (dim_unique > 0) then
        do j = 1, dim_unique
          write(unit,'(*(I0,1X))') this%uniqueA_tuples(:, j)
        end do
      end if

      if (allocated(this%map_large_to_uniqueA)) then
        dim_map = size(this%map_large_to_uniqueA)
      else
        dim_map = 0
      end if
      write(unit,'(A,1X,I0)') 'dim_map_large', dim_map
      if (dim_map > 0) then
        write(unit,'(*(I0,1X))') this%map_large_to_uniqueA
      end if

      close(unit)
    end subroutine dump_uniqueA


end module module_base_cnlm
!>------end_module_base_cnlm-------------------------------------------------------------


subroutine print_ace_descriptor
   use mld_mpi, only: mld_rank
   use module_ace_desc, only: base_params, ace_numax
   use module_base_cnlm, only: base_cnlm
   !use module_chemical_species, only: fix_ch_elements
   use module_tuples_for_given_lbold, only :array_to_string
   
   implicit none 
   integer :: ii, iib,  dim_mub, dim_nb, dim_l0b, dim_Lb, dim_Li, dim_nl0, dim_B, jj, tmpdim_Lb, idx_lb0, mumax, idx_mub 
   integer, allocatable, dimension(:) :: dim_basis_per_element

   integer :: tmp_n, tmp_l, tmp_n_max, tmp_n_min, tmp_l_max, tmp_l_min, ll, unique_mubold, iunit, idx_nb, icount, degree1
   integer :: iunit2, idx_cg  

   ! Printing header
   if (mld_rank == 0) then
   print '("------------------------------------ACE descriptor review------------------------------------------------------")'
   print '("----------------------------------------general ACE info-------------------------------------------------------")'
       
   print '("Order nmax  lmax  mumax  dim_mu  dim_n  dim_l0  dim_Lb  dim_Li  dim_nl0 dim_B     r_cut_out  r_cut_in   lambda")'
   ! Loop through each ACE descriptor and print
   do ii = 1, ace_numax
       
       if (base_params(ii)%active) then
          if (ii == 1) then 
              dim_mub = base_params(ii)%mumax
              dim_nb = base_params(ii)%nmax + 1 
              dim_l0b =  1
              dim_Lb = dim_l0b * dim_nb * dim_mub 
              dim_Li = dim_l0b * dim_nb * dim_mub 
              dim_nl0 = dim_l0b * dim_nb * dim_mub 
              dim_B =   base_params(ii)%mumax * (base_params(ii)%nmax + 1)               
          else
              dim_mub = size(base_cnlm(ii)%dico_mu)
              dim_nb = size(base_cnlm(ii)%dico_nb)
              dim_l0b = size(base_cnlm(ii)%dico_lb0)
              dim_B = size(base_cnlm(ii)%dicoB)

              dim_Lb = 0
              do jj = 1, size(base_cnlm(ii)%dico_mub_nb_lb0)
                  idx_lb0 = base_cnlm(ii)%dico_mub_nb_lb0(jj)%idx_lb0
                  dim_Lb = dim_Lb + size(base_cnlm(ii)%dico_lb0(idx_lb0)%cg, 2)
              end do     

              dim_Li = 0
              do jj = 1, size(base_cnlm(ii)%dico_mub_nb_lb0)
                  idx_cg = base_cnlm(ii)%dico_mub_nb_lb0(jj)%idx_cg
                  tmpdim_Lb = size(base_cnlm(ii)%dicoCG(idx_cg)%cg, 2)
                  dim_Li = dim_Li + tmpdim_Lb
              end do
              dim_nl0 = size(base_cnlm(ii)%dico_mub_nb_lb0)
          end if 

           print '(I5, 3I5,  7I8, 3F12.3)', ii, base_params(ii)%nmax, base_params(ii)%lmax, &
                 base_params(ii)%mumax, &
                 dim_mub, dim_nb, dim_l0b, &
                 dim_Lb, dim_Li, dim_nl0, dim_B, &
                 base_params(ii)%r_cut_out, &
                 base_params(ii)%r_cut_in, base_params(ii)%lambda
       else
           print '(I5, 3A5, 7A8, 3A12)', ii, " -", " -", " -", " -", " -", " -", " -", " -", " -", " -", " -", " -"
       end if
   end do

   print '("------------------------------------non zero ACE basis per element------------------------------------------------------")'
  
  
   print '("Order  dim_B  non-zero  uni_nb_max  uni_nb_min  uni_nb_mean  uni_mub_max  uni_mub_min  uni_mub_mean")'
  
  do ii = 2, ace_numax
     if (base_params(ii)%active) then
        mumax=base_params(ii)%mumax
        ! dim_basis_per_element = 0
        if (allocated(dim_basis_per_element)) deallocate(dim_basis_per_element)      
        allocate(dim_basis_per_element(mumax))
      end if 
   end do 

   do ii = 2, ace_numax
     if (base_params(ii)%active) then
       dim_basis_per_element(:) = 0
       do jj = 1, size(base_cnlm(ii)%dico_mub_nb_lb0)
          idx_mub = base_cnlm(ii)%dico_mub_nb_lb0(jj)%idx_mub
          !debug write(*,*) 'mustamp!!!!!!!!!!!!!!!!!!!!!!!!!!!', base_cnlm(ii)%dico_mu(idx_mub)%mustamp(:)
          do iib = 1, mumax 
            if (base_cnlm(ii)%dico_mu(idx_mub)%mustamp(iib) == 1) then
              dim_basis_per_element(iib) = dim_basis_per_element(iib) + 1
            end if 
          end do 
       end do 

       dim_B =    size(base_cnlm(ii)%dicoB) 
       tmp_n = 0 
       tmp_l = 0 
       tmp_n_max = 0
       tmp_n_min = 1000000000
       tmp_l_max = 0
       tmp_l_min = 1000000000
       do ll = 1, size(base_cnlm(ii)%dico_lb0)
          unique_mubold = base_cnlm(ii)%dico_lb0(ll)%max_storage_nmu / base_cnlm(ii)%dico_lb0(ll)%unique_nbold
          ! unique_mubold = size(base_cnlm(ii)%dico_lb0(ll)%storage_nmu,2) / base_cnlm(ii)%dico_lb0(ll)%unique_nbold
          tmp_n = tmp_n + base_cnlm(ii)%dico_lb0(ll)%unique_nbold
          tmp_l = tmp_l + unique_mubold
          if (unique_mubold > tmp_n_max) tmp_n_max = unique_mubold
          if (unique_mubold < tmp_n_min) tmp_n_min = unique_mubold
          if (base_cnlm(ii)%dico_lb0(ll)%unique_nbold > tmp_l_max) tmp_l_max = base_cnlm(ii)%dico_lb0(ll)%unique_nbold
          if (base_cnlm(ii)%dico_lb0(ll)%unique_nbold < tmp_l_min) tmp_l_min = base_cnlm(ii)%dico_lb0(ll)%unique_nbold
       end do   
       tmp_n = int(tmp_n / size(base_cnlm(ii)%dico_lb0))
       tmp_l = int(tmp_l / size(base_cnlm(ii)%dico_lb0))
       write(*,'(i4, i7, i8, 3i12, 3i13)') ii, dim_B, sum(dim_basis_per_element)/mumax, tmp_n_max, tmp_n_min, tmp_n, tmp_l_max, tmp_l_min, tmp_l
       !$! do iib =1, mumax
       !$!  write(*,'(I5, "  ", a, "  ", I10)') iib, fix_ch_elements(iib), dim_basis_per_element(iib) 
       !$! end do 
    else 
        print '(I4, A7, A8, 3A12, 3A13)', ii, " -", " -", " -", " -", " -", " -", " -", " -"

    end if 
   end do 


   print '("--------------------------------------------------------------------------------------------------------------")'

   open(newunit=iunit, file='base_B_ace.info', status='unknown')
   open(newunit=iunit2, file='base_munl_ace.info', status='unknown')

    do ii = 2, ace_numax
       
      if (base_params(ii)%active) then
      
        icount = 0
        do jj = 1, size(base_cnlm(ii)%dico_mub_nb_lb0)
          idx_lb0 = base_cnlm(ii)%dico_mub_nb_lb0(jj)%idx_lb0
          idx_mub = base_cnlm(ii)%dico_mub_nb_lb0(jj)%idx_mub
          idx_nb = base_cnlm(ii)%dico_mub_nb_lb0(jj)%idx_nb
          degree1 = sum(base_cnlm(ii)%dico_lb0(idx_lb0)%lbold(:)) + sum(base_cnlm(ii)%dico_nb(idx_nb)%nbold(:))
          write(iunit2,'(I6, I8, A)') jj,  degree1, " l=" // trim(array_to_string(base_cnlm(ii)%dico_lb0(idx_lb0)%lbold(:))) // &
                    " n=" // trim(array_to_string(base_cnlm(ii)%dico_nb(idx_nb)%nbold(:))) // &
                    " mu=" // trim(array_to_string(base_cnlm(ii)%dico_mu(idx_mub)%mubold(:))) 
                    
 
          do ll = 1, size(base_cnlm(ii)%dico_lb0(idx_lb0)%Lmat, 1)
            icount = icount + 1
            !degree1 = sum(base_cnlm(ii)%dico_lb0(idx_lb0)%lbold(:)) + sum(base_cnlm(ii)%dico_nb(idx_nb)%nbold(:))
            write(iunit,'(3I6, I8, A)') icount, jj, ll, degree1, " l=" // trim(array_to_string(base_cnlm(ii)%dico_lb0(idx_lb0)%lbold(:))) // &
                                " n=" // trim(array_to_string(base_cnlm(ii)%dico_nb(idx_nb)%nbold(:))) // &
                                " mu=" // trim(array_to_string(base_cnlm(ii)%dico_mu(idx_mub)%mubold(:))) // & 
                                " Lb=" // trim(array_to_string(base_cnlm(ii)%dico_lb0(idx_lb0)%Lmat(ll, :)))
          end do                       

        end do !! jj
      end if 
    end do
    close(iunit, status='keep') 
   end if ! mld_rank 
end subroutine print_ace_descriptor



