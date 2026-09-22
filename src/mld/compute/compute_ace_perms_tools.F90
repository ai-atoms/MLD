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


!>------begin_tuple_module---------------------------------------------------------------
module tuple_module
  implicit none
  type vector_tuple
      integer, allocatable :: values(:)
  end type vector_tuple

  contains

  subroutine generate_combinations(vectors, all_combinations)
    type(vector_tuple), allocatable, intent(in) :: vectors(:)
    integer, allocatable, intent(out) :: all_combinations(:,:)
    integer, allocatable :: indices(:)
    integer :: total_combinations, i, j, nvec
    nvec = size(vectors)
    total_combinations = 1
    allocate(indices(nvec))

    ! Calculate total combinations and initialize indices
    do i = 1, nvec
        total_combinations = total_combinations * size(vectors(i)%values)
        indices(i) = 1
    end do

    allocate(all_combinations(total_combinations, nvec))

    ! Iteratively generate all combinations
    do i = 1, total_combinations
        ! Fill the current combination
        do j = 1, nvec
            all_combinations(i, j) = vectors(j)%values(indices(j))
        end do
        ! Update indices for next combination
        call update_indices(indices, vectors)
    end do

  end subroutine generate_combinations


  subroutine update_indices(indices, vectors)
    integer, intent(inout) :: indices(:)
    type(vector_tuple), intent(in) :: vectors(:)
    integer :: i

    ! Iterate from the last index backwards to update
    do i = size(indices), 1, -1
      ! Increment the current index
      indices(i) = indices(i) + 1
      ! Check if the current index exceeds its maximum value
      if (indices(i) > size(vectors(i)%values)) then
          ! Reset the current index and move to the previous one
          indices(i) = 1
          ! If this is the first index, then all combinations have been generated
          if (i == 1) exit
      else
          ! Correct index found, no need to adjust previous indices
          exit
      end if
    end do
  end subroutine update_indices

end module tuple_module
!>------end_tuple_module-----------------------------------------------------------------


!>------begin_permutation_module---------------------------------------------------------
module permutation_module
  use mld_logger
  implicit none  
  type permutation_object
    integer, allocatable :: permutations(:,:)
    integer :: nu  ! Number of elements in the permutation
    integer :: factorial_nu 
    logical :: lalready = .false. ! Added to track if the object already
    contains
    procedure :: initialize
    procedure :: generate_permutations
    procedure :: apply_permutation
    procedure :: extract_permutations
    procedure :: check
    procedure :: extract_common_permutations
  end type permutation_object

  contains
  integer function factorial(n) result(f)
    integer, intent(in) :: n
    integer :: i  
    f = 1
    do i = 1, n
        f = f * i
    end do
  end function factorial

  subroutine initialize(this, nu)
    class(permutation_object), intent(inout) :: this
    integer, intent(in) :: nu
    integer, dimension(nu) :: vector
    integer :: perm_index
    integer, dimension(nu) :: vector_copy
    integer :: ii         
    do ii = 1, nu
        vector(ii) = ii
    end do        
    this%nu = nu
    this%factorial_nu = factorial(nu)
    allocate(this%permutations(this%factorial_nu, nu))
    perm_index = 0
    vector_copy = vector  ! Make a copy to pass to the recursive subroutine
    this%lalready = .true.
    call generate_heap_permutations(this, nu, vector_copy, perm_index)
  end subroutine initialize

  subroutine check(this, status)
    class(permutation_object), intent(in) :: this
    logical, intent(inout) :: status
    if (this%nu == -1) then
        !call log_info("Permutation object not initialized")
        status = .false.
    else 
        call log_info("Permutation object already initialized")
        status = .true.
    end if

    if (this%lalready) then
        !call log_warning("Permutation object already initialized")
        status = .false. 
    else 
        !call log_info("Permutation object not initialized")
        status  = .true.
    end if 
  end subroutine check


  subroutine generate_permutations(this, vector, full_permutations)
    class(permutation_object), intent(inout) :: this
    integer, dimension(:,:), allocatable, intent(out) :: full_permutations
    integer :: perm_index
    integer, dimension(:), intent(in) :: vector
    integer :: i 
    integer, dimension(size(vector)) :: vector_out 

    perm_index = 0
    if (allocated(full_permutations)) deallocate(full_permutations)
    allocate(full_permutations(factorial(this%nu), this%nu))
    do i = 1, size(this%permutations,1)
      call apply_permutation(this, vector, vector_out, i)
      full_permutations(i, :) = vector_out 
    end do
  end subroutine generate_permutations

  subroutine apply_permutation(this, tuple_in, tuple_out, index)
    class(permutation_object), intent(inout) :: this
    integer, dimension(:), intent(in) :: tuple_in
    integer, dimension(:), intent(out) :: tuple_out
    integer, intent(in) :: index
    integer :: i
    _NAMECURRENT_("apply_permutation")
    
    if (index > size(this%permutations,1)) then
        call log_critical("Error: index out of range in "//NAMECURRENT//" subroutine ")
        call log_critical("Error: in module permutation_module")
        stop
    end if
    ! Reorder the elements of tuple based on the selected permutation
    do i = 1, size(tuple_in)
        tuple_out(i) = tuple_in (this%permutations(index, i))
    end do
  end subroutine apply_permutation

  recursive subroutine generate_heap_permutations(this, n, a, perm_index)
    ! from here: Heap, B. R. (1963). "Permutations by Interchanges". The Computer Journal. 6 (3): 293–4. 
    !            doi:10.1093/comjnl/6.3.293 
    !            and here  https://en.wikipedia.org/wiki/Heap%27s_algorithm
    ! the main outcome the lexyographic order of permutations
    class(permutation_object), intent(inout) :: this
    integer, intent(in) :: n
    integer, dimension(this%nu), intent(inout) :: a
    integer, intent(inout) :: perm_index
    integer :: i

    if (n == 1) then
      perm_index = perm_index + 1
      this%permutations(perm_index, :) = a
    else
      do i = 1, n
          call generate_heap_permutations(this, n - 1, a, perm_index)
          if (mod(n, 2) == 0) then
              call swap(a(i), a(n))
          else
              call swap(a(1), a(n))
          end if
      end do
    end if
  end subroutine generate_heap_permutations


  !old! recursive subroutine generate_heap_permutations(this, n, a, perm_index)
  !old!   class(permutation_object), intent(inout) :: this
  !old!   integer, intent(in) :: n
  !old!   integer, dimension(this%nu), intent(inout) :: a
  !old!   integer, intent(inout) :: perm_index
  !old!   integer :: i
  !old!   integer, dimension(this%nu) :: temp_a
  !old! 
  !old!   if (n == 1) then
  !old!     perm_index = perm_index + 1
  !old!     this%permutations(perm_index, :) = a
  !old!   else
  !old!     do i = 1, n
  !old!         call generate_heap_permutations(this, n - 1, a, perm_index)
  !old!         if (n > 2) then
  !old!           temp_a = a
  !old!           if (mod(n, 2) == 0) then
  !old!               call swap(temp_a(1), temp_a(n))
  !old!           else
  !old!               call swap(temp_a(i), temp_a(n))
  !old!           end if
  !old!           a(this%nu-n+1:this%nu) = temp_a(this%nu-n+1:this%nu)
  !old!         else
  !old!           if (mod(n, 2) == 0) then
  !old!               call swap(a(1), a(2))
  !old!           else
  !old!               call swap(a(i), a(n))
  !old!           end if
  !old!         end if
  !old!     end do
  !old!   end if
  !old! end subroutine generate_heap_permutations

  subroutine swap(x, y)
    integer, intent(inout) :: x, y
    integer :: temp  
    temp = x
    x = y
    y = temp
  end subroutine swap

  subroutine extract_permutations(this, vector, vect_list)
    class(permutation_object), intent(in) :: this
    integer, dimension(:), intent(in) :: vector
    integer, dimension(:), allocatable, intent(out) :: vect_list
    integer :: i, dim_list 
    integer, dimension(:), allocatable  :: tmp_vect_list


    allocate(tmp_vect_list(size(this%permutations, 1)))
    tmp_vect_list(:) = 0

    do i = 1, size(this%permutations, 1)
        if (all(vector == vector(this%permutations(i, :)))) then
            tmp_vect_list(i) = 1
        end if
    end do

    dim_list = sum(tmp_vect_list)
    if (allocated(vect_list)) deallocate(vect_list)
    allocate(vect_list(dim_list))
    ! Use pack to extract indices where vect_list is 1
    vect_list = pack([(i, i=1, size(tmp_vect_list))], tmp_vect_list == 1)
  
  end subroutine extract_permutations

  subroutine extract_common_permutations(this, vect_list1, vect_list2, vect_list_out)
    class(permutation_object), intent(in) :: this
    integer, dimension(:), allocatable, intent(out) :: vect_list_out
    integer, dimension(:), intent(in) :: vect_list1, vect_list2
    !local
    integer :: s1, s2, dim_max, icount, is1, is2  
    dim_max = this%factorial_nu
    s1 = size(vect_list1)
    s2 = size(vect_list2)

    ! Special case where vectors are identical
    if ( (s1 == s2) .and. (dim_max == s1) ) then 
      if (allocated(vect_list_out)) deallocate(vect_list_out)
      allocate(vect_list_out(dim_max))
      vect_list_out = vect_list1
      return 
    end if 

    ! Use two-pointer technique for general case
    icount = 0
    is1 = 1
    is2 = 1

    ! First pass to count the common elements
    do while (is1 <= s1 .and. is2 <= s2)
      if (vect_list1(is1) < vect_list2(is2)) then
        is1 = is1 + 1
      else if (vect_list1(is1) > vect_list2(is2)) then
        is2 = is2 + 1
      else
        icount = icount + 1
        is1 = is1 + 1
        is2 = is2 + 1
      end if
    end do

    if (allocated(vect_list_out)) deallocate(vect_list_out)
    allocate(vect_list_out(icount))

    ! Reset pointers for second pass to fill the output vector
    icount = 0
    is1 = 1
    is2 = 1

    do while (is1 <= s1 .and. is2 <= s2)
      if (vect_list1(is1) < vect_list2(is2)) then
        is1 = is1 + 1
      else if (vect_list1(is1) > vect_list2(is2)) then
        is2 = is2 + 1
      else
        icount = icount + 1
        vect_list_out(icount) = vect_list1(is1)
        is1 = is1 + 1
        is2 = is2 + 1
      end if
    end do
end subroutine extract_common_permutations

!old! subroutine extract_common_permutations_old(this, vect_list1, vect_list2, vect_list_out)
!old!   class(permutation_object), intent(in) :: this
!old!   integer, dimension(:), allocatable, intent(out) :: vect_list_out
!old!   integer, dimension(:), intent(in) :: vect_list1, vect_list2
!old!   !local 
!old!   integer :: s1, s2, dim_max, icount, is1, is2  
!old!   dim_max = this%factorial_nu
!old!   s1 = size(vect_list1)
!old!   s2 = size(vect_list2)
!old!   if ( (s1 == s2 ) .and. (dim_max == s1 ) ) then 
!old!     ! if the list is full 
!old!     if (allocated(vect_list_out)) deallocate(vect_list_out) ; allocate(vect_list_out(dim_max))
!old!     vect_list_out = vect_list1
!old!     return 
!old!   end if 
!old!   icount = 0 ! s1*s2 operations ... 
!old!   do is1 = 1, s1 
!old!     do is2 = 1, s2
!old!       if (vect_list1(is1) == vect_list2(is2)) then
!old!         icount = icount + 1  
!old!       end if 
!old!     end do
!old!   end do   
!old!   if (allocated(vect_list_out)) deallocate(vect_list_out) ; allocate(vect_list_out(icount))
!old!   icount = 0 
!old!   do is1 = 1, s1 
!old!     do is2 = 1, s2
!old!       if (vect_list1(is1) == vect_list2(is2)) then
!old!         icount = icount + 1  
!old!         vect_list_out(icount) = vect_list1(is1)
!old!       end if 
!old!     end do
!old!   end do   
!old! end subroutine extract_common_permutations_old

  integer function lehmer_rank_1based(n, perm)
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: perm(n)      ! values 1..n
    integer :: used(n), i, j, idx, rank0
    do i = 1, n
      used(i) = i
    end do
    rank0 = 0
    do i = 1, n
      idx = 0
      do j = 1, n - i + 1
        if (used(j) == perm(i)) then
          idx = j - 1
          exit
        end if
      end do
      rank0 = rank0 * (n - i + 1) + idx
      do j = idx + 1, n - i
        used(j) = used(j + 1)
      end do
    end do
    lehmer_rank_1based = rank0 + 1
  end function lehmer_rank_1based
end module permutation_module
!>------end_permutation_module-----------------------------------------------------------


!>------begin_module_l_space_first-------------------------------------------------------
module module_l_space_first
    use mld_logger
    implicit none 

contains 

    subroutine generate_l_space_ini(nu, lmax, l_space, dim_l)
        integer, intent(in) :: nu, lmax
        integer, dimension(:,:), allocatable, intent(out) :: l_space 
        integer, intent(out) :: dim_l
        integer, dimension(nu) :: current_combination
        integer :: index

        ! Initialize the index and allocate space for the combinations
        index = 0
        allocate(l_space(((lmax + 1)**nu), nu))

        ! Generate combinations recursively
        call generate_combinations(nu, lmax, current_combination, 1, l_space, index)

        ! Update the dimension based on the combinations generated
        dim_l = index
        if (dim_l < size(l_space, 1)) then
            l_space = l_space(1:dim_l, :)
        end if
    end subroutine generate_l_space_ini

    recursive subroutine generate_combinations(nu, lmax, comb, depth, l_space, index)
        integer, intent(in) :: nu, lmax, depth
        integer, dimension(nu), intent(inout) :: comb
        integer, dimension(:, :), intent(inout) :: l_space
        integer, intent(inout) :: index
        integer :: i, sum_l

        if (depth > nu) then
            ! Check if the sum of the components is even
            sum_l = sum(comb)
            if (mod(sum_l, 2) == 0) then
                index = index + 1
                l_space(index, :) = comb
            end if
        else
            do i = 0, lmax
                comb(depth) = i
                call generate_combinations(nu, lmax, comb, depth + 1, l_space, index)
            end do
        end if
    end subroutine generate_combinations

end module module_l_space_first
!>------end_module_l_space_first---------------------------------------------------------

!>------begin_module_n_space_first-------------------------------------------------------
module module_n_space_first
  implicit none
  type type_base_n0 
    integer, dimension(:), allocatable :: nbold
    integer, dimension(:), allocatable :: list_permutations
  end type type_base_n0

  type(type_base_n0), dimension(:), allocatable :: base_n0
  contains 


    subroutine generate_n_space_ini(nu, nvals, n_space, dim_n)
      use permutation_module, only: permutation_object
      integer, intent(in) :: nu, nvals
      integer, dimension(:,:), allocatable, intent(out) :: n_space 
      integer, intent(out) :: dim_n
      integer, dimension(nu) :: current_combination
      integer :: index, ii 
      !integer, dimension(:), allocatable :: list_permutations, another_vector
      !type(permutation_object) :: perm_obj
  
      ! Initialize the index and allocate space for the combinations
      index = 0
      allocate(n_space(((nvals)**nu), nu))  ! Adjusted allocation size
  
      ! Generate combinations recursively
      call generate_combinations(nu, nvals, current_combination, 1, n_space, index)
  
      ! Update the dimension based on the combinations generated
      dim_n = index
      if (dim_n < size(n_space, 1)) then
          n_space = n_space(1:dim_n, :)
      end if
      if (allocated(base_n0)) deallocate(base_n0) ; allocate(base_n0(dim_n))
      do ii = 1, dim_n
          if (allocated(base_n0(ii)%nbold)) deallocate(base_n0(ii)%nbold) ; allocate(base_n0(ii)%nbold(nu))
          base_n0(ii)%nbold = n_space(ii, :)
      end do
      !$! call perm_obj%initialize(nu)
      !$! allocate(another_vector(nu))
      !$! do ii = 1, dim_n 
      !$!     another_vector(:) = base_n0(ii)%nbold(:)
      !$!     call perm_obj%extract_permutations(another_vector, list_permutations)
      !$!     if (allocated(base_n0(ii)%list_permutations))  deallocate(base_n0(ii)%list_permutations)
      !$!     allocate(base_n0(ii)%list_permutations(size(list_permutations)))
      !$!     base_n0(ii)%list_permutations = list_permutations
      !$! end do 

  end subroutine generate_n_space_ini

    subroutine generate_n_space_ini_acek(nu, nvals, n_space, dim_n)
      use permutation_module, only: permutation_object
      integer, intent(in) :: nu, nvals
      integer, dimension(:,:), allocatable, intent(out) :: n_space 
      integer, intent(in) :: dim_n
      integer, dimension(nu) :: current_combination
      integer :: index, ii 
      !integer, dimension(:), allocatable :: list_permutations, another_vector
      !type(permutation_object) :: perm_obj
  
      ! Initialize the index and allocate space for the combinations
      index = 0
      allocate(n_space(dim_n, nu))  ! Adjusted allocation size
  
      if (allocated(base_n0)) deallocate(base_n0) ; allocate(base_n0(dim_n))
      do ii = 1, dim_n
          n_space(ii,:) = ii
          if (allocated(base_n0(ii)%nbold)) deallocate(base_n0(ii)%nbold) ; allocate(base_n0(ii)%nbold(nu))
          base_n0(ii)%nbold = n_space(ii, :)
      end do

  end subroutine generate_n_space_ini_acek


  subroutine find_permutations_for_n(nu, n_space)
    use permutation_module, only: permutation_object
    integer, intent(in) :: nu
    integer, dimension(:,:), intent(in) :: n_space
    integer :: ii 
    !$! character(len=80) :: CHFMT
    integer, dimension(:), allocatable :: list_permutations, another_vector
    type(permutation_object) :: perm_obj
    call perm_obj%initialize(nu)
    allocate(another_vector(nu))
    do ii = 1, size(n_space, 1)
      another_vector = n_space(ii, :)
      call perm_obj%extract_permutations(another_vector, list_permutations)
      if (allocated(base_n0(ii)%list_permutations))  deallocate(base_n0(ii)%list_permutations)
      allocate(base_n0(ii)%list_permutations(size(list_permutations)))
      base_n0(ii)%list_permutations = list_permutations
      !debug print '("Permutations for lbold ", I12, ": ", *(I5, ", "))', ii, list_permutations(:)
    end do
  end subroutine find_permutations_for_n

  recursive subroutine generate_combinations(nu, nvals, comb, depth, n_space, index)
    integer, intent(in) :: nu, nvals, depth
    integer, dimension(nu), intent(inout) :: comb
    integer, dimension(:, :), intent(inout) :: n_space
    integer, intent(inout) :: index
    integer :: i
    if (depth > nu) then
        ! Check if the sum of the components is even
        !sum_n = sum(comb)
        !if (mod(sum_n, 2) == 0) then
            index = index + 1
            n_space(index, :) = comb
        !end if
    else
        do i = 1, nvals  ! Adjusted loop range from 0 to nmax embedded into 1 to nvals= nmax + 1 
            comb(depth) = i
            call generate_combinations(nu, nvals, comb, depth + 1, n_space, index)
        end do
    end if
  end subroutine generate_combinations
end module module_n_space_first
!>------end_module_n_space_first---------------------------------------------------------


!>------begin_module_mu_space_first-------------------------------------------------------
module module_mu_space_first
  use mld_logger
  implicit none
  type type_base_mu 
    integer, dimension(:), allocatable :: mubold
    integer, dimension(:), allocatable :: list_permutations
  end type type_base_mu

  type(type_base_mu), dimension(:), allocatable :: base_mu
  contains 


    subroutine generate_mu_space_ini(nu, mumax, mu_space, dim_mu)
      use permutation_module, only: permutation_object
      integer, intent(in) :: nu, mumax
      integer, dimension(:,:), allocatable, intent(out) :: mu_space 
      integer, intent(out) :: dim_mu
      integer, dimension(nu) :: current_combination
      integer :: index, ii 
  
      ! Initialize the index and allocate space for the combinations
      index = 0
      allocate(mu_space(((mumax)**nu), nu))  ! Adjusted allocation size
  
      ! Generate combinations recursively
      call generate_combinations(nu, mumax, current_combination, 1, mu_space, index)
  
      ! Update the dimension based on the combinations generated
      dim_mu = index
      if (dim_mu < size(mu_space, 1)) then
          mu_space = mu_space(1:dim_mu, :)
      end if
      if (allocated(base_mu)) deallocate(base_mu) ; allocate(base_mu(dim_mu))
      do ii = 1, dim_mu
          if (allocated(base_mu(ii)%mubold)) deallocate(base_mu(ii)%mubold) ; allocate(base_mu(ii)%mubold(nu))
          base_mu(ii)%mubold = mu_space(ii, :)
      end do

  end subroutine generate_mu_space_ini

  subroutine generate_mu_space_ini_acek(nu, mumax, mu_space, dim_mu)
      use permutation_module, only: permutation_object
      integer, intent(in) :: nu, mumax
      integer, dimension(:,:), allocatable, intent(out) :: mu_space 
      integer, intent(in) :: dim_mu
      integer, dimension(nu) :: current_combination
      integer :: index, ii 

      if (dim_mu /= 1)  then
        call log_critical("Error: dim_mu is not 1 in generate_mu_space_ini_acek")
        stop
      end if

  
      ! Initialize the index and allocate space for the combinations
      index = 0
      allocate(mu_space(dim_mu, nu))  ! Adjusted allocation size
  
      if (allocated(base_mu)) deallocate(base_mu) ; allocate(base_mu(dim_mu))
      do ii = 1, dim_mu
          mu_space(ii, :) = ii
          if (allocated(base_mu(ii)%mubold)) deallocate(base_mu(ii)%mubold) ; allocate(base_mu(ii)%mubold(nu))
          base_mu(ii)%mubold = mu_space(ii, :)
      end do

  end subroutine generate_mu_space_ini_acek 

  subroutine generate_mu_space_ini_acek_block(nu, mumax, mu_space, dim_mu)
    integer, intent(in) :: nu, mumax
    integer, dimension(:,:), allocatable, intent(out) :: mu_space
    integer, intent(out) :: dim_mu
    integer :: ii

    ! BLOCK_HSVD: diagonal mu-bold, one channel per species, only diagonal terms
    dim_mu = mumax
    allocate(mu_space(dim_mu, nu))

    if (allocated(base_mu)) deallocate(base_mu) ; allocate(base_mu(dim_mu))
    do ii = 1, dim_mu
      mu_space(ii, :) = ii
      if (allocated(base_mu(ii)%mubold)) deallocate(base_mu(ii)%mubold) ; allocate(base_mu(ii)%mubold(nu))
      base_mu(ii)%mubold = mu_space(ii, :)
    end do
  end subroutine generate_mu_space_ini_acek_block

  subroutine find_permutations_for_mu(nu, mu_space)
    use permutation_module, only: permutation_object
    integer, intent(in) :: nu
    integer, dimension(:,:), intent(in) :: mu_space
    integer :: ii 
    !$! character(len=80) :: CHFMT
    integer, dimension(:), allocatable :: list_permutations, another_vector
    type(permutation_object) :: perm_obj
    call perm_obj%initialize(nu)
    allocate(another_vector(nu))
    do ii = 1, size(mu_space, 1)
      another_vector = mu_space(ii, :)
      call perm_obj%extract_permutations(another_vector, list_permutations)
      if (allocated(base_mu(ii)%list_permutations))  deallocate(base_mu(ii)%list_permutations)
      allocate(base_mu(ii)%list_permutations(size(list_permutations)))
      base_mu(ii)%list_permutations = list_permutations
    end do
  end subroutine find_permutations_for_mu

  recursive subroutine generate_combinations(nu, nmax, comb, depth, mu_space, index)
    integer, intent(in) :: nu, nmax, depth
    integer, dimension(nu), intent(inout) :: comb
    integer, dimension(:, :), intent(inout) :: mu_space
    integer, intent(inout) :: index
    integer :: i
    if (depth > nu) then
        ! Check if the sum of the components is even
        !sum_n = sum(comb)
        !if (mod(sum_n, 2) == 0) then
            index = index + 1
            mu_space(index, :) = comb
        !end if
    else
        do i = 1, nmax  ! Adjusted loop range from 1 to lmax
            comb(depth) = i
            call generate_combinations(nu, nmax, comb, depth + 1, mu_space, index)
        end do
    end if
  end subroutine generate_combinations

  subroutine stamp_permutations_for_mu(nu, mumax, mu_space, stamp_mu_space)
    implicit none 
    integer, intent(in) :: nu, mumax 
    integer, dimension(:,:), intent(in) :: mu_space
    integer, dimension(:,:), allocatable, intent(inout) :: stamp_mu_space
    integer :: ii, jj, imuu 
    !$! character(len=80) :: CHFMT
    integer, dimension(:), allocatable :: another_vector
    logical :: found
    if (allocated(stamp_mu_space)) deallocate(stamp_mu_space) ; allocate(stamp_mu_space(size(mu_space, 1), mumax))

    allocate(another_vector(nu))
    do ii = 1, size(mu_space, 1)
      another_vector = mu_space(ii, :)
      
      stamp_mu_space(ii, :) = 0
      ! Loop over each value of mubold
      do imuu = 1, mumax
        found = .false.
        do jj = 1, nu
            if (another_vector(jj) == imuu) then
                found = .true.
                exit
            end if
        end do
        if (found) stamp_mu_space(ii,imuu) = 1
      end do
    end do !ii 
  end subroutine stamp_permutations_for_mu


end module module_mu_space_first
!>------end_module_mu_space_first---------------------------------------------------------

!>------begin_module_m_space_first--------------------------------------------------------
module module_m_space_first
    implicit none
    integer, allocatable :: M0(:,:)
    integer :: count = 0
    !integer, parameter :: nu = 3 ! Change nu according to your specific problem size
contains

    subroutine generate_M0_space_for_l (lbold, M0sort)
      !----------------------------------------------------------------
      ! input  - lbold : array of l values in R^nu                    !
      ! output - M0sort : array of M0 values in R^(no_of_mbold X nu)  !
      !        - mbold has the restriction \sum m_i =  0              !
      !----------------------------------------------------------------
      integer, dimension(:), allocatable, intent(in) :: lbold
      integer, dimension(:,:), allocatable, intent(inout) :: M0sort
      integer, dimension(size(lbold)) :: currentComb
      integer :: totalCombinations, nu 
      nu = size(lbold)
      ! Calculate total number of combinations and allocate M0
      totalCombinations = product((2*lbold + 1))
      allocate(M0(totalCombinations, nu))
      ! Start recursive combination generation
      call recursiveGenerate(lbold, currentComb, 1, 0)
      
      ! Resize M0 to the number of valid combinations
      if (count > 0) then
            M0 = M0(1:count, :)
        else
            deallocate(M0)
      endif

      if (allocated(M0sort)) deallocate(M0sort) ; allocate(M0sort(count, nu))
      M0sort(:,:) = M0(1:count, :)
      deallocate(M0)
      count = 0

    end subroutine generate_M0_space_for_l

    recursive subroutine recursiveGenerate(lbold, currentComb, index, currentSum)
      integer, dimension(:), intent(in) :: lbold 
      integer, intent(in) :: index, currentSum
      integer, dimension(size(lbold)), intent(inout) :: currentComb
      integer :: i, nu 

      nu = size(lbold)
      if (index > nu) then
          ! Check if the sum is zero
          if (currentSum == sum(currentComb)) then
              count = count + 1
              M0(count, :) = currentComb
          endif
      else
          do i = -lbold(index), lbold(index)
              currentComb(index) = i
              call recursiveGenerate(lbold, currentComb, index + 1, currentSum)
          enddo
      endif
    end subroutine recursiveGenerate
end module module_m_space_first
!>------end_module_m_space_first---------------------------------------------------------

!>------begin_module_tuples_for_given_lbold----------------------------------------------
module  module_tuples_for_given_lbold
  use mld_logger
  implicit none
  contains 
  subroutine  generate_tuples_for_given_lbold(nu, n_max, mu_max, lbold, unique_nbold, storage)
    ! input - nu : body ACE order 
    !       - n_max : maximum value of n
    !       - mu_max : maximum value of mu (number of elements - we consider that starts at 1 
    !       - lbold : array of lbold  
    !--------------------------
    ! output - storage : array of tuples (n, mu) 
    !                     We'll store 2*nu integers per column: (n_1,...,n_nu, mu_1,...,mu_nu)^T
    !--------------------------
    ! 1) Generate all (n, mu) tuples
    ! 2) Compute the canonical representative of (n, mu) given l
    ! 3) Store canonical form in the next column of 'storage'
    ! 4) Sort all columns in 'storage' lexicographically and remove duplicates
    ! The algo has a complexity  (n_max*mu_max)^nu * nu * log(n_max * mu_max) 
    !                instead  of (n_max*mu_max)^(2*nu) the brute force approach
    !--------------------------
    implicit none
    integer, intent(in) :: nu, n_max,  mu_max
    integer, dimension(:), intent(in) :: lbold
    integer, dimension(:,:), allocatable, intent(inout) :: storage
    integer, intent(out) :: unique_nbold

    integer :: i, j
    integer :: total_combinations, idx
    integer, dimension(:), allocatable ::  nbold, mubold
    integer, dimension(:), allocatable ::  new_n, new_mu
    integer :: brute_count, max_storage 
    integer, dimension(:,:), allocatable :: tmp_storage 

    ! Example l:
    allocate(nbold(nu), mubold(nu))
    !debug l = (/0,1,1,2,1,1/)   ! For instance, indices with same l-value can be permuted among themselves

    ! total number of (n, mu) tuples before symmetry/dedup
    total_combinations = (n_max**nu)*(mu_max**nu)

    ! We'll store 2*nu integers per column: (n_1,...,n_nu, mu_1,...,mu_nu)
    allocate(storage(2*nu, total_combinations))

    idx = 0

    !----------------------------------------------------------------------
    ! 1) Generate all (n, mu) tuples
    !----------------------------------------------------------------------
    do i = 1, n_max**nu
        call decode_tuple(i, n_max, nu, nbold)   ! decode i-th combination into n-array
        do j = 1, mu_max**nu
            call decode_tuple(j, mu_max, nu, mubold) ! decode j-th combination into mu-array

            ! 2) Compute the canonical representative of (n, mu) given l
            call canonical_form(lbold, nbold, mubold, new_n, new_mu)

            ! 3) Store canonical form in the next column of 'storage'
            idx = idx + 1
            storage(1:nu, idx)     = new_n
            storage(nu+1:2*nu,idx) = new_mu
        end do
    end do

    brute_count = total_combinations
    !debug print *, "Finished generating all tuples => ", brute_count, " total."

    !----------------------------------------------------------------------
    ! 4) Sort all columns in 'storage' lexicographically and remove duplicates
    !----------------------------------------------------------------------
    call sort_tuples(storage, idx, 2*nu)
    call remove_duplicates(storage, idx, 2*nu)

    ! idx now is the count of unique canonical tuples
    !debug call log_info("ML: ACE body "//vtoa(nu)//" unique canonical (l, n, mu) vs brute:  "//vtoa(idx)//" "//vtoa(brute_count))
    max_storage = idx
    allocate(tmp_storage(2*nu, idx))
    tmp_storage(:,:) = storage(:,1:idx)
    deallocate(storage)
    allocate(storage(2*nu, idx))
    storage(:,:) = tmp_storage(:,:)
    deallocate(tmp_storage)
    unique_nbold = count_distinct_n(storage, nu)
    !----------------------------------------------------------------------
    ! 5) (Optional) Print results to file unit=777
    !----------------------------------------------------------------------
    !debug! do i = 1, idx
    !debug!     write(777,'(I12, A)') i, " l=" // trim(array_to_string(lbold)) // &
    !debug!                           " n=" // trim(array_to_string(storage(1:nu,i))) // &
    !debug!                           " mu=" // trim(array_to_string(storage(nu+1:2*nu,i)))
    !debug! end do

    deallocate(nbold, mubold, new_n, new_mu)
  end subroutine generate_tuples_for_given_lbold

  !======================================================================
  subroutine decode_tuple(num, base, length, arr)
      ! Converts integer 'num' into an array of length 'length'
      ! in base 'base', with each component in [1..base].
      implicit none
      integer, intent(in)    :: num, base, length
      integer, intent(out)   :: arr(:)
      integer               :: tmp, i
      tmp = num - 1
      do i = length, 1, -1
          arr(i) = mod(tmp, base) + 1
          tmp    = tmp / base
      end do
  end subroutine decode_tuple
  !======================================================================
  subroutine canonical_form(l, n, mu, new_n, new_mu)
      ! Produce the canonical form of (n, mu) given a fixed l.
      ! 1) Identify subsets of indices that share the same l-value
      ! 2) Sort the corresponding (n_i, mu_i) pairs lexicographically within each subset
      ! 3) Rebuild (n, mu) in that canonical order
      implicit none
      integer, intent(in)               :: l(:), n(:), mu(:)
      integer, allocatable, intent(out) :: new_n(:), new_mu(:)
      integer                           :: nu, i, j, k, count, unique_count
      integer, allocatable             :: unique_lvals(:), subset_indices(:)
      logical                           :: found
      nu = size(l)
      allocate(new_n(nu), new_mu(nu))
      ! (A) Collect unique l-values
      allocate(unique_lvals(nu))
      unique_count = 0
      do i = 1, nu
          found = .false.
          do j = 1, unique_count
              if (l(i) == unique_lvals(j)) then
                  found = .true.
                  exit
              end if
          end do
          if (.not. found) then
              unique_count = unique_count + 1
              unique_lvals(unique_count) = l(i)
          end if
      end do
      if (unique_count < nu) then
          unique_lvals = unique_lvals(1:unique_count)
      end if
      ! (B) Sort the distinct l-values (smallest to largest)
      call sort_int_array(unique_lvals, unique_count)
      ! (C) For each distinct l-value, gather indices, sort them by (n_i, mu_i), and build
      k = 1
      do i = 1, unique_count
          ! Gather indices with current l-value
          count = 0
          do j = 1, nu
              if (l(j) == unique_lvals(i)) count = count + 1
          end do
          allocate(subset_indices(count))
          count = 0
          do j = 1, nu
              if (l(j) == unique_lvals(i)) then
                  count = count + 1
                  subset_indices(count) = j
              end if
          end do
          ! Sort these indices by (n, mu) lexicographically
          call sort_by_nm(n, mu, subset_indices, count)
          ! Place them into (new_n, new_mu)
          do j = 1, count
              new_n(k)  = n(subset_indices(j))
              new_mu(k) = mu(subset_indices(j))
              k = k + 1
          end do
          deallocate(subset_indices)
      end do
      deallocate(unique_lvals)
  end subroutine canonical_form
  !======================================================================
  subroutine sort_int_array(arr, size)
      ! Simple ascending insertion sort for 1D array of integers
      implicit none
      integer, intent(inout) :: arr(:)
      integer, intent(in)    :: size
      integer                :: i, j, key
      do i = 2, size
          key = arr(i)
          j   = i - 1
          do while (j >= 1 .and. arr(j) > key)
              arr(j+1) = arr(j)
              j = j - 1
          end do
          arr(j+1) = key
      end do
  end subroutine sort_int_array
  !======================================================================
  subroutine sort_by_nm(n, mu, indices, m)
      ! Sort subset 'indices' of length m by (n(indices), mu(indices)) lexicographically
      ! Currently uses a simple insertion sort for demonstration.
      implicit none
      integer, intent(in)    :: n(:), mu(:), m
      integer, intent(inout) :: indices(:)
      integer                :: i, j, key_idx
      integer                :: key_n, key_mu
      do i = 2, m
          key_idx = indices(i)
          key_n   = n(key_idx)
          key_mu  = mu(key_idx)
          j       = i - 1
          do while (j >= 1)
              if ( (n(indices(j))  > key_n) .or. &
                   ((n(indices(j)) == key_n) .and. (mu(indices(j)) > key_mu)) ) then
                  indices(j+1) = indices(j)
                  j = j - 1
              else
                  exit
              end if
          end do
          indices(j+1) = key_idx
      end do
  end subroutine sort_by_nm
  !======================================================================
  subroutine sort_tuples(storage, count, length)
      ! Public routine to sort columns [1..count] in 'storage' using Quicksort
      ! for a lexicographic ordering of each column (size 'length').
      implicit none
      integer, intent(inout) :: storage(:,:)
      integer, intent(in)    :: count, length
      if (count <= 1) return
      call quicksort_columns(storage, 1, count, length)
  end subroutine sort_tuples
  !----------------------------------------------------------------------
  recursive subroutine quicksort_columns(storage, left, right, length)
      ! Recursive quicksort on columns from 'left' to 'right'.
      implicit none
      integer, intent(inout) :: storage(:,:)
      integer, intent(in)    :: left, right, length
      integer                :: i, j
      if (left < right) then
          call partition_columns(storage, left, right, length, i, j)
          call quicksort_columns(storage, left, j, length)
          call quicksort_columns(storage, i, right, length)
      end if
  end subroutine quicksort_columns
  !----------------------------------------------------------------------
  subroutine partition_columns(storage, left, right, length, i, j)
      ! Standard partition step for quicksort. We pick a pivot column
      ! from the middle, then move i and j inward until we need a swap.
      implicit none
      integer, intent(inout) :: storage(:,:)
      integer, intent(in)    :: left, right, length
      integer, intent(out)   :: i, j
      integer, allocatable   :: pivot(:)
      integer                :: mid
      integer, allocatable   :: tmp_col(:)
      allocate(pivot(length), tmp_col(length))
      i   = left
      j   = right
      mid = (left + right) / 2
      ! Copy pivot
      pivot = storage(:, mid)
      do
          ! Move i to the right while the column is less than pivot
          do while (i <= right)
              if (.not. tuple_less(storage(:, i), pivot, length)) exit
              i = i + 1
          end do
          ! Move j to the left while the column is greater than pivot
          do while (j >= left)
              if (.not. tuple_less(pivot, storage(:, j), length)) exit
              j = j - 1
          end do
          if (i <= j) then
              tmp_col         = storage(:, i)
              storage(:, i)   = storage(:, j)
              storage(:, j)   = tmp_col
              i = i + 1
              j = j - 1
          else
              exit
          end if
      end do
      deallocate(pivot, tmp_col)
  end subroutine partition_columns
  !----------------------------------------------------------------------
  logical function tuple_less(t1, t2, length)
      ! Returns TRUE if t1 < t2 in lexicographic order, else FALSE.
      implicit none
      integer, intent(in) :: length
      integer, intent(in) :: t1(:), t2(:)
      integer             :: k
      do k = 1, length
          if (t1(k) < t2(k)) then
              tuple_less = .true.
              return
          else if (t1(k) > t2(k)) then
              tuple_less = .false.
              return
          end if
      end do
      ! If they are completely equal, t1 is not "less" than t2
      tuple_less = .false.
  end function tuple_less
  !======================================================================
  subroutine remove_duplicates(storage, count, length)
      ! Removes consecutive duplicates after sorting. 
      ! 'count' is adjusted to the new number of unique columns.
      implicit none
      integer, intent(inout) :: storage(:,:), count
      integer, intent(in)    :: length
      integer                :: i, j
      if (count <= 1) return
      j = 1
      do i = 2, count
          if (.not. are_equal(storage(:, j), storage(:, i), length)) then
              j = j + 1
              storage(:, j) = storage(:, i)
          end if
      end do
      count = j
  end subroutine remove_duplicates
  !----------------------------------------------------------------------
  logical function are_equal(t1, t2, length)
      ! Checks if t1 == t2 (element-wise).
      implicit none
      integer, intent(in) :: length
      integer, intent(in) :: t1(:), t2(:)
      integer             :: k
      do k = 1, length
          if (t1(k) /= t2(k)) then
              are_equal = .false.
              return
          end if
      end do
      are_equal = .true.
  end function are_equal
  !======================================================================
  function array_to_string(arr) result(str)
      ! Utility function to convert integer array to a string like "(1,2,3,...)".
      implicit none
      integer, intent(in) :: arr(:)
      character(len=:), allocatable :: str
      character(len=32) :: tmp_str
      integer :: i, n
      n = size(arr)
      str = "("
      do i = 1, n
          write(tmp_str, '(I0)') arr(i)
          str = trim(str)//tmp_str
          if (i < n) str = trim(str)//","
      end do
      str = trim(str)//")"
  end function array_to_string
  !======================================================================
  function count_distinct_n(storage, nu) result(ndistinct)
    ! Counts the number of distinct n = (n_1,...,n_nu) in the sorted 'storage'.
    ! We assume that 'storage' is sorted lexicographically by (n, mu).
    ! storage(1..nu,  i ) = n^i
    ! storage(nu+1..2nu,i) = mu^i
    !
    ! num_cols = total columns in 'storage' after duplicates are removed.
    ! nu       = dimension of n and mu.

    implicit none
    integer, intent(in) :: nu
    integer, intent(in) :: storage(:,:)
    integer :: ndistinct
    integer :: i
    integer, allocatable :: last_n(:), current_n(:)
    
    if (size(storage, 2) <= 1) then
       ndistinct = 1
       return
    end if

    allocate(last_n(nu), current_n(nu))

    ! The first column is definitely a new n
    ndistinct = 1
    last_n = storage(1:nu, 1)

    do i = 2, size(storage, 2)
        current_n = storage(1:nu, i)
        ! Compare current_n with last_n
        if (.not. are_equal(current_n, last_n, nu)) then
            ndistinct = ndistinct + 1
            last_n = current_n
        end if
    end do

    deallocate(last_n, current_n)
end function count_distinct_n



end module  module_tuples_for_given_lbold
!>------end_module_tuples_for_given_lbold------------------------------------------------


module module_perm_basis_remove 
  implicit none
  contains 

  subroutine sort_basis (nu, Mb, allBases, index_of_sorting)
   ! -----------------------------------------------------------------------
   ! Data structures
   ! -----------------------------------------------------------------------
   ! allBases holds Mb bases, each basis is shape (3,nu):
   !   * row 1 => mu-values
   !   * row 2 => n-values
   !   * row 3 => l-values
   integer, intent(in) :: nu, Mb 
   integer, dimension(3, nu, Mb), intent(in)  :: allBases
   integer, dimension(:), allocatable :: index_of_sorting 
   !integer, parameter :: nu = 5         ! Max # of orbitals in each basis, for demo
   !integer, parameter :: Mb = 100        ! # of total basis sets, for demo
   integer, parameter :: maxDictionary = 100000  ! capacity of dictionary
   integer, parameter :: maxStringLen  = 200     ! max length of the "key" string

   ! We'll store the "canonical forms" (string keys) in canonicalForms
   character(len=maxStringLen), dimension(maxDictionary) :: canonicalForms
   integer :: dictionarySize  ! how many unique entries so far
   
   ! Temporary variables
   integer :: iBasis
   integer, dimension(3, nu) :: Btemp, Bsorted
   character(len=maxStringLen) :: keyString
   integer, dimension(:), allocatable :: tmp_sorting

   ! -----------------------------------------------------------------------
   ! Initialize dictionary
   ! -----------------------------------------------------------------------
   dictionarySize = 0
   if (allocated(index_of_sorting)) deallocate(index_of_sorting) ; allocate(index_of_sorting(Mb))

   ! -----------------------------------------------------------------------
   ! Main loop: Process each basis
   ! -----------------------------------------------------------------------
   do iBasis = 1, Mb

      ! 1) Extract the iBasis-th basis (3 x nu)
      Btemp(1:3, 1:nu) = allBases(1:3, 1:nu, iBasis)

      ! 2) Sort the orbital tuples by (mu, n, l)
      call sort_tuples(Btemp, Bsorted, nu)

      ! 3) Flatten & stringify into a single "key"
      call flatten_and_stringify(Bsorted, nu, keyString)
      !debug write(*,*) keyString
      ! 4) Check if this key already exists in our dictionary
      if (.not. is_already_in_dictionary(keyString, canonicalForms, dictionarySize)) then
         ! 5) If not, add it
         call add_to_dictionary(keyString, canonicalForms, dictionarySize)
         index_of_sorting(dictionarySize) = iBasis
      end if
   end do

   if (allocated(tmp_sorting)) deallocate(tmp_sorting) ; allocate(tmp_sorting(dictionarySize))
   tmp_sorting = index_of_sorting(1:dictionarySize)
   if (allocated(index_of_sorting)) deallocate(index_of_sorting) ; allocate(index_of_sorting(dictionarySize))
   index_of_sorting = tmp_sorting(1:dictionarySize)
   deallocate(tmp_sorting)

   ! -----------------------------------------------------------------------
   ! Print how many distinct bases we found
   ! -----------------------------------------------------------------------
   ! print *, "Number of distinct bases (up to permutation) =", dictionarySize

  end subroutine sort_basis

  ! =======================================================================
  !  Subroutine: sort_tuples
  !
  !  Sorts the input array B(3,nu) by ascending (mu, n, l),
  !  placing the result into Bsorted(3,nu).
  ! =======================================================================
  subroutine sort_tuples(B, Bsorted, nu)
    implicit none
    integer, intent(in)    :: nu
    integer, dimension(3, nu), intent(in)  :: B
    integer, dimension(3, nu), intent(out) :: Bsorted

    integer :: j, i
    integer, dimension(3) :: tmp

    ! First, copy the input B -> Bsorted
    Bsorted = B

    ! We'll do a simple insertion sort (fine for small nu)
    do j = 2, nu
       tmp = Bsorted(:, j)
       i = j - 1

       do while (i >= 1)
          if ( compare(Bsorted(:, i), tmp) > 0 ) then
             Bsorted(:, i+1) = Bsorted(:, i)
             i = i - 1
          else
             exit
          end if
       end do
       Bsorted(:, i+1) = tmp
    end do

  end subroutine sort_tuples

  ! =======================================================================
  !  Function: compare
  !
  !  Compares two triplets (mu1,n1,l1) and (mu2,n2,l2) lexicographically:
  !     returns -1 if (r1) < (r2)
  !     returns  0 if (r1) = (r2)
  !     returns +1 if (r1) > (r2)
  ! =======================================================================
  integer function compare(r1, r2)
    implicit none
    integer, dimension(3), intent(in) :: r1, r2

    ! Compare mu
    if (r1(1) < r2(1)) then
       compare = -1
       return
    else if (r1(1) > r2(1)) then
       compare = +1
       return
    end if

    ! Mu are equal, compare n
    if (r1(2) < r2(2)) then
       compare = -1
       return
    else if (r1(2) > r2(2)) then
       compare = +1
       return
    end if

    ! Mu & n are equal, compare l
    if (r1(3) < r2(3)) then
       compare = -1
    else if (r1(3) > r2(3)) then
       compare = +1
    else
       compare = 0
    end if

  end function compare

  ! =======================================================================
  !  Subroutine: flatten_and_stringify
  !
  !  Converts Bsorted(3, nu) into a single string key, e.g.:
  !    "(4,1,0):(5,1,1):(5,1,1) ..."
  ! =======================================================================
  subroutine flatten_and_stringify(Bsorted, nu, key)
    implicit none
    integer, intent(in)             :: nu
    integer, dimension(3, nu), intent(in) :: Bsorted
    character(len=*), intent(out)   :: key

    integer i
    character(len=30) :: tupleString  ! temporary for one tuple
    key = ''  ! start empty

    do i = 1, nu
       write(tupleString, '( "(",I0,",",I0,",",I0,")" )') Bsorted(1,i), Bsorted(2,i), Bsorted(3,i)
       if (i == 1) then
          key = trim(tupleString)
       else
          key = trim(key) // ':' // trim(tupleString)
       end if
    end do

  end subroutine flatten_and_stringify

  ! =======================================================================
  !  Function: is_already_in_dictionary
  !
  !  Returns .true. if key is found in canonicalForms(1..dictionarySize),
  !  else .false.
  ! =======================================================================
  logical function is_already_in_dictionary(key, canonicalForms, dictionarySize)
    implicit none
    character(len=*), intent(in) :: key
    character(len=*), dimension(:), intent(in) :: canonicalForms
    integer, intent(in)            :: dictionarySize

    integer i
    is_already_in_dictionary = .false.

    do i = 1, dictionarySize
       if (trim(canonicalForms(i)) == trim(key)) then
          is_already_in_dictionary = .true.
          return
       end if
    end do

  end function is_already_in_dictionary

  ! =======================================================================
  !  Subroutine: add_to_dictionary
  !
  !  Appends 'key' to the end of canonicalForms and increments dictionarySize.
  ! =======================================================================
  subroutine add_to_dictionary(key, canonicalForms, dictionarySize)
    implicit none
    character(len=*), intent(in)  :: key
    character(len=*), dimension(:), intent(inout) :: canonicalForms
    integer, intent(inout)        :: dictionarySize

    dictionarySize = dictionarySize + 1
    canonicalForms(dictionarySize) = key
  end subroutine add_to_dictionary
 

end module module_perm_basis_remove 
