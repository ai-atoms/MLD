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


subroutine init_kernel_zetabody 
  use module_kind_variables, only: kind_double
  use module_chemical_species, only: fix_no_of_elements
  use module_neigh_local, only :  r_cut, r_cut_in,  build_cos_and_fcut_ja 
  use module_body_desc, only: l_body_order, dim_desc_body
  use module_kernel_zetabody, only: zetabody_order, zetabody_dim, dim_grid_zetabody,  & 
                                    time_for_desc_kernel_nb, tnn_nb, tuu_nb, tkk_nb, & 
                                    ispec2, ispec3, &
                                    dim_length_zetabody, length_zetabody2, length_zetabody3, &
                                    delta_zetabody2, delta_zetabody3, dim_delta_zetabody
  use module_ftnbody, only: dim_qbody
  use grid_and_permutations_nbody, only: grid_nbody

  use mld_logger
  implicit none 
  integer :: ii
  logical, parameter :: desc_forces_local = .false.  
  integer :: itmp, i1, i2, i3  
  character(LEN=80) :: chlog
  logical :: auto_zpoints
  integer :: nitems2, ntest 

  _NAMECURRENT_("init_kernel_zetabody")

  ! this is complicated object until 3-body: All the orthers body comes with zeta-body power   
  allocate(grid_nbody(3))
  
  auto_zpoints = .true.
   
  !$! write(*,*) 'debug --> fix_no_of_elements', fix_no_of_elements
  !TESTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
  !fix_no_of_elements = 3
  !zeta_order = 3
  !l_body_order(2) = .true.
  !debug! write(*,*) 'debug --> l_body_order(2)', l_body_order(2)
  !debug! write(*,*) 'debug --> l_body_order(3)', l_body_order(3)
  do ii = 4, size(l_body_order,1)
    if (l_body_order(ii)) then
     call log_warning("init_kernel_zetabody:---------------ATTENTION---------------------------")
     call log_warning("init_kernel_zetabody: l_body_order("//trim(vtoa(ii))//") is detected true")
     call log_warning("init_kernel_zetabody: l_body_order( "//trim(vtoa(ii))//" will be set to false") 
     l_body_order(ii) = .false.
    end if 
  end do
  !TESTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

  if (l_body_order(2)) then 
     !call log_critical("init_kernel_nbody:----------2-BODY NOT YET --------------------------")
     !stop 
     dim_qbody(2) = 1

     read(dim_length_zetabody(2), *) length_zetabody2
     if (allocated(delta_zetabody2)) deallocate(delta_zetabody2)
     allocate(delta_zetabody2(zetabody_order))
     ntest = nitems2(dim_delta_zetabody(2))
     if (ntest .lt. zetabody_order) then 
       call log_critical("init_kernel_zetabody: the number of items in dim_delta_zetabody(2) is not equal to zeta_order")
       call log_critical("--- change and relaunch  ---")
       stop "items in dim_delta_zetabody(2)"
     end if
     read(dim_delta_zetabody(2), *) delta_zetabody2

     call grid_nbody(2)%init_grid_nbody(2, r_cut_in, r_cut)
     call grid_nbody(2)%init_tuples_nbody(fix_no_of_elements)
     !$! write(*,*) size(grid_nbody(2)%unique_classes,1), size(grid_nbody(2)%unique_classes,2)
     write(chlog, '(i6)') size(grid_nbody(2)%unique_classes,2) 
     call log_info("ML:  zeta-body --> for 2-body there are "//trim(chlog)//" unique classes")  
     ! .... number_of_zpoints 
     !$! write(*,*) 'debug --> dim_desc_body(3)', dim_grid_nbody(3)
     call grid_nbody(2)%init_zpoints(auto_zpoints, dim_grid_zetabody(2)) 
     dim_desc_body(2) = size(grid_nbody(2)%zpoints,1)

     write(chlog, '(i20)') dim_desc_body(2) 
     call log_info("ML:  zeta-body --> the 2-body dimension by z-points "//trim(chlog))
     if (allocated(ispec2)) deallocate(ispec2) 
     allocate(ispec2(fix_no_of_elements, fix_no_of_elements))
     
     do ii = 1, size(grid_nbody(2)%mu_space_ini, 1)
        i1 = grid_nbody(2)%mu_space_ini(ii,1)
        i2 = grid_nbody(2)%mu_space_ini(ii,2)
        ispec2(i1,i2) = grid_nbody(2)%map_into_unique_classes(ii)
     end do 
  else 
    call log_warning("init_kernel_nbody:---------------ATTENTION---------------------------")
    call log_warning("init_kernel_nbody: 2-body order is excluded from zeta-body descriptor")
    call log_warning("init_kernel_nbody: to recover add --->         l_body_order(2)=.true.")
  end if 

  if (l_body_order(3)) then 
     dim_qbody(3) = 3
     read(dim_length_zetabody(3), *) length_zetabody3(1:3)
     if (allocated(delta_zetabody3)) deallocate(delta_zetabody3)
     allocate(delta_zetabody3(zetabody_order))
     ntest = nitems2(dim_delta_zetabody(3))
     if (ntest .lt.  zetabody_order) then 
       call log_critical("init_kernel_zetabody: the number of items in dim_delta_zetabody(3) is not equal to zeta_order")
       call log_critical("--- change and relaunch  ---")
       stop "items in dim_delta_zetabody(3)"
     end if
     read(dim_delta_zetabody(3), *) delta_zetabody3

     !auto_zpoints = .true.
     call grid_nbody(3)%init_grid_nbody(3, r_cut_in, r_cut)
     call grid_nbody(3)%init_tuples_nbody(fix_no_of_elements)
     !$! write(*,*) size(grid_nbody(3)%unique_classes,1), size(grid_nbody(3)%unique_classes,2)
     write(chlog, '(i6)') size(grid_nbody(3)%unique_classes,2) 
     call log_info("ML:  zeta-body --> for 3-body there are "//trim(chlog)//" unique classes")
     ! .... number_of_zpoints 
     !$! write(*,*) 'debug --> dim_desc_body(3)', dim_grid_nbody(3)
     call grid_nbody(3)%init_zpoints(auto_zpoints, dim_grid_zetabody(3)) 
     dim_desc_body(3) = size(grid_nbody(3)%zpoints,1)

     write(chlog, '(i20)') dim_desc_body(3) 
     call log_info("ML:  zeta-body --> the 3-body dimension by z-points "//trim(chlog))
     if (allocated(ispec3)) deallocate(ispec3) 
     allocate(ispec3(fix_no_of_elements, fix_no_of_elements, fix_no_of_elements))
     do ii = 1, size(grid_nbody(3)%mu_space_ini, 1)
        i1 = grid_nbody(3)%mu_space_ini(ii,1)
        i2 = grid_nbody(3)%mu_space_ini(ii,2)
        i3 = grid_nbody(3)%mu_space_ini(ii,3)
        ispec3(i1,i2,i3) = grid_nbody(3)%map_into_unique_classes(ii)
     end do 


  else 
    call log_warning("init_kernel_nbody:---------------ATTENTION---------------------------")
    call log_warning("init_kernel_nbody: 3-body order is excluded from zeta-body descriptor")
    call log_warning("init_kernel_nbody: to recover add --->         l_body_order(3)=.true.")
  end if  

  call log_info("ML:  zeta-body --> the order "//trim(vtoa(zetabody_order))//" is included in the descriptor")

  itmp = 0
  do ii = 1, zetabody_order
    if (l_body_order(2)) then
      itmp = itmp + dim_desc_body(2) 
    end if 
    if (l_body_order(3)) then
      itmp = itmp + dim_desc_body(3) 
    end if
  end do

  write(chlog, '(i20)') itmp 
  zetabody_dim = itmp 
  call log_info("ML:  zeta-body --> the total draft dimension is "//trim(chlog))
  
!$! if (l_body_order(3)) then
!$!   dim_desc_body(3) = dim_rff(3)  
!$!   !TODOspec for multispecies: 
!$!   count_3b = 0.d0
!$!   if (allocated(covar_matrix_3b)) deallocate (covar_matrix_3b); allocate (covar_matrix_3b(dim_qbody(3), dim_qbody(3)))
!$!   if (allocated(length_rff_3b)) deallocate (length_rff_3b); allocate (length_rff_3b(dim_qbody(3)))
!$!   if (allocated(mean_3b)) deallocate (mean_3b); allocate (mean_3b(dim_qbody(3)))
!$!   covar_matrix_3b(:,:) = 0.d0
!$!   mean_3b(:)=0.d0 
!$!   if (allocated(ftnbody_omega_3b)) deallocate (ftnbody_omega_3b); allocate (ftnbody_omega_3b(dim_qbody(3), dim_desc_body(3)))
!$!   if (allocated(ftnbody_phase_random_3b)) deallocate (ftnbody_phase_random_3b); allocate (ftnbody_phase_random_3b(dim_desc_body(3)))
!$!   if (.not.(find_best_length)) then
!$!     call  get_rff_val_sigma(krff_type, mu, dim_qbody(3), dim_desc_body(3), ftnbody_omega_3b, ftnbody_phase_random_3b)
!$!     ftnbody_omega_3b = ftnbody_omega_3b / (length_rff(3)*sqrt(2.d0))
!$!   end if
!$! end if   


  time_for_desc_kernel_nb = 0.d0 
  tnn_nb = 0.d0 
  tuu_nb = 0.d0 
  tkk_nb = 0.d0 
  _MLD_END_
end subroutine init_kernel_zetabody 


subroutine init_zetabody 
  use module_kind_variables, only: kind_double
  use module_body_desc, only: l_body_order, dim_desc_body 
  use module_ftnbody, only: ftnbody_dim, dim_qbody, dim_rff, length_rff, ftnbody_phase_random_2b, &
            ftnbody_phase_random_3b, &
            ftnbody_omega_2b, ftnbody_omega_3b,  &
            covar_matrix_2b, covar_matrix_3b,   &
            length_rff_2b, length_rff_3b, count_2b, count_3b, &
            find_best_length, mean_2b, mean_3b
  ! use module_kernel, only: np_kernel_full
  use module_kernel, only:  krff_type
  use mld_logger
  use module_sample_rand_ker, only: get_rff_val_sigma
  implicit none
  real(kind_double)    :: mu
  integer :: ii, icnt 
    
  _NAMECURRENT_("init_zetabody")
  _MLD_BEGIN_

  dim_desc_body(:)=0

  ! find_best_length = .true.
  if (.not.(find_best_length)) then
    call log_info("ML: FT-nBody warning: without finding best length!!!")
  end if
  if (l_body_order(2)) then
    ! set the F 
    ! dim_rff(2) =  np_kernel_full
    ! size of the q vectors ... 
    dim_qbody(2) = 1
    ! true only for one type of atoms  ... 
    dim_desc_body(2) = dim_rff(2)  
    !TODOspec for multispecies: 

    count_2b = 0.d0
    if (allocated(covar_matrix_2b)) deallocate (covar_matrix_2b); allocate (covar_matrix_2b(dim_qbody(2), dim_qbody(2)))
    if (allocated(length_rff_2b)) deallocate (length_rff_2b); allocate (length_rff_2b(dim_qbody(2)))
    if (allocated(mean_2b)) deallocate (mean_2b); allocate (mean_2b(dim_qbody(2)))
    covar_matrix_2b(:,:) = 0.d0
    mean_2b(:)=0.d0 

    if (allocated(ftnbody_omega_2b)) deallocate (ftnbody_omega_2b); allocate (ftnbody_omega_2b(dim_qbody(2), dim_desc_body(2)))
    if (allocated(ftnbody_phase_random_2b)) deallocate (ftnbody_phase_random_2b); allocate (ftnbody_phase_random_2b(dim_desc_body(2)))
    if (.not.(find_best_length)) then
      call  get_rff_val_sigma(krff_type, mu, dim_qbody(2), dim_desc_body(2), ftnbody_omega_2b, ftnbody_phase_random_2b)
      ftnbody_omega_2b = ftnbody_omega_2b / (length_rff(2)*sqrt(2.d0))
    end if 
  end if   

  if (l_body_order(3)) then
    ! set the F 
    ! dim_rff(3) =  np_kernel_full
    ! size of the q vectors ... 
    dim_qbody(3) = 3
    ! true only for one typer of atoms  ... 
    dim_desc_body(3) = dim_rff(3)  
    !TODOspec for multispecies: 

    count_3b = 0.d0
    if (allocated(covar_matrix_3b)) deallocate (covar_matrix_3b); allocate (covar_matrix_3b(dim_qbody(3), dim_qbody(3)))
    if (allocated(length_rff_3b)) deallocate (length_rff_3b); allocate (length_rff_3b(dim_qbody(3)))
    if (allocated(mean_3b)) deallocate (mean_3b); allocate (mean_3b(dim_qbody(3)))
    covar_matrix_3b(:,:) = 0.d0
    mean_3b(:)=0.d0 

    if (allocated(ftnbody_omega_3b)) deallocate (ftnbody_omega_3b); allocate (ftnbody_omega_3b(dim_qbody(3), dim_desc_body(3)))
    if (allocated(ftnbody_phase_random_3b)) deallocate (ftnbody_phase_random_3b); allocate (ftnbody_phase_random_3b(dim_desc_body(3)))
    if (.not.(find_best_length)) then
      call  get_rff_val_sigma(krff_type, mu, dim_qbody(3), dim_desc_body(3), ftnbody_omega_3b, ftnbody_phase_random_3b)
      ftnbody_omega_3b = ftnbody_omega_3b / (length_rff(3)*sqrt(2.d0))
    end if
  end if   

  icnt = 0 
  do ii = 1, size(l_body_order,1)
    if (l_body_order(ii)) icnt = icnt + dim_desc_body(ii)
  end do   

  ftnbody_dim = icnt 

  do ii = 1, size(dim_desc_body,1)
    if (l_body_order(ii)) then 
    call log_info("ML: FT-nBody with "//vtoa(ii)//"-body has the size "//vtoa(dim_desc_body(ii)))
    end if 
  end do   
  call log_info("ML: FT-nBody with total size "// vtoa(ftnbody_dim))

  _MLD_END_ 
end subroutine init_zetabody

subroutine init_length_zetabody
  use module_kind_variables, only: kind_double
  use module_body_desc, only: l_body_order
  use module_body_desc, only: dim_desc_body
  use module_ftnbody, only: dim_qbody, init_mode_ftnbody, covar_matrix_2b, covar_matrix_3b,  count_2b, count_3b, &
                              length_rff_2b, length_rff_3b,  length_order,   &
                              ftnbody_omega_2b, ftnbody_omega_3b,  &
                              ftnbody_phase_random_2b, ftnbody_phase_random_3b, &
                              covar_matrix_2b, number_2b_ftnb, mean_2b, & 
                              covar_matrix_3b, number_3b_ftnb, mean_3b, &
                              ltmp_mean_ftnbody
  use module_kernel, only:  krff_type
  use module_sample_rand_ker, only: get_rff_diag_sigma
  use mld_logger
  use mld_mpi
  use mld_string, only: vtoa
  use my_mpi_subroutines, only: subworlds_allreduce_from_evrywhere_double, &
                                  subworlds_allreduce_from_evrywhere_vect_double, & 
                                  subworlds_allreduce_from_evrywhere_matrix_double
  implicit none
  integer :: iiii, ii 
  real(kind_double) :: mu
  real(kind_double), allocatable, dimension(:) :: total_mean 
  real(kind_double), allocatable, dimension(:,:) :: total_cov 
  character(len=80) :: chstring, string  

  _NAMECURRENT_("init_length_ftnbody_")
  _MLD_BEGIN_

  ! First step computing the means ...........................
  init_mode_ftnbody = .true.
  ltmp_mean_ftnbody = .true. 
  call log_info("ML: ... ftnbody begins , in "//NAMECURRENT//" estimation of length, step 1. the mean " )
  call main_compute_descriptors()
  if (l_body_order(2)) then
    number_2b_ftnb = 0.d0 
 
    call subworlds_allreduce_from_evrywhere_double(count_2b, number_2b_ftnb)
      
    if (allocated(total_mean)) deallocate(total_mean) 
    allocate(total_mean(size(mean_2b,1)))
      
    call subworlds_allreduce_from_evrywhere_vect_double(mean_2b, total_mean)

    mean_2b = total_mean / number_2b_ftnb 

    write(chstring, '(i10)') int(number_2b_ftnb)
    call log_info("ML: ... 2b number ftnbody  : "//trim(chstring))
    string=''
    do ii = 1, size(mean_2b,1)
      write(chstring, '(es20.9)') mean_2b(ii)
      string = trim(string)//trim(chstring)//" "
    end do   
    call log_info("ML: ... 2b mean   ftnbody  : "//trim(string))
    count_2b = 0 
    covar_matrix_2b(:,:)=0.d0 
  end if

  if (l_body_order(3)) then
    number_3b_ftnb = 0.d0 
    
    call subworlds_allreduce_from_evrywhere_double(count_3b, number_3b_ftnb)
    if (allocated(total_mean)) deallocate(total_mean) 
    allocate(total_mean(size(mean_3b,1)))
    call subworlds_allreduce_from_evrywhere_vect_double(mean_3b, total_mean)
    mean_3b = total_mean / number_3b_ftnb 
    write(chstring, '(i10)') int(number_3b_ftnb)
    call log_info("ML: ... 3b number ftnbody  : "//trim(chstring))
    string=''
    do ii = 1, size(mean_3b,1)
      write(chstring, '(es20.9)') mean_3b(ii)
      string = trim(string)//trim(chstring)//" "
    end do  
    call log_info("ML: ... 3b mean   ftnbody  : "//trim(string))

    count_3b = 0 
    covar_matrix_3b(:,:)=0.d0 
  end if 


  !  step computing the covariance matrix  ...........................
  call log_info("ML: ... ftnbody mean estimated , in "//NAMECURRENT//" estimation of length, step 2. the covariance " )
  init_mode_ftnbody = .true.
  ltmp_mean_ftnbody = .false. 
  call main_compute_descriptors()    
  init_mode_ftnbody = .false. 
  ltmp_mean_ftnbody = .false. 

  call log_info("ML: ... ftnbody completed the covariance " )

  if (l_body_order(2)) then
    if (allocated(total_cov)) deallocate(total_cov) 
    allocate(total_cov(size(covar_matrix_2b,1), size(covar_matrix_2b,2)))
    call subworlds_allreduce_from_evrywhere_matrix_double(covar_matrix_2b, total_cov)
    covar_matrix_2b(:,:) = total_cov(:,:) / number_2b_ftnb
    length_rff_2b(1) = SQRT(covar_matrix_2b(1,1))
    length_rff_2b(:) =  length_rff_2b(:)/2
    string=''
    do ii = 1, size(length_rff_2b,1)
      write(chstring, '(es20.9)') length_rff_2b(ii)
      string = trim(string)//trim(chstring)//" "
    end do  
    call log_info("ML: ... 2b length   ftnbody  : "//trim(string))
    length_order = 2
    call  get_rff_diag_sigma(krff_type, mu, length_rff_2b, dim_qbody(2), dim_desc_body(2), ftnbody_omega_2b, ftnbody_phase_random_2b)
  end if

  if (l_body_order(3)) then
    if (allocated(total_cov)) deallocate(total_cov) 
    allocate(total_cov(size(covar_matrix_3b,1), size(covar_matrix_3b,2)))

    call subworlds_allreduce_from_evrywhere_matrix_double(covar_matrix_3b, total_cov)
    covar_matrix_3b(:,:) = total_cov(:,:) / number_3b_ftnb
    do iiii = 1, dim_qbody(3)
      length_rff_3b(iiii) = SQRT(covar_matrix_3b(iiii,iiii))
    end do
    length_rff_3b(:) = length_rff_3b(:)/2
    string=''
    do ii = 1, size(length_rff_3b,1)
      write(chstring, '(es20.9)') length_rff_3b(ii)
      string = trim(string)//trim(chstring)//" "
    end do  
    call log_info("ML: ... 3b length   ftnbody  : "//trim(string))
    length_order = 3
    call  get_rff_diag_sigma(krff_type, mu, length_rff_3b, dim_qbody(3), dim_desc_body(3), ftnbody_omega_3b, ftnbody_phase_random_3b)
  end if 

  _MLD_END_

end subroutine init_length_zetabody



module module_compute_zetabody_order
    use module_kind_variables, ONLY: kind_double
    implicit none 
    contains 
    subroutine zetabody_order_2(type_db_ja, max_neigh_local, i_type_db, i_central, r_central, tmp_dxp, ur_central, d_ur_central, desc_forces_local)
                                
      use module_body_desc, only : dim_desc_body
      use module_kernel_zetabody, only: r_cut_z2b, r_cut_width_z2b, & 
                                        ispec2, delta_zetabody2,  & 
                                        length_zetabody2, tmp_zeta_real, tmp_zeta_dreal, tmp_zeta_norm 
      use module_neigh_local, only: r_cut_in, r_cut_width_in, r_cut_pair_in, fcut_rij_inout
      use grid_and_permutations_nbody, only: grid_nbody
      
      implicit none 
      integer, intent(in) ::  type_db_ja,  max_neigh_local
      integer, dimension(0:max_neigh_local), intent(in)  :: i_central
      integer, dimension(0:max_neigh_local), intent(in)  :: i_type_db 
      real(kind_double), dimension(:), intent(in) :: r_central, ur_central
      real(kind_double), dimension(:,:), intent(in) :: tmp_dxp
      real(kind_double), dimension(:,:), intent(in) :: d_ur_central
      logical, intent(in) :: desc_forces_local 
      integer :: ia, ix, ii 
      real(kind_double), dimension(:), allocatable :: tmpf 
      real(kind_double) :: norm 
      real(kind_double), dimension(:), allocatable :: fcut_all, dfcut_all 
      real(kind_double) :: fcut, dfcut
      real(kind_double) :: local_r_cut_in
      integer :: type_fcut_in, type_fcut_out, type_db_ia
      !TODOftnbody...
      integer :: izozo, izz, itt  
      integer, dimension(size(i_type_db)) :: iv1_zozo  
      integer, dimension(size(i_central)) :: iv2_zozo  
      real(kind_double) :: zz_k2,  d_zz_k2(3)
      real(kind_double), dimension(:,:), allocatable  :: zz_tmp 
      real(kind_double), dimension(dim_desc_body(2)) :: tmp_ker2
      real(kind_double), dimension(dim_desc_body(2), 3) :: d_tmp_ker2
      real(kind_double) :: x1
      real(kind_double), dimension(3) :: d_x1
 
      izozo =  type_db_ja
      iv1_zozo = i_type_db
      iv2_zozo = i_central


      allocate(tmpf(dim_desc_body(2)))
      allocate(fcut_all(max_neigh_local))
      allocate(dfcut_all(max_neigh_local))

      !TO_DO_ZETABODY ... 
      type_fcut_in = 3 
      type_fcut_out = 2 
      norm = sqrt(1.d0)/sqrt(dble(dim_desc_body(2)))*delta_zetabody2(1)
      !$! write(*,*) 'debug --> norm', norm, delta_zetabody2(1), dble(dim_desc_body(2))
      tmp_zeta_norm = norm  
      ! compute_qja_a_partir_de_rcentral 
      tmp_zeta_real(:) = 0.d0 
      do ii = 1, max_neigh_local
        ! Pair-specific inner cutoff
        local_r_cut_in = r_cut_pair_in(type_db_ja, i_type_db(ii))
        call fcut_rij_inout(type_fcut_in, type_fcut_out, r_central(ii), r_cut_z2b, r_cut_width_z2b, local_r_cut_in, r_cut_width_in, desc_forces_local, fcut, dfcut)
        fcut_all(ii) = fcut
        dfcut_all(ii) = dfcut
      end do
      !write(*,*) r_cut_z2b, r_cut_width_z2b

      do ia = 1, max_neigh_local
          type_db_ia =  i_type_db(ia) 
          itt = ispec2(type_db_ja, type_db_ia)
          !tmpf(:) = r_central(ia) - zpoints_by_tuple(iss)%

          !$! if (init_mode_ftnbody) then
          !$!   count_2b = count_2b + 1
          !$!   if (ltmp_mean_ftnbody) then 
          !$!     mean_2b(:) = mean_2b(:) + r_central(ia)
          !$!   else   
          !$!     covar_matrix_2b(:,:) = covar_matrix_2b(:,:) + (r_central(ia) - mean_2b(1) ) **2 
          !$!   end if   
          !$!   cycle
          !$! end if
          fcut = fcut_all(ia)
          dfcut = dfcut_all(ia)
          x1 = ur_central(ia)
          if (desc_forces_local) then 
            d_x1(:) = d_ur_central(:,ia)
          end if 
          
          if (allocated(zz_tmp)) deallocate(zz_tmp)
          allocate(zz_tmp(size(grid_nbody(2)%zpoints_by_tuple(itt)%zz,1), size(grid_nbody(2)%zpoints_by_tuple(itt)%zz,2)))
          zz_tmp(:,:) = grid_nbody(2)%zpoints_by_tuple(itt)%zz(:,:)

          tmp_ker2(:) = 0.d0 
          d_tmp_ker2(:, :) = 0.d0 
          do izz = grid_nbody(2)%zpoints_by_tuple(itt)%idx_min, grid_nbody(2)%zpoints_by_tuple(itt)%idx_max
              !$! zz_k2 = dexp(-0.5d0 * (r_central(ia) - zz_tmp(1, izz))**2  / length_zetabody2**2) 
              zz_k2 = dexp(-0.5d0 * (x1 - zz_tmp(1, izz))**2  / length_zetabody2**2) * tmp_zeta_norm
              tmp_ker2(izz) = zz_k2 * fcut 

              if (desc_forces_local) then 
                do ix=1,3 
                  !$! d_zz_k2(ix)  = - (r_central(ia) - zz_tmp(1, izz)) / length_zetabody2**2 * zz_k2 * fcut * tmp_dxp(ix, ia)/r_central(ia) + & 
                  !$!                  zz_k2 * tmp_dxp(ix, ia)/r_central(ia) * dfcut
                  d_zz_k2(ix)  = - (x1 - zz_tmp(1, izz)) / length_zetabody2**2 * zz_k2 * fcut * d_x1(ix) + & 
                                   zz_k2 * tmp_dxp(ix, ia)/r_central(ia) * dfcut
                end do
                d_tmp_ker2(izz, : ) = d_zz_k2(:)
              end if               
          end do 

          tmp_zeta_real(:) = tmp_zeta_real(:) + tmp_ker2(:)
          if (desc_forces_local) then
            do ix=1,3 
              tmp_zeta_dreal(:, ia, ix) =  d_tmp_ker2(:, ix)
            end do
          end if

      
          !old! tmpf(:) = r_central(ia)*ftnbody_omega_2b(1,:) + ftnbody_phase_random_2b(:)
          !old! fcut = fcut_all(ia)
          !old! dfcut = dfcut_all(ia)
          !old! !tmp_real(:) = tmp_real(:)  + cos (r_central(ia)*ftnbody_omega_2b(1,:) + ftnbody_phase_random_2b(:))
          !old! tmp_real(:) = tmp_real(:)  + cos (tmpf(:)) * norm * fcut 
          !old! if (desc_forces_local) then
          !old!   do ix =1,3 
          !old!   !tmp_dreal(:, ia, ix) =  - sin (r_central(ia)*ftnbody_omega_2b(1,:) + ftnbody_phase_random_2b(:)) * ftnbody_omega_2b(1,:) * tmp_dxp(ix, ia)/r_central(ia)*norm
          !old!   tmp_dreal(:, ia, ix) =  - sin (tmpf(:)) * ftnbody_omega_2b(1,:) * tmp_dxp(ix, ia)/r_central(ia)*norm * fcut + cos (tmpf(:)) * tmp_dxp(ix, ia)/r_central(ia) * norm * dfcut
          !old!   end do 
          !old! end if  
      end do      

      deallocate(tmpf) 
    end subroutine zetabody_order_2 


    subroutine zetabody_order_3(type_db_ja, max_neigh_local, i_type_db, i_central, r_central, tmp_dxp, ur_central, d_ur_central, desc_forces_local)
    ! compute all the body triangle with the following structure. Please pay attention that the distances
     ! are already r_cut -ed.
     !
     !     i2
     !    /
     !  x1
     ! / x3
     ! j--x2--i3
      use module_body_desc, only : dim_desc_body
      use module_neigh_local, only: r_cut_in, r_cut_width_in, r_cut_pair_in, fcut_rij_inout
      use grid_and_permutations_nbody, only: grid_nbody
      use module_kernel_zetabody, only: r_cut_z3b, r_cut_width_z3b,  & 
                                        ispec3, delta_zetabody3,  & 
                                        length_zetabody3, tmp_zeta_real, tmp_zeta_dreal, tmp_zeta_norm
      !use module_kernel, only:  length_kernel, sigma_kernel
      implicit none 
      integer, intent(in) ::  type_db_ja,  max_neigh_local
      integer, dimension(0:max_neigh_local), intent(in)  :: i_central
      integer, dimension(0:max_neigh_local), intent(in)  :: i_type_db 
      real(kind_double), dimension(:), intent(in) :: r_central, ur_central 
      real(kind_double), dimension(:,:), intent(in) :: tmp_dxp
      real(kind_double), dimension(:,:), intent(in) :: d_ur_central
      logical, intent(in) :: desc_forces_local 
      integer :: i2, i3, ix, ii, izz 
      real(kind_double), dimension(:), allocatable :: tmpf
      integer, parameter :: dimq = 3
      real(kind_double), dimension(dimq)   :: qsym
      real(kind_double), dimension(3, dimq)      :: d2_qsym, d3_qsym
      real(kind_double), dimension(dim_desc_body(3)) :: d2t, d3t
      real(kind_double), dimension(dim_desc_body(3)) :: tmp_ker3

      real(kind_double) :: cos_2j3, x1, x2, x3, norm, d_theta_d_cos !, one_to_r_central_i2=0.d0
      real(kind_double), dimension(3)   :: d2_x1, d2_x2, d2_x3, &
                                           d3_x1, d3_x2, d3_x3
      real(kind_double), dimension(:), allocatable :: fcut_all, dfcut_all 
      real(kind_double) :: fcut, dfcut, fcut_i2, dfcut_i2, fcut_i3, dfcut_i3
      real(kind_double) :: local_r_cut_in
      real(kind_double), dimension(:,:), allocatable  :: zz_tmp 
      integer :: type_fcut_in, type_fcut_out 
      real(kind_double) :: zz_k3
      real(kind_double) :: qz1, qz2, qz3, lz1, lz2, lz3 
      !TODOftnbody...
      integer :: izozo, type_db_i2, type_db_i3, itt
      integer, dimension(size(i_type_db)) :: iv1_zozo  
      integer, dimension(size(i_central)) :: iv2_zozo  
 
      izozo =  type_db_ja
      iv1_zozo = i_type_db 
      iv2_zozo = i_central
      
      allocate(tmpf(dim_desc_body(3)))
      allocate(fcut_all(max_neigh_local))
      allocate(dfcut_all(max_neigh_local))
      ! allocate(ttmp(dim_desc_body(3),3))

      !TODO_ZETABODY 
      type_fcut_in=3
      type_fcut_out=2
      ! norm = sqrt(1.d0)/sqrt(dble(dim_desc_body(3)))*delta_zetabody3(1)
      norm = delta_zetabody3(1)
      tmp_zeta_norm = norm 

      ! compute_qja_a_partir_de_rcentral 
      !dimf = dim_desc_body(3)
      !tmpf(:) = ftnbody_phase_random_3b

      if (max_neigh_local <= 1) return
      tmp_zeta_real(:) = 0.d0 
      if (desc_forces_local) then
        tmp_zeta_dreal(:, :, :) = 0.d0
        d2t(:) = 0.d0
        d3t(:) = 0.d0
    
        d2_x1(:) = 0.d0
        d2_x2(:) = 0.d0
        d2_x3(:) = 0.d0
        !
        d3_x1(:) = 0.d0
        d3_x2(:) = 0.d0
        d3_x3(:) = 0.d0
      end if

      do ii = 1, max_neigh_local
        ! Pair-specific inner cutoff
        local_r_cut_in = r_cut_pair_in(type_db_ja, i_type_db(ii))
        call fcut_rij_inout(type_fcut_in, type_fcut_out, r_central(ii), r_cut_z3b, r_cut_width_z3b, local_r_cut_in, r_cut_width_in, desc_forces_local, fcut, dfcut)
        fcut_all(ii) = fcut
        dfcut_all(ii) = dfcut
      end do

      do i2 = 1, max_neigh_local
        type_db_i2 = i_type_db(i2)    
        
#ifdef MLD_NDM
        ! ttmp(:,:)= 0.d0 
        !one_to_r_central_i2 = 1.d0/r_central(i2)
#endif
        fcut_i2 = fcut_all(i2)
        if (fcut_i2 == 0.d0) cycle
        dfcut_i2 = dfcut_all(i2)

        x1 = ur_central(i2)
        if (desc_forces_local) then
          d2_x1(1:3) = d_ur_central(1:3, i2)
        end if
        do i3 = i2 + 1, max_neigh_local
            type_db_i3 = i_type_db(i3)
            fcut_i3 = fcut_all(i3)
            if (fcut_i3 == 0.d0) cycle
            dfcut_i3 = dfcut_all(i3)
          ! if (i2/=i3) then  
            x2 = ur_central(i3)
            cos_2j3 = dot_product(tmp_dxp(1:3, i2), tmp_dxp(1:3, i3))/(r_central(i2)*r_central(i3))
            
            ! Fix for angular coordinate: use acos(cos_theta)
            if (cos_2j3 > 1.d0) cos_2j3 = 1.d0
            if (cos_2j3 < -1.d0) cos_2j3 = -1.d0
            x3 = acos(cos_2j3)

            if (desc_forces_local) then
              d3_x2(1:3) = d_ur_central(1:3, i3)
              
              ! Chain rule: d(theta)/d(cos) = -1/sin(theta) = -1/sqrt(1-cos^2)
              if (abs(cos_2j3) < 1.0d0 - 1.0d-12) then
                 d_theta_d_cos = -1.d0 / sqrt(1.d0 - cos_2j3**2)
              else
                 d_theta_d_cos = 0.d0 
              end if

              d2_x3(1:3) = (tmp_dxp(1:3, i3)/(r_central(i2)*r_central(i3)) - tmp_dxp(1:3, i2)*cos_2j3/r_central(i2)**2) * d_theta_d_cos
              d3_x3(1:3) = (tmp_dxp(1:3, i2)/(r_central(i2)*r_central(i3)) - tmp_dxp(1:3, i3)*cos_2j3/r_central(i3)**2) * d_theta_d_cos
            end if

            qsym(1) = x1 + x2
            qsym(2) = x1 * x2
            qsym(3) = x3 


            if (desc_forces_local) then
              ! 2 is in 1 - dist ;  3 - angles
              d2_qsym(1:3, 1) = d2_x1(1:3)
              d2_qsym(1:3, 2) = x2 * d2_x1(1:3)
              d2_qsym(1:3, 3) = d2_x3(1:3)
    
              ! 3 is in 2 - dist ;  3 - angles
              d3_qsym(1:3, 1) = d3_x2(1:3)
              d3_qsym(1:3, 2) = x1 * d3_x2(1:3)
              d3_qsym(1:3, 3) = d3_x3(1:3)
            end if

            itt = ispec3(type_db_ja, type_db_i2, type_db_i3)
            if (allocated(zz_tmp)) deallocate(zz_tmp)
            allocate(zz_tmp(size(grid_nbody(3)%zpoints_by_tuple(itt)%zz,1), size(grid_nbody(3)%zpoints_by_tuple(itt)%zz,2)))
            zz_tmp(:,:) = grid_nbody(3)%zpoints_by_tuple(itt)%zz(:,:)

            lz1 = 1.d0/length_zetabody3(1)**2
            lz2 = 1.d0/length_zetabody3(2)**2
            lz3 = 1.d0/length_zetabody3(3)**2
            tmp_ker3(:) = 0.d0
            do izz = grid_nbody(3)%zpoints_by_tuple(itt)%idx_min, grid_nbody(3)%zpoints_by_tuple(itt)%idx_max
              qz1 = qsym(1) - zz_tmp(1, izz)
              qz2 = qsym(2) - zz_tmp(2, izz)
              qz3 = qsym(3) - zz_tmp(3, izz)
              !$! write(*,*) 'debug --> qz1, qz2, qz3', qz1, qz2, qz3, lz1, lz2, lz3, tmp_zeta_norm 
              zz_k3 =  dexp(-0.5d0 *  (qz1**2  * lz1  + &
                                       qz2**2  * lz2  + &
                                       qz3**2  * lz3 )) * tmp_zeta_norm
              tmp_ker3(izz) =  zz_k3 * fcut_i2 * fcut_i3
              if (desc_forces_local) then
                do ix=1,3 
                  d2t(izz)  =      - qz1 * lz1  * zz_k3 *  d2_qsym(ix,1)  &
                                   - qz2 * lz2  * zz_k3 *  d2_qsym(ix,2)  &
                                   - qz3 * lz3  * zz_k3 *  d2_qsym(ix,3)

                  d3t(izz)  =      - qz1 * lz1  * zz_k3 *  d3_qsym(ix,1)  &
                                   - qz2 * lz2  * zz_k3 *  d3_qsym(ix,2)  &
                                   - qz3 * lz3  * zz_k3 *  d3_qsym(ix,3)
                 
                  tmp_zeta_dreal(izz, i2,ix) =  tmp_zeta_dreal(izz, i2,ix) + d2t(izz) * fcut_i2 * fcut_i3 + zz_k3 * fcut_i3 * dfcut_i2 * tmp_dxp(ix, i2)/r_central(i2)
                  tmp_zeta_dreal(izz, i3,ix) =  tmp_zeta_dreal(izz, i3,ix) + d3t(izz) * fcut_i2 * fcut_i3 + zz_k3 * fcut_i2 * dfcut_i3 * tmp_dxp(ix, i3)/r_central(i3)
                end do
              end if
            end do ! izz 
            tmp_zeta_real(:) = tmp_zeta_real(:) + tmp_ker3(:)
            !$! if (init_mode_ftnbody) then
            !$!   count_3b = count_3b + 1
            !$!   if (ltmp_mean_ftnbody) then 
            !$!     mean_3b(:) = mean_3b(:)  + qsym(:) 
            !$!   else 
            !$!     do iiii = 1, dim_qbody(3)
            !$!       do jjjj = 1, dim_qbody(3)
            !$!         covar_matrix_3b(iiii,jjjj) = covar_matrix_3b(iiii,jjjj) + ( qsym(iiii)  & 
            !$!                                      - mean_3b(iiii))* (qsym(jjjj) - mean_3b(jjjj)) 
            !$!       end do
            !$!     end do
            !$!   end if 
            !$!   cycle
            !$! end if


            !$! tmpf(:) = ftnbody_phase_random_3b(:)
            !$! call dgemv('T', dimq, dimf, one, ftnbody_omega_3b, dimq, qsym, incx, one, tmpf, incy)
            !$! !!! tmpf(:) = qsym(1:3)*ftnbody_omega_3b(1:3,:) + ftnbody_phase_random_3b(:)
            !$! tmp_real(:) = tmp_real(:) + cos(tmpf(:)) * fcut_i2 * fcut_i3
            !$! do ii = 1,dimq  
            !$!   sin_tmp(:,ii) =   -sin(tmpf(:)) * ftnbody_omega_3b(ii,:)
            !$! end do   

            !$! if (desc_forces_local) then
            !$!   do ix =1,3
            !$!     d2t(:) =                    sin_tmp(:,1) * d2_qsym(ix,1) & 
            !$!                               + sin_tmp(:,2) * d2_qsym(ix,2) &
            !$!                               + sin_tmp(:,3) * d2_qsym(ix,3)
            !$!     
            !$!     d3t(:) =                    sin_tmp(:,1) * d3_qsym(ix,1) & 
            !$!                               + sin_tmp(:,2) * d3_qsym(ix,2) &
            !$!                               + sin_tmp(:,3) * d3_qsym(ix,3)
            !$!     tmp_dreal(:,i2,ix) = tmp_dreal(:,i2,ix) + d2t(:) * fcut_i2 * fcut_i3 + cos(tmpf(:)) * fcut_i3 * dfcut_i2 * tmp_dxp(ix, i2)/r_central(i2)                  
            !$!     tmp_dreal(:,i3,ix) = tmp_dreal(:,i3,ix) + d3t(:) * fcut_i2 * fcut_i3 + cos(tmpf(:)) * fcut_i2 * dfcut_i3 * tmp_dxp(ix, i3)/r_central(i3)
            !$!   end do 
            !$! end if
          ! end if 
        end do
        ! do ix =1,3
        !    tmp_dreal(:,i2,ix) =  ttmp(:,ix)*norm
        ! end do 
      end do
      ! norm = norm!/2.d0      
      !$! tmp_real(:) = tmp_real(:)*norm
      !$! if (desc_forces_local) then 
      !$!   do kk = 1, max_neigh_local
      !$!     do ix =1, 3
      !$!       tmp_dreal(:,kk,ix) = tmp_dreal(:,kk,ix)*norm  
      !$!     end do       
      !$!   end do   ! end kk 
      !$! end if 

      deallocate(tmpf) 
    end subroutine zetabody_order_3


end module module_compute_zetabody_order



module module_compute_zetabody 
  contains 
  subroutine compute_zetabody(i_start_at, i_final_at, d_n_neigh, d_kind_neigh, iconf)
    use module_kind_variables, ONLY: kind_double
#ifdef MLD_NDM
    use gen_com_m, only : lperiod
    use gen_com_m_ml, only : imm,at,bg
    use tab_imm_m_ml, only : xp
#else
    use ondm_gen_com_m, only : imm, lperiod
    use ondm_tab_imm_m, only : iwmax2, xp
#endif
    use ml_in_ndm_module, only: imm_neigh, desc_forces, rangml
    use derived_types, only: config_real, config_desc
    use module_neigh_local, only: r_cut, & 
                                  max_neigh_local, i_central, i_type, i_type_db, r_central, iw2, &
                                  preallocate_neigh_ja, build_local_neighbours_ja, reallocate_neigh_ja, & 
                                  compute_ur_transformed_distances
    use time_check_general, only: debug_time, MY_MPI_WTIME
    use module_ftnbody, only: tnn_ftbd, t2b_ftbd, t3b_ftbd, init_mode_ftnbody
    use module_body_desc, only: l_body_order, dim_desc_body
    use module_compute_body_order, only: ftnbody_order_2, ftnbody_order_3, ftnbody_order_4, ftnbody_order_5
    use module_compute_zetabody_order, only: zetabody_order_2, zetabody_order_3
    use module_kernel_zetabody, only: zetabody_order, tmp_zeta_dreal, tmp_zeta_real, tmp_zeta_norm, delta_zetabody2, delta_zetabody3 
#ifdef MLD_NDM
    use notperiod_mod
#else
    use ondm_transform_coord, only: ondm_notperiod
#endif

    use mld_logger
    use ondm_transform_coord, only: ondm_notperiod
    
    implicit none 
    integer, intent(in)  :: i_start_at, i_final_at
    integer, dimension(imm), intent(out)   :: d_n_neigh
    integer, dimension(imm, imm_neigh), intent(out)    :: d_kind_neigh
    integer, optional    :: iconf

    integer  :: max_neigh
    logical :: desc_forces_local, small 
    real(kind_double), dimension(:, :), allocatable   :: xpnp, tmp_dxp, tmp_xp
    integer :: ja, ja_atom, type_db_ja, icnt, ii, ix, ia, iz  
    real(kind_double) :: t00, t11, t22, t33, temp_dja 
    ! real(kind_double), dimension(:), allocatable :: fcut_all, dfcut_all 
    ! real(kind_double) :: fcut, dfcut
    real(kind_double), dimension(:), allocatable  :: ur_central 
    real(kind_double), dimension(:,:), allocatable :: d_ur_central, cos_dxp   
  
    _NAMECURRENT_("compute_zetabody")

    desc_forces_local = desc_forces .and. (config_real(iconf)%has_force .or. config_real(iconf)%has_stress)
    if (init_mode_ftnbody) desc_forces_local = .false.
    max_neigh = imm_neigh

    if ((i_start_at == 0) .and. (i_final_at == 0)) then
      d_n_neigh(:) = 0
      d_kind_neigh(:, :) = 0
      config_desc(iconf)%energy(:, :) = 0.d0
      return
    end if

    small = .false.
    if (present(iconf)) then
      small = config_real(iconf)%small
    end if
  
    allocate (xpnp(3, imm))
    if (lperiod) then
      xpnp(:, :) = xp(:, :)
    else
#ifdef MLD_NDM
         call notperiod(imm,xp, xpnp,at,bg,.false.)
#else
      call ondm_notperiod(xp, xpnp)
#endif
    end if
    ! call cryst_to_cart (imm, xpnp, bg, -1)
  
    d_n_neigh(:) = 0
    d_kind_neigh(:, :) = 0
    config_desc(iconf)%energy(:, :) = 0.d0
    if (desc_forces_local) config_desc(iconf)%force(:, :, :, :) = 0.d0

    call preallocate_neigh_ja(r_central, i_type, i_central, tmp_dxp, tmp_xp, imm_neigh)

#ifdef MLD_NDM
!JPC    if (i_start_at == 1)  iw2 = 0
!JPC    if (i_start_at >  1)  iw2 = iwmax2(i_start_at - 1)
    iw2=0
#else
    if (i_start_at == 1)  iw2 = 0
    if (i_start_at >  1)  iw2 = iwmax2(i_start_at - 1)  
#endif
    do ja = i_start_at, i_final_at

      ja_atom = ja
      type_db_ja = config_real(iconf)%itype_db(ja)
      if (debug_time) t00 = MY_MPI_WTIME()
      call build_local_neighbours_ja( iconf, ja, imm, xpnp, r_cut, d_n_neigh, d_kind_neigh, &
                                      r_central, i_type, i_type_db, i_central, tmp_dxp, tmp_xp, &
                                      max_neigh_local, iw2)

      if (max_neigh_local == 0) cycle
      if (max_neigh_local > max_neigh) then
        if (rangml == 0) then
          call log_warning("ML: Fatal error in compute_cluster_bonds_and_angles the number of atoms. ")
          call log_warning("ML: The number of max_neig_local is bigger than max_neigh:  r_cut "//vtoa(r_cut)//&
                           " max_neigh_local "// vtoa(max_neigh)//" max_neigh_local "//vtoa(max_neigh_local))
          call log_warning("ML:possible solutions: decrease r_cut or increase max_neigh")
        end if
        call log_critical("ML: decrease rcut for that descriptor")
      end if

      call reallocate_neigh_ja (max_neigh_local, ja_atom, i_central, i_type, i_type_db, r_central, tmp_dxp)

      ! if (allocated(ur_central)) deallocate (ur_central); allocate (ur_central(max_neigh_local))
      ! if (desc_forces_local) then
      !   if (allocated(d_ur_central)) deallocate (d_ur_central); allocate (d_ur_central(3, max_neigh_local))
      ! end if

      call compute_ur_transformed_distances(desc_forces_local, max_neigh_local, r_central, tmp_dxp, cos_dxp, ur_central, d_ur_central)  
      
      ! do ii = 1, max_neigh_local
      !   call fcut_rij_inout(r_central(ii), r_cut, r_cut_width, r_cut_in, r_cut_width_in, desc_forces_local, fcut, dfcut)
      !   fcut_all(ii) = fcut
      !   dfcut_all(ii) = dfcut
      ! end do

      if (debug_time) then 
        t11 = MY_MPI_WTIME()
        tnn_ftbd = tnn_ftbd + t11 - t00    
      end if
      
      icnt = 0 
      if (l_body_order(2)) then 
        if (allocated(tmp_zeta_real)) deallocate (tmp_zeta_real); allocate (tmp_zeta_real(dim_desc_body(2)))
        if (desc_forces_local) then
          if (allocated(tmp_zeta_dreal)) deallocate (tmp_zeta_dreal); allocate (tmp_zeta_dreal(dim_desc_body(2), max_neigh_local, 3))
        end if 

        !call ftnbody_order_2(type_db_ja, max_neigh_local, i_type_db, i_central, r_central, tmp_dxp,   desc_forces_local)
        call zetabody_order_2(type_db_ja, max_neigh_local, i_type_db, i_central, r_central, tmp_dxp, ur_central, d_ur_central, desc_forces_local)
                
        do ii = 1, dim_desc_body(2)
          icnt = icnt  + 1
          config_desc(iconf)%energy(icnt, ja) = tmp_zeta_real(ii)
          if (desc_forces_local) then 
            do ix = 1, 3
              config_desc(iconf)%force(icnt, ja, 1:max_neigh_local, ix) =  tmp_zeta_dreal(ii, 1:max_neigh_local, ix)
            end do
          end if         
        end do


        do iz = 2, zetabody_order
          do ii = 1, dim_desc_body(2)
            icnt = icnt  + 1
            temp_dja = ( tmp_zeta_real(ii) / tmp_zeta_norm )**(iz-1)
            config_desc(iconf)%energy(icnt, ja) =  temp_dja * tmp_zeta_real(ii) / tmp_zeta_norm * delta_zetabody2(iz)
            
            if (desc_forces_local) then 
              do ix = 1, 3
                config_desc(iconf)%force(icnt, ja, 1:max_neigh_local, ix) =  dble(iz) * delta_zetabody2(iz) * &
                                     tmp_zeta_dreal(ii, 1:max_neigh_local, ix) * temp_dja / tmp_zeta_norm
              end do
            end if         
          end do 

        end do 
      end if   
      if (debug_time) then 
        t22 = MY_MPI_WTIME()
        t2b_ftbd = t2b_ftbd + t22 - t11    
      end if

      if (l_body_order(3)) then 
        if (allocated(tmp_zeta_real)) deallocate (tmp_zeta_real); allocate (tmp_zeta_real(dim_desc_body(3)))
        if (desc_forces_local) then
          if (allocated(tmp_zeta_dreal)) deallocate (tmp_zeta_dreal); allocate (tmp_zeta_dreal(dim_desc_body(3), max_neigh_local, 3))
        end if 
        call zetabody_order_3(type_db_ja, max_neigh_local, i_type_db, i_central, r_central, tmp_dxp, ur_central, d_ur_central, desc_forces_local)

        do ii = 1, dim_desc_body(3)
          icnt = icnt  + 1
          config_desc(iconf)%energy(icnt, ja) = tmp_zeta_real(ii)
          if (desc_forces_local) then 
            do ix = 1, 3
              config_desc(iconf)%force(icnt, ja, 1:max_neigh_local, ix) =  tmp_zeta_dreal(ii, 1:max_neigh_local, ix)
            end do
          end if         
        end do


        do iz = 2, zetabody_order
          do ii = 1, dim_desc_body(3)
            icnt = icnt  + 1
            temp_dja = ( tmp_zeta_real(ii) / tmp_zeta_norm )**(iz-1)
            config_desc(iconf)%energy(icnt, ja) =  temp_dja * tmp_zeta_real(ii) / tmp_zeta_norm * delta_zetabody3(iz)/dble(dim_desc_body(3))
            
            if (desc_forces_local) then 
              do ix = 1, 3
                config_desc(iconf)%force(icnt, ja, 1:max_neigh_local, ix) =  dble(iz) * delta_zetabody3(iz) /dble(dim_desc_body(3)) *&
                                     tmp_zeta_dreal(ii, 1:max_neigh_local, ix) * temp_dja / tmp_zeta_norm
              end do
            end if  
          end do          
        end do 

      end if   
      if (debug_time) then 
        t33 = MY_MPI_WTIME()
        t3b_ftbd = t3b_ftbd + t33 - t22    
      end if

      !$! if (l_body_order(4)) then 
      !$!   if (allocated(tmp_real)) deallocate (tmp_real); allocate (tmp_real(dim_desc_body(4)))
      !$!   if (desc_forces_local) then
      !$!     if (allocated(tmp_dreal)) deallocate (tmp_dreal); allocate (tmp_dreal(dim_desc_body(4), max_neigh_local, 3))
      !$!   end if 
      !$!   call ftnbody_order_4(type_db_ja, max_neigh_local, i_type_db, i_central, r_central, tmp_dxp, ur_central, d_ur_central, desc_forces_local)
      !$!   do ii = 1, dim_desc_body(4)
      !$!     icnt = icnt  + 1
      !$!     config_desc(iconf)%energy(icnt, ja) = tmp_real(ii)
      !$!     if (desc_forces_local) then 
      !$!       do ix = 1, 3
      !$!         config_desc(iconf)%force(icnt, ja, 1:max_neigh_local, ix) =  tmp_dreal(ii, 1:max_neigh_local, ix)
      !$!       end do
      !$!     end if         
      !$!   end do
      !$! end if   
      !$! if (debug_time) then 
      !$!   t44 = MY_MPI_WTIME()
      !$!   t4b_ftbd = t4b_ftbd + t44 - t33    
      !$! end if
  
      !$! if (l_body_order(5)) then 
      !$!   if (allocated(tmp_real)) deallocate (tmp_real); allocate (tmp_real(dim_desc_body(5)))
      !$!   if (desc_forces_local) then
      !$!     if (allocated(tmp_dreal)) deallocate (tmp_dreal); allocate (tmp_dreal(dim_desc_body(5), max_neigh_local, 3))
      !$!   end if 
      !$!   call ftnbody_order_5(type_db_ja, max_neigh_local, i_type_db, i_central, r_central, tmp_dxp, ur_central, d_ur_central,  desc_forces_local)
      !$!   do ii = 1, dim_desc_body(5)
      !$!     icnt = icnt  + 1
      !$!     config_desc(iconf)%energy(icnt, ja) = tmp_real(ii)
      !$!     if (desc_forces_local) then 
      !$!       do ix = 1, 3
      !$!         config_desc(iconf)%force(icnt, ja, 1:max_neigh_local, ix) =  tmp_dreal(ii, 1:max_neigh_local, ix)
      !$!       end do
      !$!     end if         
      !$!   end do
      !$! end if 
      !$! if (debug_time) then 
      !$!   t55 = MY_MPI_WTIME()
      !$!   t5b_ftbd = t5b_ftbd + t55 - t44    
      !$! end if

      if (desc_forces_local) then
        do ia = 1, max_neigh_local
          config_desc(iconf)%force(:, ja, 0, :) = config_desc(iconf)%force(:, ja, 0, :) - config_desc(iconf)%force(:, ja, ia, :)
        end do
      end if
      
    end do ! end_ja   

    _MLD_BEGIN_

    _MLD_END_ 

  end subroutine compute_zetabody
end module module_compute_zetabody 





