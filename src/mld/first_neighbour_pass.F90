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

#include "../MLD_MACROS.INC"

module module_first_neighbour_pass
   implicit none
contains

!> Compute the actual max neighbour count for configuration icount
!! by looping over atoms and calling build_local_neighbours_ja (same
!! approach as compute_ace).  Must be called BEFORE descriptor allocations.
subroutine compute_max_n_neigh(icount, i_start_at, i_final_at)
   use iso_fortran_env, only: dp => real64
   use mpi
   use mld_mpi, only: codeml
#ifdef MLD_NDM
   use gen_com_m, ONLY: A2cm, lperiod
   use gen_com_m_ml, ONLY: imm, bg, at
   use tab_imm_m_ml, ONLY: xp
#else
   use ondm_gen_com_m, ONLY: imm, lperiod
   use ondm_tab_imm_m, ONLY: iwmax2, xp
#endif
   use derived_types, only: config_real
   use ml_in_ndm_module, ONLY: imm_neigh
   use module_neigh_local, only: r_central, i_central, i_type, i_type_db, &
                                 tmp_dxp, tmp_xp, max_neigh_local, iw2, &
                                 build_local_neighbours_ja
   use module_neigh_local, only: r_cut
   use mld_subworld, only: subworld
   use mld_logger
   use mld_string
#ifdef MLD_NDM
   use notperiod_mod
#else
   use ondm_transform_coord, only: ondm_notperiod
#endif

   implicit none

   integer, intent(in) :: icount, i_start_at, i_final_at
   integer :: ja, overall_max
   real(dp), dimension(:,:), allocatable :: xpnp
   integer, dimension(:), allocatable :: d_n_neigh
   integer, dimension(:,:), allocatable :: d_kind_neigh
   _NAMECURRENT_("compute_max_n_neigh")
   _MLD_BEGIN_

   overall_max = 0

   ! Procs with no atoms skip the neighbour loop but MUST still
   ! participate in the MPI_ALLREDUCE below (otherwise: deadlock).
   if (.not. ((i_start_at == 0) .and. (i_final_at == 0))) then

      ! Setup positions (same pattern as compute_ace)
      allocate(xpnp(3, imm))
      if (lperiod) then
         xpnp(:,:) = xp(:,:)
      else
#ifdef MLD_NDM
         call notperiod(imm, xp, xpnp, at, bg, .false.)
#else
         call ondm_notperiod(xp, xpnp)
#endif
      end if

      ! Temporary neighbour arrays (same dimensions as in compute_ace)
      allocate(d_n_neigh(imm))
      allocate(d_kind_neigh(imm, imm_neigh))
      d_n_neigh(:) = 0
      d_kind_neigh(:,:) = 0

      if (i_start_at == 1) iw2 = 0
#ifdef MLD_NDM
      iw2 = 0
#else
      if (i_start_at > 1) iw2 = iwmax2(i_start_at - 1)
#endif

      do ja = i_start_at, i_final_at
         call build_local_neighbours_ja(icount, ja, imm, xpnp, r_cut, &
                                        d_n_neigh, d_kind_neigh, &
                                        r_central, i_type, i_type_db, i_central, &
                                        tmp_dxp, tmp_xp, &
                                        max_neigh_local, iw2)
         if (max_neigh_local > overall_max) overall_max = max_neigh_local
      end do

      deallocate(xpnp, d_n_neigh, d_kind_neigh)

   end if

#if(PARA)
   call MPI_ALLREDUCE(MPI_IN_PLACE, overall_max, 1, MPI_INTEGER, MPI_MAX, subworld, codeml)
#endif
   config_real(icount)%max_n_neigh = overall_max
!   call log_info('max_n_neigh for config '//vtoa(icount)//' = '//vtoa(overall_max))

   _MLD_END_
end subroutine compute_max_n_neigh

end module module_first_neighbour_pass
