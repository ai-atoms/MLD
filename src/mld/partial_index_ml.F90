! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   MiLaDy (Machine Learning Dynamics)
! HND X
! HND X
! HND X   Milady was written and designed by Mihai-Cosmin Marinica, Alexandra M. Goryaeva
! HND X   Copyright 2015-2024.
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

module set_limits
   use ml_in_ndm_module, only: debug
   use mld_logger

   implicit none

   integer, dimension(:), allocatable     :: i_start_on_proc, i_end_on_proc, i_size_on_proc

contains

   !old! subroutine set_unlimit_for_atoms(rangml, im, i_start_at, i_final_at)
   !old!    implicit none
   !old!    integer, intent(in)  :: im, rangml
   !old!    integer, intent(out) :: i_start_at, i_final_at
   !old!    _NAMECURRENT_("set_unlimit_for_atoms")
   !old!    _MLD_BEGIN_
   !old!    i_start_at = 1
   !old!    i_final_at = im
   !old!    _MLD_END_
   !old! end subroutine set_unlimit_for_atoms

   !old! subroutine set_unlimit_for_cov(rangml, dim_train, i_final_cov, i_start_cov)
   !old!    implicit none
   !old!    integer, intent(in)  :: dim_train, rangml
   !old!    integer, intent(out) :: i_final_cov, i_start_cov
   !old!    _NAMECURRENT_("set_unlimit_for_cov")
   !old!    _MLD_BEGIN_
   !old!    i_start_cov = 1
   !old!    i_final_cov = dim_train
   !old!    _MLD_END_
   !old! end subroutine set_unlimit_for_cov


   subroutine set_limit_for_atoms_with_MPI_grid(im, i_start_at, i_final_at)

      use mld_subworld
      !use mpi

      implicit none
      integer, intent(in)  :: im                       ! numbers of atoms to be distributed on procs
      integer, intent(out) :: i_start_at, i_final_at
#if(PARA)
      integer  :: nratio1, nratio2, nrest
#endif
      integer  :: ip, ierr

      _NAMECURRENT_("set_limit_for_atoms_with_MPI_grid")


      _MLD_BEGIN_
#if(PARA)
      if (im >= subworld_size) then
         nrest = mod(im, subworld_size)
         nratio1 = im/subworld_size + 1
         nratio2 = im/subworld_size

         if (subrank < nrest) then
            i_start_at = subrank*nratio1 + 1
            if (subrank /= (subworld_size - 1)) i_final_at = (subrank + 1)*nratio1
         end if

         if (subrank >= nrest) then
            i_start_at = nrest*nratio1 + (subrank - nrest)*nratio2 + 1
            if (subrank /= (subworld_size - 1)) i_final_at = nrest*nratio1 + (subrank - nrest + 1)*nratio2
         end if

         if (subrank == (subworld_size - 1)) i_final_at = im
      else
         if (im == 1) then
            if (subrank == 0) then
               i_start_at = 1
               i_final_at = 1
            else
               i_start_at = 0
               i_final_at = 0
            end if
         else
            if (subrank + 1 <= im) then
               i_start_at = subrank + 1
               i_final_at = subrank + 1
            else
               i_start_at = 0
               i_final_at = 0
            end if
         end if
      end if


      if (allocated(i_start_on_proc)) deallocate (i_start_on_proc); allocate (i_start_on_proc(0:subworld_size - 1))
      if (allocated(i_end_on_proc)) deallocate (i_end_on_proc); allocate (i_end_on_proc(0:subworld_size - 1))
      if (allocated(i_size_on_proc)) deallocate (i_size_on_proc); allocate (i_size_on_proc(0:subworld_size - 1))
      i_start_on_proc(:) = 0
      i_end_on_proc(:) = 0
      i_start_on_proc(subrank) = i_start_at
      i_end_on_proc(subrank) = i_final_at
      call MPI_ALLREDUCE(MPI_IN_PLACE, i_start_on_proc, subworld_size, MPI_INTEGER, MPI_SUM, subworld, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, i_end_on_proc,   subworld_size, MPI_INTEGER, MPI_SUM, subworld, ierr)
      do ip = 0, subworld_size - 1                        ! TODO nb_procsml?
         i_size_on_proc(ip) = i_end_on_proc(ip) - i_start_on_proc(ip) + 1
      end do
#endif

      _MLD_END_
   end subroutine set_limit_for_atoms_with_MPI_grid

   subroutine set_limit_for_atoms(rangml, im, i_start_at, i_final_at)
#if(PARA)
      use mld_mpi, ONLY: nb_procsml
      implicit none
#else
      implicit none
#endif
      integer, intent(in)  :: im                       ! numbers of atoms to be distributed on procs
      integer, intent(inout)     :: rangml             ! rang of the proc for MPI
      integer, intent(out) :: i_start_at, i_final_at
#if(PARA)
      integer  :: nratio1, nratio2, nrest
#else
      integer  :: nb_procsml
#endif

      _NAMECURRENT_("set_limit_for_atoms")


      _MLD_BEGIN_
#if(PARA)
      if (im >= nb_procsml) then
         nrest = mod(im, nb_procsml)
         nratio1 = im/nb_procsml + 1
         nratio2 = im/nb_procsml

         if (rangml < nrest) then
            i_start_at = rangml*nratio1 + 1
            if (rangml /= (nb_procsml - 1)) i_final_at = (rangml + 1)*nratio1
         end if

         if (rangml >= nrest) then
            i_start_at = nrest*nratio1 + (rangml - nrest)*nratio2 + 1
            if (rangml /= (nb_procsml - 1)) i_final_at = nrest*nratio1 + (rangml - nrest + 1)*nratio2
         end if

         if (rangml == (nb_procsml - 1)) i_final_at = im
      else
         if (im == 1) then
            if (rangml == 0) then
               i_start_at = 1
               i_final_at = 1
            else
               i_start_at = 0
               i_final_at = 0
            end if
         else
            if (rangml + 1 <= im) then
               i_start_at = rangml + 1
               i_final_at = rangml + 1
            else
               i_start_at = 0
               i_final_at = 0
            end if
         end if
      end if

#else

      rangml = 0
      nb_procsml = 1
      i_start_at = 1
      i_final_at = im

#endif


      if (allocated(i_start_on_proc)) deallocate (i_start_on_proc); allocate (i_start_on_proc(0:nb_procsml - 1))
      if (allocated(i_end_on_proc)) deallocate (i_end_on_proc); allocate (i_end_on_proc(0:nb_procsml - 1))
      if (allocated(i_size_on_proc)) deallocate (i_size_on_proc); allocate (i_size_on_proc(0:nb_procsml - 1))
      i_start_on_proc(:) = 0
      i_end_on_proc(:) = 0
      i_start_on_proc(rangml) = i_start_at
      i_end_on_proc(rangml) = i_final_at

      _MLD_END_

   end subroutine set_limit_for_atoms


   subroutine set_limit_for_configs(rangml, im, i_start_conf, i_final_conf, iconf_start_on_proc, &
      iconf_end_on_proc, iconf_size_on_proc, iconf_to_proc)
      use mld_mpi, ONLY: nb_procsml,  comm_mld
      !use mpi
      implicit none
      integer, intent(in)  :: im           ! numbers of atoms to be distributed: im atoms over nb_procsml procs
      integer, intent(in)  :: rangml       ! rang of the proc for MPI
      integer, intent(out) :: i_start_conf, i_final_conf
      integer, dimension(:), allocatable, intent(inout) :: iconf_start_on_proc, iconf_end_on_proc, &
         iconf_size_on_proc, iconf_to_proc
      integer  :: nratio1, nratio2, nrest

      integer  :: ii, rr

      _NAMECURRENT_("set_limit_for_configs")
      _MLD_BEGIN_

      if (im >= nb_procsml) then
         nrest = mod(im, nb_procsml)
         nratio1 = im/nb_procsml + 1
         nratio2 = im/nb_procsml

         if (rangml < nrest) then
            i_start_conf = rangml*nratio1 + 1
            if (rangml /= (nb_procsml - 1)) i_final_conf = (rangml + 1)*nratio1
         end if

         if (rangml >= nrest) then
            i_start_conf = nrest*nratio1 + (rangml - nrest)*nratio2 + 1
            if (rangml /= (nb_procsml - 1)) i_final_conf = nrest*nratio1 + (rangml - nrest + 1)*nratio2
         end if

         if (rangml == (nb_procsml - 1)) i_final_conf = im
      else
         if (im == 1) then
            if (rangml == 0) then
               i_start_conf = 1
               i_final_conf = 1
            else
               i_start_conf = 0
               i_final_conf = 0
            end if
         else
            if (rangml + 1 <= im) then
               i_start_conf = rangml + 1
               i_final_conf = rangml + 1
            else
               i_start_conf = 0
               i_final_conf = 0
            end if
         end if
      end if

      if (allocated(iconf_start_on_proc)) deallocate (iconf_start_on_proc); allocate (iconf_start_on_proc(0:nb_procsml - 1))
      if (allocated(iconf_end_on_proc)) deallocate (iconf_end_on_proc); allocate (iconf_end_on_proc(0:nb_procsml - 1))
      if (allocated(iconf_size_on_proc)) deallocate (iconf_size_on_proc); allocate (iconf_size_on_proc(0:nb_procsml - 1))
      iconf_start_on_proc(:) = 0
      iconf_end_on_proc(:) = 0
      iconf_size_on_proc(:)=0
      iconf_start_on_proc(rangml) = i_start_conf
      iconf_end_on_proc(rangml) = i_final_conf
      if ((i_start_conf == i_final_conf).and.(i_start_conf==0)) then
         iconf_size_on_proc(rangml) = 0
      else
         iconf_size_on_proc(rangml) = i_final_conf - i_start_conf + 1
      end if
      !TORC! call MPI_ALLREDUCE(MPI_IN_PLACE, iconf_start_on_proc, nb_procsml, MPI_INTEGER, MPI_SUM, , codeml)
      call comm_mld%sum(iconf_start_on_proc)
      !TORC! call MPI_ALLREDUCE(MPI_IN_PLACE, iconf_end_on_proc, nb_procsml, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, codeml)
      call comm_mld%sum(iconf_end_on_proc)
      !TORC! call MPI_ALLREDUCE(MPI_IN_PLACE, iconf_size_on_proc, nb_procsml, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, codeml)
      call comm_mld%sum(iconf_size_on_proc)

      if (allocated(iconf_to_proc)) deallocate (iconf_to_proc); allocate (iconf_to_proc(im))
      !$! do ii = 1, im
      !$!   if ((ii <= i_final_conf).and.(ii>=i_start_conf)) then
      !$!     iconf_to_proc(ii) = rangml
      !$!   end if
      !$!
      !$! end do

      do ii = 1 , nb_procsml
         do rr = iconf_start_on_proc(ii-1), iconf_end_on_proc(ii-1)
            iconf_to_proc(rr) = ii-1
         end do
      end do

      _MLD_END_
   end subroutine set_limit_for_configs

   subroutine set_limit_for_configs_with_MPI_grid(im, i_start_conf, i_final_conf)
      use my_mpi_subroutines, only : subworlds_allreduce_int, subworlds_allreduce_vect_int
      use mld_subworld
      !use mpi

      implicit none
      integer, intent(in)  :: im                    ! numbers of atoms to be distributed on procs
      integer, intent(out) :: i_start_conf, i_final_conf
      integer  :: nratio1, nratio2, nrest


      _NAMECURRENT_("set_limit_for_configs_with_MPI_grid")


      _MLD_BEGIN_

      if (im >= nb_subworlds) then
         nrest = mod(im, nb_subworlds)
         nratio1 = im/nb_subworlds + 1
         nratio2 = im/nb_subworlds

         if (id_subworld < nrest) then
            i_start_conf = id_subworld*nratio1 + 1
            if (id_subworld /= (nb_subworlds - 1)) i_final_conf = (id_subworld + 1)*nratio1
         end if

         if (id_subworld >= nrest) then
            i_start_conf = nrest*nratio1 + (id_subworld - nrest)*nratio2 + 1
            if (id_subworld /= (nb_subworlds - 1)) i_final_conf = nrest*nratio1 + (id_subworld - nrest + 1)*nratio2
         end if

         if (id_subworld == (nb_subworlds - 1)) i_final_conf = im
      else
         if (im == 1) then
            if (id_subworld == 0) then
               i_start_conf = 1
               i_final_conf = 1
            else
               i_start_conf = 0
               i_final_conf = 0
            end if
         else
            if (id_subworld + 1 <= im) then
               i_start_conf = id_subworld + 1
               i_final_conf = id_subworld + 1
            else
               i_start_conf = 0
               i_final_conf = 0
            end if
         end if
      end if

      _MLD_END_
   end subroutine set_limit_for_configs_with_MPI_grid

   subroutine set_limit_for_cov(rangml, dim_train, i_final_cov, i_start_cov)
      ! set the limit for covariance matrix. all atoms are defined on all procs
      use mld_mpi, ONLY: nb_procsml
      implicit none
      integer, intent(in)  :: dim_train                ! dimension of the covariance matrix to be distributed
      integer, intent(inout)     :: rangml             ! the rang of the proc for MPI
      integer, intent(out) :: i_final_cov, i_start_cov ! local limits defined on each procs.
      integer  :: nratio1, nratio2, nrest

      _NAMECURRENT_("set_limit_for_cov")
      _MLD_BEGIN_

      nrest = mod(dim_train, nb_procsml)
      nratio1 = dim_train/nb_procsml + 1
      nratio2 = dim_train/nb_procsml
      if (rangml < nrest) then
         i_start_cov = rangml*nratio1 + 1
         if (rangml /= (nb_procsml - 1)) i_final_cov = (rangml + 1)*nratio1
      end if
      if (rangml >= nrest) then
         i_start_cov = nrest*nratio1 + (rangml - nrest)*nratio2 + 1
         if (rangml /= (nb_procsml - 1)) i_final_cov = nrest*nratio1 + (rangml - nrest + 1)*nratio2
      end if
      if (rangml == (nb_procsml - 1)) i_final_cov = dim_train

      _MLD_END_
   end subroutine set_limit_for_cov

end module set_limits
