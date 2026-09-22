! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   MiLaDy (Machine Learning Dynamics)
! HND X
! HND X
! HND X   Milady was written and designed by Mihai-Cosmin Marinica, Alexandra M. Goryaeva
! HND X   Copyright 2015-2023.
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

!TORC! module mld_subworld
!TORC! 
!TORC!    use mpi
!TORC!    implicit none
!TORC! 
!TORC!    integer  :: subworld, subrank, subworld_size
!TORC!    integer  :: masters_world, nb_subworlds, id_subworld
!TORC! 
!TORC! contains
!TORC! 
!TORC!    subroutine init_subworld(procs_per_file)
!TORC!       ! initialisation de la grille MPI
!TORC! 
!TORC!       use mld_mpi, ONLY : nb_procsml
!TORC!       use ondm_gen_com_m, ONLY: rangml
!TORC!       implicit none
!TORC! 
!TORC!       integer, intent(in) :: procs_per_file
!TORC!       integer :: color, key
!TORC!       integer :: ierr
!TORC! 
!TORC! 
!TORC!       ! the default value is the number of all procs in the world
!TORC!       color = rangml / procs_per_file
!TORC!       key = modulo(rangml,procs_per_file)
!TORC! 
!TORC!       !An initial communicator subworld is created by splitting MPI_COMM_WORLD using the color and key. 
!TORC!       !This means that all processes with the same color value will end up in the same subworld communicator, 
!TORC!       ! and within that communicator, the ranks will be ordered by key.
!TORC!       call MPI_Comm_split(MPI_COMM_WORLD, color, key, subworld, ierr)
!TORC!       call MPI_Comm_rank(subworld, subrank, ierr)
!TORC!       call MPI_Comm_size(subworld, subworld_size, ierr)
!TORC! 
!TORC!       !After the creation of subworld, the code redefines the color. 
!TORC!       !If the key is 0 (meaning the process is the first in its subworld), color is set to 0. 
!TORC!       !Otherwise, color is set to MPI_UNDEFINED. 
!TORC!       !This effectively separates out the first process in each subworld.
!TORC!       if(key == 0) then
!TORC!          color = 0
!TORC!       else
!TORC!          color = MPI_UNDEFINED
!TORC!       end if
!TORC! 
!TORC!       !Another call to MPI_Comm_split is made with the new color value. 
!TORC!       !This creates a new communicator masters_world that only includes the 
!TORC!       !processes that had color set to 0 (which are the first processes in their 
!TORC!       !respective subworlds, based on the previous key check). 
!TORC!       !If color is MPI_UNDEFINED, the process will not be part of the new masters_world communicator.
!TORC!       call MPI_Comm_split(MPI_COMM_WORLD, color, rangml, masters_world, ierr)
!TORC!       if (MPI_COMM_NULL /= masters_world) then
!TORC!          call MPI_Comm_rank(masters_world, id_subworld, ierr)
!TORC!          call MPI_Comm_size(masters_world, nb_subworlds, ierr)
!TORC!       end if
!TORC!       call MPI_BCAST(id_subworld, 1, MPI_INTEGER, 0, subworld, ierr)
!TORC!       call MPI_BCAST(nb_subworlds, 1, MPI_INTEGER, 0, subworld, ierr)
!TORC!       !So, masters_world is a communicator that groups the "master" (first) processes of 
!TORC!       !each division created by the procs_per_file scheme. 
!TORC!       !These "masters" can then communicate among themselves, which is often used for 
!TORC!       !tasks like aggregating data from each subworld or coordinating actions across subworlds.
!TORC! 
!TORC!    end subroutine init_subworld
!TORC! 
!TORC!    subroutine close_subworld()
!TORC!       implicit none
!TORC!       integer :: ierr
!TORC!       call MPI_Comm_free(subworld, ierr)
!TORC!       if (MPI_COMM_NULL /= masters_world) call MPI_Comm_free(masters_world, ierr)
!TORC!    end subroutine close_subworld
!TORC! 
!TORC!
!TORC! 
!TORC! end module mld_subworld
