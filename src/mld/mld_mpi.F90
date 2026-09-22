
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

#include "../MLD_MACROS.INC"

module mpi_tools_by_FB
  !-------------------------------------------------------------
  !
  ! MPI TOOLS after an original ideea of Dr. Fabien Bruneval
  !
  ! How to CITE? MolGW 2.0: A new software for quantum chemistry
  !------------------------------------------------------------
   use mpi
   use mld_logger
   use, intrinsic :: iso_fortran_env, dp=>real64
   implicit none

   type, public :: mpi_communicator
     integer    :: comm       ! MPI communicator
     integer    :: nproc      ! number of procs in the communicator comm
     integer    :: rank       ! index           in the communicator comm
     contains
     !init 
     procedure :: init => mpic_init
     procedure :: barrier => mpic_barrier
     !procedure :: finalize => mpic_finalize
     ! broadcast
     generic :: bcast  => mpic_bcast_dp
     generic :: bcast  => mpic_bcast_cdp
     generic :: bcast  => mpic_bcast_i
     generic :: bcast  => mpic_bcast_vali
     generic :: bcast  => mpic_bcast_log
     generic :: bcast  => mpic_bcast_char
     procedure :: mpic_bcast_dp
     procedure :: mpic_bcast_cdp
     procedure :: mpic_bcast_i
     procedure :: mpic_bcast_vali
     procedure :: mpic_bcast_log
     procedure :: mpic_bcast_char
     ! sum 
     generic :: sum  => mpic_sum_dp
     generic :: sum  => mpic_sum_cdp
     generic :: sum  => mpic_sum_i
     procedure :: mpic_sum_dp
     procedure :: mpic_sum_cdp
     procedure :: mpic_sum_i

   end type mpi_communicator

   contains

   !=========================================================================
   subroutine mpic_init(this, comm_in)
     implicit none
   
     class(mpi_communicator),intent(inout) :: this
     integer,intent(in)                    :: comm_in
     !=====
     integer :: ierror
     !=====
     _NAMECURRENT_('mpic_init')
   
     ! initialize locals so sanitizers are happy
     ierror     = MPI_SUCCESS
     this%comm  = comm_in
     this%nproc = -1
     this%rank  = -1
   
     this%comm = comm_in
     
     call MPI_COMM_SIZE(this%comm,this%nproc,ierror)
     if (ierror /= MPI_SUCCESS) then
       call log_critical('error in MPI_COMM_SIZE in '//NAMECURRENT)
     end if   
     
     call MPI_COMM_RANK(this%comm,this%rank,ierror)
     if (ierror /= MPI_SUCCESS) then
       call log_critical('error in MPI_COMM_RANK in '//NAMECURRENT)
     end if

   end subroutine mpic_init

   !=========================================================================
   subroutine mpic_barrier(this)
     implicit none
   
     class(mpi_communicator),intent(in) :: this
     !=====
     integer  :: ierror = 0 
     !=====
      _NAMECURRENT_('mpic_barrier')
   
     call MPI_BARRIER(this%comm, ierror)
     if (ierror /= 0) then
       call log_critical('error in MPI_BARRIER in '//NAMECURRENT)
     end if 

   end subroutine mpic_barrier

   !=========================================================================
   subroutine mpic_bcast_i(this, rank, array)
     implicit none
     class(mpi_communicator), intent(in) :: this
     integer,intent(in)    :: rank
     integer,intent(inout) :: array(:)
     !=====
     integer :: nsize
     integer :: ierror=0
     !=====
     _NAMECURRENT_('mpic_bcast_i')
   
     if( this%nproc == 1 ) return
   
     nsize = SIZE(array)
     call MPI_BCAST(array, nsize, MPI_INTEGER, rank, this%comm, ierror)
     if( ierror /= 0 ) then
       call log_critical('error in MPI_BCAST in '//NAMECURRENT)
     endif
   
   end subroutine mpic_bcast_i


   !=========================================================================
   subroutine mpic_bcast_vali(this, rank, val)
     implicit none
     class(mpi_communicator), intent(in) :: this
     integer,intent(in)    :: rank
     integer,intent(inout) :: val
     !=====
     integer :: nsize
     integer :: ierror=0
     !=====
     _NAMECURRENT_('mpic_bcast_vali')
   
     if( this%nproc == 1 ) return
   
     nsize = 1
     call MPI_BCAST(val, nsize, MPI_INTEGER, rank, this%comm, ierror)
     if( ierror /= 0 ) then
       call log_critical('error in MPI_BCAST in '//NAMECURRENT)
     endif
   
   end subroutine mpic_bcast_vali

   !=========================================================================
   subroutine mpic_bcast_dp(this, rank, array)
     implicit none
     class(mpi_communicator),intent(in) :: this
     integer,intent(in)     :: rank
     real(dp),intent(inout) :: array(..)
     !=====
     integer :: nsize
     integer :: ierror=0
     !=====
     _NAMECURRENT_('mpic_bcast_dp')
   
     if( this%nproc == 1 ) return
   
     nsize = SIZE(array)
   
     call MPI_BCAST(array, nsize, MPI_DOUBLE_PRECISION, rank, this%comm, ierror)
   
     if( ierror /= 0 ) then
       call log_critical('error in MPI_BCAST in '//NAMECURRENT)
     endif
   
   end subroutine mpic_bcast_dp

   !=========================================================================
   subroutine mpic_bcast_cdp(this, rank, array)
     implicit none
     class(mpi_communicator),intent(in) :: this
     integer,intent(in)     :: rank
     complex(dp),intent(inout) :: array(..)
     !=====
     integer :: nsize
     integer :: ierror=0
     !=====
     _NAMECURRENT_('mpic_bcast_cdp')
   
     if( this%nproc == 1 ) return
   
     nsize = SIZE(array)
   
     call MPI_BCAST(array,nsize, MPI_DOUBLE_COMPLEX, rank, this%comm, ierror)
   
     if( ierror /= 0 ) then
       call log_critical('error in MPI_BCAST in '//NAMECURRENT)
     endif
   
   end subroutine mpic_bcast_cdp

      !=========================================================================
   subroutine mpic_bcast_log(this, rank, larray)
     implicit none
     class(mpi_communicator),intent(in) :: this
     integer, intent(in)     :: rank
     logical, intent(inout) :: larray(..)
     !=====
     integer :: nsize
     integer :: ierror=0
     !=====
     _NAMECURRENT_('mpic_bcast_log')
   
     if( this%nproc == 1 ) return
   
     nsize = SIZE(larray)
   
     call MPI_BCAST(larray, nsize, MPI_LOGICAL, rank, this%comm, ierror)
   
     if( ierror /= 0 ) then
       call log_critical('error in MPI_BCAST in '//NAMECURRENT)
     endif
   
   end subroutine mpic_bcast_log

      subroutine mpic_bcast_char(this, rank, carray)
     implicit none
     class(mpi_communicator),intent(in) :: this
     integer, intent(in)     :: rank
     character(len=*), intent(inout) :: carray
      
     !=====
     integer :: nsize
     integer :: ierror=0
     !=====
     _NAMECURRENT_('mpic_bcast_char')
   
     if( this%nproc == 1 ) return
   
     nsize = len(carray)
   
     call MPI_BCAST(carray, nsize, MPI_CHARACTER, rank, this%comm, ierror)
   
     if( ierror /= 0 ) then
       call log_critical('error in MPI_BCAST in '//NAMECURRENT)
     endif
   
   end subroutine mpic_bcast_char

   !=========================================================================
   subroutine mpic_sum_dp(this, array)
     implicit none
     class(mpi_communicator),intent(in) :: this
     real(dp),intent(inout) :: array(..)
     !=====
     integer :: nsize
     integer :: ierror=0
     !=====
     _NAMECURRENT_('mpic_sum_dp')
   
     if( this%nproc == 1 ) return
   
     nsize = SIZE(array)

     call MPI_ALLREDUCE( MPI_IN_PLACE, array, nsize, MPI_DOUBLE_PRECISION, MPI_SUM, this%comm, ierror)
     if (ierror /= 0) then
       call log_critical('error in MPI_ALLREDUCE in '//NAMECURRENT)
     end if   
   
   end subroutine mpic_sum_dp


   !=========================================================================
   subroutine mpic_sum_cdp(this, array)
     implicit none
     class(mpi_communicator),intent(in) :: this
     complex(dp),intent(inout) :: array(..)
     !=====
     integer :: nsize
     integer :: ierror=0
     !=====
     _NAMECURRENT_('mpic_sum_cdp')
   
     if( this%nproc == 1 ) return
   
     nsize = SIZE(array)
   
     call MPI_ALLREDUCE( MPI_IN_PLACE, array, nsize, MPI_DOUBLE_COMPLEX, MPI_SUM, this%comm, ierror)
     if (ierror /= 0) then
       call log_critical('error in MPI_ALLREDUCE in '//NAMECURRENT)
     end if
   
   end subroutine mpic_sum_cdp

   !=========================================================================
   subroutine mpic_sum_i(this, array)
     implicit none
     class(mpi_communicator),intent(in) :: this
     integer, intent(inout) :: array(..)
     !=====
     integer :: nsize
     integer :: ierror=0
     !=====
     _NAMECURRENT_('mpic_sum_i')
     
     if( this%nproc == 1 ) return
   
     nsize = SIZE(array)
   
     call MPI_ALLREDUCE( MPI_IN_PLACE, array, nsize, MPI_INTEGER, MPI_SUM, this%comm, ierror)

     if( ierror /= 0 ) then
        call log_critical('error in MPI_ALLREDUCE in '//NAMECURRENT)
     endif
   
   end subroutine mpic_sum_i
   
   
   end module mpi_tools_by_FB


!-----------------------------!
!       MODULE mld_mpi        !
!-----------------------------!

module mld_mpi
  use module_kind_variables, only: kind_double
  use mpi
  use mld_string
  use mld_logger
#ifdef MLD_NDM
  use gen_com_m, ONLY: rangml
#else
  use ondm_gen_com_m, ONLY: rangml
#endif
  use mpi_tools_by_FB, ONLY: mpi_communicator

  implicit none

  !old mod_mpi_ml stuff + mpi_comm_mld 
  integer :: mpi_comm_mld = MPI_COMM_NULL
  integer :: nb_procsml = -1 , codeml = -1  
  integer, dimension(MPI_STATUS_SIZE)    :: statut_ml = 0  

  integer  :: &
    mld_rank = -1, &
    mld_size = -1, &
    mld_exit_failure = 1, &
    mld_ierror = 0 , &
    mld_status(mpi_status_size) = 0 

  logical, private     :: is_mpi_init_done = .false.

  real(kind_double) :: temps_deb = 0.0_kind_double 
  integer :: myid = -1 , nprocs = -1  
  type(mpi_communicator) :: comm_mld 

contains

  subroutine mld_mpi_init()
  use mpi
  implicit none
  logical :: flag = .false.        ! ✅ initialize to a defined value
  integer :: i2, len
  character(len=MPI_MAX_ERROR_STRING) :: es
  character(len=32) :: s_rank, s_size

  if (is_mpi_init_done) return

  call MPI_Initialized(flag, mld_ierror)
  if (mld_ierror /= MPI_SUCCESS) then
     call MPI_Error_string(mld_ierror, es, len, i2)
     write(*,*) 'MPI_Initialized failed: ', trim(es(1:len))
     call MPI_Abort(MPI_COMM_WORLD, 1, i2)
  end if

#ifdef MLD_NDM  
  if (.not. flag) then
     write(*,*) 'ERROR: MPI not initialized in NDM path; expected external init.'
     call MPI_Abort(MPI_COMM_WORLD, 1, mld_ierror)
  end if
#else 
  if (.not. flag) then
     call MPI_Init(mld_ierror)
     if (mld_ierror /= MPI_SUCCESS) then
        call MPI_Error_string(mld_ierror, es, len, i2)
        write(*,*) 'MPI_Init failed: ', trim(es(1:len))
        call MPI_Abort(MPI_COMM_WORLD, 1, i2)
     end if
  end if
#endif       


    ! Rank/size are needed in BOTH branches
    call MPI_Comm_size(MPI_COMM_WORLD, mld_size, mld_ierror)
    if (mld_ierror /= MPI_SUCCESS) then
       call MPI_Error_string(mld_ierror, es, len, i2)
       write(*,*) 'MPI_Comm_size failed: ', trim(es(1:len))
       call MPI_Abort(MPI_COMM_WORLD, 1, i2)
    end if

    call MPI_Comm_rank(MPI_COMM_WORLD, mld_rank, mld_ierror)
    if (mld_ierror /= MPI_SUCCESS) then
       call MPI_Error_string(mld_ierror, es, len, i2)
       write(*,*) 'MPI_Comm_rank failed: ', trim(es(1:len))
       call MPI_Abort(MPI_COMM_WORLD, 1, i2)
    end if

    ! You mentioned the missing DUP — ensure it’s done in BOTH paths
    call MPI_Comm_dup(MPI_COMM_WORLD, mpi_comm_mld, mld_ierror)
    if (mld_ierror /= MPI_SUCCESS) then
       call MPI_Error_string(mld_ierror, es, len, i2)
       write(*,*) 'MPI_Comm_dup failed: ', trim(es(1:len))
       call MPI_Abort(MPI_COMM_WORLD, 1, i2)
    end if

    call comm_mld%init(mpi_comm_mld)

    ! Mirror legacy fields
    logger_mpi_rank = mld_rank
    nb_procsml      = mld_size
    rangml          = mld_rank
    myid            = mld_rank
    nprocs          = mld_size
    temps_deb       = MPI_Wtime()
  
    ! Safe logging (avoid vtoa trailing garbage)
    write(s_rank,'(I0)') mld_rank
    write(s_size,'(I0)') mld_size
    call log_info('mld_mpi_init ! hello world ! from rank '//trim(s_rank)//' of '//trim(s_size))
    is_mpi_init_done = .true.
end subroutine mld_mpi_init


!#ifdef MLD_NDM  
!      call mpi_comm_dup(mpi_comm_world, mpi_comm_mld, mld_ierror)
!      call comm_mld%init(mpi_comm_mld)
!#else 
!      call mpi_init(mld_ierror)
!      call mpi_comm_size(mpi_comm_world, mld_size, mld_ierror)
!      call mpi_comm_rank(mpi_comm_world, mld_rank, mld_ierror)
!      call mpi_comm_dup(mpi_comm_world, mpi_comm_mld, mld_ierror)
!      call comm_mld%init(mpi_comm_mld)
!#endif       
!      logger_mpi_rank = mld_rank                       ! avoid circular dependency
!      nb_procsml = mld_size   ! TODO remove that
!      rangml = mld_rank       ! TODO remove that
!      myid = mld_rank         ! TODO remove that
!      nprocs = mld_size       ! TODO remove that
!      temps_deb = MPI_Wtime() ! TODO remove that
!      call log_info("mld_mpi_init ! hello world ! from rank "//vtoa(mld_rank)//"of "//vtoa(mld_size))
!    end if
!  end subroutine

  subroutine mld_mpi_finalize(mess)
    ! finish mpi stuff on regular situation
    character(*), optional     :: mess
    call log_info("mld_mpi_finalize ! goodbye world ! from rank "//vtoa(mld_rank)//"of "//vtoa(mld_size))
    if (present(mess)) call log_info(mess)
    call mpi_finalize(mld_ierror)
  end subroutine

  function mld_is_mpi_init_done() result(res)
    ! TODO assume if not mpi compiled, but milady is ever mpi
    logical  :: res
    res = is_mpi_init_done
  end function


  subroutine mld_mpi_abort(mess)
    ! finish mpi stuff on problem for log(s)
    ! never mind finish stuff milady other opened files etc
    character(*)   :: mess
    call log_critical("mld_mpi_abort:"//nwl//mess)
    call log_flush()
    call mpi_barrier(mpi_comm_mld, mld_ierror)     ! may be useful, or not ?
    call mpi_abort(mpi_comm_mld, mld_exit_failure, mld_ierror)
    ! stop 1    ! bash error ko useful TODO ?
  end subroutine

  subroutine mld_critical_abort(mess)
    ! general finish mpi stuff on problem
    ! log a final critical message,
    ! try to finish stuff on opened files etc TODO
    ! then stop 1 as bash error TODO ?
    ! TODO use mld_unit, only: str_opened_files
    use mpi
    character(*)   :: mess
    integer  :: exit_failure = 1, ierror

    call log_critical(mess)
    !call log_info('Stay opened files'//nwl//str_opened_files())
    ! try to finish stuff on opened files etc TODO
    call mpi_abort(mpi_comm_mld, exit_failure, ierror)
    ! stop 1    ! ko TODO?
  end subroutine

end module mld_mpi 

!-----------------------------!
!       MODULE mld_subworld   !
!-----------------------------!


module mld_subworld

   use mpi
   implicit none

   integer  :: subworld, subrank, subworld_size
   integer  :: masters_world, nb_subworlds, id_subworld
   integer, dimension(:), allocatable :: list_of_masters_in_mld

contains

   subroutine init_subworld(procs_per_file)
      ! initialisation de la grille MPI

      use mld_mpi, ONLY : mpi_comm_mld
#ifdef MLD_NDM
      use gen_com_m, ONLY: rangml
#else
      use ondm_gen_com_m, ONLY: rangml
#endif
      implicit none

      integer, intent(in) :: procs_per_file
      integer :: color, key
      integer :: ierr, tmp_id 


      ! the default value is the number of all procs in the world
      color = rangml / procs_per_file
      key = modulo(rangml,procs_per_file)

      !An initial communicator subworld is created by splitting MPI_COMM_WORLD using the color and key. 
      !This means that all processes with the same color value will end up in the same subworld communicator, 
      ! and within that communicator, the ranks will be ordered by key.
      !TORC! call MPI_Comm_split(MPI_COMM_WORLD, color, key, subworld, ierr)
      call MPI_Comm_split(mpi_comm_mld, color, key, subworld, ierr)
      call MPI_Comm_rank(subworld, subrank, ierr)
      call MPI_Comm_size(subworld, subworld_size, ierr)

      !After the creation of subworld, the code redefines the color. 
      !If the key is 0 (meaning the process is the first in its subworld), color is set to 0. 
      !Otherwise, color is set to MPI_UNDEFINED. 
      !This effectively separates out the first process in each subworld.
      if(key == 0) then
         color = 0
      else
         color = MPI_UNDEFINED
      end if

      !Another call to MPI_Comm_split is made with the new color value. 
      !This creates a new communicator masters_world that only includes the 
      !processes that had color set to 0 (which are the first processes in their 
      !respective subworlds, based on the previous key check). 
      !If color is MPI_UNDEFINED, the process will not be part of the new masters_world communicator.

      !TORC! call MPI_Comm_split(MPI_COMM_WORLD, color, rangml, masters_world, ierr)
      call MPI_Comm_split(mpi_comm_mld, color, rangml, masters_world, ierr)
      if (MPI_COMM_NULL /= masters_world) then
         call MPI_Comm_rank(masters_world, id_subworld, ierr)
         call MPI_Comm_size(masters_world, nb_subworlds, ierr)
      end if
      call MPI_BCAST(id_subworld, 1, MPI_INTEGER, 0, subworld, ierr)
      call MPI_BCAST(nb_subworlds, 1, MPI_INTEGER, 0, subworld, ierr)
      !So, masters_world is a communicator that groups the "master" (first) processes of 
      !each division created by the procs_per_file scheme. 
      !These "masters" can then communicate among themselves, which is often used for 
      !tasks like aggregating data from each subworld or coordinating actions across subworlds.
      if (allocated(list_of_masters_in_mld)) deallocate(list_of_masters_in_mld)
      allocate(list_of_masters_in_mld(0:nb_subworlds-1))
      
      tmp_id = 0 
      if (key == 0) then
          ! Master processes send their global rank
          list_of_masters_in_mld(id_subworld) = rangml
          tmp_id = rangml 
          !call MPI_GATHER(rangml, 1, MPI_INTEGER, list_of_masters_in_mld, 1, MPI_INTEGER, 0, mpi_comm_mld, ierr)
      end if 
      list_of_masters_in_mld(:) = 0
      list_of_masters_in_mld(id_subworld) = tmp_id
      !write(*,*) 'id', rangml, id_subworld, "ll", list_of_masters_in_mld(:) 
      if (MPI_COMM_NULL /= masters_world) call MPI_ALLREDUCE(MPI_IN_PLACE, list_of_masters_in_mld, nb_subworlds, &
                                                            MPI_INTEGER, MPI_SUM, masters_world, ierr)
      call MPI_BCAST(list_of_masters_in_mld, nb_subworlds, MPI_INTEGER, 0, subworld, ierr)


   end subroutine init_subworld

   subroutine close_subworld()
      implicit none
      integer :: ierr
      call MPI_Comm_free(subworld, ierr)
      if (MPI_COMM_NULL /= masters_world) call MPI_Comm_free(masters_world, ierr)
   end subroutine close_subworld

end module mld_subworld


!----------------------------------------!
!       MODULE my_mpi_subroutines        !
!----------------------------------------!


module my_mpi_subroutines
   use module_kind_variables, only: kind_double
   use mpi
contains

   subroutine my_barrier_subworld(ierr)
      ! barrier for mld
      use mld_subworld, ONLY: subworld
      implicit none
      integer, intent(inout) :: ierr
      call MPI_BARRIER(subworld, ierr)
   end subroutine my_barrier_subworld

   subroutine my_broadcast_val_int (val, no_proc, mpi_comm, mpi_grid)
      ! no_proc  -> proc from where the broadcast is made
      ! mpi_comm -> mpi communicator
      ! mpi_gris -> mpi grid
      implicit none
      integer, intent(inout)  :: val
      integer :: no_proc, mpi_comm, mpi_grid
      integer, parameter :: one=1

      call MPI_BCAST(val, one, MPI_INTEGER, no_proc, mpi_comm, mpi_grid)
   end subroutine my_broadcast_val_int

   subroutine my_broadcast_val_real (val, no_proc, mpi_comm, mpi_grid)
      ! no_proc  -> proc from where the broadcast is made
      ! mpi_comm -> mpi communicator
      ! mpi_gris -> mpi grid
      implicit none
      real(kind_double), intent(inout)  :: val
      integer :: no_proc, mpi_comm, mpi_grid
      integer, parameter :: one=1

      call MPI_BCAST(val, one, MPI_DOUBLE_PRECISION, no_proc, mpi_comm, mpi_grid)
   end subroutine my_broadcast_val_real

   subroutine my_broadcast_vect_real (val, no_proc, mpi_comm, mpi_grid)
      ! no_proc  -> proc from where the broadcast is made
      ! mpi_comm -> mpi communicator
      ! mpi_gris -> mpi grid
      implicit none
      real(kind_double), dimension(:), intent(inout)  :: val
      integer :: no_proc, mpi_comm, mpi_grid

      call MPI_BCAST(val, size(val,1), MPI_DOUBLE_PRECISION, no_proc, mpi_comm, mpi_grid)
   end subroutine my_broadcast_vect_real

   subroutine my_broadcast_vect_int (val, no_proc, mpi_comm, mpi_grid)
      ! no_proc  -> proc from where the broadcast is made
      ! mpi_comm -> mpi communicator
      ! mpi_gris -> mpi grid
      implicit none
      integer, dimension(:), intent(inout)  :: val
      integer :: no_proc, mpi_comm, mpi_grid

      call MPI_BCAST(val, size(val,1), MPI_INTEGER, no_proc, mpi_comm, mpi_grid)
   end subroutine my_broadcast_vect_int


   subroutine my_broadcast_char (val, no_proc, mpi_comm, mpi_grid)
      ! no_proc  -> proc from where the broadcast is made
      ! mpi_comm -> mpi communicator
      ! mpi_gris -> mpi grid
      implicit none
      character(len=* ), intent(inout)  :: val
      integer :: no_proc, mpi_comm, mpi_grid
      integer  :: length

      length = len(val)
      call MPI_BCAST(val, length, MPI_CHARACTER, no_proc, mpi_comm, mpi_grid)

   end subroutine my_broadcast_char

   subroutine subworlds_allreduce_matrix_double(namemat, val)
      ! reduction between all the masters of subworlds, then spreading info in each subworld (for a matrix of double)
!old #ifdef MLD_NDM
!old       use gen_com_m, ONLY: rangml
!old #else
!old       use ondm_gen_com_m, ONLY: rangml
!old #endif
      use mld_subworld
      implicit none
      real(kind=kind(0.d0)), dimension(:,:) :: val
      character(len=*), intent(in)  :: namemat 
      integer :: ierr
      character(len=60) :: czozo 
      czozo = trim(namemat)
       
      if (MPI_COMM_NULL /= masters_world) then 
         !$! ! Check for floating underflow
         !$! if (any(val < 1e-10)) then
         !$!    do i = 1, size(val, 1)
         !$!       do j = 1, size(val, 2)
         !$!          if ((abs(val(i,j)) < 1e-40).and.(abs(val(i,j)) > 0.d0)) then
         !$!             write(*,*) trim(namemat), " Underflow detected at index (", i, ",", j, ")", val(i,j)
         !$!          end if
         !$!       end do
         !$!    end do
         !$! end if
         call MPI_ALLREDUCE(MPI_IN_PLACE, val, size(val,1)*size(val,2), MPI_DOUBLE_PRECISION, MPI_SUM, masters_world, ierr)
      end if 

      call MPI_BCAST(val, size(val,1)*size(val,2), MPI_DOUBLE_PRECISION, 0, subworld, ierr)
   end subroutine subworlds_allreduce_matrix_double

   subroutine subworlds_allreduce_vect_double(val)
      ! reduction between all the masters of subworlds, then spreading info in each subworld (for an array of double)
      use mld_subworld
      implicit none
      real(kind=kind(0.d0)), dimension(:) :: val
      integer :: ierr
      if (MPI_COMM_NULL /= masters_world) call MPI_ALLREDUCE(MPI_IN_PLACE, val, size(val,1), MPI_DOUBLE_PRECISION, MPI_SUM, masters_world, ierr)
      call MPI_BCAST(val, size(val,1), MPI_DOUBLE_PRECISION, 0, subworld, ierr)
   end subroutine subworlds_allreduce_vect_double

   subroutine subworlds_allreduce_vect_logical(val)
   ! Reduction between all the masters of subworlds, then spreading info in each subworld (for an array of logical)
   use mld_subworld
   use mpi
   implicit none
   logical, dimension(:), intent(inout) :: val
   integer, dimension(size(val)) :: int_val
   integer :: ierr

   ! Convert logical array to integer array (1 for .true., 0 for .false.)
   int_val = merge(1, 0, val)
   ! Perform the MPI reduction on the integer representation
   if (MPI_COMM_NULL /= masters_world) call MPI_ALLREDUCE(MPI_IN_PLACE, int_val, size(val), MPI_INTEGER, MPI_SUM, masters_world, ierr)
   ! Broadcast the result to all subworlds
   call MPI_BCAST(int_val, size(val), MPI_INTEGER, 0, subworld, ierr)
   ! Convert the integer array result back to logical array (any non-zero value is .true., 0 is .false.)
   val = (int_val /= 0)
end subroutine subworlds_allreduce_vect_logical


   subroutine subworlds_allreduce_vect_int(val)
      ! reduction between all the masters of subworlds, then spreading info in each subworld (for an array of int)
      use mld_subworld
      implicit none
      integer, dimension(:), intent(inout)  :: val
      integer :: ierr
      if (MPI_COMM_NULL /= masters_world) call MPI_ALLREDUCE(MPI_IN_PLACE, val, size(val,1), MPI_INTEGER, MPI_SUM, masters_world, ierr)
      call MPI_BCAST(val, size(val,1), MPI_INTEGER, 0, subworld, ierr)
   end subroutine subworlds_allreduce_vect_int

   subroutine subworlds_allreduce_int(val)
      ! reduction between all the masters of subworlds, then spreading info in each subworld (for an int)
      use mld_subworld
      implicit none
      integer, intent(inout)  :: val
      integer :: ierr
      if (MPI_COMM_NULL /= masters_world) call MPI_ALLREDUCE(MPI_IN_PLACE, val, 1, MPI_INTEGER, MPI_SUM, masters_world, ierr)
      call MPI_BCAST(val, 1, MPI_INTEGER, 0, subworld, ierr)
   end subroutine subworlds_allreduce_int

   subroutine subworlds_allreduce_from_evrywhere_int(val_in, val_reduce)
      ! reduction between all the procs of subworlds into masters of subworlds, 
      ! then spreading info in each subworld 
      ! (for an int)
      use mld_subworld
      implicit none
      integer, intent(in)  :: val_in
      integer, intent(inout)  :: val_reduce
      integer :: ierr

      call MPI_Reduce(val_in, val_reduce, 1, MPI_INTEGER, MPI_SUM, 0, subworld, ierr)
      if (MPI_COMM_NULL /= masters_world) call MPI_AllReduce(MPI_IN_PLACE, val_reduce, 1, MPI_INTEGER, MPI_SUM, masters_world, ierr)
      call MPI_Bcast(val_reduce , 1, MPI_INTEGER, 0, subworld, ierr)
   end subroutine subworlds_allreduce_from_evrywhere_int

   subroutine subworlds_allreduce_from_evrywhere_double(val_in, val_reduce)
      ! reduction between all the procs of subworlds into masters of subworlds, 
      ! then spreading info in each subworld 
      ! (for an double)
      use mld_subworld
      use module_kind_variables, only: kind_double
      implicit none
      real(kind_double), intent(in)  :: val_in   
      real(kind_double), intent(inout)  :: val_reduce
      integer :: ierr

      call MPI_Reduce(val_in, val_reduce, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, subworld, ierr)
      if (MPI_COMM_NULL /= masters_world) call MPI_AllReduce(MPI_IN_PLACE, val_reduce, 1, MPI_DOUBLE_PRECISION, MPI_SUM, masters_world, ierr)
      call MPI_Bcast(val_reduce , 1, MPI_DOUBLE_PRECISION, 0, subworld, ierr)
   end subroutine subworlds_allreduce_from_evrywhere_double


   subroutine subworlds_allreduce_from_evrywhere_vect_double(val_in, val_reduce)
      ! reduction between all the procs of subworlds into masters of subworlds, 
      ! then spreading info in each subworld 
      ! (for an vector double)
      use mld_subworld
      use module_kind_variables, only: kind_double
      implicit none
      real(kind_double), dimension(:), intent(in)  :: val_in   
      real(kind_double), dimension(:), intent(inout)  :: val_reduce
      integer :: ierr

      call MPI_Reduce(val_in, val_reduce, size(val_in,1), MPI_DOUBLE_PRECISION, MPI_SUM, 0, subworld, ierr)
      if (MPI_COMM_NULL /= masters_world) call MPI_AllReduce(MPI_IN_PLACE, val_reduce, size(val_reduce,1), MPI_DOUBLE_PRECISION, MPI_SUM, masters_world, ierr)
      call MPI_Bcast(val_reduce , size(val_reduce,1), MPI_DOUBLE_PRECISION, 0, subworld, ierr)
   end subroutine subworlds_allreduce_from_evrywhere_vect_double


   subroutine subworlds_get_max_double(val_in, val_max)
      ! get the max value of a double from all the procs of subworlds into masters of subworlds, 
      ! then spreading info in each subworld 
      ! (for an double)
      use mld_subworld
      use module_kind_variables, only: kind_double
      implicit none
      real(kind_double), intent(in)  :: val_in   
      real(kind_double), intent(inout)  :: val_max
      integer :: ierr
      !call MPI_Allreduce(local_val, max_val, 1, MPI_DOUBLE_PRECISION, MPI_MAX, masters_world, ierr)
      if (MPI_COMM_NULL /= masters_world) call MPI_AllReduce(val_in, val_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, masters_world, ierr)
      call MPI_Bcast(val_max , 1, MPI_DOUBLE_PRECISION, 0, subworld, ierr)
   end subroutine subworlds_get_max_double 


   subroutine subworlds_allreduce_from_evrywhere_matrix_double(val_in, val_reduce)
      ! reduction between all the procs of subworlds into masters of subworlds, 
      ! then spreading info in each subworld 
      ! (for an vector double)
      use mld_logger
      use mld_subworld
      use module_kind_variables, only: kind_double
      implicit none
      real(kind_double), dimension(:,:), intent(in)  :: val_in   
      real(kind_double), dimension(:,:), intent(inout)  :: val_reduce
      integer :: ierr, size1, size2
      _NAMECURRENT_('subworlds_allreduce_from_evrywhere_matrix_double')
      _MLD_BEGIN_

      size1 = size(val_in,1)
      size2 = size(val_in,2)
      if (size1 /= size(val_reduce,1)) then 
        call log_critical('critical error: mismatch size  size1 /= size(val_reduce,1) in '//NAMECURRENT)
        stop 'mismatch size1'
      end if 

      if (size2 /= size(val_reduce,2)) then 
        call log_critical('critical error: mismatch size  size2 /= size(val_reduce,2) in '//NAMECURRENT)
        stop 'mismatch size2'
      end if 

      call MPI_Reduce(val_in, val_reduce, size1*size2, MPI_DOUBLE_PRECISION, MPI_SUM, 0, subworld, ierr)
      if (MPI_COMM_NULL /= masters_world) call MPI_AllReduce(MPI_IN_PLACE, val_reduce, size1*size2, MPI_DOUBLE_PRECISION, MPI_SUM, masters_world, ierr)
      call MPI_Bcast(val_reduce , size1*size2, MPI_DOUBLE_PRECISION, 0, subworld, ierr)

      _MLD_END_ 
   end subroutine subworlds_allreduce_from_evrywhere_matrix_double

end module my_mpi_subroutines







