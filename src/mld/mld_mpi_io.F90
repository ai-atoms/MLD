
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

module mld_mpi_io

  use mpi
  use mld_mpi
  use mld_string
  use mld_logger

  implicit none

  integer, private     :: &
    ioerr, &
    ioerrbis, &
    error_len, &
    fh

  integer(kind=MPI_OFFSET_KIND), private :: offset
  integer, dimension(MPI_STATUS_SIZE), private :: iostat
  character(len=mpi_max_error_string)    :: error_text


contains


  subroutine mpi_io_abort_test(mess, ioerror)
    character(*)   :: mess  ! a message from caller mandatory
    integer  :: ioerror

    if (ioerror /= mpi_success) then
      call mpi_error_string(ioerr, error_text, error_len, ioerrbis)
      call log_error(mess//" error"//nwl//'"'//error_text(1:error_len)//'"')
      call mld_mpi_abort(mess//" error")
    end if
  end subroutine


  function mpi_io_read_asci(afile) result(res)
    ! use mpi_file_read_at_all
    ! to get ALL file (as C do not stop on new lines)
    ! mandatory expected all file with only ascii bytes

    _NAMECURRENT_("mpi_io_read_asci")


    character(*)   :: afile
    character(:), allocatable  :: res
    ! local
    logical  :: ok
    integer  :: i, count
    integer(KIND=MPI_OFFSET_KIND)    :: siz
    character, dimension(:), allocatable   :: buffer_byte

    _MLD_BEGIN_
    ok = .false.
    offset = 0
    res = ""
    call mpi_file_open(mpi_comm_mld, afile, mpi_mode_rdonly, mpi_info_null, fh, ioerr)
    call mpi_io_abort_test(NAMECURRENT, ioerr)
    call mpi_file_get_size(fh, siz, ioerr)
    call mpi_io_abort_test(NAMECURRENT, ioerr)
    count = int(siz)
    allocate (buffer_byte(count))
    call mpi_file_read_at_all(fh, offset, buffer_byte, count, mpi_byte, iostat, ioerr)
    call mpi_io_abort_test(NAMECURRENT, ioerr)
    call log_debug('mpi_io_read_asci iostat ioerr siz '//vtoa([iostat, ioerr, count]))
    call mpi_file_close(fh, ioerr)
    call mpi_io_abort_test(NAMECURRENT, ioerr)

    res = repeat('x', siz)  ! dimensioning res (as allocate)
    ! standard mode casting, suppose all asci bytes, could protect TODO
    do i = 1, int(siz)
      res(i:i) = buffer_byte(i)                        ! res is preallocated
    end do

    ok = .true.             ! if we are here, it is ok, supposedly
    call log_debug('res from '//afile//nwl//'->'//res//'<-')                  ! TODO log only if short file.
    deallocate (buffer_byte)
    _MLD_END_
  end function

  !-----------------------------------------------
  ! namelist utilities (are with mpi_io)
  !-----------------------------------------------
  function remove_comments_namelist(mess) result(res)
    ! remove fortran comments in as-namelist string from mpi_io_read_asci
    ! because  'dmtype=18 !17 ABF, 7 PHONDY' cause problem (as 'end-of-file during read')
    ! when read namelist in memory (but not with read namelist file old way)

    character(*)   :: mess  ! a content string (as a namelist) with '! etc' to remove
    character(:), allocatable  :: res
    integer  :: i, j, imax
    character(1)   :: next_char
    logical  :: ok

    res = mess
    imax = len(mess)
    res = repeat(' ', imax)
    ok = .true.
    j = 1
    do i = 1, imax
      next_char = mess(i:i)
      if (next_char == '!') ok = .false.
      if (next_char == nwl) ok = .true.
      if (ok) then
        res(j:j) = next_char
        j = j + 1
      end if
    end do
    j = get_useful_len(res)
    res = res(1:j)
    ! call log_debug('remove_comments_namelist'//nwl//'->'//res//'<-')
  end function


  function tabulate_namelist(mess) result(res)
    ! format namelist as remove whitespaces/tabulations before name = value
    ! to get pretty print

    character(*)   :: mess  ! a content string (as a namelist) with '! etc' to remove
    character(:), allocatable  :: res
    integer  :: i, j, imax
    character(1)   :: next_char
    logical  :: ok, begin_line

    res = mess
    imax = len(mess)
    res = repeat(' ', imax)
    ok = .false.
    j = 1
    begin_line = .true.
    do i = 1, imax
      next_char = mess(i:i)
      if (next_char == tbl) next_char = ' '
      if (next_char == '=') begin_line = .false.
      if (next_char /= ' ') begin_line = .false.
      if ((next_char == ' ') .and. (begin_line)) then
        ok = .false.
      else
        ok = .true.
        begin_line = .false.
      end if
      if (next_char == nwl) then
        ok = .true.
        begin_line = .true.
      end if
      if (ok) then
        res(j:j) = next_char
        j = j + 1
      end if
    end do
    res = res(1:j)
    ! call log_debug('tabulate_namelist'//nwl//'->'//res//'<-')
  end function


  function replace_comma_namelist(mess) result(res)
    ! insert new line in place of ','
    ! to get pretty print

    character(*)   :: mess  ! a content string (as a namelist) with '! etc' to remove
    character(:), allocatable  :: res
    integer  :: i, j, imax
    character(1)   :: next_char
    logical  :: ok, begin_line

    res = mess
    imax = len(mess)
    res = repeat(' ', imax)
    ok = .false.
    j = 1
    begin_line = .true.
    do i = 1, imax
      next_char = mess(i:i)
      if (next_char == ',') next_char = nwl
      res(j:j) = next_char
      j = j + 1
    end do
    res = res(1:j)
    ! call log_debug('replace_comma_namelist'//nwl//'->'//res//'<-')
  end function


  function reformat_namelist(mess) result(res)
    character(*)   :: mess  ! a content string (as a namelist) with '! etc' to remove
    character(:), allocatable  :: res

    res = remove_comments_namelist(mess)
    res = tabulate_namelist(res)
    res = replace_comma_namelist(res)
    call log_debug('reformat_namelist'//nwl//'->'//res//'<-')
  end function


end module
