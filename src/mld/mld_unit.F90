
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

module mld_unit

  ! This is a simple module to dispatch files unit in MILADY code.
  ! If no units are available, critical message (stop code) is returned.

  use mld_string
  use mld_logger

  implicit none


  ! unit of files range for milady from 110 to 150, as checked no more simultaneous opened files
  integer, parameter, private      :: &
    unit_min = 110, &
    unit_max = 150, &
    nbmax = unit_max - unit_min

  type     :: unit_type
    integer  :: unit = -1
    character(len=100)   :: name = ''
  end type

  ! singleton common program milady used opened file units
  type(unit_type), dimension(:), private :: com_mld_unit(0:nbmax) = unit_type(-1, '')                ! first init
  logical, private     :: init_done = .false.      ! overprecaution

contains


  subroutine mld_unit_init()
    ! seems overprecautions if first init exist
    com_mld_unit(:)%unit = -1
    com_mld_unit(:)%name = ''
    init_done = .true.
  end subroutine


  function get_unit(name) result(res)
    integer  :: res, i
    character(*), intent(in)   :: name

    if (.not. init_done) call mld_unit_init()
    res = -1
    do i = 0, nbmax
      if (com_mld_unit(i)%unit == -1) then
        res = unit_min + i
        com_mld_unit(i)%unit = res
        com_mld_unit(i)%name = trim(name)
        return
      end if
    end do
    !call mld_critical_abort('maximum simultaneous opened files reached, fix it, opened files are'//nwl//str_opened_files())
    stop "maximum simultaneous opened files reached, fix it" ! opened files are'//nwl//str_opened_files()
    
  end function


  function open_file(name, form, action) result(res)
    ! new or old file as unknown status
    ! res is file unit that mld_unit_mod decides
    integer  :: res, iostat
    character(*), intent(in)   :: name
    character(*), intent(in), optional     :: form, action
    character(len=:), allocatable    :: formval, actionval

    formval = 'formatted'
    if (present(form)) formval = trim(form)
    actionval = 'readwrite'
    if (present(action)) actionval = trim(action)

    res = get_unit(name)
    open (res, file=name, status='unknown', form=formval, action=actionval, iostat=iostat)
    if (iostat /= 0) then
      call log_critical('problem open file "'//name//'" unit '//vtoa(res)//' io_status '//vtoa(iostat))
      res = -1
    end if
    call log_debug('Open file "'//name//'" unit '//vtoa(res))
  end function


  function open_new_file(name, form, action) result(res)
    ! res is file unit that mld_unit_mod decides
    integer  :: res, iostat
    character(*), intent(in)   :: name
    character(*), intent(in), optional     :: form, action
    character(len=:), allocatable    :: formval, actionval

    formval = 'formatted'
    if (present(form)) formval = trim(form)
    actionval = 'readwrite'
    if (present(action)) actionval = trim(action)

    res = get_unit(name)
    open (res, file=name, status='new', form=formval, action=actionval, iostat=iostat)
    if (iostat /= 0) then
      !call mld_critical_abort('problem open new file "'//name//'" unit '//vtoa(res)//' io_status '//vtoa(iostat))
      stop 'problem open new file "'//name//'" unit '//vtoa(res)//' io_status '//vtoa(iostat)
    end if
    call log_debug('open new file "'//name//'" unit '//vtoa(res))
  end function


  function open_old_file(name, form, action) result(res)
    ! res is file unit that mld_unit_mod decides
    integer  :: res, iostat
    character(*), intent(in)   :: name
    character(*), intent(in), optional     :: form, action
    character(len=:), allocatable    :: formval, actionval

    formval = 'formatted'
    if (present(form)) formval = trim(form)
    actionval = 'readwrite'
    if (present(action)) actionval = trim(action)

    res = get_unit(name)
    open (res, file=name, status='old', form=formval, action=actionval, iostat=iostat)
    if (iostat /= 0) then
      call log_critical('problem open old file "'//name//'" unit '//vtoa(res)//' io_status '//vtoa(iostat))
      ! call mld_critical_abort('problem open old file "'//name//'" unit '//vtoa(res)//' io_status '//vtoa(iostat))
      res = -1
      ! never mind that it will comes bad, critical message done...
    else
      call log_debug('open old file "'//name//'" unit '//vtoa(res))
    end if
  end function


  function close_unit(unit) result(res)
    logical  :: res
    integer, intent(in)  :: unit
    integer  :: i

    if (.not. init_done) call mld_unit_init()
    res = .false.
    if ((unit < unit_min) .or. (unit >= unit_max)) then
      call log_warning('close unit '//vtoa(unit)// &
                       'not in allowed range '//str_range(unit_min, unit_max))
    end if

    i = unit - unit_min
    if (com_mld_unit(i)%unit /= -1) then
      call log_debug('close file "'//trim(com_mld_unit(i)%name)//'" unit '//vtoa(unit))
      close (unit)
      com_mld_unit(i)%unit = -1
      com_mld_unit(i)%name = ''
      res = .true.
    else
      call log_warning('close unit '//vtoa(unit)// &
                       'not in known open file'//nwl//str_opened_files())
    end if
  end function


  function str_opened_files() result(res)
    ! returns string resume of currently opened files
    character(len=:), allocatable    :: res
    integer  :: i

    if (.not. init_done) call mld_unit_init()
    res = ''
    do i = 0, nbmax
      if (com_mld_unit(i)%name /= '') then
        res = res//vtoa(com_mld_unit(i)%unit)//trim(com_mld_unit(i)%name)//nwl
      end if
    end do
    if (res == '') then
      res = 'no opened file known'
    end if
  end function


  function file_read_stream(afile, lg_max) result(res)
    ! classical no-mpi read io stream all file (as C do not stop on new lines)
    ! mandatory expected all ascii bytes
    ! lg_max bytes in files expected

    character(*)   :: afile
    integer  :: lg_max
    logical  :: ok
    character(lg_max)    :: buffer
    character(:), allocatable  :: res
    integer  :: iunit, iostat

    ok = .false.
    res = ""
    buffer = repeat(' ', len(buffer))                ! mandatory precaution
    iunit = get_unit(afile)
    open (iunit, file=afile, status="old", access="stream")
    read (iunit, iostat=iostat) buffer
    ok = close_unit(iunit)

    ! iostat < 0 then end of the input has reached
    call log_info('read_stream iostat '//vtoa(iostat)//'(-1 as EOF is normal as overdimensioned buffer is expected)')
    if (iostat >= 0) then
      call log_warning('ead_stream buffer_nml too small or other problem, unexpected '//vtoa([iostat, len(buffer)]))
      ! call mld_mpi_abort('read_stream error')
    end if

    res = buffer(1:get_useful_len(buffer))           ! overdimensioned
    call log_debug('buffer from '//afile//nwl//'->'//res//'<-')
  end function


  function read_first_line(name) result(res)
    ! return line trimmed without line feed, max 1024 characters
    character(len=:), allocatable    :: res
    character(*), intent(in)   :: name
    character(1024)      :: tmp
    integer  :: iunit
    logical  :: ok

    tmp = repeat("x", len(tmp))
    iunit = open_old_file(name, 'formatted')
    read (iunit, '(a)') tmp
    ok = close_unit(iunit)
    res = trim(tmp)
    call log_debug('read_first_line of file '//name//nwl//'"'//res//'"')
  end function


  subroutine test_mld_unit()
    integer  :: i1, i2, i3

    i1 = get_unit('first file')
    i2 = get_unit('second file')
    i3 = get_unit('third file')
    call log_info('Unit files are'//nwl//str_opened_files())
  end subroutine

end module
