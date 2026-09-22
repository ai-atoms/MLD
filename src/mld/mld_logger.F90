
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

module mld_logger

  use iso_fortran_env, only: error_unit, output_unit, input_unit
  use mld_string
  use mld_colors_ansi

#if(IFORT)
  use ifport, only: isatty                         ! isatty gnu is intrinsic https://gcc.gnu.org/onlinedocs/gfortran/ISATTY.html
#endif

  implicit none


  integer, parameter   :: &
    log_lev_debug = 10, &
    log_lev_info = 20, &
    log_lev_test = 25, &
    log_lev_warning = 30, &
    log_lev_error = 40, &
    log_lev_critical = 50

  character(8), parameter    :: &
    tag_unknown = "UNKNOWN", &
    tag_debug = "DEBUG", &
    tag_test = "TEST", &
    tag_info = "INFO", &
    tag_warning = "WARNING", &
    tag_error = "ERROR", &
    tag_critical = "CRITICAL"

  ! STDIN_ FILENO 0 Standard C input  value stdin,  5 Standard fortran
  ! STDOUT_FILENO 1 Standard C output value stdout, 6 Standard fortran
  ! STDERR_FILENO 2 Standard C error  value stderr, 0 Standard fortran
  integer, parameter, private      :: &
    cte_stdout = output_unit, &                      ! from iso_fortran_env
    cte_stderr = error_unit, &
    cte_stdin = input_unit

  integer, private     :: &
    ! SetDebug 
    !current_level = log_lev_debug                     ! debug     ! TODO log_lev_info
    current_level = log_lev_info                     ! debug     ! TODO log_lev_info

  logical, private     :: &
    verbose = .true., &     ! set .false. implies avoid execute _MLD_BEGIN_ & _MLD_END_, see ../MLD_MACRO.INC
    init_log_done = .false.

  integer, public      :: &
    io_log = cte_stdout, &  ! io_log is stdout for all log_xxxx()
    io_fic = -1, &          ! io_fic is iunit file milady.log
    logger_mpi_rank = 0     ! have to be set by mld_mpi_init

  logical, public      :: &
    ! SetColor 
    is_colorized_io_log = .true., &                  ! set to 1 implies color as tty as default but sometimes not
    is_colorized_io_fic = .false., &                 ! set to 0 implies no color as file, ever as file, never as colorized
    is_log_only_rank_0 = .true.                      ! set implies every log else rank 0 are invalid



contains

  pure function mld_verbose() result(res)
    logical  :: res
    res = verbose
  end function

  function get_tag_level(level) result(res)
    integer  :: level
    character(len=:), allocatable    :: res

    select case (level)
    case (log_lev_debug)
      res = tag_debug
    case (log_lev_test)
      res = tag_test
    case (log_lev_info)
      res = tag_info
    case (log_lev_warning)
      res = tag_warning
    case (log_lev_error)
      res = tag_error
    case (log_lev_critical)
      res = tag_critical
    case default
      res = tag_unknown
    end select
  end function get_tag_level


  function get_tag_level_color(level) result(res)
    integer  :: level
    character(len=:), allocatable    :: res

    select case (level)
    case (log_lev_debug)
      res = '<cyan>'
    case (log_lev_test)
      res = ''
    case (log_lev_info)
      res = '<green>'
    case (log_lev_warning)  ! theorically continue
      res = '<red>'
    case (log_lev_error)    ! theorically try to continue
      res = '<red><bold>'
    case (log_lev_critical) ! theorically better to stop
      res = '<red><bold>'
    case default
      res = ''
    end select
  end function get_tag_level_color


  subroutine log_example()
    call log_debug('this is log_debug')
    call log_test('this is log_test')
    call log_info('this is log_info')
    call log_warning('this is log_warning')
    call log_error('this is log_error')
    call log_critical('this is log_critical')
  end subroutine


  function log_to_ansi(mess, is_colorized) result(res)
    ! convert '<red>' to ansi code etc.
    ! only if is_colorized, else remove '<red>' etc.
    character(*)   :: mess
    character(len=:), allocatable    :: res
    logical  :: is_colorized

    if (is_colorized) then
      res = ansi_set(mess)
    else
      res = ansi_reset(mess)
    end if
  end function log_to_ansi

  subroutine log_set_all_ranks()
    ! to get log from all mpi processus
    ! usage for debug on small number of mpi processus
    is_log_only_rank_0 = .false.
  end subroutine

  subroutine log_set_only_rank_0()
    ! to get log from only mpi processus of rank 0
    is_log_only_rank_0 = .true.
  end subroutine

  subroutine log_flush()
    ! to flush log(s) in case of mpi abort for example
    ! flush is standard statement in Fortran 2003
    if (io_log /= -1) flush (io_log)
    if (io_fic /= io_log) flush (io_fic)
  end subroutine

  subroutine log_set_nocolor()
    is_colorized_io_log = .false.
    call log_debug('set no colorized')
  end subroutine log_set_nocolor


  subroutine log_set_color()
    is_colorized_io_log = .true.
    call log_debug('set colorized')
  end subroutine log_set_color


  subroutine log_set_level(level)
    integer  :: level
    if ((level < 0) .or. (level > log_lev_critical)) then
      call log_error('log_set_level: problem level value '//str_convi(level))
      return
    end if
    current_level = level
  end subroutine log_set_level


  subroutine log_init(mld_log)
    character(*), optional     :: mld_log
    character(len=:), allocatable    :: mld_log_val
    integer  :: iostat

    if (init_log_done) then
      call log_error('log_init done yet, skipped.')
      return
    end if

    mld_log_val = "milady.log"
    if (present(mld_log)) mld_log_val = trim(mld_log)

    ! io_fic = mld_open_new_file(mld_log_val) ! remove circular dependency
    io_fic = 99             ! as unusued in mld_unit
    open (io_fic, file=mld_log_val, status='new', form='formatted', iostat=iostat)
    if (iostat /= 0) then
      !call mld_critical_abort('Problem open new file "'//mld_log_val//'" unit '//vtoa(io_fic)//' io_status '//vtoa(iostat))
      stop 'Problem open new file "'//mld_log_val//'" unit '//vtoa(io_fic)//' io_status '//vtoa(iostat)
    end if
    call log_debug('Open new file for logger "'//mld_log_val//'" unit '//vtoa(io_fic))

    ! as redirection tee pipe etc ...
    ! TODOcvw invalid gfortran 10 is_colorized_io_log = isatty(io_log)
    is_colorized_io_log = .true.

    ! unconditionnaly
    is_colorized_io_fic = .false.

    init_log_done = .true.
  end subroutine


  subroutine log_write(level, mess)
    integer  :: level
    character(*)   :: mess
    integer  :: ilog
    character(len=:), allocatable    :: &
      mess_color, level_color, tag_color, tag_clear, mess_indented

    if (current_level > level) return                ! no message
    if ((logger_mpi_rank .ne. 0) .and. (is_log_only_rank_0)) return           ! no message

    ! output log as stdout if not redirected
    if (io_log /= io_fic) then                       ! only if stdout
      ilog = io_log
      mess_color = log_to_ansi(mess, is_colorized_io_log)
      level_color = log_to_ansi(get_tag_level(level), is_colorized_io_log)
      tag_color = log_to_ansi(get_tag_level_color(level), is_colorized_io_log)
      tag_clear = log_to_ansi('<clear>', is_colorized_io_log)
#if (openMP)
      mess_indented = str_repl(mess_color, nwl, nwl//"             ")
      ! Suspect write is not thread-safe ceinture et bretelles
!$OMP CRITICAL(the_tty)
      write (ilog, '(a,a1,a8,i2,a2,a,a)') tag_color, "[", level_color, omp_get_thread_num(), "] ", mess_indented, tag_clear
!$OMP END CRITICAL(the_tty)
#else
      mess_indented = str_repl(mess_color, nwl, nwl//"           ")
      write (ilog, '(a,a1,a8,a2,a,a)') tag_color, "[", level_color, "] ", mess_indented, tag_clear
#endif
    end if

    ! output fic as file unconditionaly, get log message also, AND old write (io_fic, *)
    if (io_fic /= -1) then  ! only if io_fic opened
      ilog = io_fic
      mess_color = log_to_ansi(mess, is_colorized_io_fic)
      if (len(mess_color) == 0) mess_color = ''
      level_color = log_to_ansi(get_tag_level(level), is_colorized_io_fic)
      tag_color = log_to_ansi(get_tag_level_color(level), is_colorized_io_fic)
      tag_clear = log_to_ansi('<clear>', is_colorized_io_fic)
#if (openMP)
      mess_indented = str_repl(mess_color, nwl, nwl//"             ")
      ! Suspect write is not thread-safe ceinture et bretelles
!$OMP CRITICAL(the_file)
      write (ilog, '(a,a1,a8,i2,a2,a,a)') &
        tag_color, "[", level_color, omp_get_thread_num(), "] ", mess_indented, tag_clear
!$OMP END CRITICAL(the_ file)
#else
      mess_indented = str_repl(mess_color, nwl, nwl//"           ")
      write (ilog, '(a,a1,a8,a2,a,a)') &
        tag_color, "[", level_color, "] ", tag_clear, mess_indented, tag_clear
#endif
    end if

  end subroutine log_write


  subroutine log_debug(mess)
    character(*)   :: mess
    call log_write(log_lev_debug, mess)
  end subroutine log_debug


  subroutine log_test(mess)
    ! used for asserts in unittest
    ! avoid problem comparing expected/obtained test output strings
    character(*)   :: mess
    call log_write(log_lev_test, mess)
  end subroutine log_test


  subroutine log_info(mess)
    character(*)   :: mess
    call log_write(log_lev_info, mess)
  end subroutine log_info


  subroutine log_warning(mess)
    character(*)   :: mess
    call log_write(log_lev_warning, mess)
  end subroutine log_warning


  subroutine log_error(mess)
    character(*)   :: mess
    call log_write(log_lev_error, mess)
  end subroutine log_error


  subroutine log_critical(mess)
    character(*)   :: mess
    call log_write(log_lev_critical, mess)
  end subroutine log_critical

end module
