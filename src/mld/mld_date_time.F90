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

module mld_date_time

  use module_kind_variables, only: kind_double
  use cc_utils_f2c_mod

  implicit none

  ! kind=8 -> count_max 9223372036854775807 as 106751 days
  ! kind=4 -> count_max 2147483647 as 24 days
  integer, parameter, private      :: kind_i = 8

  character(*), parameter, private :: &
    fma = '(a)', &
    nl = achar(10), &       ! NEW_LINE('A')
    tab = achar(9)


  type time_counter
    ! SYSTEM CLOCK TIME
    integer(kind_i)      :: &
      cpt_init = 0, &
      cpt_fin = 0, &
      cpt_max = 0, &
      freq = 0                ! nb clock periods per sec
    ! CPU TIME
    real(kind_double)   :: &
      t_init = 0., &
      t_fin = 0.
  end type time_counter

contains


  subroutine init_timer(counter)
    type(time_counter)   :: counter
    ! initialisations
    call system_clock(count_rate=counter%freq, count_max=counter%cpt_max, count=counter%cpt_init)
    call cpu_time(time=counter%t_init)
  end subroutine


  subroutine write_timer(counter, iunit)
    type(time_counter)   :: counter
    integer  :: iunit
    write (iunit, fma) get_timer(counter)
  end subroutine


  function get_timer(counter) result(res)
    character(len=:), allocatable    :: res
    character(*), parameter    :: fmt = '(a,f10.2,a,f10.2,a)'
    character(100) :: tmp
    type(time_counter)   :: counter
    real(kind_double)   :: time_elapsed, t_cpu

    time_elapsed = get_time_elapsed(counter)
    if (counter%freq == 0) then                      ! init_timer not called
      t_cpu = -1
    else
      t_cpu = counter%t_fin - counter%t_init
    end if
    ! write (iunit, *) 'COUNT_MAX : ', real(counter%cpt_max) / counter%freq / 3600 / 24, ' days'
    write (tmp, fmt) 'Time elapsed   : ', time_elapsed, 's'//nl//'Time cpu user  : ', t_cpu, 's'
    res = trim(tmp)
  end function

  function get_usage(counter) result(res)
    character(len=:), allocatable    :: res
    character(*), parameter    :: fmt = '(a,f10.2,a,i10,a)'
    character(100) :: tmp
    type(time_counter)   :: counter

    integer(kind=c_int)  :: memory_usage
    real(kind=c_double)  :: time_user, time_sys


    res = get_timer(counter)                         ! Time elapsed and Time cpu user
    call c_getrusage(time_user, time_sys, memory_usage)                       ! Time system and max memory
    write (tmp, fmt) 'Time system    : ', time_sys, 's'//nl//'Max memory     : ', memory_usage, 'k'
    res = res//nl//trim(tmp)

  end function



  function get_time_elapsed(counter) result(time_elapsed)
    real(kind_double)   :: time_elapsed
    type(time_counter)   :: counter
    integer(kind_i)      :: cpt

    call cpu_time(time=counter%t_fin)
    call system_clock(count=counter%cpt_fin)
    cpt = counter%cpt_fin - counter%cpt_init
    if (counter%cpt_fin < counter%cpt_init) then     ! if huge reached
      cpt = cpt + counter%cpt_max
    end if
    if (counter%freq == 0) then                      ! init_timer not
      time_elapsed = -1       ! as not initiated
    else
      time_elapsed = real(cpt)/counter%freq
    end if
  end function


  function get_date() result(res)
    character, parameter :: sl = '/'
    character(len=:), allocatable    :: res
    integer, dimension(8)      :: values             ! 2020    3   24   60   15   33   28  191
    character(30)  :: tmp

    call date_and_time(values=values)
    write (tmp, '(i2.2,a,i2.2,a,i4.4,a,i2.2,a,i2.2,a,i2.2,a)') &
      values(3), sl, values(2), sl, values(1), ' ', &
      values(5), 'h', values(6), 'm', values(7), 's'
    res = trim(tmp)
  end function

end module mld_date_time
