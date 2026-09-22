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

module mld_environ

  use mld_string
  use mld_logger

  implicit none


contains


  function get_env(name) result(res)
    ! copy in code memory value
    character(len=:), allocatable    :: res
    character(*)   :: name
    character(len=255)   :: tmp

    call get_environment_variable(name, tmp)
    res = trim(tmp)
  end function

  function get_env_s(name, default) result(res)
    ! copy in code memory value
    character(len=:), allocatable    :: res
    character(*)   :: name, default
    character(len=:), allocatable    :: strvalue

    strvalue = get_env(name)
    if (strvalue == "") then
      res = trim(default)
      call log_info('Get value from environ '//name//' as <red>default '//trim(default)//'<default>')
    else
      res = trim(strvalue)
    end if
  end function

  function get_env_i(name, default, range) result(res)
    ! get integer value
    integer  :: res, default
    integer, optional    :: range(2)

    character(*)   :: name
    character(len=:), allocatable    :: strvalue
    logical  :: ko

    strvalue = get_env(name)
    if (strvalue == "") then
      res = default
      call log_info('Get value from environ '//name//' as <red>default '//itoa(default)//'<default>')
    else
      res = atoi(strvalue)
      call log_info('Get value from environ '//name//' as '//strvalue)
      ko = .false.
      if (present(range)) then
        if (res < range(1)) then
          res = range(1)
          ko = .true.
        elseif (res > range(2)) then
          res = range(2)
          ko = .true.
        end if
        if (ko) call log_warning('Get value from environ '//name//' outside range '// &
                                 str_range(range(1), range(2))//', set as '//vtoa(res))
      end if
    end if
  end function

  function get_hostname() result(res)
    character(len=:), allocatable    :: res
    res = get_env('HOSTNAME')
  end function

  function get_pwd() result(res)
    character(len=:), allocatable    :: res
    res = get_env('PWD')
  end function

  function get_user() result(res)
    character(len=:), allocatable    :: res
    res = get_env('USER')
    if (res == "") then
      res = get_env('USERNAME')                        ! sometimes windows
    end if
    if (res == "") res = 'UNKNOWN'
  end function

end module
