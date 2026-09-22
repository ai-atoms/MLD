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

module mld_colors_ansi
  ! console ansi colors for ascii character strings
  use mld_string
  implicit none
  ! private

  character(1), parameter    :: &
    c_esc = achar(27), &    ! Escape
    c_end = 'm'             ! End ansi code, "m"

  character(2), parameter    :: &
    c_sta = c_esc//'['      ! Start ansi code, "\["

  character(4), parameter    :: &
    c_clear = c_sta//'0'//c_end, &                   ! Clear all styles, "\[0m"
    c_BOLD = c_sta//'1'//c_end, &
    c_DIM = c_sta//'2'//c_end

  character(5), parameter    :: &
    c_BLACK = c_sta//'30'//c_end, &
    c_RED = c_sta//'31'//c_end, &
    c_GREEN = c_sta//'32'//c_end, &
    c_YELLOW = c_sta//'33'//c_end, &
    c_BLUE = c_sta//'34'//c_end, &
    c_MAGENTA = c_sta//'35'//c_end, &
    c_CYAN = c_sta//'36'//c_end, &
    c_WHITE = c_sta//'37'//c_end, &
    c_DEFAULT = c_sta//'39'//c_end


contains

  function ansi_set(astr) result(res)
    ! replace '<red>' as c_RED in astr, and etc. for other tags
    character(*)   :: astr
    character(len=:), allocatable    :: res

    res = aStr
    res = str_repl(res, '<clear>', c_clear)
    res = str_repl(res, '<bold>', c_bold)
    res = str_repl(res, '<dim>', c_dim)
    res = str_repl(res, '<black>', c_black)
    res = str_repl(res, '<red>', c_red)
    res = str_repl(res, '<green>', c_green)
    res = str_repl(res, '<yellow>', c_yellow)
    res = str_repl(res, '<blue>', c_blue)
    res = str_repl(res, '<magenta>', c_magenta)
    res = str_repl(res, '<cyan>', c_cyan)
    res = str_repl(res, '<white>', c_white)
    res = str_repl(res, '<default>', c_default)
  end function ansi_set

  function ansi_reset(astr) result(res)
    ! replace '<red>' as nothing in astr, and etc. for other tags
    character(*)   :: astr
    character(len=:), allocatable    :: res

    res = astr
    if (len(astr) == 0) return
    ! if compile -fcheck=all avoid
    ! Fortran runtime error: Allocatable actual argument 'res' is not allocated
    res = aStr//'x'         ! do not like empty res
    res = str_repl(res, '<clear>', '')
    res = str_repl(res, '<bold>', '')
    res = str_repl(res, '<dim>', '')
    res = str_repl(res, '<black>', '')
    res = str_repl(res, '<red>', '')
    res = str_repl(res, '<green>', '')
    res = str_repl(res, '<yellow>', '')
    res = str_repl(res, '<blue>', '')
    res = str_repl(res, '<magenta>', '')
    res = str_repl(res, '<cyan>', '')
    res = str_repl(res, '<white>', '')
    res = str_repl(res, '<default>', '')
    ! write (0,*) 'ansi_reset in  "'//astr//'"'//nl//'ansi_reset out "'//res(:)//'"', len(res)
    res = res(1:len(res) - 1)
  end function ansi_reset


end module


