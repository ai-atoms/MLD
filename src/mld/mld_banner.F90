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

module mld_banner

  use iso_fortran_env
  use mld_string
  use mld_environ

  implicit none

  character(*), parameter, private :: &
    fma = '(a)', &
    fmesl = '(es23.15)', &  ! float long
    fmess = '(es12.5)', &   ! float short
    nl = achar(10), &       ! NEW_LINE('A')
    tab = achar(9), &
    NULL = achar(0)         ! used as non visible character, non empty character(len=:), allocatable , and NULL as end of string C

  character(*), parameter    :: banner = &
                                "            _ _           _       "//nl// &
                                "           (_) |         | |      "//nl// &
                                "   ___ ___  _| | __ _  __| |_   _ "//nl// &
                                " | '_ ` _ \| | |/ _` |/ _` | | | |"//nl// &
                                " | | | | | | | | (_| | (_| | |_| |"//nl// &
                                " |_| |_| |_|_|_|\__,_|\__,_|\__, |"//nl// &
                                "                             __/ |"//nl// &
                                "                            |___/ "//nl

  character(*), parameter    :: banner_cosmin = &
                                " ,---.    ,---.-./`)   .---.       ____    ______        ____     __  "//nl// &
                                " |    \  /    \ .-.')  | ,_|     .'  __ `.|    _ `''.    \   \   /  / "//nl// &
                                " |  ,  \/  ,  / `-' \,-./  )    /   '  \  \ _ | ) _  \    \  _. /  '  "//nl// &
                                " |  |\_   /|  |`-'`'`\  '_ '`)  |___|  /  |( ''_'  ) |     _( )_ .'   "//nl// &
                                " |  _( )_/ |  |.---.  > (_)  )     _.-`   | . (_) `. | ___(_ o _)'    "//nl// &
                                " | (_ o _) |  ||   | (  .  .-'  .'   _    |(_    ._) '|   |(_,_)'     "//nl// &
                                " |  (_,_)  |  ||   |  `-'`-'|___|  _( )_  |  (_.\.' / |   `-'  /      "//nl// &
                                " |  |      |  ||   |   |        \ (_ o _) /       .'   \      /       "//nl// &
                                " '--'      '--''---'   `--------`'.(_,_).''-----'`      `-..-'        "//nl

  character(*), parameter    :: copyright = &
                                "This code is under copyright & licence."//nl// &
                                "The distribution of the package or parts of package is not allowed."//nl// &
                                "Contributions A. M. Goryaeva, W. Unn-Toc, J. Baima, C. Lapointe, C. Van Wambeke"//nl// &
                                "For more details, please contact: mihai-cosmin.marinica@cea.fr"

contains

  function get_banner(c_opt) result(res)
    character(len=:), allocatable    :: res          ! provide res as long as you want
    character(*)   :: c_opt
    character(len=:), allocatable    :: c_ver, s_dat, s_time, c_opt_nl


    c_ver = compiler_version()
    c_ver = str_repl(c_ver, "Version", nl//"  Version")
    c_opt_nl = "  "//str_repl(c_opt, " -", " "//nl//"  -")
    ! c_opt = compiler_options()
    ! TODO touch banner.F90 after make to get real compile time
    s_dat = __DATE__
    s_time = __TIME__

    res = nl//banner//nl//copyright//nl// &
          ! nl//banner_cosmin//nl// &
          nl//"Compiled by"//nl//c_ver//nl// &
          nl//"Compilation at "//__DATE__//" "//__TIME__//nl// &
          nl//"Compiler options was"//nl//c_opt_nl//nl// &
          nl//"Execution user is "//get_env_s("USER", "UNKNOWN")// &
          nl//"Execution hostname is "//get_env_s("HOSTNAME", "UNKNOWN")//nl
  end function

end module
