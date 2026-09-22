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

module mld_string

  implicit none

  character(*), parameter, private :: &
    fma = '(a)', &
    nl = achar(10), &       ! NEW_LINE('A') !is to short to be public, use
    tab = achar(9), &
    NULL = achar(0)         ! used as non visible character, non empty character(len=:), allocatable , and NULL as end of string C

  integer, parameter   :: dbl = 8

  ! user format for readability
  integer, parameter   :: indent_fmt = 15          ! as next 'a15'
  character(*), parameter    :: &
    nwl = nl, &             ! NEW_LINE('A')
    tbl = tab, &
    fmt_indent = '                ', &               ! as next 'a15' + 1x
    fmt_l = '(a15, l10)', &
    fmt_t = '(a15, 1x, a)', &
    fmt_i = '(a15, i10)', &
    fmt_f = '(a15, es%es%)', &                       ! %es% is something as es12.4
    fmt_f2d = '(a15, a)', &
    fmt_f_vtoa = '10es%es%', &
    ! useless fmt_i_vtoa = 'i%es%', &
    fmt_f1d = '(a15, a)'    ! , &



  ! conversions prefix international system
  real(dbl), parameter :: &
    tofemto = 1d15, &
    topico = 1d12, &
    tonano = 1d9, &
    tomicro = 1d6, &
    tomilli = 1d3, &
    tokilo = 1d-3, &
    tomega = 1d-6, &
    togiga = 1d-9, &
    totera = 1d-12, &
    topeta = 1d-15

  integer  :: cmdc_precision = 0
  integer  :: tnr_precision = 0

  ! as value numeric to asci with one trailing whitespace
  interface vtoa
    module procedure &
      ltoa, &                 ! logical
      itoa, &                 ! integer
      iltoa, &                ! long integer
      array_ltoa, &           ! logical
      array_itoa, &
      array_iltoa, &
      array_2d_itoa, &
      array_2d_ftoa, &
      ftoa, &
      array_ftoa
  end interface vtoa

  ! as value numeric to asci without whitespace
  interface str_conv
    module procedure &
      str_convl, &            ! logical
      str_convi, &            ! integer
      str_convf
  end interface str_conv

  ! as user format for readability
  interface get_out
    module procedure &
      get_out_t, &
      get_out_comment, &
      get_out_l, &
      get_out_i, &
      get_out_f, &
      get_out_i2d, &
      get_out_i1d, &
      get_out_f2d, &
      get_out_f1d
  end interface get_out


contains

  function ltoi(aVal) result(res)
    ! 1 as true
    ! 0 as false
    logical  :: aVal
    integer  :: res

    if (aVal) then
      res = 1
    else
      res = 0
    end if
  end function


  !-----------------------------------------------
  ! asci to binary
  !-----------------------------------------------

  function atoi(astr) result(res)
    ! conversion from asci
    character(*), intent(in)   :: astr
    integer  :: res

    read (astr, *) res
  end function

  function atof(astr) result(res)
    ! conversion from asci
    character(*), intent(in)   :: astr
    real(dbl)      :: res

    read (astr, *) res
  end function


  function str_repl(astr, ori, rep) result(res)
    ! replace in astr
    character(*), intent(in)   :: astr, ori, rep
    character(len=:), allocatable    :: res          ! provide res as long as you want
    integer  :: i, j, lori

    res = ''
    if (len(astr) == 0) return
    lori = len_trim(ori) - 1
    i = 1                   ! start index 1
    do
      j = index(astr(i:), ori)
      if (j == 0) then        ! complete res
        res = res//astr(i:)
#if (CHECK == 1)
        ! res = '' implies in caller (example in color_ansi ansi_reset line 60)
        ! Fortran runtime error: Allocatable actual argument 'res' is not allocated
        ! if compile -fcheck=all
        ! write (0,*) 'str_repl "'//res//'"', len(res)
        ! if (len(res) == 0) res = ''  ! '£' or ' ' if problem stay, NULL is risky
#endif
        exit
      end if
      res = res//astr(i:i + j - 2)//rep
      i = i + j + lori
    end do
  end function


  !------------------------------------------------------
  ! binary to asci
  ! avoid '****' on leak formats, and digit truncatures
  !------------------------------------------------------

  function trunc_if(astr, lg) result(res)
    ! trunc left characters from astr only if they are whitespaces,
    ! target final length lg characters
    ! WITH appended trailing ONE whitespace
    character(*), intent(in)   :: astr
    character(len=:), allocatable    :: res
    integer  :: lg, lgcou, icou

    lgcou = len(astr)
    if (lg > 1) then        ! coherent (and if /= -1 of course)
      icou = lgcou - lg
      if ((icou < 1) .or. (icou > lgcou)) then         ! stupid lg ?
        res = trim(adjustl(astr))//' '
        return
      end if
      if (astr(icou:icou) == ' ') then                 ! ok
        res = astr(icou + 1:lgcou)//' '
        return
      end if
    end if
    res = trim(adjustl(astr))//' '
  end function


  function to_len(astr, lg) result(res)
    ! append leading whitespaces
    ! target final length lg characters,
    ! do not trunc if len astr greater lg yet
    character(*), intent(in)   :: astr
    character(len=:), allocatable    :: res
    integer  :: lg
    character(*), parameter    :: spaces = "                                               "

    if (lg > 1) then        ! coherent (and if /= -1 of course)
      if (len(astr) < lg) then
        res = spaces(1:lg - len(astr))//astr
        return
      end if
    end if
    res = astr              ! no shorter, no trunc digits
  end function


  function ltoa(aVal, lg) result(res)
    ! conversion from logical
    character(len=:), allocatable    :: res
    integer, optional    :: lg
    character(7)   :: tmp   ! '.true.' '.false.'
    logical  :: aVal
    integer  :: lgval

    lgval = -1
    if (present(lg)) lgval = lg

    write (tmp, '(L7)') aVal
    res = trunc_if(tmp, lgval)
  end function


  function itoa(aVal, lg) result(res)
    ! conversion from int
    character(len=:), allocatable    :: res
    integer, optional    :: lg
    character(15)  :: tmp
    integer  :: aVal, lgval

    lgval = -1
    if (present(lg)) lgval = lg

    write (tmp, '(i15)') aVal
    res = trunc_if(tmp, lgval)
  end function


  function iltoa(aVal, lg) result(res)
    ! conversion from int
    character(len=:), allocatable    :: res
    integer, optional    :: lg
    character(20)  :: tmp
    integer(8)     :: aVal
    integer  :: lgval

    lgval = -1
    if (present(lg)) lgval = lg

    write (tmp, '(i20)') aVal
    res = trunc_if(tmp, lgval)
  end function


  function array_ltoa(aVal, lg) result(res)
    ! conversion from int
    character(len=:), allocatable    :: res
    integer, optional    :: lg
    logical  :: aVal(:)
    integer  :: i, lgval

    lgval = -1
    if (present(lg)) lgval = lg

    res = ''
    do i = 1, size(aVal)
      res = res//ltoa(aVal(i), lgval)
    end do
  end function


  function array_itoa(aVal, lg) result(res)
    ! conversion from int
    character(len=:), allocatable    :: res
    integer, optional    :: lg
    integer  :: aVal(:)
    integer  :: i, lgval, ii, inl = 20

    lgval = -1
    if (present(lg)) lgval = lg

    res = ''
    ii = 0
    do i = 1, size(aVal)
      res = res//itoa(aVal(i), lgval)
      ii = ii + 1
      if (modulo(ii, inl) == 0) res = res//nl
    end do
  end function

  function array_ftoa(aVal, fm) result(res)
    ! conversion from float array
    character(*), optional     :: fm
    character(len=:), allocatable    :: res
    real(dbl)      :: aVal(:)
    integer  :: i, ii, inl = 10

    res = ''
    ii = 0
    do i = 1, size(aVal)
      if (present(fm)) then
        res = res//to_len(ftoa(aVal(i), fm), 10)
      else
        res = res//to_len(ftoa(aVal(i)), 10)
      end if
      ii = ii + 1
      if (modulo(ii, inl) == 0) res = res//nl
    end do
  end function

  function array_2d_itoa(aVal, lg) result(res)
    ! conversion from int
    character(len=:), allocatable    :: res
    integer, optional    :: lg
    integer  :: aVal(:, :)
    integer  :: i, j, lgval

    lgval = -1
    if (present(lg)) lgval = lg

    res = ''
    do i = lbound(aVal, 1), ubound(aVal, 1)
      do j = lbound(aVal, 2), ubound(aVal, 2)
        res = res//itoa(aVal(i, j), lgval)
      end do
      res = res//nl
    end do
  end function

  function array_2d_ftoa(aVal, fm) result(res)
    ! conversion from int
    character(len=:), allocatable    :: res
    character(*), optional     :: fm
    real(dbl)      :: aVal(:, :)
    integer  :: i, j


    res = ''
    do i = lbound(aVal, 1), ubound(aVal, 1)
      do j = lbound(aVal, 2), ubound(aVal, 2)
        if (present(fm)) then
          res = res//to_len(ftoa(aVal(i, j), fm), 10)
        else
          res = res//to_len(ftoa(aVal(i, j)), 10)
        end if
      end do
      res = res//nl
    end do
  end function


  function array_iltoa(aVal, lg) result(res)
    ! conversion from int long
    character(len=:), allocatable    :: res
    integer, optional    :: lg
    integer(8)     :: aVal(:)
    integer  :: i, lgval

    lgval = -1
    if (present(lg)) lgval = lg

    res = ''
    do i = 1, size(aVal)
      res = res//iltoa(aVal(i), lgval)
    end do
  end function


  function ftoa(aVal, fm) result(res)
    ! conversion from float real passe partout
    character(*), optional     :: fm
    real(dbl)      :: aVal
    character(len=:), allocatable    :: res

    if (present(fm)) then
      res = ftoaf(aVal, fm)
      return
    end if

    if (aval == 0d0) then
      res = '0.0 '
      return
    end if

    ! write (0, *) 'ftoa',aVal
    if ((aVal < 1e-2) .and. (aVal > -1e-2)) then
      res = ftoae(aVal)
      return
    end if
    res = ftoaf(aVal)
  end function


  function ftoae(aVal, fm) result(res)
    ! conversion from float real es12.3...
    character(*), optional     :: fm
    character(len=:), allocatable    :: fmval, res
    character(30)  :: tmp
    real(dbl)      :: aVal

    fmval = 'es10.2'
    if (present(fm)) fmval = trim(fm)

    write (tmp, '('//fmval//')') aVal
    res = trim(adjustl(tmp))//' '
  end function


  function ftoaf(aVal, fm) result(res)
    ! conversion from float real f10.2...
    character(*), optional     :: fm
    character(len=:), allocatable    :: fmval, res
    character(30)  :: tmp
    real(dbl)      :: aVal

    fmval = 'f10.2'
    if (present(fm)) fmval = trim(fm)

    write (tmp, '('//fmval//')') aVal
    if (tmp(1:1) == '*') then                        ! debordement format insuffisant
      res = ftoae(aVal)
    else
      if (aVal < 0d0) then    ! negative add leading '-' so add space if positive to align
        res = trim(adjustl(tmp))//' '
      else
        res = ' '//trim(adjustl(tmp))//' '
      end if
    end if
  end function


  function ftoag(aVal, fm) result(res)
    ! conversion from float real g10.4...
    character(*), optional     :: fm
    character(len=:), allocatable    :: fmval, res
    ! character(13) :: res
    character(20)  :: tmp
    real(dbl)      :: aVal

    fmval = 'g10.4'
    if (present(fm)) fmval = trim(fm)

    write (tmp, '('//fmval//')') aVal
    if (tmp(1:1) == '*') then                        ! debordement format insuffisant
      res = ftoae(aVal)
    else
      res = trim(adjustl(tmp))//' '
    end if
  end function



  function str_convi(aVal, fm) result(res)
    ! conversion from int, no leading and trailing whitespace
    character(*), optional     :: fm
    character(len=:), allocatable    :: fmval, res
    integer  :: aVal
    character(len=20)    :: conv

    fmval = 'i10'
    if (present(fm)) fmval = trim(fm)

    write (conv, '('//fmval//')') aVal
    res = trim(adjustl(conv))
  end function

  function str_0_convi(aVal, lg) result(res)
    ! returns 00123 for str_0_convi(123, 5)
    ! call log_info("begin"//str_0_convi(123, 5)//"end")
    character(len=:), allocatable    :: fmval, res
    integer  :: aVal, lg, lgcou

    lgcou = lg
    if ((lg > 10) .or. (lg < 2)) lgcou = 10
    fmval = 'i'//vtoa(lg)//'.'//vtoa(lg)
    res = str_convi(aVal, fmval)
  end function


  function str_convf(aVal, fm) result(res)
    ! conversion from real
    character(*), optional     :: fm
    character(len=:), allocatable    :: fmval, res
    real(dbl)      :: aVal
    character(len=30)    :: conv

    fmval = 'es12.5'
    if (present(fm)) fmval = trim(fm)

    write (conv, '('//fmval//')') aVal
    res = trim(adjustl(conv))
  end function


  function str_convl(aVal) result(res)
    ! conversion from logical
    character(len=:), allocatable    :: res
    logical  :: aVal
    character(len=2)     :: conv

    write (conv, "(l1)") aVal
    res = trim(adjustl(conv))
  end function

  function str_range(vmin, vmax) result(res)
    ! conversion from 2 integer values
    character(len=:), allocatable    :: res
    integer  :: vmin, vmax

    res = '['//trim(vtoa(vmin))//', '//trim(vtoa(vmax))//']'
  end function


  function str_trim(astr) result(res)
    ! remove TAB and NEW_LINE
    character(*)   :: astr
    character(len=:), allocatable    :: res          ! provide res as long as you want
    integer  :: i, lori

    res = ""
    lori = len_trim(astr)
    do i = 1, lori
      if ((astr(i:i) == nl) .or. (astr(i:i) == tab) .or. (astr(i:i) == ' ')) cycle
      res = res//astr(i:i)
    end do
  end function


  function toline(tit, val) result(res)
    ! format single 'data = value' for cmdc output
    ! align at position 15 the '=' (for readability)
    character(*)   :: tit, val
    character(len=:), allocatable    :: res
    character(15)  :: tmp

    write (tmp, '(a15)') tit
    res = trim(tmp)//' = '//trim(adjustl(val))
  end function


  function tosize(lbou, ubou) result(res)
    ! return 'lbou:ubou' as '-3:6' for example
    ! as format lbound:ubound
    character(len=:), allocatable    :: res
    integer  :: lbou, ubou

    res = str_convi(lbou)//':'//str_convi(ubou)
  end function


  !-----------------------------------------------
  ! interface get_out
  ! user formats for easy readable writes
  ! but not too long
  !-----------------------------------------------

  function get_fmt(fmt) result(res)
    character(*)   :: fmt
    character(len=:), allocatable    :: res
    integer  :: lmax

    if (tnr_precision < 1) then
      cmdc_precision = 2      ! get_env_i('tnr_precision', 5)
    end if
    if (cmdc_precision < 1) then
      cmdc_precision = 2      ! get_env_i('cmdc_precision', 5)
    end if

    lmax = cmdc_precision + 8                        ! ' +1.xxxxe+23' as len 7 (x not count)
    ! replace %es% as 12.5 for example
    res = str_repl(fmt, '%es%', trim(itoa(lmax))//'.'//trim(itoa(cmdc_precision)))
  end function


  function get_out_t(name) result(res)
    ! t as title
    character(len=:), allocatable    :: res
    character(*)   :: name
    character(1000)      :: tmp

    write (tmp, get_fmt(fmt_t)) '#########', ' '//name
    res = trim(tmp)
  end function

  function get_out_comment(name, comment) result(res)
    ! t as title
    character(len=:), allocatable    :: res
    character(*)   :: name, comment
    character(1000)      :: tmp

    write (tmp, get_fmt(fmt_t)) name, comment
    res = trim(tmp)
  end function

  function get_out_l(name, value) result(res)
    character(len=:), allocatable    :: res
    character(*)   :: name
    logical  :: value
    character(1000)      :: tmp

    write (tmp, get_fmt(fmt_l)) name, value
    res = trim(tmp)
  end function

  function get_out_i(name, value) result(res)
    character(len=:), allocatable    :: res
    character(*)   :: name
    integer  :: value
    character(1000)      :: tmp

    write (tmp, get_fmt(fmt_i)) name, value
    res = trim(tmp)
  end function

  function get_out_i2d(name, value) result(res)
    character(len=:), allocatable    :: res
    character(*)   :: name
    integer, dimension(:, :)   :: value
    character(len=:), allocatable    :: tmp

    tmp = array_2d_itoa(value, 3)
    tmp = str_repl(tmp, nl, nl//fmt_indent)          ! indent_fmt
    res = to_len(name, indent_fmt)//' '//tmp
  end function

  function get_out_i1d(name, value) result(res)
    character(len=:), allocatable    :: res
    character(*)   :: name
    integer, dimension(:)      :: value
    character(len=:), allocatable    :: tmp

    tmp = array_itoa(value, 3)
    tmp = str_repl(tmp, nl, nl//fmt_indent)          ! indent_fmt
    res = to_len(name, indent_fmt)//' '//tmp
  end function

  function get_out_f(name, value) result(res)
    character(len=:), allocatable    :: res
    character(*)   :: name
    real(dbl)      :: value
    character(1000)      :: tmp

    write (tmp, get_fmt(fmt_f)) name, value
    res = trim(tmp)
  end function

  function get_out_f2d(name, value) result(res)
    character(len=:), allocatable    :: res
    character(*)   :: name
    real(dbl), dimension(:, :) :: value
    character(len=:), allocatable    :: tmp

    tmp = array_2d_ftoa(value, get_fmt(fmt_f_vtoa))
    tmp = str_repl(tmp, nl, nl//fmt_indent)          ! indent
    res = to_len(name, indent_fmt)//' '//tmp
  end function

  function get_out_f1d(name, value) result(res)
    character(len=:), allocatable    :: res
    character(*)   :: name
    real(dbl), dimension(:)    :: value
    character(len=:), allocatable    :: tmp

    tmp = array_ftoa(value, get_fmt(fmt_f_vtoa))
    tmp = str_repl(tmp, nl, nl//fmt_indent)          ! indent
    res = to_len(name, indent_fmt)//' '//tmp
  end function


  function get_useful_len(mess) result(res)
    character(*)   :: mess
    integer  :: res, i

    do i = len(mess), 1, -1
      if (mess(i:i) == ' ') cycle
      if (mess(i:i) == nwl) cycle
      if (mess(i:i) == tbl) cycle
      exit
    end do
    res = i + 1
  end function


  function count_integers_in_string(str) result(num_integers)
    implicit none
    character(len=*) :: str
    integer :: num_integers
    integer :: i, len_str, start, end_pos, ios
    character(len=10) :: substr
    
    ! Get the length of the string
    len_str = len_trim(str)
    
    ! Initialize count
    num_integers = 0
    i = 1
    
    ! Loop through the string
    do while (i <= len_str)
        ! Find the first non-space character
        start = scan(str(i:), '0123456789')
        
        if (start == 0) exit  ! Exit if no digits are found
        
        ! Find the next space or the end of the string
        end_pos = index(str(i+start-1:), ' ')
        if (end_pos == 0) then
            end_pos = len_str - i + 2
        end if
        
        ! Extract the substring containing the number
        substr = str(i+start-1:i+start+end_pos-2)
        
        ! Try to read it as an integer
        read(substr, *, iostat=ios)
        if (ios == 0) then
            num_integers = num_integers + 1
        end if
        
        ! Move the index forward
        i = i + start + end_pos - 2
    end do
    
  end function count_integers_in_string


  !-----------------------------------------------
  ! test
  !-----------------------------------------------
  subroutine string_test()
    integer, parameter   :: im = 2, jm = 3
    integer, allocatable :: t1(:)
    integer, allocatable :: t2(:, :)
    real(dbl), allocatable     :: t3(:, :)
    integer  :: &
      i, j, &
      imin = -im, &
      imax = im, &
      jmin = -jm, &
      jmax = jm

    allocate (t1(imin:imax))
    allocate (t2(imin:imax, jmin:jmax))
    allocate (t3(imin:imax, jmin:jmax))

    do i = lbound(t1, 1), ubound(t1, 1)
      t1(i) = i
    end do

    do i = lbound(t2, 1), ubound(t2, 1)
      do j = lbound(t2, 2), ubound(t2, 2)
        t2(i, j) = i*100 + j
        if (t2(i, j) == 1) t2(i, j) = -10100000
      end do
    end do

    t3 = t2

    write (0, *) 't1('//tosize(imin, imax)//') = '//vtoa(t1)
    write (0, *) 't2('//tosize(imin, imax)//', '//tosize(jmin, jmax)//') = '//nl//vtoa(t2)
    write (0, *) 't2('//tosize(imin, imax)//', '//tosize(jmin, jmax)//') = '//nl//vtoa(t2, 6)
    write (0, *) 't3('//tosize(imin, imax)//', '//tosize(jmin, jmax)//') = '//nl//vtoa(t3)
    write (0, *) 't3('//tosize(imin, imax)//', '//tosize(jmin, jmax)//') = '//nl//vtoa(t3, 'es10.2')
  end subroutine


end module mld_string

