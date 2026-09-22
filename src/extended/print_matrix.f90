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

function r8mat_norm_li(m, n, a)
  ! R8MAT_NORM_LI returns the matrix L-infinity norm of an R8MAT.
  !    An R8MAT is a two dimensional matrix of double precision real values.
  implicit none

  integer(kind=4) m
  integer(kind=4) n

  real(kind=8) a(m, n)
  integer(kind=4) i
  real(kind=8) r8mat_norm_li

  r8mat_norm_li = 0.0D+00
  do i = 1, m
    r8mat_norm_li = max(r8mat_norm_li, sum(abs(a(i, 1:n))))
  end do
end function r8mat_norm_li



subroutine r8mat_print(m, n, a, title)
  ! R8MAT_PRINT prints an R8MAT.
  ! An R8MAT is a two dimensional matrix of double precision real values.
  implicit none

  integer(kind=4) m
  integer(kind=4) n
  real(kind=8) a(m, n)
  character(len=*) title

  call r8mat_print_some(m, n, a, 1, 1, m, n, title)
end subroutine


subroutine r8mat_print_some(m, n, a, ilo, jlo, ihi, jhi, title)
  ! R8MAT_PRINT_SOME prints some of an R8MAT.
  ! An R8MAT is a two dimensional matrix of double precision real values.
  implicit none

  integer(kind=4), parameter :: incx = 5
  integer(kind=4) m
  integer(kind=4) n

  real(kind=8) a(m, n)
  character(len=14) ctemp(incx)
  integer(kind=4) i
  integer(kind=4) i2hi
  integer(kind=4) i2lo
  integer(kind=4) ihi
  integer(kind=4) ilo
  integer(kind=4) inc
  integer(kind=4) j
  integer(kind=4) j2
  integer(kind=4) j2hi
  integer(kind=4) j2lo
  integer(kind=4) jhi
  integer(kind=4) jlo
  character(len=*) title

  if (0 < len_trim(title)) then
    write (*, '(a)') ' '
    write (*, '(a)') trim(title)
  end if

  do j2lo = max(jlo, 1), min(jhi, n), incx

    j2hi = j2lo + incx - 1
    j2hi = min(j2hi, n)
    j2hi = min(j2hi, jhi)

    inc = j2hi + 1 - j2lo

    write (*, '(a)') ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write (ctemp(j2), '(i8,6x)') j
    end do

    write (*, '(''  Col   '',5a14)') ctemp(1:inc)
    write (*, '(a)') '  Row'
    write (*, '(a)') ' '

    i2lo = max(ilo, 1)
    i2hi = min(ihi, m)

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if (a(i, j) == real(int(a(i, j)), kind=8)) then
          write (ctemp(j2), '(f8.0,6x)') a(i, j)
        else
          write (ctemp(j2), '(g14.6)') a(i, j)
        end if

      end do

      write (*, '(i5,1x,5a14)') i, (ctemp(j), j=1, inc)

    end do
  end do
end subroutine


subroutine r8mat_uniform_01(m, n, seed, r)
  ! R8MAT_UNIFORM_01 fills an R8MAT with unit pseudorandom numbers.
  ! An R8MAT is a two dimensional matrix of double precision real values.
  implicit none

  integer(kind=4) m
  integer(kind=4) n

  integer(kind=4) i
  integer(kind=4) j
  integer(kind=4) k
  integer(kind=4) seed
  real(kind=8) r(m, n)

  do j = 1, n
    do i = 1, m
      k = seed/127773
      seed = 16807*(seed - k*127773) - k*2836
      if (seed < 0) then
        seed = seed + huge(seed)
      end if
      r(i, j) = real(seed, kind=8)*4.656612875D-10
    end do
  end do
end subroutine


subroutine r8vec_print(n, a, title)
  ! R8VEC_PRINT prints an R8VEC.
  !  Discussion:
  !    An R8VEC is an array of double precision real values.
  implicit none

  integer(kind=4) n

  real(kind=8) a(n)
  integer(kind=4) i
  character(len=*) title

  if (0 < len_trim(title)) then
    write (*, '(a)') ' '
    write (*, '(a)') trim(title)
  end if

  write (*, '(a)') ' '
  do i = 1, n
    write (*, '(2x,i8,2x,es14.5)') i, a(i)
  end do

end subroutine


subroutine r8vec_print_some(n, a, i_lo, i_hi, title)
  ! R8VEC_PRINT_SOME prints "some" of an R8VEC.
  !    An R8VEC is a vector of R8 values.
  implicit none

  integer(kind=4) n

  real(kind=8) a(n)
  integer(kind=4) i
  integer(kind=4) i_hi
  integer(kind=4) i_lo
  character(len=*) title

  if (0 < len_trim(title)) then
    write (*, '(a)') ' '
    write (*, '(a)') trim(title)
  end if

  write (*, '(a)') ' '
  do i = max(i_lo, 1), min(i_hi, n)
    write (*, '(2x,i8,2x,es14.5)') i, a(i)
  end do

end subroutine


subroutine timestamp()
  ! TIMESTAMP prints the current YMDHMS date as a time stamp.
  !  Example:
  !    31 May 2001   9:45:54.872 AM
  implicit none

  character(len=8) ampm
  integer(kind=4) d
  integer(kind=4) h
  integer(kind=4) m
  integer(kind=4) mm
  character(len=9), parameter, dimension(12)   :: month = (/ &
                                                  'January  ', 'February ', 'March    ', 'April    ', &
                                                  'May      ', 'June     ', 'July     ', 'August   ', &
                                                  'September', 'October  ', 'November ', 'December '/)
  integer(kind=4) n
  integer(kind=4) s
  integer(kind=4) values(8)
  integer(kind=4) y

  call date_and_time(values=values)

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if (h < 12) then
    ampm = 'AM'
  else if (h == 12) then
    if (n == 0 .and. s == 0) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if (h < 12) then
      ampm = 'PM'
    else if (h == 12) then
      if (n == 0 .and. s == 0) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write (*, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)') &
    d, trim(month(m)), y, h, ':', n, ':', s, '.', mm, trim(ampm)

end subroutine

module matrix_printer
    implicit none
contains
    subroutine int_printMatrix(matrix, title)
        integer, intent(in) :: matrix(:,:)
        character(len=*), intent(in) :: title
        integer :: ii, numRows, numCols
        integer, dimension(:), allocatable :: pcols
        integer, dimension(:), allocatable :: vect  
        character(len=128) :: CHFMT,  CHFMT_c
        
        numRows = size(matrix, 1)
        numCols = size(matrix, 2)
        ! Determine the field width based on the largest number
        if (allocated(vect)) deallocate(vect) ; allocate(vect(numCols))
        if (allocated(pcols)) deallocate(pcols) ; allocate(pcols(numCols))
        do ii = 1, numCols
          pcols(ii) = ii 
        end do 
        write (CHFMT_c, *) '("   -    ", ', numCols, '(i6), "  - " )'

        write (CHFMT, *) '(i5, " | ", ', numCols, '(i6), "  | " )'
      

        write(*,*) 
        write(*,'(a)') trim(title)
        write(*,'(a, i0, a, i0)') 'Matrix dimensions: ', numRows, 'x', numCols
        write(*,*) 


        ! Print the matrix with aligned columns
        write(*, CHFMT_c) pcols 
        do ii = 1, numRows
            ! Print row number
            vect(:) = matrix(ii, :)
            write(*, CHFMT) ii, vect(:)
        end do
    end subroutine int_printMatrix
end module matrix_printer