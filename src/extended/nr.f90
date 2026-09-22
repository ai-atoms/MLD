
! (C) Copr. 1986-92 Numerical Recipes Software 'k'1k30m,t+W.
SUBROUTINE indexx(n, arr, indx)
  implicit none

  INTEGER n, indx(n), M, NSTACK
  DOUBLE PRECISION arr(n)
  PARAMETER(M=7, NSTACK=50)
  INTEGER i, indxt, ir, itemp, j, jstack, k, l, istack(NSTACK)
  DOUBLE PRECISION a


  do j = 1, n             ! 11
    indx(j) = j
  end do                  ! 11
  jstack = 0
  l = 1
  ir = n

1 if (ir - l .lt. M) then
    do j = l + 1, ir        ! 13
      indxt = indx(j)
      a = arr(indxt)
      do i = j - 1, 1, -1     ! 12
        if (arr(indx(i)) .le. a) go to 2
        indx(i + 1) = indx(i)
      end do                  ! 12
      i = 0
2     indx(i + 1) = indxt
    end do                  ! 13
    if (jstack .eq. 0) return
    ir = istack(jstack)
    l = istack(jstack - 1)
    jstack = jstack - 2
  else
    k = (l + ir)/2
    itemp = indx(k)
    indx(k) = indx(l + 1)
    indx(l + 1) = itemp
    if (arr(indx(l + 1)) .gt. arr(indx(ir))) then
      itemp = indx(l + 1)
      indx(l + 1) = indx(ir)
      indx(ir) = itemp
    end if
    if (arr(indx(l)) .gt. arr(indx(ir))) then
      itemp = indx(l)
      indx(l) = indx(ir)
      indx(ir) = itemp
    end if
    if (arr(indx(l + 1)) .gt. arr(indx(l))) then
      itemp = indx(l + 1)
      indx(l + 1) = indx(l)
      indx(l) = itemp
    end if
    i = l + 1
    j = ir
    indxt = indx(l)
    a = arr(indxt)
3   continue
    i = i + 1
    if (arr(indx(i)) .lt. a) go to 3
4   continue
    j = j - 1
    if (arr(indx(j)) .gt. a) go to 4
    if (j .lt. i) go to 5
    itemp = indx(i)
    indx(i) = indx(j)
    indx(j) = itemp
    go to 3
5   indx(l) = indx(j)
    indx(j) = indxt
    jstack = jstack + 2
    if (jstack .gt. NSTACK) write (6, *) 'NSTACK too small in indexx'
    if (ir - i + 1 .ge. j - l) then
      istack(jstack) = ir
      istack(jstack - 1) = i
      ir = j - 1
    else
      istack(jstack) = j - 1
      istack(jstack - 1) = l
      l = i
    end if
  end if
  go to 1
END subroutine




