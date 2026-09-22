
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

subroutine diago_scalapack(N, W)
#ifdef MLD_NDM
  use gen_com_m, ONLY: angst
  use gen_com_m_ml, ONLY: im, imm, at
#else  
  use ondm_gen_com_m, ONLY: im, imm, angst, at
#endif
  !use mpi
  !use cholesky

  implicit none

  logical  :: leigenvectors = .true.
  integer  :: istat, info, i
  real(kind=8), dimension(:, :), allocatable   :: a
  real(kind=8), dimension(:, :), allocatable   :: c
  real(kind=8), dimension(:, :), allocatable   :: z
  integer  :: N
  real(kind=8)   :: W(N)
  ! real(kind=8),dimension(:),allocatable  :: solution
  integer, dimension(:), allocatable     :: ipiv
  real(kind=8)   :: temps1, temps2
  integer  :: root
  real(kind=8), dimension(:), allocatable      :: memory_a

  integer  :: context, iam
  integer  :: mycol, myrow, nb
  integer  :: npcol, nprocs, nprow
  integer  :: l_nrowsa, l_ncolsa
  integer  :: l_nrowsc, l_ncolsc

  ! PDSYEVX
  integer, allocatable, dimension(:, :)  :: ra, ca
  integer  :: imin, imax
  integer  :: nvectors, nvalues
  integer  :: nns, npp, npq, lwork, liwork         ! ,numroc
  integer, allocatable, dimension(:)     :: ifail, iclustr, iwork
  real(kind=8), allocatable, dimension(:)      :: work, gap
  real(kind=8)   :: orfac, vl, vu, abstol
  real(kind=8), external     :: pslamch

  integer, parameter   :: descriptor_len = 9
  integer  :: desca(descriptor_len)
  integer  :: descz(descriptor_len)
  integer  :: descc(descriptor_len)

  integer  :: numroc, ITEM
  integer  :: neig, nn, np0, mq0, anb, sqnpc, nps, nsytrd_lwopt
  integer, external    :: pjlaenv, iceil
  integer  :: itest


  root = 0

  ! Initialize the blacs. Note: processors are counted starting at 0.
  call blacs_pinfo(iam, nprocs)
  ITEM = 1000 + iam
  if (iam == 0) write (6, '(" nprocs = ",i4)') nprocs
  call gridsetup_ml(nprocs, nprow, npcol)
  call blacs_get(-1, 0, context)
  call BLACS_GRIDINIT(context, 'r', nprow, npcol)
  call BLACS_GRIDINFO(context, nprow, npcol, myrow, mycol)
  call blockset_ml(nb, 64, n, nprow, npcol)
  if (iam == 0) write (*, *) 'The block was set to ....', nb

  ! Distributed matrices: get num. local rows/cols. Create description.
  ! INTEGER FUNCTION NUMROC( N, NB, IPROC, ISRCPROC, NPROCS )
  ! NUMROC computes the NUMber of Rows Or Columns of a distributed
  ! matrix owned by the process indicated by IPROC.

  l_nrowsa = numroc(n, nb, myrow, 0, nprow)
  l_ncolsa = numroc(n, nb, mycol, 0, npcol)
  ! write (*,'("1 proc, n,nb,l_nrowsa,l_ncolsa  ",5i7)') iam, n,nb,l_nrowsa,l_ncolsa

  call descinit(desca, n, n, nb, nb, 0, 0, context, l_nrowsa, info)
  call descinit(descz, n, n, nb, nb, 0, 0, context, l_nrowsa, info)

  l_nrowsc = numroc(n, nb, myrow, 0, nprow)
  l_ncolsc = numroc(1, nb, mycol, 0, npcol)
  call descinit(descc, n, 1, nb, nb, 0, 0, context, l_nrowsc, info)

  ! Allocate LHS, RHS, pivot, and solution -----
  allocate (a(l_nrowsa, l_ncolsa), stat=istat)
  if (istat /= 0) stop "ERR:ALLOCATE FAILS for A"
  allocate (z(l_nrowsa, l_ncolsa), stat=istat)
  if (istat /= 0) stop "ERR:ALLOCATE FAILS for Z"
  allocate (c(l_nrowsc, l_ncolsc), stat=istat)
  if (istat /= 0) stop "ERR:ALLOCATE FAILS for C"
  allocate (ipiv(n), stat=istat)
  if (istat /= 0) stop "ERR:ALLOCATE FAILS for IPIV"
  allocate (memory_a(nprocs), stat=istat)
  if (istat /= 0) stop "ERR:ALLOCATE FAILS for memory_a"

  ! PDSYEVX
  allocate (gap(nprocs), stat=istat)
  if (istat /= 0) stop "ERR:ALLOCATE FAILS for gap"
  allocate (ifail(n), stat=istat)
  if (istat /= 0) stop "ERR:ALLOCATE FAILS for ifail"
  allocate (iclustr(2*nprow*npcol), stat=istat)
  if (istat /= 0) stop "ERR:ALLOCATE FAILS for iclustr"
  allocate (ra(nb, nb), stat=istat)
  if (istat /= 0) stop "ERR:ALLOCATE FAILS for ra"
  allocate (ca(nb, nb), stat=istat)
  if (istat /= 0) stop "ERR:ALLOCATE FAILS for ca"

  ! Initialize LHS and RHS
  ! call  read_my_matrix (context,iam,nprocs,n,a,nprow,npcol,myrow,mycol,desca)
  call init_my_matrix_ml(context, iam, nprocs, N, A, nprow, npcol, myrow, mycol, desca)

  ! symmetrization of the hessian matrix A = 1/2(A + A^T)
  Z(:, :) = A(:, :)
  call pdtran(N, N, 0.5d0, A, 1, 1, desca, 0.5d0, Z, 1, 1, descz)
  A(:, :) = Z(:, :)

  ! Show how arrays distributed
  ! PDSYEVX
  vl = -15.               ! eig min. if  'V' instead 'A'
  vu = 15.                ! eig max  if 'V' instead 'A'
  imin = -13              ! Ind  min.if  'I'  'A'
  imax = 13               ! Ind max. if 'I instead 'A'
  abstol = -1.            ! if  abstol=pslamch(context,'U') but dob't TOUCH
  abstol = pslamch(context, 'U')                   ! If  JOBZ='V',  setting  ABSTOL to PDLAMCH( CONTEXT, 'U') yields
  ! the most orthogonal eigenvectors
  orfac = 0.              ! Ortho

  temps1 = MPI_Wtime()

  nns = MAX(n, nb, 2)
  npp = numroc(nns, nb, 0, 0, nprow)
  npq = max(n, nprow*npcol + 1, 4)
  if (leigenvectors) then
    neig = n                ! number of requested eigenvectors
    nn = max(n, nb, 2)
    np0 = numroc(nn, nb, 0, 0, nprow)
    mq0 = numroc(max(neig, nb, 2), nb, 0, 0, npcol)
    lwork = 5*n + max(5*nn, np0*mq0 + 2*NB*NB) + iceil(neig, nprow*npcol)*nn

    anb = pjlaenv(desca, 3, 'pdsyttrd', 'L', 0, 0, 0, 0)
    sqnpc = int(sqrt(dble(NPROW*NPCOL)))
    nps = max(numroc(n, 1, 0, 0, sqnpc), 2*anb)
    nsytrd_lwopt = n + 2*(anb + 1)*(4*nps + 2) + (nps + 3)*nps

    lwork = max(2*lwork, 5*n + nsytrd_lwopt)
    ! write (*,*) 'setup', lwork, 7*n + nsytrd_lwopt
  else
    lwork = max(nb*(npp + 1), 5*nns) + 5*n
  end if

  liwork = 6*npq
  allocate (iwork(liwork), stat=istat)
  if (istat /= 0) stop "ERR:ALLOCATE FAILS for iwork"
  allocate (work(lwork), stat=istat)
  if (istat /= 0) stop "ERR:ALLOCATE FAILS for work"

  ! V/N   = eigenvalue and eigenvetore/eigenvalues
  ! A/V/I = all eignvalues / aleigenvalues in the interval [vl,vu] / eigenvalues index from  imin to imax
  ! U/L = upper / lower triangular
  ! n  = matrix dimension
  ! A  = the matrix
  ! 1,1 = the beginning of the matrix
  ! desca = descriptor
  ! vl,vu,imin,imax - described above
  ! abstol - if jobz (the first argument is set to V
  ! nvalues (output) - total number of eigenvalues
  ! nvectore (output) - total number of eigenvectors which fill the matrix Z

  if (leigenvectors) then
    if (iam == 0) print *, 'Eigenvalues and eigenvectors are computed'
    call PDSYEVX('V', 'A', 'L', n, A, 1, 1, desca, vl, vu, imin, imax, &
                 abstol, nvalues, nvectors, W, orfac, Z, 1, 1, descz, work, &
                 lwork, iwork, liwork, ifail, iclustr, gap, info)

    if (iam == 0) then
      print *, 'The number of eigenvalues  computed', nvalues
      print *, 'The number of eigenvectors computed', nvectors
      itest = 0
      do i = 1, n
        if (ifail(i) /= 0) then
          itest = 1
        end if
      end do
      if (itest == 1) then
        print *, 'WARNING .... some eigenvalues failed to converge'
      end if
    end if
    ! call get_solution (N,3*504,Z,DESCZ,eigvect)
    ! do i= 1, n
    !   call pdelget('A',' ',eigvect(i),z,i,3*504,descz)
    ! enddo

    ! if (iam==0) write (*,*)  eigvect (:)
  else
    call PDSYEVX('N', 'A', 'L', n, A, 1, 1, desca, vl, vu, imin, imax, &
                 abstol, nvalues, nvectors, W, orfac, Z, 1, 1, descz, work, &
                 lwork, iwork, liwork, ifail, iclustr, gap, info)
  end if

  if (info /= 0) then
    write (*, *) 'Diago not OK'
    write (*, *) 'info', info
  end if

  temps2 = MPI_Wtime()
  if (iam == 0) print *, 'That diago time ', temps2 - temps1, ' seconds'
  if (iam == 0) then
    do i = 1, n
      ! print '(25X,"W(",I5,") : ",1PE12.4)', i, W(i)
      ! write (10,'(i6,1D18.10)') i, W(i)
      write (10, *) i, w(i)
    end do
  end if

  deallocate (a, ipiv, memory_a, Z)
  ! Exit BLACS cleanly
  call blacs_gridexit(context)
  ! call blacs_exit( 0 )
end subroutine diago_scalapack

