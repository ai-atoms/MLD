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

#include "../MLD_MACROS.INC"

module module_scalapack_tools
   use module_kind_variables, only: kind_double
   implicit none

contains

   subroutine gridsetup_ml(nproc, nprow, npcol)
      ! factorizes the number of processors (nproc) into nprow and npcol
      ! that are the sizes of the 2d processors mesh.
      !implicit none
      integer, intent(in)  ::  nproc
      integer, intent(out) ::  nprow, npcol
      integer ::  sqrtnp, i

      sqrtnp = int(sqrt(dble(nproc)) + 1)
      do i = 1, sqrtnp
         if (mod(nproc, i) .eq. 0) nprow = i
      end do
      npcol = nproc/nprow
   end subroutine gridsetup_ml

   subroutine blockset_ml(nb, nbuser, n, nprow, npcol)
      ! try to choose an optimal block size
      ! for the distributd matrix.
      !implicit none
      integer, intent(in)   :: n, nbuser, nprow, npcol
      integer, intent(out)  ::  nb

      nb = min(n/nprow, n/npcol)
      if (nbuser .gt. 0) then
         nb = min(nb, nbuser)
      end if
      nb = max(nb, 1)
   end subroutine blockset_ml

   subroutine up_scale_scalpack_matrix(AA, desc_AA, dimr_AA, dimc_AA, anrm, iascl) ! , l_dimr_AA, l_dimc_AA)
      ! Scalapack utilities ... 
      integer, parameter  ::  ctxt_ = 2
      real(kind_double), parameter ::  zero = 0.0d+0, one = 1.0d+0 

      ! Scalapack input output variables
      real(kind_double), dimension(:,:),  intent(inout) :: AA
      integer, dimension(:),  intent(inout)   :: desc_AA
      integer, intent(in) :: dimr_AA, dimc_AA ! ,  l_dimr_AA, l_dimc_AA
      integer, intent(out)  :: iascl
      real(kind_double), intent(out) :: anrm
      real(kind_double) :: smlnum, bignum
      real(kind_double), dimension(1) :: rwork  
      real(kind_double) ::  pdlamch, pdlange
      integer :: context, nprow, npcol, myrow, mycol, info  
      integer, parameter :: ia=1, ja=1

       context = desc_AA( ctxt_ )
       call blacs_gridinfo( context, nprow, npcol, myrow, mycol )

      ! Get machine parameters for possible scaling
      smlnum = pdlamch( context, 'S' )
      smlnum = smlnum / pdlamch( context, 'P' )
      bignum = one / smlnum
      call pdlabad( context, smlnum, bignum )

      !scale A, B matrix inthe range [SMLNUM,BIGNUM]
      anrm = pdlange( 'M', dimr_AA, dimc_AA, AA, ia, ja, desc_AA, rwork )
      iascl=0 
      if ((anrm > zero).and.(anrm < smlnum)) then 
        ! ...zero...anrm....smlnum.......bignum............
        ! scale matrix norm up to smlnum 
        call pdlascl( 'G', anrm, smlnum, dimr_AA, dimc_AA, AA, ia, ja, desc_AA, info )
        iascl=1 
      else if (anrm > bignum )  then 
        ! ...zero..........smlnum.......bignum.....anrm....
        ! scale matrix norm down to bignum
         call pdlascl( 'G', anrm, bignum, dimr_AA, dimc_AA, AA, ia, ja, desc_AA, info )
         iascl=2 
      else if (anrm==zero) then 
         !matrix all zero. Return zero solution ... 
         continue 
      end if    



!$!  *     Get machine parameters
!$!  *
!$!        smlnum = pdlamch( ictxt, 'S' )
!$!        smlnum = smlnum / pdlamch( ictxt, 'P' )
!$!        bignum = one / smlnum
!$!        CALL pdlabad( ictxt, smlnum, bignum )
!$!  *
!$!  *     Scale A, B if max entry outside range [SMLNUM,BIGNUM]
!$!  *
!$!        anrm = pdlange( 'M', m, n, a, ia, ja, desca, rwork )
!$!        iascl = 0
!$!        IF( anrm.GT.zero .AND. anrm.LT.smlnum ) THEN
!$!  *
!$!  *        Scale matrix norm up to SMLNUM
!$!  *
!$!           CALL pdlascl( 'G', anrm, smlnum, m, n, a, ia, ja, desca,
!$!       $                 info )
!$!           iascl = 1
!$!        ELSE IF( anrm.GT.bignum ) THEN
!$!  *
!$!  *        Scale matrix norm down to BIGNUM
!$!  *
!$!           CALL pdlascl( 'G', anrm, bignum, m, n, a, ia, ja, desca,
!$!       $                 info )
!$!           iascl = 2
!$!        ELSE IF( anrm.EQ.zero ) THEN
!$!  *
!$!  *        Matrix all zero. Return zero solution.
!$!  *
!$!           CALL pdlaset( 'F', max( m, n ), nrhs, zero, zero, b, ib, jb,
!$!       $                 descb )
!$!           GO TO 10
!$!        END IF
!$!  *
!$!        brow = m
!$!        IF( tpsd )
!$!       $   brow = n
!$!  *
!$!        bnrm = pdlange( 'M', brow, nrhs, b, ib, jb, descb, rwork )
!$!  *
!$!        ibscl = 0
!$!        IF( bnrm.GT.zero .AND. bnrm.LT.smlnum ) THEN
!$!  *
!$!  *        Scale matrix norm up to SMLNUM
!$!  *
!$!           CALL pdlascl( 'G', bnrm, smlnum, brow, nrhs, b, ib, jb,
!$!       $                 descb, info )
!$!           ibscl = 1
!$!        ELSE IF( bnrm.GT.bignum ) THEN
!$!  *
!$!  *        Scale matrix norm down to BIGNUM
!$!  *
!$!           CALL pdlascl( 'G', bnrm, bignum, brow, nrhs, b, ib, jb,
!$!       $                 descb, info )
!$!           ibscl = 2
!$!        END IF

   end subroutine up_scale_scalpack_matrix   

   subroutine up_scale_scalpack_systemAy(AA, desc_AA, dimr_AA, dimc_AA, yy, desc_yy, dimr_yy, dimc_yy, &
                                         anrm, bnrm, iascl, ibscl) 
      ! Scalapack utilities ... 
      integer, parameter  :: ctxt_ = 2
      real(kind_double), parameter ::  zero = 0.0d+0, one = 1.0d+0 

      ! Scalapack input output variables
      real(kind_double), dimension(:,:),  intent(inout) :: AA, yy 
      integer, dimension(:),  intent(inout)   :: desc_AA, desc_yy 
      integer, intent(in) :: dimr_AA, dimc_AA, dimr_yy, dimc_yy  ! ,  l_dimr_AA, l_dimc_AA
      integer, intent(out)  :: iascl, ibscl
      real(kind_double), intent(out) :: anrm, bnrm
      real(kind_double) :: smlnum, bignum 
      real(kind_double), dimension(1) :: rwork  
      real(kind_double) ::  pdlamch, pdlange
      integer :: context, nprow, npcol, myrow, mycol, info  
      integer, parameter :: ia=1, ja=1, ib=1, jb=1 


       context = desc_AA( ctxt_ )
       call blacs_gridinfo( context, nprow, npcol, myrow, mycol )

      ! Get machine parameters for possible scaling
      smlnum = pdlamch( context, 'S' )
      smlnum = smlnum / pdlamch( context, 'P' )
      bignum = one / smlnum
      call pdlabad( context, smlnum, bignum )

      !scale A, B matrix inthe range [SMLNUM,BIGNUM]
      anrm = pdlange( 'M', dimr_AA, dimc_AA, AA, ia, ja, desc_AA, rwork )
      iascl=0 
      if ((anrm > zero).and.(anrm < smlnum)) then 
        ! ...zero...anrm....smlnum.......bignum............
        ! scale matrix norm up to smlnum 
        call pdlascl( 'G', anrm, smlnum, dimr_AA, dimc_AA, AA, ia, ja, desc_AA, info )
        iascl=1 
      else if (anrm > bignum )  then 
        ! ...zero..........smlnum.......bignum.....anrm....
        ! scale matrix norm down to bignum
         call pdlascl( 'G', anrm, bignum, dimr_AA, dimc_AA, AA, ia, ja, desc_AA, info )
         iascl=2 
      else if (anrm==zero) then 
         !matrix all zero. Return zero solution ... 
         continue 
         call pdlaset('F', dimr_yy, dimc_yy, zero, zero, yy, ia, ja, desc_yy)
      end if    

! Now for yy ... 

      bnrm = pdlange( 'M', dimr_yy, dimc_yy, yy, ib, jb, desc_yy, rwork )
      ibscl = 0 
      if ((bnrm > zero).and.(bnrm < smlnum)) then 
        ! ...zero...anrm....smlnum.......bignum............
        call pdlascl( 'G', bnrm, smlnum, dimr_yy, dimc_yy, yy, ib, jb, desc_yy, info )
        ibscl=1 
      else if (bnrm > bignum) then 
        ! ...zero..........smlnum.......bignum.....anrm....
        call pdlascl( 'G', bnrm, bignum, dimr_yy, dimc_yy, yy, ib, jb, desc_yy, info )
        ibscl=2 
      end if  

   end subroutine up_scale_scalpack_systemAy   

   subroutine allocation_scalapack_matrix (subname, AA, desc_AA, dimr_AA, dimc_AA, &
      l_dimr_AA, l_dimc_AA, &
      nbr_AA, nbc_AA, &
      myrow, mycol, nprow, npcol, context, iam)

      use module_ml_scalapack, only: debug_scalapack, nbr_predefined, nbc_predefined
      use mld_logger
      !use module_kind_variables, only: kind_double
      integer, intent(in) :: dimr_AA, dimc_AA, myrow, mycol, nprow, npcol, iam, context
      integer, intent(inout) :: nbr_AA, nbc_AA, l_dimr_AA, l_dimc_AA
      real(kind_double), dimension(:,:), allocatable, intent(inout) :: AA
      integer, dimension(:), allocatable, intent(inout)   :: desc_AA
      character(len=*) :: subname
      integer, parameter :: descriptor_len = 9
      integer :: info, istat, lld
      integer :: numroc

      if (allocated(desc_AA)) deallocate (desc_AA); allocate (desc_AA(descriptor_len))
      if (nbr_AA < 0 ) then
         call blockset_ml(nbr_AA, nbr_predefined, dimr_AA, nprow, npcol)
      end if
      if (nbc_AA < 0 ) then
         call blockset_ml(nbc_AA, nbc_predefined, dimc_AA, nprow, npcol)
      end if

      if ((iam == 0).and.debug_scalapack)  write (6, '("ML sca: the block was set for row and column to ...........:", (a), 2i7)') trim(subname), nbr_AA, nbc_AA
      l_dimr_AA = numroc(dimr_AA, nbr_AA, myrow, 0, nprow)
      !dsca
      !l_dimc_AA = max(1, numroc(dimc_AA, nbc_AA, mycol, 0, npcol))
      l_dimc_AA = numroc(dimc_AA, nbc_AA, mycol, 0, npcol)
      if (l_dimr_AA < 0) then 
        call log_warning('ML sca: l_dimr_AA < 0')
      end if 
      if (l_dimc_AA < 0) then 
        call log_warning('ML sca: l_dimc_AA < 0')
      end if
      lld = max (1, l_dimr_AA)
      call descinit(desc_AA, dimr_AA, dimc_AA, nbr_AA, nbc_AA, 0, 0, context, lld, info)

      if (info < 0) then
         write (6, '( "ML sca: Illegal argument in descinit for ", i7,i6, (a) )') info, iam, trim(subname)
         write(6,'(6i8)') dimr_AA, dimc_AA, nbr_AA, nbc_AA, l_dimr_AA, l_dimc_AA
      end if
      if (info > 0) then
         write(*,*) "Error in descinit: On process ", iam, " local dimension is not large enough for the distributed matrix."
         write(*,*) "This often means l_dimr_AA or l_dimc_AA is too small."
         write(*,'(6i8)') iam, dimr_AA, dimc_AA, nbr_AA, nbc_AA, l_dimr_AA, l_dimc_AA
         stop "descinit failure: insufficient local array size"
      end if
      if (allocated(AA)) deallocate (AA); allocate (AA(l_dimr_AA, l_dimc_AA), stat=istat)

      if (istat /= 0) then
         write(6,*) 'error:  allocate fails for  ', trim(subname)
         stop "error allocate in allocation_scalapack_matrix with stop "
      end if
      call blacs_barrier(context, 'A')
   end subroutine allocation_scalapack_matrix

   subroutine pdgemm_nn(m, n, k, alpha, AA, desc_AA, BB, desc_BB, beta, CC, desc_CC)
      ! C(mxn) = alpha * A(mxk) * B(kxn) + beta * C, all matrices distributed, offsets (1,1).
      ! Replaces pdgemm('N','N',...): with MKL BLACS for OpenMPI (tested MKL 2025.1 + OpenMPI 5.0.7)
      ! pdgemm('N','N') segfaults on 2 or more MPI ranks, whereas the other transpose flags work.
      !  n == 1 : pdgemv('N'), no extra memory.
      !  n  > 1 : A^T is formed with pdtran, then pdgemm('T','N') gives A*B.
      integer, intent(in) :: m, n, k
      real(kind_double), intent(in) :: alpha, beta
      real(kind_double), dimension(:,:), intent(in) :: AA, BB
      real(kind_double), dimension(:,:), intent(inout) :: CC
      integer, dimension(:), intent(in) :: desc_AA, desc_BB, desc_CC

      real(kind_double), dimension(:,:), allocatable :: AT
      integer :: desc_AT(9)
      integer :: ctxt, nprow, npcol, myrow, mycol, info, mbAT, nbAT, lldAT, lcAT
      integer :: numroc
      real(kind_double), parameter :: zero = 0.d0, one = 1.d0

      if (n == 1) then
         call pdgemv('N', m, k, alpha, AA, 1, 1, desc_AA, BB, 1, 1, desc_BB, 1, &
            beta, CC, 1, 1, desc_CC, 1)
         return
      end if

      ! A^T (kxm): row blocks = column blocks of A, column blocks = row blocks of A
      ctxt = desc_AA(2)
      call blacs_gridinfo(ctxt, nprow, npcol, myrow, mycol)
      mbAT = desc_AA(6)
      nbAT = desc_AA(5)
      lldAT = max(1, numroc(k, mbAT, myrow, 0, nprow))
      lcAT  = max(1, numroc(m, nbAT, mycol, 0, npcol))
      call descinit(desc_AT, k, m, mbAT, nbAT, 0, 0, ctxt, lldAT, info)
      allocate(AT(lldAT, lcAT))

      call pdtran(k, m, one, AA, 1, 1, desc_AA, zero, AT, 1, 1, desc_AT)
      call pdgemm('T', 'N', m, n, k, alpha, AT, 1, 1, desc_AT, BB, 1, 1, desc_BB, &
         beta, CC, 1, 1, desc_CC)
      deallocate(AT)
   end subroutine pdgemm_nn


end module module_scalapack_tools


module module_init_ScaMatrix

   use module_kind_variables, only: kind_double
   implicit none

contains

   subroutine init_ScaMatrix_with_LocalVector_from_iproc(iproc, AA,  descAA, Iini, Jini, mlocal, dim_mlocal)
      ! Here we have a generic matrix AA (nprow, npcol) with ScaLapack desc integer descA(:) and integer context
      ! the vector mlocal on the the proc iproc (of pblas_grid) is distributed in AA
      ! integer : myrow, mycol are the local row and col index in the pblas_grid
      ! iproc : is the id of local proc on pblas_grid (not MPI d'accord ? )
      ! nproc_ml_sca : are the total number of pblas_grid.
      !
      !use module_kind_variables, only: kind_double
      use ml_in_ndm_module, only: rangml
      use time_check_general, only: MY_MPI_WTIME
      use module_ml_scalapack, only: nprocs_ml_sca, &
         npcol, myrow, mycol, context
      use mld_logger
      !implicit none

      integer, intent(in)   :: dim_mlocal, Iini, Jini, iproc
      integer, intent(in)   :: descAA(:)
      real(kind_double), intent(in)   :: AA(:, :), mlocal(:)

      real(kind_double)  :: AAIJ
      real(kind_double), dimension(:), allocatable :: mlocali

      integer  :: I, J, jptest, ilocal
      integer  :: prank, qrank
      !real(kind_double)  :: temps_i1, temps_i2, temps_i3, temps_f, temps_b, &
      !            temps_i4

      _NAMECURRENT_("init_ScaMatrix_with_LocalVector_from_iproc")


      _MLD_BEGIN_
      call BLACS_BARRIER(context, 'A')
      !temps_i1 = MY_MPI_Wtime()
      ! the case when all the matrices are written of the HDD iread_ml==1:
      ! write(*,*) 'INNNN', iproc, size(mlocal), size(AA, dim=1), size(AA, dim=2)
      if (allocated(mlocali)) deallocate(mlocali) ; allocate(mlocali(dim_mlocal))
      do jptest = 1, nprocs_ml_sca
         !temps_i3 = MY_MPI_Wtime()
         ! the second broadcast with the matrices
         if (iproc == (jptest - 1)) then
            prank = int(iproc / npcol)
            qrank = iproc - prank*npcol
            if (prank /= myrow) then
               if (rangml == 0) write (6, '("Problems in prank vs myrow  ", 2i7)') prank, myrow
            end if
            if (qrank /= mycol) then
               if (rangml == 0) write (6, '("Problems in qrank vs myrol  ", 2i7)') qrank, mycol
            end if
            !if (debug) then
            !  write (*, *) 'ML: proc, size mlocal mlocali', iproc, size(mlocali), size(mlocal)
            !end if
            mlocali(:) = mlocal(:)
            call DGEBS2D(context, 'All', ' ', dim_mlocal, 1, mlocal, dim_mlocal)
         end if

         if (iproc /= (jptest - 1)) then
            prank = int((jptest - 1)/npcol)
            qrank = (jptest - 1) - prank*npcol
            call DGEBR2D(context, 'All', ' ', dim_mlocal, 1, mlocali, dim_mlocal, prank, qrank)
         end if

         call BLACS_BARRIER(context, 'A')

         !temps_i4 = MY_MPI_Wtime()
         !temps_b = temps_b + (temps_i4 - temps_i3)
         do ilocal = 1, dim_mlocal
            AAIJ = mlocali(ilocal)
            I = Iini + ilocal - 1
            J = Jini
            CALL PDELSET(AA, I, J, descAA, AAIJ)
         end do
      end do

      _MLD_END_
   end subroutine  init_ScaMatrix_with_LocalVector_from_iproc

   subroutine init_distributed_ScaMatrix_with_LocalMatrix_Anida(AA,  descAA, Iini, Jini, mlocal, dimI_mlocal, dimJ_mlocal)

      ! Here we have a generic matrix AA (nprow, npcol) with ScaLapack desc integer descA(:) and integer context
      ! the matrix mlocal on the the proc iproc (of pblas_grid) is distributed in AA
      ! integer : myrow, mycol are the local row and col index in the pblas_grid
      ! nproc_ml_sca : are the total number of pblas_grid.

      use time_check_general, only: MY_MPI_WTIME
      use module_ml_scalapack, only: nprocs_ml_sca, iproc_sca, &
         npcol, myrow, mycol, context
      use mld_logger
      !implicit none

      integer, dimension(:), intent(in)   :: dimI_mlocal, dimJ_mlocal, Iini, Jini
      integer, intent(in)   :: descAA(:)
      real(kind_double), intent(in)   :: AA(:, :), mlocal(:,:)

      real(kind_double)  :: AAIJ
      real(kind_double), dimension(:,:), allocatable :: mlocali

      integer  :: I, J, proc, ilocal, jlocal
      integer  :: prank, qrank
      integer  :: dimI, dimJ

      _NAMECURRENT_("init_distributed_ScaMatrix_with_LocalMatrix_Anida")


      _MLD_BEGIN_
      call BLACS_BARRIER(context, 'A')
      do proc = 1, nprocs_ml_sca
         dimI = dimI_mlocal(proc); dimJ = dimJ_mlocal(proc)
         allocate(mlocali(dimI, dimJ))
         if (iproc_sca == (proc - 1)) then
            prank = int(iproc_sca / npcol)
            qrank = iproc_sca - prank*npcol
            if (prank /= myrow) then
               write (6, '("Problems in prank vs myrow  for proc", 3i7)') prank, myrow, iproc_sca
            end if
            if (qrank /= mycol) then
               write (6, '("Problems in qrank vs myrol  for proc", 3i7)') qrank, mycol, iproc_sca
            end if
            mlocali(:,:) = mlocal(:,:)
            call DGEBS2D(context, 'All', ' ', dimI, dimJ, mlocal, dimI)
         else
            prank = int((proc - 1)/npcol)
            qrank = (proc - 1) - prank*npcol
            call DGEBR2D(context, 'All', ' ', dimI, dimJ, mlocali, dimI, prank, qrank)
         end if

         !AnidaSCA call BLACS_BARRIER(context, 'A')
         do ilocal = 1, dimI
            do jlocal = 1, dimJ
               AAIJ = mlocali(ilocal, jlocal)
               I = Iini(proc) + ilocal
               J = Jini(proc) + jlocal
               CALL PDELSET(AA, I, J, descAA, AAIJ)
            end do
         end do
         deallocate(mlocali)
      end do

      _MLD_END_
   end subroutine init_distributed_ScaMatrix_with_LocalMatrix_Anida


subroutine init_distributed_ScaMatrix_with_LocalMatrix_mine(AA, descAA, Iini, Jini, mlocal, dimI_mlocal, dimJ_mlocal)
    use time_check_general, only: MY_MPI_WTIME
    use module_ml_scalapack, only: nprocs_ml_sca, iproc_sca, npcol, myrow, mycol, context
    use mld_logger
    implicit none

    integer, dimension(:), intent(in) :: dimI_mlocal, dimJ_mlocal, Iini, Jini
    integer, dimension(:), intent(in) :: descAA
    real(kind_double), dimension(:,:), intent(in) :: AA, mlocal

    integer :: I, J, proc, ilocal, jlocal, prank, qrank, dimI, dimJ
    real(kind_double), dimension(:,:), allocatable :: mlocali
    integer, dimension(2) :: iminmaxvec
    integer :: c_dimI, c_dimJ 
    real(kind_double)  :: AAIJ
    _NAMECURRENT_("init_distributed_ScaMatrix_with_LocalMatrix_mine")

    _MLD_BEGIN_
    call BLACS_BARRIER(context, 'A')

    ! Calculate the maximum size needed for the vectors

   

      _MLD_BEGIN_
      call BLACS_BARRIER(context, 'A')

      dimI = dimI_mlocal(iproc_sca+1); dimJ = dimJ_mlocal(iproc_sca+1)


      do proc = 1, nprocs_ml_sca
         if (iproc_sca == (proc - 1)) then
            c_dimI = dimI_mlocal(iproc_sca+1); c_dimJ = dimJ_mlocal(iproc_sca+1)
            iminmaxvec(1) = c_dimI 
            iminmaxvec(2) = c_dimJ
            prank = int(iproc_sca / npcol)
            qrank = iproc_sca - prank*npcol
            call IGEBS2D(context,'All',' ', 2, 1,iminmaxvec,2)
         else
            prank = int((proc-1) / npcol)
            qrank = (proc-1) - prank*npcol
            call IGEBR2D(context,'All',' ',2,1,iminmaxvec,2,prank,qrank)
            c_dimI = iminmaxvec(1)
            c_dimJ = iminmaxvec(2)
         end if   

         call BLACS_BARRIER( context, 'A' )

         allocate(mlocali(c_dimI, c_dimJ))

         if (iproc_sca == (proc - 1)) then

            prank = int(iproc_sca / npcol)
            qrank = iproc_sca - prank*npcol
            if (prank /= myrow) then
               write (6, '("Problems in prank vs myrow  for proc", 3i7)') prank, myrow, iproc_sca
            end if
            if (qrank /= mycol) then
               write (6, '("Problems in qrank vs myrol  for proc", 3i7)') qrank, mycol, iproc_sca
            end if
            mlocali(:,:) = mlocal(:,:)
            call DGEBS2D(context, 'All', ' ', c_dimI, c_dimJ, mlocal, c_dimI)
         else
            prank = int((proc - 1)/npcol)
            qrank = (proc - 1) - prank*npcol
            call DGEBR2D(context, 'All', ' ', c_dimI, c_dimJ, mlocali, c_dimI, prank, qrank)
         end if

         do ilocal = 1, c_dimI
            do jlocal = 1, c_dimJ
               AAIJ = mlocali(ilocal, jlocal)
               I = Iini(proc) + ilocal
               J = Jini(proc) + jlocal
               CALL PDELSET(AA, I, J, descAA, AAIJ)
            end do
         end do
         deallocate(mlocali)
      end do

      _MLD_END_


    _MLD_END_
end subroutine init_distributed_ScaMatrix_with_LocalMatrix_mine

subroutine init_distributed_ScaMatrix_with_LocalMatrix_from_proc(AA, descAA, Iini, Jini, mlocal, dimI_mlocal, dimJ_mlocal)
   use time_check_general, only: MY_MPI_WTIME
   use module_db_poscar, only: procs_per_file
   use module_ml_scalapack, only: nprocs_ml_sca, iproc_sca, context, nbr_predefined, nbc_predefined
   use module_scalapack_tools, only:  blockset_ml
   use module_scalapack_interfaces, only: my_pdgemr2d
   use mld_mpi, only: mld_rank, comm_mld
   use mld_logger
   implicit none
   integer, dimension(:), intent(in) :: dimI_mlocal, dimJ_mlocal, Iini, Jini
   integer, dimension(:), intent(in) :: descAA
   real(kind_double), dimension(:,:), intent(inout) :: AA
   real(kind_double), dimension(:,:), intent(in) :: mlocal
   integer :: iproc_mpi, group, root_proc, local_rows, local_cols, info
   integer, dimension(9) :: descLocal
   integer :: context_src,  isca, nbr_local, nbc_local  
   !$! integer :: ii, jj 
   integer, dimension(1) :: umap_dummy
   _NAMECURRENT_("init_distributed_ScaMatrix_with_LocalMatrix_from_proc")

   _MLD_BEGIN_
   iproc_mpi = mld_rank
   ! Define the group based on MPI rank
   group = iproc_mpi / procs_per_file
   root_proc = group * procs_per_file  ! The root process of each group

   call BLACS_BARRIER(context, 'A')

   do isca = 0, nprocs_ml_sca-1
     call blacs_get( 0, 0, context_src)
     umap_dummy(1) = isca
     call blacs_gridmap(context_src, umap_dummy, 1, 1, 1)
     if (iproc_sca .eq. isca) then
       local_rows = size(mlocal, 1)
       local_cols = size(mlocal, 2)
       nbr_local = -1 
       call blockset_ml(nbr_local, nbr_predefined, local_rows, 1, 1)
       nbc_local = -1 
       call blockset_ml(nbc_local, nbc_predefined, local_cols, 1, 1)
       !call descinit(descLocal, local_rows, local_cols, local_rows, local_cols, 0, 0, context_src, max(1,local_rows), info)
       call descinit(descLocal, local_rows, local_cols, nbr_local, nbc_local, 0, 0, context_src, max(1,local_rows), info)
     else
       descLocal(2) = -1
     end if
     call comm_mld%barrier()
     call my_pdgemr2d(dimI_mlocal(isca+1), dimJ_mlocal(isca+1), mlocal, 1, 1, descLocal, AA, Iini(isca+1)+1, Jini(isca+1)+1, descAA, context)                     
     call comm_mld%barrier()
     if (iproc_sca .eq. isca) then
     !if (iproc_sca .eq. isca .and. iproc_mpi == root_proc) then
       call blacs_gridexit(context_src)
     end if                
   end do 

   !$! local_rows = size(AA, 1)
   !$! local_cols = size(AA, 2)
   !$! do ii = 1, local_rows
   !$!   do jj = 1, local_cols
   !$!     if (dabs(AA(ii,jj)) < 1.d-50) AA(ii, jj) = 0.0d0
   !$!   end do
   !$! end do

   call BLACS_BARRIER(context, 'A')

   _MLD_END_
end subroutine init_distributed_ScaMatrix_with_LocalMatrix_from_proc


subroutine init_distributed_ScaMatrix_with_LocalMatrix_from_group(AA, descAA, Iini, Jini, mlocal, dimI_mlocal, dimJ_mlocal)
   use time_check_general, only: MY_MPI_WTIME
   use module_db_poscar, only: procs_per_file
   use module_ml_scalapack, only: nprocs_ml_sca, iproc_sca, context, nbr_predefined, nbc_predefined
   use module_scalapack_tools, only:  blockset_ml
   use mld_mpi, only: mld_rank, comm_mld
   use mld_logger
   implicit none

   integer, dimension(:), intent(in) :: dimI_mlocal, dimJ_mlocal, Iini, Jini
   integer, dimension(:), intent(in) :: descAA
   real(kind_double), dimension(:,:), intent(inout) :: AA
   real(kind_double), dimension(:,:), intent(in) :: mlocal

   integer :: iproc_mpi, group, root_proc, local_rows, local_cols, info
   integer, dimension(9) :: descLocal
   integer :: context_src
   integer, dimension(1) :: umap_dummy
   integer :: n_groups, gg, nbr_local, nbc_local

   _NAMECURRENT_("init_distributed_ScaMatrix_with_LocalMatrix_from_group")

   _MLD_BEGIN_
   iproc_mpi = mld_rank
   n_groups = nprocs_ml_sca / procs_per_file

   call BLACS_BARRIER(context, 'A')

   ! Iterate over groups
   do gg = 0, n_groups - 1
      root_proc = gg * procs_per_file  ! The root process of each group
      group = gg
      call blacs_get(0, 0, context_src)
      umap_dummy(1) = iproc_sca
      call blacs_gridmap(context_src, umap_dummy, 1, 1, 1)
      if (iproc_mpi == root_proc) then
         local_rows = size(mlocal, 1)
         local_cols = size(mlocal, 2)
         nbr_local = -1 
         call blockset_ml(nbr_local, nbr_predefined, local_rows, 1, 1)
         nbc_local = -1 
         call blockset_ml(nbc_local, nbc_predefined, local_cols, 1, 1)
         !call descinit(descLocal, local_rows, local_cols, local_rows, local_cols, 0, 0, context_src, max(1,local_rows), info)
         call descinit(descLocal, local_rows, local_cols, nbr_local, nbc_local, 0, 0, context_src, max(1,local_rows), info)
      else 
         descLocal(2) = -1
      end if

      call comm_mld%barrier()
      call pdgemr2d(dimI_mlocal(root_proc+1), dimJ_mlocal(root_proc+1), mlocal, 1, 1, descLocal, &
                    AA, Iini(root_proc+1)+1, Jini(root_proc+1)+1, descAA, context, info)
      call comm_mld%barrier()

      if (iproc_mpi == root_proc) then
        call blacs_gridexit(context_src)
      end if 

      call comm_mld%barrier()
    end do

    call BLACS_BARRIER(context, 'A')

    _MLD_END_
end subroutine init_distributed_ScaMatrix_with_LocalMatrix_from_group


   subroutine init_ScaMatrix_with_LocalVector(AA,  descAA, Iini, Jini, mlocal, dim_mlocal)
      ! Here we have a generic matrix AA (nprow, npcol) with ScaLapack desc integer descA(:) and integer context
      ! the vector mlocal with the dimension dim_mlocal is distributed over AA
      !
      !use module_kind_variables, only: kind_double
      use time_check_general, only: MY_MPI_WTIME
      use module_ml_scalapack, only: context
      use mld_logger
      !implicit none

      integer, intent(in)   :: dim_mlocal, Iini, Jini
      integer, intent(in)   :: descAA(:)
      real(kind_double), intent(in)   :: AA(:, :), mlocal(:)

      real(kind_double)  :: AAIJ

      integer  :: I, J, ilocal
      !real(kind_double)  :: temps_i1, temps_i2, temps_i3, temps_f, temps_b, &
      !            temps_i4

      _NAMECURRENT_("init_ScaMatrix_with_LocalVector")


      _MLD_BEGIN_

      do ilocal = 1, dim_mlocal
         AAIJ = mlocal(ilocal)
         I = Iini + ilocal - 1
         J = Jini
         CALL PDELSET(AA, I, J, descAA, AAIJ)
      end do

      call BLACS_BARRIER(context, 'A')
      _MLD_END_
   end subroutine  init_ScaMatrix_with_LocalVector


   subroutine init_ScaMatrix_with_LocalMatrix(AA,  descAA, Iini, Jini, matlocal, dim_in_I, dim_in_J)
      ! Here we have a generic matrix AA (nprow, npcol) with ScaLapack desc integer descA(:) and integer context
      ! the matrix matlocal available on all procs is distributed in AA pblas_grid.
      !
      !use time_check_general, only: MY_MPI_WTIME, debug_time
      !use module_kind_variables, only: kind_double
      use mld_logger
      use module_ml_scalapack, only: context
      !implicit none

      integer, intent(in)   :: Iini, Jini, dim_in_I, dim_in_J
      integer, intent(in)   :: descAA(:)
      real(kind_double), intent(in)   :: AA(:, :), matlocal(:,:)

      real(kind_double)  :: AAIJ
      integer  :: I, J, ilocal, jlocal

      _NAMECURRENT_("init_ScaMatrix_with_LocalMatrix")


      _MLD_BEGIN_

      do ilocal = 1, dim_in_I
         do jlocal = 1, dim_in_J
            AAIJ = matlocal(ilocal, jlocal)
            I = Iini + ilocal - 1
            J = Jini + jlocal - 1
            call PDELSET(AA, I, J, descAA, AAIJ)
         end do
      end do

      call BLACS_BARRIER(context, 'A')
      _MLD_END_
   end subroutine  init_ScaMatrix_with_LocalMatrix


end module module_init_ScaMatrix


module module_fit_ScaMatrix
   implicit none
contains

   subroutine prepare_sca_phi
      use module_kind_variables, only:  kind_double
      use module_ml_scalapack, only: sca_Amat, desc_sca_Amat, l_dimr_sca_Amat, l_dimc_sca_Amat, &
         nbr_Amat, nbc_Amat,  &
         dimr_sca_Amat, dimc_sca_Amat, &
         dimr_sca_Cmat, dimc_sca_Cmat, sca_Cmat, desc_sca_Cmat, &
         l_dimr_sca_Cmat, l_dimc_sca_Cmat, nbr_Cmat, nbc_Cmat, &
         dimr_sca_phi, dimc_sca_phi, sca_phi, desc_sca_phi, &
         l_dimr_sca_phi, l_dimc_sca_phi, nbr_phi, nbc_phi, &
         myrow, mycol, nprow, npcol, iam, context, &
         sca_ymat, desc_sca_ymat, &
         !dimr_sca_ymat, 
         dimc_sca_ymat, &
         sca_ymat_qr_svd, desc_sca_ymat_qr_svd, &
         dimr_sca_ymat_qr_svd, dimc_sca_ymat_qr_svd, &
         l_dimr_sca_ymat_qr_svd, l_dimc_sca_ymat_qr_svd, &
         nbr_ymat_qr_svd, nbc_ymat_qr_svd

      use module_scalapack_tools, only: allocation_scalapack_matrix, pdgemm_nn
      use ml_in_ndm_module, only : lambda_krr
      use snap, only: weights_snap
      use math, only: dreal_matmul
      use mld_logger
      ! implicit none
      integer :: ilocal, jglobal, iglobal, jlocal
      integer :: INDXL2G
      real(kind_double), parameter :: one=1.d0, zero=0.d0
      real(kind_double) :: alpha
      _NAMECURRENT_("prepare_sca_phi")


      _MLD_BEGIN_
      !sca_Cmat(:, :) = sca_Amat(:, :)
      ! the solution is (A^T w A)^-1 A^T W ymat or (Amat w Amat^T)^-1 Amat w ymat
      dimr_sca_phi = dimr_sca_Amat
      dimc_sca_phi = dimr_sca_Amat
      nbr_phi = nbr_Amat
      nbc_phi = nbr_Amat
      call allocation_scalapack_matrix (" sca_phi ", sca_phi, desc_sca_phi, dimr_sca_phi, dimc_sca_phi, &
         l_dimr_sca_phi, l_dimc_sca_phi, nbr_phi, nbc_phi, &
         myrow, mycol, nprow, npcol, context, iam  )

      ! allocate desc for Cmat generic version.
      nbr_Cmat = nbr_Amat
      nbc_Cmat = nbc_Amat
      dimr_sca_Cmat = dimr_sca_Amat
      dimc_sca_Cmat = dimc_sca_Amat
      call allocation_scalapack_matrix (" sca_Cmat ", sca_Cmat, desc_sca_Cmat, &
         dimr_sca_Cmat, dimc_sca_Cmat, &
         l_dimr_sca_Cmat, l_dimc_sca_Cmat, nbr_Cmat, nbc_Cmat, &
         myrow, mycol, nprow, npcol, context, iam  )


      ! ymat_qr_svd vector distribution: generic version.
      !         |      ...   |
      !         |      ...   |
      !  ymat = |  f(D) x 1  |  f(D) is line numbers of Amat
      !         |      ...   |
      !         |      ...   |
      nbr_ymat_qr_svd = nbr_Amat
      nbc_ymat_qr_svd = 1
      dimr_sca_ymat_qr_svd = dimr_sca_Amat
      dimc_sca_ymat_qr_svd = 1
      call allocation_scalapack_matrix (" sca_ymat_qr_svd ", sca_ymat_qr_svd, desc_sca_ymat_qr_svd, &
         dimr_sca_ymat_qr_svd, dimc_sca_ymat_qr_svd, &
         l_dimr_sca_ymat_qr_svd, l_dimc_sca_ymat_qr_svd, nbr_ymat_qr_svd, nbc_ymat_qr_svd, &
         myrow, mycol, nprow, npcol, context, iam  )


      !sca_Cmat(:, :) = sca_Amat(:, :)
      ! the solution is (A^T w A)^-1 A^T W ymat or (Amat w Amat^T)^-1 Amat w ymat
      ! Cmat = Amat x w

      !serial operation
      !do i = 1, size(Cmat, 2)
      !  Cmat(:, i) = Amat(:, i)*weights_snap(i)
      !end do
      !write(6,'("before parallel version  ", 5i8)')  rangml, dimr_sca_Amat, dimc_sca_Amat,l_dimr_sca_Amat, l_dimc_sca_Amat
      call BLACS_BARRIER(context, 'A')

      !parallel version
      do ilocal = 1, l_dimr_sca_Amat
         !iglobal = INDXL2G(ilocal, nb_r, myrow, 0, nprow)
         do jlocal = 1, l_dimc_sca_Amat
            jglobal = INDXL2G(jlocal, nbc_Amat, mycol, 0, npcol)
            sca_Cmat(ilocal, jlocal) = sca_Amat(ilocal, jlocal)*weights_snap(jglobal)
         end do
      end do

      !if (rangml==0) write(6,*) 'before BLACS_BARRIER'


      call BLACS_BARRIER(context, 'A')

      !! ToScaRemove
      !! Compare Amat to sca_Amat
      !write(fname_pdwrite, '(A)') 'sca_Amat_matrix.dat'
      !call PDLAWRITE(trim(fname_pdwrite), dimr_sca_Amat, dimc_sca_Amat, sca_Amat, 1, 1, desc_sca_Amat, 0, 0, prnwork)
      !if (rangml == 0) then
      !open(unit=unitA, name='ser_Amat_matrix.dat', status='unknown')
      !  do jlocal = 1, dimc_sca_Amat
      !    do ilocal = 1, dimr_sca_Amat
      !      write(unitA, '(E30.18)') Amat(ilocal, jlocal)
      !  end do
      !end do
      !end if
      !close(unitA, status='keep')
      ! end ToScaRemove

      !$! ! ToScaRemove
      !$! M_obsv = size(Amat, 2)
      !$! P_desc = size(Amat, 1)
      !$!
      !$! if (allocated(phi_diag)) deallocate (phi_diag); allocate (phi_diag(P_desc, P_desc))
      !$! if (allocated(phi)) deallocate (phi); allocate (phi(P_desc, P_desc))
      !$!
      !$! if (allocated(Cmat)) deallocate (Cmat); allocate (Cmat(P_desc, M_obsv))
      !$! if (allocated(Cmat_transpose)) deallocate (Cmat_transpose); allocate (Cmat_transpose(M_obsv, P_desc))
      !$! ! end ToScaRemove
      !$!
      !$! ! ToScaRemove
      !$! ! compare Cmat to sca_Cmat
      !$! write(fname_pdwrite, '(A)') 'sca_Cmat_matrix.dat'
      !$! call PDLAWRITE(trim(fname_pdwrite), dimr_sca_Amat, dimc_sca_Amat, sca_Cmat, 1, 1, desc_sca_Cmat, 0, 0, prnwork)
      !$! ! Cmat = Amat x w
      !$! do iglobal = 1, size(Cmat, 2)
      !$!   Cmat(:, iglobal) = Amat(:, iglobal)*weights_snap(iglobal)
      !$! end do
      !$! if (rangml == 0) then
      !$!   open(unit=unitA, name='ser_Cmat_matrix.dat', status='unknown')
      !$!     do jlocal = 1, dimc_sca_Amat
      !$!       do ilocal = 1, dimr_sca_Amat
      !$!         write(unitA, '(E30.18)') Cmat(ilocal, jlocal)
      !$!     end do
      !$!   end do
      !$! end if
      !$! close(unitA, status='keep')
      !$! ! end ToScaRemove


      ! serial operation
      ! Cmat_transpose(:, :) = transpose(Cmat)
      ! phi = matmul ( Amat, Cmat_transpose)
      !call dreal_matmul(Amat, size(Amat, 1), size(Amat, 2), Cmat_transpose, size(Cmat_transpose, 1), size(Cmat_transpose, 2), &
      !  phi, size(phi, 1), size(phi, 2))
      ! phi = Amat * Amat^T
      ! c =a * b   mxn = mxk * kxn
      ! call pdgemm(transa, transb, m, n, k, alpha, a, ia, ja, desca, b, ib, jb, descb, beta, c, ic, jc, descc)
      !                                             b, ib, jb, descb,
      !                                      beta,  c, ic, jc, descc)
      call pdgemm ('N', 'T', dimr_sca_Amat, dimr_sca_Amat, dimc_sca_Amat, one, &
         sca_Amat, 1, 1, desc_sca_Amat, &
         sca_Cmat, 1, 1, desc_sca_Cmat, &
         zero, sca_phi, 1, 1, desc_sca_phi)

      ! solving with regularization
      if (lambda_krr < 0 ) lambda_krr = 0.d0

      !$! ! ToScaRemove
      !$! ! Compare phi to sca_phi
      !$! Cmat_transpose(:, :) = transpose(Cmat)
      !$! call dreal_matmul(Amat, size(Amat, 1), size(Amat, 2), Cmat_transpose, size(Cmat_transpose, 1), size(Cmat_transpose, 2), &
      !$!               phi, size(phi, 1), size(phi, 2))
      !$! if (lambda_krr < 0) lambda_krr = 0.d0
      !$! do ilocal = 1, size(Amat, 1)
      !$!   phi(ilocal,ilocal) = phi(ilocal,ilocal) + lambda_krr
      !$! end do
      !$! ! end ToScaRemove


      do iglobal = 1, dimr_sca_phi
         call pdelget('A', ' ', alpha, sca_phi, iglobal, iglobal, desc_sca_phi)
         call pdelset(sca_phi, iglobal, iglobal, desc_sca_phi, alpha + lambda_krr)
      end do

      !$! ! ToScaRemove
      !$! write(fname_pdwrite, '(A)') 'sca_phi_matrix.dat'
      !$! call PDLAWRITE(trim(fname_pdwrite), dimr_sca_phi, dimr_sca_phi, sca_phi, 1, 1, desc_sca_phi, 0, 0, prnwork)
      !$! if (rangml == 0) then
      !$!   open(unit=unitA, name='ser_phi_matrix.dat', status='unknown')
      !$!   do jlocal = 1, dimr_sca_phi
      !$!     do ilocal = 1, dimr_sca_phi
      !$!       write(unitA, '(E30.18)') phi(ilocal, jlocal)
      !$!     end do
      !$!   end do
      !$! end if
      !$! close(unitA, status='keep')
      !$! ! end ToScaRemove

      !Mline = size(phi, 1)
      !Ncolm = size(phi, 2)
      !!LDA = max(Mline, 1)
      !!LDB = max(Mline, max(Ncolm, 1))
      !! prepare the ymat_qr_svd = A w y ... je crois :)

      !$! ! ToScaRemove
      !$! if (allocated(ymat_qr_svd)) deallocate (ymat_qr_svd); allocate (ymat_qr_svd(size(Cmat, 1), size(ymat, 2)))
      !$! call dreal_matmul(Cmat, size(Cmat, 1), size(Cmat, 2), ymat, size(ymat, 1), size(ymat, 2), &
      !$!                 ymat_qr_svd, size(ymat_qr_svd, 1), size(ymat_qr_svd, 2))
      !$! ! end ToScaRemove

      ! c =Cmat * y   mxn = mxk * kxn
      ! call pdgemm(transa, transb, m, n, k, alpha, a, ia, ja, desca, b, ib, jb, descb, beta, c, ic, jc, descc)
      call pdgemm_nn(dimr_sca_Cmat, dimc_sca_ymat, dimc_sca_Cmat, one, sca_Cmat, desc_sca_Cmat, &
         sca_ymat, desc_sca_ymat, &
         zero, sca_ymat_qr_svd, desc_sca_ymat_qr_svd)

      !$! ! ToScaRemove
      !$! write(fname_pdwrite, '(A)') 'sca_ymat_matrix.dat'
      !$! call PDLAWRITE(trim(fname_pdwrite), dimr_sca_ymat_qr_svd, dimc_sca_ymat_qr_svd, sca_ymat_qr_svd, 1, 1, desc_sca_ymat_qr_svd, 0, 0, prnwork)
      !$! if (rangml == 0) then
      !$!   write(6,*) dimr_sca_ymat, dimc_sca_ymat, size(ymat,1),  size(ymat,2)
      !$!   open(unit=unitA, name='ser_ymat_matrix.dat', status='unknown')
      !$!   do jlocal = 1, dimc_sca_ymat_qr_svd
      !$!     do ilocal = 1, dimr_sca_ymat_qr_svd
      !$!       write(unitA, '(E30.18)') ymat_qr_svd(ilocal, jlocal)
      !$!     end do
      !$!   end do
      !$! end if
      !$! close(unitA, status='keep')
      !$! ! end ToScaRemove

      _MLD_END_
   end subroutine prepare_sca_phi

   subroutine prepare_sca_Amat_big
      use module_kind_variables, only:  kind_double
      use module_ml_scalapack, only: sca_Amat,  &
         nbr_Amat, nbc_Amat,  &
         dimr_sca_Amat, dimc_sca_Amat, &
         myrow, mycol, nprow, npcol, iam, context, &
         dimr_sca_ymat, dimc_sca_ymat, &
         dimr_sca_ymat, dimc_sca_ymat, &
         nbr_ymat, nbc_ymat, &
         sca_ymat_qr_svd, desc_sca_ymat_qr_svd, &
         dimr_sca_ymat_qr_svd, dimc_sca_ymat_qr_svd, &
         l_dimr_sca_ymat_qr_svd, l_dimc_sca_ymat_qr_svd, &
         nbr_ymat_qr_svd, nbc_ymat_qr_svd, &
         sca_Amat_big, desc_sca_Amat_big, l_dimc_sca_Amat_big, l_dimr_sca_Amat_big, &
         dimc_sca_Amat_big, dimr_sca_Amat_big, nbr_Amat_big, nbc_Amat_big

      use module_scalapack_tools, only: allocation_scalapack_matrix
      use module_scalapack_interfaces, only: my_pdgemr2d
      use ml_in_ndm_module, only : lambda_krr
      use snap, only: weights_snap, ymat
      use math, only: dreal_matmul
      use mld_logger
      ! implicit none
      integer :: ilocal, jglobal, iglobal, jlocal
      integer :: INDXL2G
      real(kind_double), parameter :: one=1.d0, zero=0.d0
      !print matrix ...
      !character(len=80) :: fname_pdwrite
      !real(kind_double) :: prnwork(500)
      real(kind_double) :: wtmp

      integer  :: nbr_tmpDiag, nbc_tmpDiag
      real(kind_double), dimension(:, :), allocatable   :: sca_tmpDiag
      integer, dimension(:), allocatable     :: desc_sca_tmpDiag
      integer  :: dimr_sca_tmpDiag, dimc_sca_tmpDiag
      integer  :: l_dimr_sca_tmpDiag, l_dimc_sca_tmpDiag

      integer  :: nbr_Amat_w, nbc_Amat_w
      real(kind_double), dimension(:, :), allocatable   :: sca_Amat_w
      integer, dimension(:), allocatable     :: desc_sca_Amat_w
      integer  :: dimr_sca_Amat_w, dimc_sca_Amat_w
      integer  :: l_dimr_sca_Amat_w, l_dimc_sca_Amat_w

      integer  :: nbr_ymat_w, nbc_ymat_w
      real(kind_double), dimension(:, :), allocatable   :: sca_ymat_w
      integer, dimension(:), allocatable     :: desc_sca_ymat_w
      integer  :: dimr_sca_ymat_w, dimc_sca_ymat_w
      integer  :: l_dimr_sca_ymat_w, l_dimc_sca_ymat_w


      integer :: ii


      _NAMECURRENT_("prepare_sca_Amat_big")

      _MLD_BEGIN_

      !  prepare sca_Amat_w and sca_ymat_w
      dimr_sca_Amat_w = dimr_sca_Amat
      dimc_sca_Amat_w = dimc_sca_Amat
      nbr_Amat_w = nbr_Amat
      nbc_Amat_w = nbc_Amat
      call allocation_scalapack_matrix (" sca_Amat_w ", sca_Amat_w, desc_sca_Amat_w, dimr_sca_Amat_w, dimc_sca_Amat_w, &
         l_dimr_sca_Amat_w, l_dimc_sca_Amat_w, nbr_Amat_w, nbc_Amat_w, &
         myrow, mycol, nprow, npcol, context, iam  )

      do ilocal = 1, l_dimr_sca_Amat_w
         !iglobal = INDXL2G(ilocal, nb_r, myrow, 0, nprow)
         do jlocal = 1, l_dimc_sca_Amat_w
            jglobal = INDXL2G(jlocal, nbc_Amat_w, mycol, 0, npcol)
            sca_Amat_w(ilocal, jlocal) = sca_Amat(ilocal, jlocal)*dsqrt(weights_snap(jglobal))
         end do
      end do

      dimr_sca_ymat_w = dimr_sca_ymat
      dimc_sca_ymat_w = dimc_sca_ymat
      nbr_ymat_w = nbr_ymat
      nbc_ymat_w = nbc_ymat
      call allocation_scalapack_matrix (" sca_ymat_w ", sca_ymat_w, desc_sca_ymat_w, dimr_sca_ymat_w, dimc_sca_ymat_w, &
         l_dimr_sca_ymat_w, l_dimc_sca_ymat_w, nbr_ymat_w, nbc_ymat_w, &
         myrow, mycol, nprow, npcol, context, iam  )

      !do ilocal = 1, l_dimr_sca_ymat_w
      !  iglobal = INDXL2G(ilocal, nbr_ymat_w, myrow, 0, nprow)
      !  write(23,'(3i9,2es20.10)') rangml, ilocal, iglobal, sca_ymat(ilocal,1), dsqrt(weights_snap(iglobal))
      !  sca_ymat_w(ilocal,1) = sca_ymat(ilocal,1)*dsqrt(weights_snap(iglobal))
      !end do

      do iglobal = 1, dimr_sca_ymat_w
         wtmp = ymat(iglobal,1)*dsqrt(weights_snap(iglobal))
         call pdelset(sca_ymat_w, iglobal, 1, desc_sca_ymat_w, wtmp)
      end do

      ! Amat_bid: generic version.
      !                 f(D)
      !           |      ...    |
      !           |      ...    |
      !         M |   Amat_w^T  |
      !           |      ...    |
      !Amat_big = |      ...    |
      !            -------------
      !           |      ...    |
      !           |      ...    |
      !       f(D)|  f(D) x f(D)|  f(D) is row  numbers of Amat
      !           |      ...    |  M    is colm numbers of Amat
      !           |      ...    |



      dimr_sca_Amat_big = dimc_sca_Amat + dimr_sca_Amat
      dimc_sca_Amat_big = dimr_sca_Amat
      nbr_Amat_big = nbc_Amat
      nbc_Amat_big = nbr_Amat
      call allocation_scalapack_matrix (" sca_Amat_big ", sca_Amat_big, desc_sca_Amat_big, dimr_sca_Amat_big, dimc_sca_Amat_big, &
         l_dimr_sca_Amat_big, l_dimc_sca_Amat_big, nbr_Amat_big, nbc_Amat_big, &
         myrow, mycol, nprow, npcol, context, iam  )

      ! allocate sca_tmpDiag:  generic version.
      nbr_tmpDiag = nbr_Amat
      nbc_tmpDiag = nbr_Amat
      dimr_sca_tmpDiag = dimr_sca_Amat
      dimc_sca_tmpDiag = dimr_sca_Amat
      call allocation_scalapack_matrix (" sca_tmpDiag ", sca_tmpDiag, desc_sca_tmpDiag, &
         dimr_sca_tmpDiag, dimc_sca_tmpDiag, &
         l_dimr_sca_tmpDiag, l_dimc_sca_tmpDiag, nbr_tmpDiag, nbc_tmpDiag, &
         myrow, mycol, nprow, npcol, context, iam  )


      do ilocal = 1, l_dimr_sca_tmpDiag
         do jlocal = 1, l_dimc_sca_tmpDiag
            sca_tmpDiag(ilocal, jlocal) = 0.d0
         end do
      end do

      ! solving with regularization
      if (lambda_krr < 0 ) lambda_krr = 0.d0

      do ii = 1, dimr_sca_tmpDiag
         call pdelset(sca_tmpDiag, ii, ii, desc_sca_tmpDiag, lambda_krr)
      end do

      call my_pdgemr2d(dimr_sca_tmpDiag, dimc_sca_tmpDiag,  sca_tmpDiag, 1, 1, desc_sca_tmpDiag, &
         sca_Amat_big, dimc_sca_Amat+1, 1, desc_sca_Amat_big, context)


      ! call pdtran(m, n, alpha, a, ia, ja, desca, beta, c, ic, jc, descc)
      ! The p?tran routines transpose a real distributed matrix. The operation is defined as
      ! sub(C):=beta*sub(C) + alpha*sub(A)',
      ! m Specifies the number of rows of the distributed matrix sub(C), m≥ 0.
      ! n Specifies the number of columns of the distributed matrix sub(C) , n≥ 0.

      call pdtran(dimc_sca_Amat_w, dimr_sca_Amat_w, one, sca_Amat_w, 1, 1, desc_sca_Amat_w, zero, sca_Amat_big, 1, 1, desc_sca_Amat_big)


      ! ymat_qr_svd vector distribution: generic version.
      !               |      ...   |
      !               |      ...   |
      !               |     M x 1  |   M is row numbers of Amat_w
      !               |      ...   |
      !               |      ...   |
      ! ymat_qr_svd = ------------
      !               |      ...   |
      !               |      ...   |
      !               |  f(D) x 1  |  f(D) is line numbers of Amat
      !               |      ...   |
      !               |      ...   |
      if (allocated(sca_ymat_qr_svd)) deallocate(sca_ymat_qr_svd)
      if (allocated(desc_sca_ymat_qr_svd)) deallocate(desc_sca_ymat_qr_svd)
      nbr_ymat_qr_svd =  nbc_Amat
      nbc_ymat_qr_svd = 1
      dimr_sca_ymat_qr_svd = dimc_sca_Amat + dimr_sca_Amat
      dimc_sca_ymat_qr_svd = 1

      call allocation_scalapack_matrix (" sca_ymat_qr_svd ", sca_ymat_qr_svd, desc_sca_ymat_qr_svd, &
         dimr_sca_ymat_qr_svd, dimc_sca_ymat_qr_svd, &
         l_dimr_sca_ymat_qr_svd, l_dimc_sca_ymat_qr_svd, nbr_ymat_qr_svd, nbc_ymat_qr_svd, &
         myrow, mycol, nprow, npcol, context, iam  )

      do iglobal = dimr_sca_ymat_w+1, dimr_sca_ymat_qr_svd
         call pdelset(sca_ymat_qr_svd, iglobal, 1, desc_sca_ymat_qr_svd, 0.d0 )
      end do

      call my_pdgemr2d(dimr_sca_ymat_w, dimc_sca_ymat_w,  sca_ymat_w, 1, 1, desc_sca_ymat_w, &
         sca_ymat_qr_svd, 1, 1, desc_sca_ymat_qr_svd, context)


      call BLACS_BARRIER(context, 'A')

      if (allocated(desc_sca_Amat_w))  deallocate (desc_sca_Amat_w)
      if (allocated(sca_Amat_w))  deallocate (sca_Amat_w)

      if (allocated(desc_sca_ymat_w))  deallocate (desc_sca_ymat_w)
      if (allocated(sca_ymat_w))  deallocate (sca_ymat_w)


      _MLD_END_
   end subroutine prepare_sca_Amat_big

   !subroutine scalapack_lsystem_by_homeLU(sca_phi, desc_sca_phi, sca_ymat_qr_svd, desc_sca_ymat_qr_svd, &
   !                                   Mline, Ncolm, w_params)
   subroutine scalapack_lsystem_by_homeLU(Mline, Ncolm, sca_AA, desc_sca_AA, nbr_AA, nbc_AA, &
      dimr_yy, dimc_yy, sca_yy, desc_sca_yy, w_params)
      !  Solve the system sca_AA * w = ymat   MxN * Nx1 = Mx1
      !  sca_AA is often sca_phi matrix  = Amat W Amat^T ( f(D) x f(D) matrix )
      !       N= f(D)          1         1
      ! |               |             |     |
      ! |               |   |     |   |     |
      ! |               |   |     |   |     |
      ! | M   sca_A     | x | N w | = | M y |
      ! |               |   |     |   |     |
      ! |               |   |     |   |     |
      ! |               |             |     |
      ! This composition requires:
      !       square block decomposition for PDGETRF
      !       square matrix  M = N
      use mld_logger
      use module_kind_variables, only: kind_double
      !implicit none
      integer, intent(in) :: Mline, Ncolm, nbr_AA, nbc_AA, dimr_yy, dimc_yy
      real(kind_double), dimension(:,:), intent(in)  :: sca_AA
      real(kind_double), dimension(:,:), intent(inout)  :: sca_yy
      integer, dimension(:), intent(in)  :: desc_sca_AA
      integer, dimension(:), intent(inout)  :: desc_sca_yy
      real(kind_double), dimension(:,:), intent(inout)  :: w_params
      integer :: ipvt(Mline + nbr_AA)
      integer :: nrhs
      integer :: info, ii, lwork 
      real(kind_double), dimension(:), allocatable :: tau, work
      character(len=90) :: cinfo 
      _NAMECURRENT_("scalapack_lsystem_by_homeLU")

      _MLD_BEGIN_

      if (allocated(tau)) deallocate(tau) ; allocate(tau(Ncolm))

      if (Mline > Ncolm) then 
        ! QR factorization of A
        lwork = -1 
        if (allocated(work)) deallocate(work) ; allocate(work(1))
        call pdgeqrf(Mline, Ncolm, sca_AA, 1, 1, desc_sca_AA, tau, work, lwork, INFO)
        lwork = int(work(1)) + 2 
        if (allocated(work)) deallocate(work) ; allocate(work(lwork))
        call pdgeqrf(Mline, Ncolm, sca_AA, 1, 1, desc_sca_AA, tau, work, lwork, INFO)
        if (info < 0 ) then
          write(cinfo, '(i8)') info
          call log_warning("ML sca: error in "//NAMECURRENT//" for pdgeqrf info: "//trim(cinfo)) 
        end if 
  
        ! Multiply Q^T with y (results stored in y)
        if (allocated(work)) deallocate(work) ; allocate(work(1))
        lwork = -1
        call pdormqr('L', 'T', Mline, 1, Ncolm, sca_AA, 1, 1, desc_sca_AA, tau, &
                               sca_yy, 1, 1, desc_sca_yy, work, lwork, info)
        lwork = int(work(1)) + 2 
        if (allocated(work)) deallocate(work) ; allocate(work(lwork))
        call pdormqr('L', 'T', Mline, 1, Ncolm, sca_AA, 1, 1, desc_sca_AA, tau, &
                               sca_yy, 1, 1, desc_sca_yy, work, lwork, info)
        if (info < 0 ) then
          write(cinfo, '(i8)') info
          call log_warning("ML sca: error in "//NAMECURRENT//" for pdormqr info: "//trim(cinfo)) 
        end if
  
        ! Solve triangular system using the R factor of A and the transformed y
        call pdtrtrs('U', 'N', 'N', Ncolm, 1, sca_AA, 1, 1, desc_sca_AA, sca_yy, 1, 1, desc_sca_yy, info)
        if (info < 0 ) then
          write(cinfo, '(i8)') info
          call log_warning("ML sca: error in "//NAMECURRENT//" for pdtrtrs info: "//trim(cinfo)) 
        end if

        deallocate(tau, work) 

      end if 

      if (Mline < Ncolm) then 

        ! LQ factorization of A
        lwork = -1
        if (allocated(work)) deallocate(work) ; allocate(work(1))
        call pdgelqf(Mline, Ncolm, sca_AA, 1, 1, desc_sca_AA, tau, work, lwork, info)
        lwork = int(work(1)) + 2
        if (allocated(work)) deallocate(work) ; allocate(work(lwork))
        call pdgelqf(Mline, Ncolm, sca_AA, 1, 1, desc_sca_AA, tau, work, lwork, info)
        if (info < 0 ) then
          write(cinfo, '(i8)') info
          call log_warning("ML sca: error in "//NAMECURRENT//" for pdgelqf info: "//trim(cinfo)) 
        end if
  
        ! Solve triangular system using the L factor of A
        lwork = -1
        if (allocated(work)) deallocate(work) ; allocate(work(1))
        call pdtrtrs('L', 'N', 'N', Mline, 1, sca_AA, 1, 1, desc_sca_AA, sca_yy, 1, 1, desc_sca_yy, info)
        lwork = int(work(1)) + 2
        if (allocated(work)) deallocate(work) ; allocate(work(lwork))
        call pdtrtrs('L', 'N', 'N', Mline, 1, sca_AA, 1, 1, desc_sca_AA, sca_yy, 1, 1, desc_sca_yy, info)
        if (info < 0 ) then
          write(cinfo, '(i8)') info
          call log_warning("ML sca: error in "//NAMECURRENT//" for pdtrtrs info: "//trim(cinfo)) 
        end if
        ! Multiply Q with y to get the final solution
        lwork = -1
        if (allocated(work)) deallocate(work) ; allocate(work(1))
        call pdormlq('L', 'N', Ncolm, 1, Mline, sca_AA, 1, 1, desc_sca_AA, tau, sca_yy, 1, 1, desc_sca_AA, work, lwork, info)
        lwork = int(work(1)) + 2
        if (allocated(work)) deallocate(work) ; allocate(work(lwork))
        call pdormlq('L', 'N', Ncolm, 1, Mline, sca_AA, 1, 1, desc_sca_AA, tau, sca_yy, 1, 1, desc_sca_AA, work, lwork, info)
        if (info < 0 ) then
          write(cinfo, '(i8)') info
          call log_warning("ML sca: error in "//NAMECURRENT//" for pdormlq info: "//trim(cinfo)) 
        end if
  
        deallocate(tau, work) 

      end if

       

      if ( Mline == NColm) then 
        ! Call LU factorization routine
        ! PDGETRF computes an LU factorization of a general M-by-N distributed
        ! using partial pivoting with row interchanges.
        call pdgetrf(Mline, Ncolm,sca_AA,1,1,desc_sca_AA,ipvt,info)
        if (info .ne. 0) write(6,'("ML sca: error for pdgetrf info: ", i8)') info
  
        ! Call LU solver routine
        !PDGETRS solves a system of distributed linear equations with a general
        !N-by-N distributed matrix sub A using the LU factorization computed by PDGETRF.
        nrhs = dimc_yy
        call pdgetrs('N',Mline,nrhs,sca_AA,1,1,desc_sca_AA,ipvt,sca_yy,1,1,desc_sca_yy,info)
        if (info .ne. 0) write(6,'("ML sca: error for pdgetrs info: ", i8)') info
  
      end if 

      do ii = 1, size(w_params,1)
         call pdelget('A', ' ', w_params(ii,1), sca_yy, ii, 1, desc_sca_yy)
      end do

      _MLD_END_ 

   end subroutine scalapack_lsystem_by_homeLU

   subroutine scalapack_square_lsystem_by_cholesky(Mline, Ncolm, sca_AA, desc_sca_AA, nbr_AA, nbc_AA, &
      dimr_yy, dimc_yy, sca_yy, desc_sca_yy, w_params)
!  Solve the system sca_AA * w = ymat   MxN * Nx1 = Mx1
!  sca_AA is often sca_phi matrix  = Amat W Amat^T ( f(D) x f(D) matrix )
!       N= f(D)          1         1
! |               |             |     |
! |               |   |     |   |     |
! |               |   |     |   |     |
! | M   sca_A     | x | N w | = | M y |
! |               |   |     |   |     |
! |               |   |     |   |     |
! |               |             |     |
! This composition requires:
!       square block decomposition for PDPOSV
!       square matrix  M = N
      use module_kind_variables, only: kind_double
      use mld_logger
!implicit none
      integer, intent(in) :: Mline, Ncolm, nbr_AA, nbc_AA, dimr_yy, dimc_yy
      real(kind_double), dimension(:,:), intent(in)  :: sca_AA
      real(kind_double), dimension(:,:), intent(inout)  :: sca_yy
      integer, dimension(:), intent(in)  :: desc_sca_AA
      integer, dimension(:), intent(inout)  :: desc_sca_yy
      real(kind_double), dimension(:,:), intent(inout)  :: w_params
      !integer :: ipvt(Mline + nbr_AA)
      integer :: nrhs
      integer :: info, ii

      if (nbr_AA /= nbc_AA) then
         call log_critical("ML sca: this solution requires PDPOSV with square block decomposition i.e. row and column block factor should be equal ")
         call log_critical("ML sca: now row and column factor are set to: "// vtoa(nbr_AA)//"  "//vtoa(nbc_AA)//"  ")
         stop 'non square block factor in scalapack_lsystem_by_cholesky'
      end if
      if (Mline /= Ncolm) then
         write(6,'("ML sca: this solution requires PDPOSV with square matrix ")')
         write(6,'("ML sca: number of rows ans columns are set to: ", 2i8)') Mline, Ncolm
         stop 'non square matrix in scalapack_square_lsystem_by_cholesky'
      end if
! PDPOSV computes the solution to a real system of linear equations using Cholesky decomposition
!The Cholesky decomposition is used to factor sub( A ) as
!
!                     sub( A ) = U**T * U,  if UPLO = 'U', or
!
!                     sub( A ) = L * L**T,  if UPLO = 'L',
!
!  where U is an upper triangular matrix and L is a lower triangular
!  matrix.  The factored form of sub( A ) is then used to solve the
!  system of equations.
!SUBROUTINE pdposv( UPLO, N, NRHS, A, IA, JA, DESCA, B, IB, JB,
!  $                   DESCB, INFO )

      nrhs = dimc_yy
      call pdposv('L', Mline, nrhs, sca_AA, 1, 1, desc_sca_AA, sca_yy, 1, 1, desc_sca_yy, info)
      if (info .ne. 0) write(6,'("ML sca: error for pdposv info: ", i8)') info

      do ii = 1, size(w_params,1)
         call pdelget('A', ' ', w_params(ii,1), sca_yy, ii, 1, desc_sca_yy)
      end do

   end subroutine scalapack_square_lsystem_by_cholesky


   subroutine scalapack_lsystem_by_QR(Mline, Ncolm, sca_AA, desc_sca_AA, nbr_AA, nbc_AA, &
      dimr_yy, dimc_yy, sca_yy, desc_sca_yy, w_params)
!  Solve the system sca_AA * w = ymat   MxN * Nx1 = Mx1
!  sca_AA is often sca_phi or Amat_big matrix  = Amat W Amat^T ( f(D) x f(D) matrix )
!       N= f(D)          1         1
! |               |             |     |
! |               |   |     |   |     |
! |               |   |     |   |     |
! | M   sca_A     | x | N w | = | M y |
! |               |   |     |   |     |
! |               |   |     |   |     |
! |               |             |     |
! This composition requires:
!       square block decomposition for PDGETRF
!       square matrix  M = N
      use module_kind_variables, only: kind_double
!implicit none
      integer, intent(in) :: Mline, Ncolm, nbr_AA, nbc_AA, dimr_yy, dimc_yy
      real(kind_double), dimension(:,:), intent(in)  :: sca_AA
      real(kind_double), dimension(:,:), intent(inout)  :: sca_yy
      integer, dimension(:), intent(in)  :: desc_sca_AA
      integer, dimension(:), intent(inout)  :: desc_sca_yy
      real(kind_double), dimension(:,:), intent(inout)  :: w_params
      real(kind_double), dimension(:), allocatable :: work
      !integer :: ipvt(Mline + nbr_AA)
      integer :: nrhs, lwork
      integer :: info, ii

!if (Mline /= Ncolm) then
!write(6,'("ML sca: this solution requires PDGETRF / PDGETRS with square matrix ")')
!write(6,'("ML sca: number of rows ans columns are set to: ", 2i8)') Mline, Ncolm
!stop 'non square matrix in scalapack_square_lsystem_by_homeLU'
!end if


!PDGELS solves overdetermined or underdetermined real linear
!systems involving an M-by-N matrix sub( A ) = A(IA:IA+M-1,JA:JA+N-1),
!or its transpose, using a QR or LQ factorization of sub( A ).
! IMPORTANT: It is assumed that sub( A ) has full rank.

!SUBROUTINE pdgels( TRANS, M, N, NRHS, A, IA, JA, DESCA, B, IB, JB,
!  $                   DESCB, WORK, LWORK, INFO )
      nrhs = dimc_yy
      if (allocated(work)) deallocate(work) ; allocate(work(1))
      lwork = -1
      call pdgels('N', Mline, Ncolm, nrhs, sca_AA, 1, 1, desc_sca_AA, sca_yy, &
         1, 1, desc_sca_yy,  work, lwork, info)
      if (info < 0) then
         write(6,*) 'ML sca error: the dimension odf work failed in scalapack_lsystem_by_QR'
      end if

      lwork = int(work(1)) + 10
      if (allocated(work)) deallocate(work) ; allocate(work(lwork))

      call pdgels('N', Mline, Ncolm, nrhs, sca_AA, 1, 1, desc_sca_AA, sca_yy, &
         1, 1,  desc_sca_yy, work, lwork, info)
      if (info < 0) then
         write(6,*) 'ML sca error: the LLS problem failed  in scalapack_lsystem_by_QR'
      end if

      do ii = 1, size(w_params,1)
         call pdelget('A', ' ', w_params(ii,1), sca_yy, ii, 1, desc_sca_yy)
      end do

   end subroutine scalapack_lsystem_by_QR


   subroutine scalapack_lsystem_by_SVD(Mline, Ncolm, sca_AA, desc_sca_AA, nbr_AA, nbc_AA, &
      dimr_yy, dimc_yy, sca_yy, desc_sca_yy, w_params, svd_rcond, rank_sca_AA, context)
      !  Solve the system sca_AA * w = ymat   MxN * Nx1 = Mx1
      !  sca_AA is often sca_phi or Amat_big matrix  = Amat W Amat^T ( f(D) x f(D) matrix )
      !       N= f(D)          1         1
      ! |               |             |     |
      ! |               |   |     |   |     |
      ! |               |   |     |   |     |
      ! | M   sca_A     | x | N w | = | M y |
      ! |               |   |     |   |     |
      ! |               |   |     |   |     |
      ! |               |             |     |
      ! This decomposition requires:
      !       iam, myrow, mycol, nprow, npcol
      !       square matrix  M = N
      use ml_in_ndm_module, only: debug
      use module_kind_variables, only: kind_double
      use module_ml_scalapack, only : iam, myrow, mycol, nprow, npcol
      use module_scalapack_tools, only: allocation_scalapack_matrix, up_scale_scalpack_systemAy
      use module_optimization, only: optimize_weights_chem, optimize_weights_db
      use mld_logger, only: log_info 
      !implicit none
      integer, intent(in) :: Mline, Ncolm, nbr_AA, nbc_AA, dimr_yy, dimc_yy
      integer, intent(in) :: context
      integer, intent(inout) :: rank_sca_AA
      real(kind_double), intent(inout) :: svd_rcond
      real(kind_double), dimension(:,:), intent(inout)  :: sca_AA
      real(kind_double), dimension(:,:), intent(inout)  :: sca_yy
      integer, dimension(:), intent(inout)  :: desc_sca_AA
      integer, dimension(:), intent(inout)  :: desc_sca_yy
      real(kind_double), dimension(:,:), intent(inout)  :: w_params
      real(kind_double), dimension(:), allocatable :: work
      integer :: lwork
      integer :: info, ii, rsize, icount
      real(kind_double) :: zero=0.d0, one=1.d0, alpha

      ! sca_UU utilities ...
      integer :: nbr_UU, nbc_UU, dimr_sca_UU, dimc_sca_UU
      integer :: l_dimr_sca_UU, l_dimc_sca_UU
      integer, dimension(:), allocatable :: desc_sca_UU
      real(kind_double), dimension(:,:), allocatable  :: sca_UU

      ! sca_VT utilities ...
      integer :: nbr_VT, nbc_VT, dimr_sca_VT, dimc_sca_VT
      integer :: l_dimr_sca_VT, l_dimc_sca_VT
      integer, dimension(:) , allocatable:: desc_sca_VT
      real(kind_double), dimension(:,:), allocatable :: sca_VT

      ! sca_vec utilities ...
      integer :: nbr_vec, nbc_vec, dimr_sca_vec, dimc_sca_vec
      integer :: l_dimr_sca_vec, l_dimc_sca_vec
      integer, dimension(:) , allocatable:: desc_sca_vec
      real(kind_double), dimension(:,:), allocatable :: sca_vec

      ! sca_sol utilities ...
      integer :: nbr_sol, nbc_sol, dimr_sca_sol, dimc_sca_sol
      integer :: l_dimr_sca_sol, l_dimc_sca_sol
      integer, dimension(:) , allocatable:: desc_sca_sol
      real(kind_double), dimension(:,:), allocatable :: sca_sol

      real(kind_double), dimension(:), allocatable :: Sigma_SVD
      character(len=100) :: chlog

      ! Local gather + dgemv workaround for MKL PB_CVMnpq bug:
      ! MKL pdgemv/pdgemm crash with integer divide-by-zero when
      ! block sizes are large relative to vector dimensions.
      ! We gather the input vector, do local dgemv on each rank's
      ! tile, then reduce-scatter the result.
      real(kind_double), dimension(:), allocatable :: yy_full, vec_full, vec_local, sol_full, sol_local
      integer :: il, jl, gi, gj, INDXL2G
      real(kind_double), external :: ddot


      !$! call up_scale_scalpack_systemAy(sca_AA, desc_sca_AA, Mline, Ncolm, &
      !$!                                 sca_yy, desc_sca_yy, dimr_yy, dimc_yy, &
      !$!                                 anrm, bnrm, iascl, ibscl)
      !$! 
      !$! if (iam == 0) write(6,*)  'Amat normalized by', anrm 
      !$! if (iam == 0) write(6,*)  'ymat normalized by', bnrm                                
      !  PDGESVD computes the singular value decomposition (SVD) of an
      !  M-by-N matrix A, optionally computing the left and/or right
      !  singular vectors. The SVD is written as
      !
      !       A = U * SIGMA * transpose(V)
      !      M x N = M x M * S * N x N -> this is the full SVD decomposition
      !      here we will choose the decomposition with rank max
      !      r = min(M, N) (the rank cannot be larger).
      !      M x N = M x r * S * r x N

      rsize = min ( Mline, Ncolm )
      ! Allocate scalapack UU:

      nbr_UU =  nbr_AA
      nbc_UU =  min (nbr_AA, nbc_AA)
      dimr_sca_UU = Mline
      dimc_sca_UU = rsize
      call allocation_scalapack_matrix (" sca_UU ", sca_UU, desc_sca_UU, dimr_sca_UU, dimc_sca_UU, &
         l_dimr_sca_UU, l_dimc_sca_UU, nbr_UU, nbc_UU, &
         myrow, mycol, nprow, npcol, context, iam  )

      ! Allocate scalapack VT:

      nbr_VT =  min(nbr_AA, nbc_AA)
      nbc_VT =  nbc_AA
      dimr_sca_VT = rsize
      dimc_sca_VT = Ncolm
      call allocation_scalapack_matrix (" sca_VT ", sca_VT, desc_sca_VT, dimr_sca_VT, dimc_sca_VT, &
         l_dimr_sca_VT, l_dimc_sca_VT, nbr_VT, nbc_VT, &
         myrow, mycol, nprow, npcol, context, iam  )

      if (allocated (Sigma_SVD)) deallocate(Sigma_SVD) ; allocate(Sigma_SVD (rsize))

      !SUBROUTINE PDGESVD(JOBU,JOBVT,M,N,A,IA,JA,DESCA,S,U,IU,JU,DESCU,
      !  +                   VT,IVT,JVT,DESCVT,WORK,LWORK,INFO)
      !nrhs = dimc_yy
      if (allocated(work)) deallocate(work) ; allocate(work(1))
      lwork = -1
      call pdgesvd('V', 'V', Mline, Ncolm, sca_AA, 1, 1, desc_sca_AA, &
         Sigma_SVD, sca_UU, 1, 1, desc_sca_UU,  sca_VT,  1, 1, desc_sca_VT, work, lwork, info)
      if (info /=  0) then
         if (iam == 0) write(6,*) 'ML sca error: pdgesvd  failed  allocation in scalapack_lsystem_by_SVD'
         call messages_pdgsvd(info)
      end if

      call blacs_barrier(context, 'A')
      !write(6,*) iam, 'before SVD ............'

      lwork = int(work(1)) + 10
      if (allocated(work)) deallocate(work) ; allocate(work(lwork))
      
      call pdgesvd('V', 'V', Mline, Ncolm,  sca_AA, 1, 1, desc_sca_AA, &
         Sigma_SVD, sca_UU, 1, 1, desc_sca_UU,  sca_VT,  1, 1, desc_sca_VT, work, lwork, info)
      if (info /= 0) then
         if (iam == 0) write(6,*) 'ML sca error: the LLS problem failed  in scalapack_lsystem_by_SVD'
         call messages_pdgsvd(info)
      end if
      if (debug) then
         if (iam == 0) write (6, '("ML: after SVD...............................")')
      end if

      call blacs_barrier(context, 'A')
      !write(6,*) iam, 'after SVD ............', svd_rcond
      if (svd_rcond < 0) then
         svd_rcond = 100.d0*epsilon(1.d0)
         write(chlog, '(es20.10)') svd_rcond
         call log_info("old svd_rcond was set to :"//trim(chlog))
         !svd_rcond =   max(Mline, Ncolm) * maxval(Sigma_SVD)  *epsilon(1.d0)
         svd_rcond =   maxval(Sigma_SVD)  *epsilon(0.d0)


         write(chlog, '(es20.10)') svd_rcond
         call log_info("new svd_rcond is set to  :"//trim(chlog))
         call log_info("ML sca: Remember that you also have the option to manually set the svd_rcond parameter to a positive value.")
         call log_info("If the results are crazy use larger values of svd_rcond") 
         call log_info("However, increasing svd_rcond can enhance robustness, but it may reduce the accuracy of the fit.") 
         if (minval(Sigma_SVD) /= 0.d0) then
            write(chlog, '(es20.10)') maxval(Sigma_SVD)/minval(Sigma_SVD)
            call log_info("ML Sca: The condition number of the matrix is : "//trim(chlog))
         else 
            call log_info("ML Sca: The condition number of the matrix is : infinite :) ")
         end if   
         
      end if

      ! get the rank ...
      icount = 0
      do ii = 1, size(Sigma_SVD, 1)
         if (dabs(Sigma_SVD(ii)) >= svd_rcond) then
            icount = icount + 1
         end if
      end do
      rank_sca_AA = icount

      if (.not.(optimize_weights_chem.or.optimize_weights_db)) then
         if (iam == 0) then
            write (6, '("ML: dgelsd SVD inversion info Mline Ncolm lwork RCOND RANK ....:", i4, i9, i9, i9, d20.10,i6)') &
               info, Mline, Ncolm, lwork, svd_rcond, rank_sca_AA
         end if
      end if


      if (rank_sca_AA == 0) then
         if (iam==0) write(6,'("ML error: serious problem in SVD: rank 0 of matrix")')
         stop
      end if

      !if (iam == 0) then
      !  do ii = 1, size(Sigma_SVD, 1)
      !  write(*,*)  ii, Sigma_SVD(ii), svd_rcond, rank_sca_AA
      !  end do
      !end if


      ! Get the SVD solutions:
      ! vec =  r x 1  = U.T (r x M) * ymat (M x 1)
      !
      ! ============================================================
      ! Workaround for MKL PB_CVMnpq integer divide-by-zero bug:
      ! MKL's pdgemv/pdgemm crash when block sizes are large
      ! relative to vector dimensions.  Instead we:
      !   1. Gather the input vector yy into a full local copy
      !   2. Each rank does local dgemv on its tile of sca_UU
      !   3. Assemble the partial results into vec_full via
      !      DGSUM2D (BLACS reduction along rows/columns)
      ! ============================================================

      ! --- Step 1: Gather yy_full(Mline) from distributed sca_yy ---
      allocate(yy_full(Mline))
      yy_full = zero
      do ii = 1, Mline
         call pdelget('A', ' ', yy_full(ii), sca_yy, ii, 1, desc_sca_yy)
      end do

      ! --- Step 2: Local dgemv on each rank's tile of sca_UU ---
      ! sca_UU is Mline x rsize, distributed with nbr_UU x nbc_UU blocks.
      ! Each rank owns l_dimr_sca_UU rows and l_dimc_sca_UU columns.
      ! We want vec = U^T * yy, i.e. vec(j) = sum_i UU(i,j) * yy(i)
      !
      ! For each local column jl, global column gj = INDXL2G(jl, nbc_UU, mycol, 0, npcol):
      !   vec_local(gj) += sum over local rows il:  sca_UU(il, jl) * yy_full(gi)
      !     where gi = INDXL2G(il, nbr_UU, myrow, 0, nprow)
      !
      ! This is a local dgemv: vec_local = sca_UU^T * yy_local_rows

      allocate(vec_full(rsize), vec_local(rsize))
      vec_full  = zero
      vec_local = zero

      if (allocated(work)) deallocate(work)

      if (l_dimr_sca_UU > 0 .and. l_dimc_sca_UU > 0) then
         ! Extract the local rows of yy that this rank owns
         allocate(work(l_dimr_sca_UU))
         do il = 1, l_dimr_sca_UU
            gi = INDXL2G(il, nbr_UU, myrow, 0, nprow)
            work(il) = yy_full(gi)
         end do
         ! Local dgemv: vec_part(l_dimc_sca_UU) = sca_UU^T * work
         ! But the result is indexed by local column → must map to global
         do jl = 1, l_dimc_sca_UU
            gj = INDXL2G(jl, nbc_UU, mycol, 0, npcol)
            ! dot product of column jl of sca_UU with work
            vec_local(gj) = vec_local(gj) + &
               ddot(l_dimr_sca_UU, sca_UU(1, jl), 1, work, 1)
         end do
         deallocate(work)
      end if

      ! --- Step 3: Reduce vec_local across all processes ---
      ! Use BLACS DGSUM2D to sum across the full grid
      call dgsum2d(context, 'A', ' ', rsize, 1, vec_local, rsize, -1, -1)
      vec_full(:) = vec_local(:)

      deallocate(yy_full, vec_local)

      ! --- Apply Sigma_SVD^{-1} to vec (truncated to rank) ---
      do ii = 1, rank_sca_AA
         vec_full(ii) = vec_full(ii) / Sigma_SVD(ii)
      end do
      if (rsize > rank_sca_AA) then
         do ii = rank_sca_AA + 1, rsize
            vec_full(ii) = zero
         end do
      end if

      ! ============================================================
      ! xsol = VT^T * vec   (Ncolm x 1 = (Ncolm x rsize) * (rsize x 1))
      ! Same workaround: local dgemv on each rank's tile of sca_VT.
      ! sca_VT is rsize x Ncolm, distributed with nbr_VT x nbc_VT blocks.
      ! We want sol = VT^T * vec, i.e. sol(j) = sum_i VT(i,j) * vec(i)
      ! With trans='T', the output has dimension Ncolm.
      ! ============================================================

      allocate(sol_full(Ncolm), sol_local(Ncolm))
      sol_full  = zero
      sol_local = zero

      if (l_dimr_sca_VT > 0 .and. l_dimc_sca_VT > 0) then
         ! Extract the local rows of vec that this rank owns (rows of VT)
         if (allocated(work)) deallocate(work)
         allocate(work(l_dimr_sca_VT))
         do il = 1, l_dimr_sca_VT
            gi = INDXL2G(il, nbr_VT, myrow, 0, nprow)
            work(il) = vec_full(gi)
         end do
         ! For each local column jl of VT → global column gj of sol:
         !   sol_local(gj) += dot(VT(:,jl)_local, vec_local_rows)
         do jl = 1, l_dimc_sca_VT
            gj = INDXL2G(jl, nbc_VT, mycol, 0, npcol)
            sol_local(gj) = sol_local(gj) + &
               ddot(l_dimr_sca_VT, sca_VT(1, jl), 1, work, 1)
         end do
         deallocate(work)
      end if

      ! Reduce sol_local across all processes
      call dgsum2d(context, 'A', ' ', Ncolm, 1, sol_local, Ncolm, -1, -1)
      sol_full(:) = sol_local(:)

      deallocate(vec_full, sol_local)

      ! Copy solution to w_params
      do ii = 1, size(w_params, 1)
         w_params(ii, 1) = sol_full(ii)
      end do

      deallocate(sol_full)

      if (allocated(sca_UU))  deallocate (sca_UU)
      if (allocated(desc_sca_UU))  deallocate (desc_sca_UU)

      if (allocated(sca_VT))  deallocate (sca_VT)
      if (allocated(desc_sca_VT))  deallocate (desc_sca_VT)

   end subroutine scalapack_lsystem_by_SVD

   subroutine messages_pdgsvd(info)
      use mld_logger
      implicit none 
      integer, intent(in) :: info
      

      if (info < 0 ) then  
         call log_warning("ML sca: error in scalapack_lsystem_by_SVD for pdgesvd info: " // vtoa(info)) 
         call log_warning("ML sca:  if INFO = -i, the i-th argument had an illegal value " // vtoa(info))
      end if 

      if (info > 0 ) then 
         call log_warning("ML sca:  DBDSQR did not converge. If INFO = MIN(M,N) + 1, then " // &
                   " PDGESVD has detected heterogeneity by finding that eigenvalues were not "//  &
                   " identical across the process grid. In this case, the " // &
                  " accuracy of the results from PDGESVD cannot be guaranteed. " // vtoa(info))
      end if     
   end subroutine messages_pdgsvd


   subroutine scalapack_SVD_decomposition(sca_AA, desc_sca_AA, nbr_AA, nbc_AA, dimr_sca_AA, dimc_sca_AA, l_dimr_sca_AA, l_dimc_sca_AA, &
      sca_UU, desc_sca_UU, nbr_UU, nbc_UU, dimr_sca_UU, dimc_sca_UU, l_dimr_sca_UU, l_dimc_sca_UU, &
      sca_VT, desc_sca_VT, nbr_VT, nbc_VT, dimr_sca_VT, dimc_sca_VT, l_dimr_sca_VT, l_dimc_sca_VT, &
      svd_rcond, rank_sca_AA)
      !  svd decomposition of sca_AA   MxN
      ! This decomposition requires:
      !       iam, myrow, mycol, nprow, npcol, context
      use mld_logger
      use module_kind_variables, only: kind_double
      use module_ml_scalapack, only : debug_scalapack, iam, myrow, mycol, nprow, npcol, context
      use module_scalapack_tools, only: allocation_scalapack_matrix
      !implicit none
      real(kind_double), dimension(:,:), intent(in)  :: sca_AA
      integer, dimension(:), intent(in)  :: desc_sca_AA
      integer, intent(in) ::  nbr_AA, nbc_AA, dimr_sca_AA, dimc_sca_AA
      integer, intent(in) ::  l_dimr_sca_AA, l_dimc_sca_AA
      ! sca_UU utilities ...
      real(kind_double), dimension(:,:), allocatable, intent(inout)   :: sca_UU
      integer, dimension(:), allocatable, intent(inout) :: desc_sca_UU
      integer, intent(inout) :: nbr_UU, nbc_UU, dimr_sca_UU, dimc_sca_UU
      integer, intent(inout) :: l_dimr_sca_UU, l_dimc_sca_UU
      ! sca_VT utilities ...
      real(kind_double), dimension(:,:), allocatable, intent(inout) :: sca_VT
      integer, intent(inout) :: nbr_VT, nbc_VT, dimr_sca_VT, dimc_sca_VT
      integer, intent(inout) :: l_dimr_sca_VT, l_dimc_sca_VT
      integer, dimension(:) , allocatable, intent(inout) :: desc_sca_VT

      integer, intent(inout) :: rank_sca_AA
      real(kind_double), intent(inout) :: svd_rcond

      !local variable, vectors, matrix  ...
      real(kind_double), dimension(:,:), allocatable  :: sca_AA_copy
      integer, dimension(:), allocatable  :: desc_sca_AA_copy
      integer ::  nbr_AA_copy, nbc_AA_copy, dimr_sca_AA_copy, dimc_sca_AA_copy
      integer ::  l_dimr_sca_AA_copy, l_dimc_sca_AA_copy


      real(kind_double), dimension(:), allocatable :: work
      integer :: lwork
      integer :: info, ii, jj, rsize, icount

      real(kind_double), dimension(:), allocatable :: Sigma_SVD
      integer :: Mline, Ncolm
      character(len=100) :: chlog


      !  PDGESVD computes the singular value decomposition (SVD) of an
      !  M-by-N matrix A, optionally computing the left and/or right
      !  singular vectors. The SVD is written as
      !
      !       A = U * SIGMA * transpose(V)
      !      M x N = M x M * S * N x N -> this is the full SVD decomposition
      !      here we will choose the decomposition with rank max
      !      r = min(M, N) (the rank cannot be larger).
      !      M x N = M x r * S * r x N
      Mline = dimr_sca_AA
      Ncolm = dimc_sca_AA
      rsize = min (Mline, Ncolm)


      ! Allocate scalapack UU:

      nbr_UU =  nbr_AA
      nbc_UU =  min (nbr_AA, nbc_AA)
      dimr_sca_UU = Mline
      dimc_sca_UU = rsize
      call allocation_scalapack_matrix (" sca_UU ", sca_UU, desc_sca_UU, dimr_sca_UU, dimc_sca_UU, &
         l_dimr_sca_UU, l_dimc_sca_UU, nbr_UU, nbc_UU, &
         myrow, mycol, nprow, npcol, context, iam  )

      ! Allocate scalapack VT:

      nbr_VT =  min(nbr_AA, nbc_AA)
      nbc_VT =  nbc_AA
      dimr_sca_VT = rsize
      dimc_sca_VT = Ncolm
      call allocation_scalapack_matrix (" sca_VT ", sca_VT, desc_sca_VT, dimr_sca_VT, dimc_sca_VT, &
         l_dimr_sca_VT, l_dimc_sca_VT, nbr_VT, nbc_VT, &
         myrow, mycol, nprow, npcol, context, iam  )

      nbr_AA_copy =  nbr_AA
      nbc_AA_copy =  nbc_AA
      dimr_sca_AA_copy = dimr_sca_AA
      dimc_sca_AA_copy = dimc_sca_AA
      call allocation_scalapack_matrix (" sca_AA_copy ", sca_AA_copy, desc_sca_AA_copy, dimr_sca_AA_copy, dimc_sca_AA_copy, &
         l_dimr_sca_AA_copy, l_dimc_sca_AA_copy, nbr_AA_copy, nbc_AA_copy, &
         myrow, mycol, nprow, npcol, context, iam  )
      do jj = 1, l_dimc_sca_AA
         do ii = 1, l_dimr_sca_AA
            sca_AA_copy(ii,jj) = sca_AA(ii,jj)
         end do
      end do


      if (allocated (Sigma_SVD)) deallocate(Sigma_SVD) ; allocate(Sigma_SVD (rsize))

      !SUBROUTINE PDGESVD(JOBU,JOBVT,M,N,A,IA,JA,DESCA,S,U,IU,JU,DESCU,
      !  +                   VT,IVT,JVT,DESCVT,WORK,LWORK,INFO)
      if (allocated(work)) deallocate(work) ; allocate(work(1))
      lwork = -1
      call pdgesvd('V', 'V', Mline, Ncolm, sca_AA_copy, 1, 1, desc_sca_AA_copy, &
         Sigma_SVD, sca_UU, 1, 1, desc_sca_UU,  sca_VT,  1, 1, desc_sca_VT, work, lwork, info)
      if (info < 0) then
         write(6,*) 'ML sca error: the dimension odf work failed in scalapack_lsystem_by_SVD'
      end if

      call blacs_barrier(context, 'A')

      lwork = int(work(1)) + 10
      if (allocated(work)) deallocate(work) ; allocate(work(lwork))
      call pdgesvd('V', 'V', Mline, Ncolm,  sca_AA_copy, 1, 1, desc_sca_AA_copy, &
         Sigma_SVD, sca_UU, 1, 1, desc_sca_UU,  sca_VT,  1, 1, desc_sca_VT, work, lwork, info)
      if (info < 0) then
         write(6,*) 'ML sca error: the SVD problem failed  in scalapack_SVD_decomposition'
      end if
      if (allocated(sca_AA_copy)) deallocate(sca_AA_copy)
      if (allocated(desc_sca_AA_copy)) deallocate(desc_sca_AA_copy)

      if ((iam == 0).and.debug_scalapack) write (6, '("ML: after SVD...............................")')

      call blacs_barrier(context, 'A')
      if (svd_rcond < 0) then
         ! this is the default value
         svd_rcond = 100.d0*epsilon(1.d0)
         ! this is the default value Lapack ... 
         write(chlog, '(es20.10)') svd_rcond
         call log_info("old svd_rcond was set to :"//trim(chlog))
         !svd_rcond =   max(Mline, Ncolm) * maxval(Sigma_SVD)  *epsilon(1.d0)
         svd_rcond =   maxval(Sigma_SVD)  *epsilon(0.d0)
         write(chlog, '(es20.10)') svd_rcond
         call log_info("new svd_rcond is set to  :"//trim(chlog))
         call log_info("ML sca: Remember that you also have the option to manually set the svd_rcond parameter to a positive value.")
         call log_info("If the results are crazy use larger values of svd_rcond") 
         call log_info("However, increasing svd_rcond can enhance robustness, but it may reduce the accuracy of the fit.") 
      end if
      ! get the rank ...
      icount = 0
      do ii = 1, size(Sigma_SVD, 1)
         if (dabs(Sigma_SVD(ii)) >= svd_rcond) then
            icount = icount + 1
         end if
      end do
      rank_sca_AA = icount

      if ((iam == 0).and.debug_scalapack) then
         write (6, '("ML: pdgesvd SVD decomposition of phia info Mline Ncolm lwork RCOND RANK ....:", i4, i9, i9, i9, d20.10,i6)') &
            info, Mline, Ncolm, lwork, svd_rcond, rank_sca_AA
      end if

      if (rank_sca_AA == 0) then
         if (iam==0) write(6,'("ML error: serious problem in CUR SVD: rank 0 of matrix")')
         stop
      end if

      call blacs_barrier(context, 'A')


   end subroutine scalapack_SVD_decomposition


   subroutine scalapack_pseudo_inverse (sca_AA, desc_sca_AA, &
      nbr_AA, nbc_AA, dimr_sca_AA, dimc_sca_AA, &
      l_dimr_sca_AA, l_dimc_sca_AA, &
      rank_sca_AA, &
      sca_pinvAA, desc_sca_pinvAA, &
      nbr_pinvAA, nbc_pinvAA, dimr_sca_pinvAA, dimc_sca_pinvAA, &
      l_dimr_sca_pinvAA, l_dimc_sca_pinvAA)

      use module_kind_variables, only: kind_double
      use ml_in_ndm_module, only: svd_rcond
      use module_scalapack_tools, only: allocation_scalapack_matrix
      use module_ml_scalapack, only: debug_scalapack, myrow, mycol, nprow, npcol, context, iam
      use module_scalapack_interfaces, only: my_pdgemr2d
      use mld_logger

      implicit none

      real(kind_double), dimension(:, :), allocatable, intent(in) :: sca_AA
      integer, dimension(:), allocatable, intent(in)  :: desc_sca_AA
      integer, intent(in)   :: nbr_AA, nbc_AA, dimr_sca_AA, dimc_sca_AA
      integer, intent(in)   :: l_dimr_sca_AA, l_dimc_sca_AA

      integer, intent(out)  :: rank_sca_AA

      real(kind_double), dimension(:, :), allocatable, intent(out) :: sca_pinvAA
      integer, dimension(:), allocatable, intent(out)  :: desc_sca_pinvAA
      integer, intent(out)   :: nbr_pinvAA, nbc_pinvAA, dimr_sca_pinvAA, dimc_sca_pinvAA
      integer, intent(out)   :: l_dimr_sca_pinvAA, l_dimc_sca_pinvAA

      real(kind_double), external :: pdlange
      !local variables ...
      real(kind_double), dimension(:,:), allocatable  :: sca_AA_copy
      integer, dimension(:), allocatable :: desc_sca_AA_copy
      integer  :: nbr_AA_copy, nbc_AA_copy, dimr_sca_AA_copy, dimc_sca_AA_copy
      integer  :: l_dimr_sca_AA_copy, l_dimc_sca_AA_copy


      real(kind_double), dimension(:,:), allocatable  :: sca_UU
      integer, dimension(:), allocatable :: desc_sca_UU
      integer  :: nbr_UU, nbc_UU, dimr_sca_UU, dimc_sca_UU
      integer  :: l_dimr_sca_UU, l_dimc_sca_UU

      real(kind_double), dimension(:,:), allocatable  :: sca_UUk
      integer, dimension(:), allocatable :: desc_sca_UUk
      integer  :: nbr_UUk, nbc_UUk, dimr_sca_UUk, dimc_sca_UUk
      integer  :: l_dimr_sca_UUk, l_dimc_sca_UUk

      ! sca_VT utilities ...
      real(kind_double), dimension(:,:), allocatable  :: sca_VT
      integer  :: nbr_VT, nbc_VT, dimr_sca_VT, dimc_sca_VT
      integer  :: l_dimr_sca_VT, l_dimc_sca_VT
      integer, dimension(:) , allocatable  :: desc_sca_VT

      real(kind_double), dimension(:,:), allocatable  :: sca_VTk
      integer  :: nbr_VTk, nbc_VTk, dimr_sca_VTk, dimc_sca_VTk
      integer  :: l_dimr_sca_VTk, l_dimc_sca_VTk
      integer, dimension(:) , allocatable  :: desc_sca_VTk

      real(kind_double), dimension(:), allocatable :: Sigma_SVD

      real(kind_double), dimension(:,:), allocatable  :: sca_Sigmak
      integer  :: nbr_Sigmak, nbc_Sigmak, dimr_sca_Sigmak, dimc_sca_Sigmak
      integer  :: l_dimr_sca_Sigmak, l_dimc_sca_Sigmak
      integer, dimension(:) , allocatable  :: desc_sca_Sigmak

      real(kind_double), dimension(:,:), allocatable  :: sca_tmp
      integer  :: nbr_tmp, nbc_tmp, dimr_sca_tmp, dimc_sca_tmp
      integer  :: l_dimr_sca_tmp, l_dimc_sca_tmp
      integer, dimension(:) , allocatable  :: desc_sca_tmp


      real(kind_double), dimension(:), allocatable :: work
      integer :: lwork
      integer :: Mline, Ncolm, rsize, info, ii, jj, icount
      real(kind_double), parameter:: zero=0.d0, one=1.d0
      real(kind_double) :: svd_rcond_local
      

      _NAMECURRENT_("scalapack_pseudo_inverse")

      _MLD_BEGIN_


      Mline = dimr_sca_AA
      Ncolm = dimc_sca_AA
      rsize = min (Mline, Ncolm)

      ! Allocate scalapack UU:

      nbr_UU =  nbr_AA
      nbc_UU =  min (nbr_AA, nbc_AA)
      dimr_sca_UU = Mline
      dimc_sca_UU = rsize
      call allocation_scalapack_matrix (" sca_UU in pi ", sca_UU, desc_sca_UU, dimr_sca_UU, dimc_sca_UU, &
         l_dimr_sca_UU, l_dimc_sca_UU, nbr_UU, nbc_UU, &
         myrow, mycol, nprow, npcol, context, iam  )

      ! Allocate scalapack VT:

      nbr_VT =  min(nbr_AA, nbc_AA)
      nbc_VT =  nbc_AA
      dimr_sca_VT = rsize
      dimc_sca_VT = Ncolm
      call allocation_scalapack_matrix (" sca_VT in pi ", sca_VT, desc_sca_VT, dimr_sca_VT, dimc_sca_VT, &
         l_dimr_sca_VT, l_dimc_sca_VT, nbr_VT, nbc_VT, &
         myrow, mycol, nprow, npcol, context, iam  )

      nbr_AA_copy = nbr_AA
      nbc_AA_copy = nbc_AA
      dimr_sca_AA_copy = dimr_sca_AA
      dimc_sca_AA_copy = dimc_sca_AA
      call allocation_scalapack_matrix (" sca_AA_copy in pi ", sca_AA_copy, desc_sca_AA_copy, dimr_sca_AA_copy, dimc_sca_AA_copy, &
         l_dimr_sca_AA_copy, l_dimc_sca_AA_copy, nbr_AA_copy, nbc_AA_copy, &
         myrow, mycol, nprow, npcol, context, iam  )

      do jj = 1, l_dimc_sca_AA
         do ii = 1, l_dimr_sca_AA
            sca_AA_copy(ii,jj) = sca_AA(ii,jj)
         end do
      end do

      if (allocated (Sigma_SVD)) deallocate(Sigma_SVD) ; allocate(Sigma_SVD (rsize))


      if (allocated(work)) deallocate(work) ; allocate(work(1))
      lwork = -1
      call pdgesvd('V', 'V', Mline, Ncolm, sca_AA_copy, 1, 1, desc_sca_AA_copy, &
         Sigma_SVD, sca_UU, 1, 1, desc_sca_UU,  sca_VT,  1, 1, desc_sca_VT, work, lwork, info)
      if (info < 0) then
         call log_critical("ML error: the dimension of work failed in "//NAMECURRENT)
      end if

      call blacs_barrier(context, 'A')

      lwork = int(work(1)) + 10
      if (allocated(work)) deallocate(work) ; allocate(work(lwork))
      call pdgesvd('V', 'V', Mline, Ncolm,  sca_AA_copy, 1, 1, desc_sca_AA_copy, &
         Sigma_SVD, sca_UU, 1, 1, desc_sca_UU,  sca_VT,  1, 1, desc_sca_VT, work, lwork, info)

      if (info < 0) then
         call log_critical("ML error: the SVD problem failed  in "//NAMECURRENT)
      end if
      if (allocated(sca_AA_copy)) deallocate(sca_AA_copy)
      if (allocated(desc_sca_AA_copy)) deallocate(desc_sca_AA_copy)


      !if (iam == 0) write (6, '("ML: after SVD...............................")')
      if (debug_scalapack)  call log_info("ML: calapack_pseudo_inverse -> after SVD...............................")

      call blacs_barrier(context, 'A')

      if (svd_rcond < 0) then
         svd_rcond_local = 100.d0*epsilon(1.d0)*dabs(Sigma_SVD(1))
      else 
         svd_rcond_local = maxval(Sigma_SVD)  * svd_rcond 
      end if
      ! get the rank ...
      icount = 0
      do ii = 1, size(Sigma_SVD, 1)
         if (Sigma_SVD(ii) >= svd_rcond_local) then
            icount = icount + 1
         end if
      end do
      rank_sca_AA = icount

      if ((iam == 0).and.debug_scalapack) then
         write (6, '("ML: pdgesvd SVD decomposition of phia info Mline Ncolm lwork RCOND RANK ....:", i4, i9, i9, i9, d20.10,i6)') &
            info, Mline, Ncolm, lwork, svd_rcond, rank_sca_AA
      end if

      if (rank_sca_AA == 0) then
         call log_critical("ML error:  serious problem in SVD: rank 0 of matrix in "//NAMECURRENT)
         stop
      end if

      ! Allocate scalapack UUk - cropped UU until the rank:

      nbr_UUk =  nbr_AA
      nbc_UUk =  min (nbr_AA, nbc_AA)
      dimr_sca_UUk = Mline
      dimc_sca_UUk = rank_sca_AA
      call allocation_scalapack_matrix (" sca_UUk in pinv ", sca_UUk, desc_sca_UUk, dimr_sca_UUk, dimc_sca_UUk, &
         l_dimr_sca_UUk, l_dimc_sca_UUk, nbr_UUk, nbc_UUk, &
         myrow, mycol, nprow, npcol, context, iam  )

      !call pdgemr2d(m, n, a, ia, ja, desca, b, ib, jb, descb, ictxt)
      call my_pdgemr2d(dimr_sca_UUk, dimc_sca_UUk, sca_UU, 1, 1, desc_sca_UU, sca_UUk, 1, 1, desc_sca_UUk, context)
      if (allocated(sca_UU)) deallocate(sca_UU)
      if (allocated(desc_sca_UU)) deallocate(desc_sca_UU)

      if (debug_scalapack) call log_info("scalapack_pseudo_inverse ->  first pdgemr2d  completed")
      ! Allocate scalapack Sigmak - cropped version until the rank:

      nbr_Sigmak = min(nbr_AA, nbc_AA)
      nbc_Sigmak = min(nbr_AA, nbc_AA)
      dimr_sca_Sigmak = rank_sca_AA
      dimc_sca_Sigmak = rank_sca_AA
      call allocation_scalapack_matrix (" sca_Sigmak in pinv ", sca_Sigmak, desc_sca_Sigmak, dimr_sca_Sigmak, dimc_sca_Sigmak, &
         l_dimr_sca_Sigmak, l_dimc_sca_Sigmak, nbr_Sigmak, nbc_Sigmak, &
         myrow, mycol, nprow, npcol, context, iam  )
      do jj = 1, l_dimc_sca_Sigmak
         do ii = 1, l_dimr_sca_Sigmak
            sca_Sigmak(ii, jj) = 0.d0
         end do
      end do

      do ii = 1, dimr_sca_Sigmak
         call pdelset(sca_Sigmak, ii, ii, desc_sca_Sigmak, 1.d0/Sigma_SVD(ii))
      end do

      
      if (debug_scalapack) then 
        call log_info("scalapack_pseudo_inverse -> dimr_sca_Sigmak, dimc_sca_Sigmak "// vtoa(dimr_sca_Sigmak) // " " // vtoa(dimc_sca_Sigmak))
        !write(6,*)iam, 1/Sigma_SVD(1:dimr_sca_Sigmak)
        call log_info("scalapack_pseudo_inverse ->  sca_Sigmak completed")
      end if 
      nbr_tmp = min(nbr_AA, nbc_AA)
      nbc_tmp = nbr_AA
      dimr_sca_tmp = rank_sca_AA
      dimc_sca_tmp = Mline
      call allocation_scalapack_matrix (" sca_tmp in pinv ", sca_tmp, desc_sca_tmp, dimr_sca_tmp, dimc_sca_tmp, &
         l_dimr_sca_tmp, l_dimc_sca_tmp, nbr_tmp, nbc_tmp, &
         myrow, mycol, nprow, npcol, context, iam  )
      !call pdgemm(transa, transb, m, n, k, alpha, a, ia, ja, desca, b, ib, jb, descb, beta, c, ic, jc, descc)
      ! C(mxn) = A(mxk) B(kxn)
      ! tmp(rxm) = Sigma(rxr)*transpose(Uk)(rxm)
      ! C = alpha A^transa*B^transb + beta C
      call blacs_barrier(context, 'A')
      !debug! if (debug_scalapack) then
      !debug!   if (iam == 0 ) then 
      !debug!   write(6,*) "dimr_sca_tmp=", dimr_sca_tmp, "dimc_sca_tmp=", dimc_sca_tmp
      !debug!   write(6,*) "dimr_sca_Sigmak=", dimr_sca_Sigmak, "dimc_sca_Sigmak=", dimc_sca_Sigmak
      !debug!   write(6,*) "dimr_sca_UUk=", dimr_sca_UUk, "dimc_sca_UUk=", dimc_sca_UUk
      !debug!   write(6,*) "desc_sca_Sigmak=", desc_sca_Sigmak
      !debug!   write(6,*) "desc_sca_UUk=", desc_sca_UUk
      !debug!   write(6,*) "desc_sca_tmp=", desc_sca_tmp
      !debug!   end if 
      !debug!   if (.not.allocated(sca_Sigmak)) call log_critical("sca_Sigmak not allocated!")
      !debug!   if (.not.allocated(sca_UUk)) call log_critical("sca_UUk not allocated!")
      !debug!   if (.not.allocated(sca_tmp)) call log_critical("sca_tmp not allocated!")
      !debug! end if

      !PDGEMM_CORR  CORRECTION HERE! 
      call pdgemm('N','T', dimr_sca_tmp, dimc_sca_tmp, dimc_sca_UUk, one, sca_Sigmak, 1, 1, desc_sca_Sigmak, &
         sca_UUk, 1, 1,    desc_sca_UUk, &
         zero,    sca_tmp, 1, 1,    desc_sca_tmp  )

      if (debug_scalapack) call log_info("scalapack_pseudo_inverse -> first pdgemm   completed")
      call blacs_barrier(context, 'A')

      if (allocated(sca_UUk)) deallocate(sca_UUk)
      if (allocated(desc_sca_UUk)) deallocate(desc_sca_UUk)

      ! Allocate scalapack VTk - cropped version until the rank:
      nbr_VTk =  min(nbr_AA, nbc_AA)
      nbc_VTk =  nbc_AA
      dimr_sca_VTk = rank_sca_AA
      dimc_sca_VTk = Ncolm
      call allocation_scalapack_matrix (" sca_VTk in pinv ", sca_VTk, desc_sca_VTk, dimr_sca_VTk, dimc_sca_VTk, &
         l_dimr_sca_VTk, l_dimc_sca_VTk, nbr_VTk, nbc_VTk, &
         myrow, mycol, nprow, npcol, context, iam  )

      call my_pdgemr2d(dimr_sca_VTk, dimc_sca_VTk, sca_VT, 1, 1, desc_sca_VT, sca_VTk, 1, 1, desc_sca_VTk, context)
      if (allocated(sca_VT)) deallocate(sca_VT)
      if (allocated(desc_sca_VT)) deallocate(desc_sca_VT)
      if (debug_scalapack) call log_info("scalapack_pseudo_inverse -> second pdgemr2d   completed")

      !allocate sca_pinvAA
      dimr_sca_pinvAA  = dimc_sca_AA
      dimc_sca_pinvAA  = dimr_sca_AA
      nbr_pinvAA = nbc_AA
      nbc_pinvAA = nbr_AA
      call allocation_scalapack_matrix (" sca_pinvAA in pinv ", sca_pinvAA, desc_sca_pinvAA, dimr_sca_pinvAA, dimc_sca_pinvAA, &
         l_dimr_sca_pinvAA, l_dimc_sca_pinvAA, nbr_pinvAA, nbc_pinvAA, &
         myrow, mycol, nprow, npcol, context, iam  )

      ! pinv(nxm) = transpose(VT)(nxr)*tmp(rxm)
      call pdgemm('T', 'N', dimr_sca_pinvAA, dimc_sca_pinvAA, dimr_sca_tmp, one,   sca_VTk,  1, 1, desc_sca_VTk, &
         sca_tmp,  1, 1, desc_sca_tmp, &
         zero, sca_pinvAA, 1, 1, desc_sca_pinvAA)
      if (debug_scalapack) call log_info("scalapack_pseudo_inverse -> second pdgemm   completed")

      if (allocated(sca_VTk)) deallocate(sca_VTk)
      if (allocated(desc_sca_VTk)) deallocate(desc_sca_VTk)
      if (allocated(sca_tmp)) deallocate(sca_tmp)
      if (allocated(desc_sca_tmp)) deallocate(desc_sca_tmp)

      !!------------------------testing sca_pinvAA----------------------------!
      !nbr_tmp =  nbc_AA
      !nbc_tmp =  nbc_AA
      !dimr_sca_tmp = Ncolm
      !dimc_sca_tmp = Ncolm
      !call allocation_scalapack_matrix (" sca_mtmp ", sca_tmp, desc_sca_tmp, dimr_sca_tmp, dimc_sca_tmp, &
      !                                  l_dimr_sca_tmp, l_dimc_sca_tmp, nbr_tmp, nbc_tmp, &
      !                                  myrow, mycol, nprow, npcol, context, iam  )
      !!C(mxn) = A(mxk)*B(kxn)
      !!mmtp(nxn) = pinvAA(nxm)*AA(mxn)
      !call pdgemm('N', 'N', dimr_sca_tmp, dimc_sca_tmp, dimc_sca_pinvAA, one,   sca_pinvAA,  1, 1, desc_sca_pinvAA, &
      !                                                                             sca_AA,  1, 1, desc_sca_AA, &
      !                                                                     zero, sca_tmp, 1, 1, desc_sca_tmp)
      !nbr_mtest =  nbr_AA
      !nbc_mtest =  nbc_AA
      !dimr_sca_mtest = Mline
      !dimc_sca_mtest = Ncolm
      !call allocation_scalapack_matrix (" sca_mtest ", sca_mtest, desc_sca_mtest, dimr_sca_mtest, dimc_sca_mtest, &
      !                                  l_dimr_sca_mtest, l_dimc_sca_mtest, nbr_mtest, nbc_mtest, &
      !                                  myrow, mycol, nprow, npcol, context, iam  )
      !!mtest(mxn) = AA(mxn)*mtmp(nxn)
      !call pdgemm('N', 'N', dimr_sca_AA, dimc_sca_AA, dimc_sca_AA,   one,   sca_AA,  1, 1, desc_sca_AA, &
      !                                                                             sca_tmp,  1, 1, desc_sca_tmp, &
      !                                                                     zero, sca_mtest, 1, 1, desc_sca_mtest)
      !! mtest(mxn) = AA(mxn) - mtest(mxn)
      !do jj=1, l_dimc_sca_AA
      !  do ii=1, l_dimr_sca_AA
      !    sca_mtest(ii,jj) = sca_AA(ii,jj) - sca_mtest(ii,jj)
      !  end do
      !end do
      !if (allocated(work)) deallocate(work)
      !allocate(work(1))
      !tmpval =  pdlange('F', Mline, Ncolm, sca_mtest, 1, 1, desc_sca_mtest, work)
      !write(*,*) 'pinvAA test.................', tmpval
      !!------------------------testing sca_pinvAA----------------------------!
      call blacs_barrier(context, 'A')

      _MLD_END_
   end subroutine scalapack_pseudo_inverse



end module module_fit_ScaMatrix



subroutine set_up_scalapack()
   !
   ! The design matrix Amat is distributed across ScaLAPACK.
   ! /------------------dimc_sca_MATRIX-----------------\
   ! |                                                  |
   ! |                                                  |
   ! dimr_sca_MATRIX                                    |
   ! |                                                  |
   ! |                                                  |
   ! \--------------------------------------------------/
   ! -> dimr ... : dim_design_line (r.d. 1 + D + D^2 etc )
   ! -> dimc ... : the number of database
   !               points energy, forces, stress etc
   !
   use module_ml_scalapack, only:   iam, nprocs_ml_sca, nprow, npcol, context, &
      myrow, mycol, iproc_sca
   use module_scalapack_tools, only: gridsetup_ml, blockset_ml, allocation_scalapack_matrix
   use mld_logger
   implicit none
   _NAMECURRENT_("set_up_scalapack")


   _MLD_BEGIN_
   ! -----    Initialize the blacs.  Note: processors are counted starting at 0.
   call blacs_pinfo(iam, nprocs_ml_sca)
   iproc_sca = iam
   if (iam == 0) write (6, '("ML sca: nprocs on scalapack grid = ",i6)') nprocs_ml_sca

   ! -----    Set the dimension of the 2d processors grid.
   !
   ! factorizes the number of processors (nproc) into nprow and npcol
   ! that are the sizes of the 2d processors mesh.
   call gridsetup_ml(nprocs_ml_sca, nprow, npcol)
   !
   ! -----    Initialize a single blacs context.  Determine which processor I
   !          am in the 2D process or grid.
   !
   call blacs_get(-1, 0, context)
   call blacs_gridinit(context, 'r', nprow, npcol)
   call blacs_gridinfo(context, nprow, npcol, myrow, mycol)
   if (iam == 0) write (6, '("ML sca: the blacs grid has number of row and col set to  ...........:", 2i7)') nprow, npcol

   _MLD_END_

end subroutine set_up_scalapack





