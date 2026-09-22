! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   MiLaDy (Machine Learning Dynamics) 
! HND X   
! HND X
! HND X   Milady was written and designed by Mihai-Cosmin Marinica, Alexandra M. Goryaeva 
! HND X   Copyright 2015-2025.
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

!-------------------------------------------------------------------------------
! module_nnls_bk  –  NNLS solver with warm start (Bunch-Kaufman backend)
!
! A nonnegative least-squares solver with optional ridge regularization,
! designed for robust constrained linear solves and repeated use inside
! ALS-style block-renormalization loops.
!
! The solver combines:
!   - a Lawson-Hanson active-set strategy,
!   - LAPACK DSYTRF / DSYTRS for robust reduced-system solves,
!   - optional warm starts,
!   - a specialized helper for the nonnegative alpha-update in ALS.
!
! Solves:   min_{x >= 0}  0.5 ||A x - b||^2  +  0.5 rho ||x||^2
!
! Public API:
!   nnls_solver (derived type) with methods:
!     solve            – cold start (x = 0)
!     solve_warm       – warm start from previous solution
!     solve_alpha_step – ALS alpha-step wrapper (forms G = A^T A internally)
!     solve_gram_warm  – warm start on pre-formed Gram system G x = c
!                        (avoids squaring when G and c are already available,
!                         e.g. the BtB / Bty alpha sub-problem in ALS)
!
! Adapted from nnls_mode/module_nnls.F90 to use MiLaDy kind_double.
!-------------------------------------------------------------------------------
module module_nnls_bk
   use module_kind_variables, only: kind_double
   implicit none
   private

   public :: nnls_solver

   integer, parameter :: dp = kind_double

   !===========================================================
   ! LAPACK interfaces
   !===========================================================
   interface
      subroutine dsytrf(uplo, n, a, lda, ipiv, work, lwork, info)
         import :: dp
         character(len=1), intent(in) :: uplo
         integer, intent(in) :: n, lda, lwork
         real(dp), intent(inout) :: a(lda,*)
         integer, intent(out) :: ipiv(*)
         real(dp), intent(inout) :: work(*)
         integer, intent(out) :: info
      end subroutine dsytrf

      subroutine dsytrs(uplo, n, nrhs, a, lda, ipiv, b, ldb, info)
         import :: dp
         character(len=1), intent(in) :: uplo
         integer, intent(in) :: n, nrhs, lda, ldb
         real(dp), intent(in) :: a(lda,*)
         integer, intent(in) :: ipiv(*)
         real(dp), intent(inout) :: b(ldb,*)
         integer, intent(out) :: info
      end subroutine dsytrs
   end interface

   !===========================================================
   ! OO solver
   !===========================================================
   type :: nnls_solver
      real(dp) :: tol = 1.0d-16
      integer  :: max_outer = -1
      integer  :: max_inner = -1
      real(dp) :: regularization = 0.0d0
      real(dp) :: jitter0 = 1.0d-14
      real(dp) :: jitter_max = 1.0d-4
   contains
      procedure :: solve              => nnls_solve_cold
      procedure :: solve_warm         => nnls_solve_warm
      procedure :: solve_alpha_step   => nnls_solve_alpha_step
      procedure :: solve_gram_warm    => nnls_solve_gram_warm
   end type nnls_solver

contains

   !===========================================================
   ! Cold start solve:
   !   x starts at 0
   !===========================================================
   subroutine nnls_solve_cold(self, A, b, x, rnorm, info, rho)
      class(nnls_solver), intent(in) :: self
      real(dp), intent(in)  :: A(:,:), b(:)
      real(dp), intent(out) :: x(:)
      real(dp), intent(out) :: rnorm
      integer,  intent(out) :: info
      real(dp), intent(in), optional :: rho

      logical, allocatable :: passive(:)

      allocate(passive(size(x)))
      passive = .false.
      x = 0.0_dp

      call nnls_core(self, A, b, x, passive, rnorm, info, rho)

      deallocate(passive)
   end subroutine nnls_solve_cold

   !===========================================================
   ! Warm start solve:
   !   x contains initial guess on entry, final solution on exit
   !   passive may be given and updated
   !===========================================================
   subroutine nnls_solve_warm(self, A, b, x, rnorm, info, passive, rho)
      class(nnls_solver), intent(in) :: self
      real(dp), intent(in)    :: A(:,:), b(:)
      real(dp), intent(inout) :: x(:)
      real(dp), intent(out)   :: rnorm
      integer,  intent(out)   :: info
      logical,  intent(inout), optional :: passive(:)
      real(dp), intent(in), optional    :: rho

      logical, allocatable :: passive_loc(:)

      if (present(passive)) then
         call sanitize_warm_start(x, passive)
         call nnls_core(self, A, b, x, passive, rnorm, info, rho)
      else
         allocate(passive_loc(size(x)))
         call infer_passive_from_x(x, passive_loc, self%tol)
         call nnls_core(self, A, b, x, passive_loc, rnorm, info, rho)
         deallocate(passive_loc)
      end if
   end subroutine nnls_solve_warm

   !===========================================================
   ! Specialized helper for ALS alpha-step:
   !
   ! solves
   !   min_{alpha >= 0} 0.5 || B alpha - y ||^2 + 0.5 rho ||alpha||^2
   !
   ! alpha can be warm-started.
   !===========================================================
   subroutine nnls_solve_alpha_step(self, Bmat, y, alpha, info, rnorm, rho, passive)
      class(nnls_solver), intent(in) :: self
      real(dp), intent(in)    :: Bmat(:,:), y(:)
      real(dp), intent(inout) :: alpha(:)
      integer,  intent(out)   :: info
      real(dp), intent(out), optional :: rnorm
      real(dp), intent(in), optional  :: rho
      logical,  intent(inout), optional :: passive(:)

      real(dp) :: rnorm_loc

      call self%solve_warm(Bmat, y, alpha, rnorm_loc, info, passive, rho)

      if (present(rnorm)) rnorm = rnorm_loc
   end subroutine nnls_solve_alpha_step

   !===========================================================
   ! Gram-system warm-start solve:
   !
   ! Solves   min_{x >= 0}  subject to KKT on pre-formed
   !   G x = c      where G = A^T A (+rho*I if desired)
   !                       c = A^T b
   !
   ! This avoids the internal A^T A formation, so G is NOT
   ! squared. Use this when G and c are already available
   ! (e.g., from B^T B and B^T y in the alpha sub-problem).
   !
   ! x is warm-started from its value on entry.
   ! passive mask is maintained across calls.
   !===========================================================
   subroutine nnls_solve_gram_warm(self, G, c, x, info, passive)
      class(nnls_solver), intent(in) :: self
      real(dp), intent(in)    :: G(:,:)   ! (n, n) pre-formed Gram matrix
      real(dp), intent(in)    :: c(:)     ! (n)    pre-formed rhs = A^T b
      real(dp), intent(inout) :: x(:)     ! (n)    solution, warm-started
      integer,  intent(out)   :: info
      logical,  intent(inout), optional :: passive(:)

      logical, allocatable :: passive_loc(:)
      integer :: n

      n = size(x)

      if (present(passive)) then
         call sanitize_warm_start(x, passive)
         call nnls_core_gram(self, G, c, n, x, passive, info)
      else
         allocate(passive_loc(n))
         call infer_passive_from_x(x, passive_loc, self%tol)
         call nnls_core_gram(self, G, c, n, x, passive_loc, info)
         deallocate(passive_loc)
      end if
   end subroutine nnls_solve_gram_warm

   !===========================================================
   ! Core Lawson-Hanson active-set NNLS on pre-formed Gram system
   !
   ! G and c are provided directly (not formed from A).
   ! The KKT dual vector is  w = c - G x.
   !===========================================================
   subroutine nnls_core_gram(self, G, c, n, x, passive, info)
      class(nnls_solver), intent(in) :: self
      integer,  intent(in)    :: n
      real(dp), intent(in)    :: G(n,n), c(n)
      real(dp), intent(inout) :: x(n)
      logical,  intent(inout) :: passive(n)
      integer,  intent(out)   :: info

      integer :: i, j, outer_it, inner_it
      integer :: max_outer_loc, max_inner_loc
      integer :: idx_add, npass, lap_info
      real(dp) :: tol_loc, alpha_step, denom
      logical :: done_outer, all_positive, any_zeroed

      real(dp), allocatable :: w(:), z(:), x_old(:)
      integer,  allocatable :: pset(:)

      tol_loc = self%tol

      max_outer_loc = self%max_outer
      if (max_outer_loc <= 0) max_outer_loc = 5*n + 10

      max_inner_loc = self%max_inner
      if (max_inner_loc <= 0) max_inner_loc = 10*n + 20

      ! Enforce nonnegative warm start
      do i = 1, n
         if (x(i) < 0.0_dp) x(i) = 0.0_dp
      end do
      call infer_passive_from_x(x, passive, tol_loc)

      allocate(w(n), z(n), x_old(n), pset(n))

!debug!      info = 0
!debug!
!debug!      ! --- Initial reduced solve on warm-started passive set ---
!debug!      ! When warm-starting (e.g. alpha = [1,...,1], all passive), the
!debug!      ! outer loop would immediately exit because all duals are zeroed
!debug!      ! for passive variables.  We must first solve the unconstrained
!debug!      ! system on the current passive set and fix any negatives.
!debug!      call get_passive_set(passive, pset, npass)
!debug!      if (npass > 0) then
!debug!         z = 0.0_dp
!debug!         call solve_reduced_sym_bk(G, c, pset, npass, z, lap_info, &
!debug!                                   self%jitter0, self%jitter_max)
!debug!         if (lap_info /= 0) then
!debug!            info = 3
!debug!            deallocate(w, z, x_old, pset)
!debug!            return
!debug!         end if
!debug!
!debug!         ! Check if reduced solution is all-positive on passive set
!debug!         all_positive = .true.
!debug!         do j = 1, npass
!debug!            i = pset(j)
!debug!            if (z(i) <= tol_loc) then
!debug!               all_positive = .false.
!debug!               exit
!debug!            end if
!debug!         end do
!debug!
!debug!         if (all_positive) then
!debug!            ! Accept unconstrained solution directly
!debug!            x = z
!debug!         else
!debug!            ! Interpolate toward z, removing variables that hit zero
!debug!            alpha_step = huge(1.0_dp)
!debug!            do j = 1, npass
!debug!               i = pset(j)
!debug!               if (z(i) <= tol_loc) then
!debug!                  denom = x(i) - z(i)
!debug!                  if (denom > 0.0_dp) alpha_step = min(alpha_step, x(i) / denom)
!debug!               end if
!debug!            end do
!debug!            if (.not. (alpha_step < huge(1.0_dp))) alpha_step = 1.0_dp
!debug!            x = x + alpha_step * (z - x)
!debug!            do j = 1, npass
!debug!               i = pset(j)
!debug!               if (x(i) <= tol_loc) then
!debug!                  x(i) = 0.0_dp
!debug!                  passive(i) = .false.
!debug!               end if
!debug!            end do
!debug!         end if
!debug!      end if

      ! Dual vector: w = c - G x
      w = c - matmul(G, x)
      do i = 1, n
         if (passive(i)) w(i) = 0.0_dp
      end do

      outer_it = 0
      done_outer = .false.

      do while (.not. done_outer)
         outer_it = outer_it + 1
         if (outer_it > max_outer_loc) then
            info = 1
            exit
         end if

         ! Add best inactive variable with positive dual gradient
         idx_add = argmax_inactive_positive(w, passive, tol_loc)
         if (idx_add == 0) then
            done_outer = .true.
            exit
         end if
         passive(idx_add) = .true.

         inner_it = 0
         do
            inner_it = inner_it + 1
            if (inner_it > max_inner_loc) then
               info = 2
               exit
            end if

            x_old = x
            z = 0.0_dp

            call get_passive_set(passive, pset, npass)

            if (npass > 0) then
               call solve_reduced_sym_bk(G, c, pset, npass, z, lap_info, &
                                         self%jitter0, self%jitter_max)
               if (lap_info /= 0) then
                  info = 3
                  exit
               end if
            end if

            ! If reduced solution is strictly positive on passive set, accept it
            all_positive = .true.
            do j = 1, npass
               i = pset(j)
               if (z(i) <= tol_loc) then
                  all_positive = .false.
                  exit
               end if
            end do

            if (all_positive) then
               x = z
               exit
            end if

            ! Step from x toward z until one passive component hits zero
            alpha_step = huge(1.0_dp)
            do j = 1, npass
               i = pset(j)
               if (z(i) <= tol_loc) then
                  denom = x(i) - z(i)
                  if (denom > 0.0_dp) alpha_step = min(alpha_step, x(i) / denom)
               end if
            end do

            if (.not. (alpha_step < huge(1.0_dp))) alpha_step = 1.0_dp
            x = x + alpha_step * (z - x)

            ! Remove variables that hit zero
            any_zeroed = .false.
            do j = 1, npass
               i = pset(j)
               if (x(i) <= tol_loc) then
                  x(i) = 0.0_dp
                  passive(i) = .false.
                  any_zeroed = .true.
               end if
            end do

            ! Safety fallback
            if (.not. any_zeroed) then
               do i = 1, n
                  x(i) = max(z(i), 0.0_dp)
               end do
               call infer_passive_from_x(x, passive, tol_loc)
               exit
            end if
         end do

         if (info /= 0) exit

         ! Update dual/KKT vector
         w = c - matmul(G, x)
         do i = 1, n
            if (passive(i)) w(i) = 0.0_dp
         end do
      end do

      deallocate(w, z, x_old, pset)
   end subroutine nnls_core_gram

   !===========================================================
   ! Core Lawson-Hanson active-set NNLS + ridge
   !===========================================================
   subroutine nnls_core(self, A, b, x, passive, rnorm, info, rho)
      class(nnls_solver), intent(in) :: self
      real(dp), intent(in)    :: A(:,:), b(:)
      real(dp), intent(inout) :: x(:)
      logical,  intent(inout) :: passive(:)
      real(dp), intent(out)   :: rnorm
      integer,  intent(out)   :: info
      real(dp), intent(in), optional :: rho

      integer :: m, n
      integer :: i, j, outer_it, inner_it
      integer :: max_outer_loc, max_inner_loc
      integer :: idx_add, npass, lap_info
      real(dp) :: rho_loc, tol_loc, alpha_step, denom
      logical :: done_outer, all_positive, any_zeroed

      real(dp), allocatable :: G(:,:), c(:), w(:), z(:), x_old(:), resid(:)
      integer,  allocatable :: pset(:)

      m = size(A,1)
      n = size(A,2)

      if (size(b) /= m) error stop "nnls_core: size(b) /= size(A,1)"
      if (size(x) /= n) error stop "nnls_core: size(x) /= size(A,2)"
      if (size(passive) /= n) error stop "nnls_core: size(passive) /= size(A,2)"

      tol_loc = self%tol
      rho_loc = self%regularization
      if (present(rho)) rho_loc = rho

      max_outer_loc = self%max_outer
      if (max_outer_loc <= 0) max_outer_loc = 5*n + 10

      max_inner_loc = self%max_inner
      if (max_inner_loc <= 0) max_inner_loc = 10*n + 20

      ! Enforce nonnegative warm start
      do i = 1, n
         if (x(i) < 0.0_dp) x(i) = 0.0_dp
      end do
      call infer_passive_from_x(x, passive, tol_loc)

      allocate(G(n,n), c(n), w(n), z(n), x_old(n), resid(m), pset(n))

      ! Normal equations
      G = matmul(transpose(A), A)
      if (rho_loc > 0.0_dp) then
         do i = 1, n
            G(i,i) = G(i,i) + rho_loc
         end do
      end if
      c = matmul(transpose(A), b)

      ! Dual vector: w = c - G x
      w = c - matmul(G, x)
      do i = 1, n
         if (passive(i)) w(i) = 0.0_dp
      end do

      info = 0
      outer_it = 0
      done_outer = .false.

      do while (.not. done_outer)
         outer_it = outer_it + 1
         if (outer_it > max_outer_loc) then
            info = 1
            exit
         end if

         ! Add best inactive variable with positive dual gradient
         idx_add = argmax_inactive_positive(w, passive, tol_loc)
         if (idx_add == 0) then
            done_outer = .true.
            exit
         end if
         passive(idx_add) = .true.

         inner_it = 0
         do
            inner_it = inner_it + 1
            if (inner_it > max_inner_loc) then
               info = 2
               exit
            end if

            x_old = x
            z = 0.0_dp

            call get_passive_set(passive, pset, npass)

            if (npass > 0) then
               call solve_reduced_sym_bk(G, c, pset, npass, z, lap_info, self%jitter0, self%jitter_max)
               if (lap_info /= 0) then
                  info = 3
                  exit
               end if
            end if

            ! If reduced solution is strictly positive on passive set, accept it
            all_positive = .true.
            do j = 1, npass
               i = pset(j)
               if (z(i) <= tol_loc) then
                  all_positive = .false.
                  exit
               end if
            end do

            if (all_positive) then
               x = z
               exit
            end if

            ! Otherwise step from x toward z until one passive component hits zero
            alpha_step = huge(1.0_dp)
            do j = 1, npass
               i = pset(j)
               if (z(i) <= tol_loc) then
                  denom = x(i) - z(i)
                  if (denom > 0.0_dp) alpha_step = min(alpha_step, x(i) / denom)
               end if
            end do

            if (.not. (alpha_step < huge(1.0_dp))) alpha_step = 1.0_dp
            x = x + alpha_step * (z - x)

            ! Remove variables that hit zero
            any_zeroed = .false.
            do j = 1, npass
               i = pset(j)
               if (x(i) <= tol_loc) then
                  x(i) = 0.0_dp
                  passive(i) = .false.
                  any_zeroed = .true.
               end if
            end do

            ! Safety fallback
            if (.not. any_zeroed) then
               do i = 1, n
                  x(i) = max(z(i), 0.0_dp)
               end do
               call infer_passive_from_x(x, passive, tol_loc)
               exit
            end if
         end do

         if (info /= 0) exit

         ! Update dual/KKT vector
         w = c - matmul(G, x)
         do i = 1, n
            if (passive(i)) w(i) = 0.0_dp
         end do
      end do

      resid = matmul(A, x) - b
      rnorm = sqrt(dot_product(resid, resid))

      deallocate(G, c, w, z, x_old, resid, pset)
   end subroutine nnls_core

   !===========================================================
   ! Reduced symmetric solve on active set:
   !   G_PP z_P = c_P
   !
   ! Uses DSYTRF / DSYTRS and adds diagonal jitter if needed.
   !===========================================================
   subroutine solve_reduced_sym_bk(G, c, pset, npass, z, info, jitter0, jitter_max)
      real(dp), intent(in) :: G(:,:), c(:)
      integer,  intent(in) :: pset(:), npass
      real(dp), intent(inout) :: z(:)
      integer,  intent(out) :: info
      real(dp), intent(in) :: jitter0, jitter_max

      integer :: i, j, ii, jj
      integer :: lwork, info_lap
      real(dp) :: work_query, jitter
      real(dp), allocatable :: Gp(:,:), rhs(:,:), work(:)
      integer,  allocatable :: ipiv(:)

      info = 0
      if (npass <= 0) return

      allocate(Gp(npass,npass), rhs(npass,1), ipiv(npass))

      jitter = 0.0_dp

      do
         ! Build reduced matrix and rhs
         do j = 1, npass
            jj = pset(j)
            rhs(j,1) = c(jj)
            do i = 1, npass
               ii = pset(i)
               Gp(i,j) = G(ii,jj)
            end do
         end do

         if (jitter > 0.0_dp) then
            do i = 1, npass
               Gp(i,i) = Gp(i,i) + jitter
            end do
         end if

         ! Workspace query
         lwork = -1
         allocate(work(1))
         call dsytrf('U', npass, Gp, npass, ipiv, work, lwork, info_lap)
         if (info_lap /= 0) then
            deallocate(work)
            info = 1
            deallocate(Gp, rhs, ipiv)
            return
         end if
         work_query = work(1)
         deallocate(work)

         lwork = max(1, int(work_query))
         allocate(work(lwork))

         call dsytrf('U', npass, Gp, npass, ipiv, work, lwork, info_lap)
         deallocate(work)

         if (info_lap == 0) exit

         ! Singular or numerically bad -> retry with more jitter
         if (jitter == 0.0_dp) then
            jitter = jitter0
         else
            jitter = 10.0_dp * jitter
         end if

         if (jitter > jitter_max) then
            info = 2
            deallocate(Gp, rhs, ipiv)
            return
         end if
      end do

      call dsytrs('U', npass, 1, Gp, npass, ipiv, rhs, npass, info_lap)
      if (info_lap /= 0) then
         info = 3
         deallocate(Gp, rhs, ipiv)
         return
      end if

      do j = 1, npass
         z(pset(j)) = rhs(j,1)
      end do

      deallocate(Gp, rhs, ipiv)
   end subroutine solve_reduced_sym_bk

   !===========================================================
   ! Helpers
   !===========================================================
   integer function argmax_inactive_positive(w, passive, tol) result(idx)
      real(dp), intent(in) :: w(:), tol
      logical,  intent(in) :: passive(:)
      integer :: i
      real(dp) :: best

      idx = 0
      best = tol
      do i = 1, size(w)
         if (.not. passive(i)) then
            if (w(i) > best) then
               best = w(i)
               idx = i
            end if
         end if
      end do
   end function argmax_inactive_positive

   subroutine get_passive_set(passive, pset, npass)
      logical, intent(in)  :: passive(:)
      integer, intent(out) :: pset(:)
      integer, intent(out) :: npass
      integer :: i

      npass = 0
      do i = 1, size(passive)
         if (passive(i)) then
            npass = npass + 1
            pset(npass) = i
         end if
      end do
   end subroutine get_passive_set

   subroutine infer_passive_from_x(x, passive, tol)
      real(dp), intent(in)  :: x(:), tol
      logical,  intent(out) :: passive(:)
      integer :: i

      do i = 1, size(x)
         passive(i) = (x(i) > tol)
      end do
   end subroutine infer_passive_from_x

   subroutine sanitize_warm_start(x, passive)
      real(dp), intent(inout) :: x(:)
      logical,  intent(inout) :: passive(:)
      integer :: i

      if (size(x) /= size(passive)) error stop "sanitize_warm_start: size mismatch"

      do i = 1, size(x)
         if (x(i) < 0.0_dp) x(i) = 0.0_dp
         passive(i) = (x(i) > 0.0_dp)
      end do
   end subroutine sanitize_warm_start

end module module_nnls_bk
