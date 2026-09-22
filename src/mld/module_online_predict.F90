!===================================================================!
!     MiLaDy - Machine Learning Dynamics                            !
!                                                                   !
!     module_online_predict.F90                                     !
!     Post-fit prediction and error evaluation for online fit       !
!     (mld_fit_type = 6).                                           !
!                                                                   !
!     STREAMING VERSION (March 2026):                               !
!     Train error statistics (R², correlation, RMSE, MAE) are       !
!     computed via O(1) scalar partial-sum accumulators instead      !
!     of O(C×N) global arrays. This eliminates the memory           !
!     bottleneck that made Phase 3 impossible for large databases.  !
!     Per-config file outputs (train_*.out) are still streaming.    !
!                                                                   !
!     See docs/README_ONLINE_FIT.md for design.                     !
!===================================================================!

#include "../MLD_MACROS.INC"

module module_online_predict
   use module_kind_variables, only: kind_double
   use mld_logger
   implicit none
   private

   public :: online_predict_train_errors

   ! Number of streaming statistics per quantity type:
   !   1=n, 2=sum_y, 3=sum_yh, 4=sum_y², 5=sum_yh²,
   !   6=sum_y·yh, 7=sum_(y-yh)², 8=sum_|y-yh|
   integer, parameter :: NSTAT = 8

contains

   !----------------------------------------------------------------
   ! Accumulate streaming statistics for a single (base, predicted)
   ! data point into partial-sum vector s(1:NSTAT).
   !----------------------------------------------------------------
   subroutine accumulate_stats(s, y_base, y_snap)
      real(kind_double), intent(inout) :: s(NSTAT)
      real(kind_double), intent(in) :: y_base, y_snap
      real(kind_double) :: diff
      diff = y_base - y_snap
      s(1) = s(1) + 1.0_kind_double
      s(2) = s(2) + y_base
      s(3) = s(3) + y_snap
      s(4) = s(4) + y_base * y_base
      s(5) = s(5) + y_snap * y_snap
      s(6) = s(6) + y_base * y_snap
      s(7) = s(7) + diff * diff
      s(8) = s(8) + dabs(diff)
   end subroutine accumulate_stats


   !----------------------------------------------------------------
   ! Compute correlation, determination coef, RMSE, MAE from
   ! streaming partial sums. Reproduces the same formulas as
   ! correlation_coef, determination_coef, rmse_mae in stats_tools.
   !
   ! Convention: y = base (measured), yh = snap (predicted).
   !----------------------------------------------------------------
   subroutine stats_from_sums(s, corr, detr, rmse, mae)
      real(kind_double), intent(in) :: s(NSTAT)
      real(kind_double), intent(out) :: corr, detr, rmse, mae

      real(kind_double) :: n, mean_y, mean_yh
      real(kind_double) :: var_y, var_yh, cov_yyh, sigma_y, sigma_yh
      real(kind_double) :: ss_model, ss_tot

      n = s(1)
      corr = 1.0_kind_double; detr = 1.0_kind_double
      rmse = 0.0_kind_double; mae = 0.0_kind_double

      if (n < 2.0_kind_double) return

      mean_y  = s(2) / n
      mean_yh = s(3) / n

      ! RMSE = sqrt( sum((y-yh)²) / n ),  MAE = sum(|y-yh|) / n
      rmse = dsqrt(s(7) / n)
      mae  = s(8) / n

      ! Population variances and covariance
      var_y   = s(4) / n - mean_y * mean_y
      var_yh  = s(5) / n - mean_yh * mean_yh
      cov_yyh = s(6) / n - mean_y * mean_yh
      sigma_y  = dsqrt(dabs(var_y))
      sigma_yh = dsqrt(dabs(var_yh))

      ! Correlation coefficient (same guard as stats_tools.F90 correlation_coef)
      if (sigma_yh <= 1.0e-20_kind_double .and. &
          dabs(n * cov_yyh) <= 1.0e-20_kind_double) then
         corr = 1.0_kind_double
      else
         corr = cov_yyh / (sigma_y * sigma_yh)
      end if

      ! Determination coefficient as in stats_tools.F90 determination_coef:
      !   R² = sum((y_predicted - mean_measured)²) / sum((y_measured - mean_measured)²)
      ss_tot   = s(4) - n * mean_y * mean_y
      ss_model = s(5) - 2.0_kind_double * mean_y * s(3) + n * mean_y * mean_y
      if (ss_tot > 0.0_kind_double) then
         detr = ss_model / ss_tot
      else
         detr = 1.0_kind_double
      end if
   end subroutine stats_from_sums


   !----------------------------------------------------------------
   ! Ordered, low-memory merge of per-subworld scratch files into a
   ! single output file. Each subworld master (subrank==0) buffers
   ! its own configs' lines in a scratch unit during the streaming
   ! loop. Here the masters append their scratch content to
   ! <final_name> one after another, in id_subworld order, using a
   ! token passed over masters_world. The result is a single file
   ! containing every subworld's configs, without ever gathering the
   ! full O(C x N) data onto a single rank.
   !----------------------------------------------------------------
   subroutine merge_master_scratch(scratch_unit, final_name)
      use mld_subworld, only: subrank, id_subworld, nb_subworlds, masters_world
#if(PARA)
      use mpi
#endif
      integer, intent(in) :: scratch_unit
      character(len=*), intent(in) :: final_name
      integer :: u_final, ios, ierr
      integer :: token
      character(len=1024) :: line
#if(PARA)
      integer :: mpi_stat(MPI_STATUS_SIZE)
#endif

      ! Only subworld masters hold scratch data.
      if (subrank /= 0) return

      token = 0
#if(PARA)
      ! Wait until the previous master has finished writing the file.
      if (id_subworld > 0) &
         call MPI_Recv(token, 1, MPI_INTEGER, id_subworld - 1, 99, masters_world, mpi_stat, ierr)
#endif

      if (id_subworld == 0) then
         open(newunit=u_final, file=final_name, status='replace', action='write')
      else
         open(newunit=u_final, file=final_name, status='unknown', position='append', action='write')
      end if

      rewind(scratch_unit)
      do
         read(scratch_unit, '(A)', iostat=ios) line
         if (ios /= 0) exit
         write(u_final, '(A)') trim(line)
      end do
      close(u_final)

#if(PARA)
      ! Hand the token to the next master.
      if (id_subworld < nb_subworlds - 1) &
         call MPI_Send(token, 1, MPI_INTEGER, id_subworld + 1, 99, masters_world, ierr)
#endif
   end subroutine merge_master_scratch


   !----------------------------------------------------------------
   ! Recompute train predictions and evaluate train RMSE/MAE
   ! using STREAMING partial-sum accumulators — O(1) memory in C.
   !
   ! Instead of allocating global arrays of size O(C×N) for all
   ! predictions, we accumulate 8 scalar sums per quantity type
   ! (energy, energy/atom, force, stress) during the config loop,
   ! then reduce across subworlds with a single MPI allreduce of
   ! 32 doubles. File output (train_*.out) is still per-config.
   !----------------------------------------------------------------
   subroutine online_predict_train_errors(w_params)
      use ml_in_ndm_module, only: rangml, desc_forces
      use derived_types, only: config_real
      use module_db_poscar, only: iread_energy, i_start_conf, i_final_conf
      use snap, only: ene_snap, stress_snap, fp_snap, &
         dim_ene_train_lml, dim_force_train_lml, dim_stress_train_lml, &
         train_rmse_energy, train_rmse_force, train_rmse_stress, &
         train_mae_energy, train_mae_force, train_mae_stress
      use main_mld_mod, only: md_mld_compute_energy, md_mld_compute_force, md_mld_compute_stress
      use module_md_mld, only: md_allocate_mld_desc
      use my_mpi_subroutines, only: subworlds_allreduce_vect_double
      use time_check_general, only: MY_MPI_WTIME
      use mld_subworld, only: subrank

      real(kind_double), intent(in) :: w_params(:,:)

      integer :: i, ix, ik
      real(kind_double) :: local_ref_energy, yb, yh
      real(kind_double) :: corr, detr, rmse, mae
      real(kind_double) :: t0, t1
      logical :: post_desc
      integer :: einp, finp, sinp

      ! Streaming accumulators: NSTAT sums × 4 quantity types
      real(kind_double) :: sums_e(NSTAT), sums_epa(NSTAT)
      real(kind_double) :: sums_f(NSTAT), sums_s(NSTAT)
      ! Single reduction buffer for one MPI allreduce (32 doubles)
      real(kind_double) :: reduce_buf(4 * NSTAT)

      _NAMECURRENT_("online_predict_train_errors")

      _MLD_BEGIN_

      sums_e   = 0.0_kind_double
      sums_epa = 0.0_kind_double
      sums_f   = 0.0_kind_double
      sums_s   = 0.0_kind_double

      ! Open per-subworld scratch buffers. Every subworld master
      ! (subrank==0) streams its own configs here; the buffers are
      ! merged into the final single output files after the loop.
      einp = -1; finp = -1; sinp = -1
      if (dim_ene_train_lml > 0 .and. subrank == 0) &
         open(newunit=einp, status='scratch', action='readwrite', form='formatted')
      if (dim_force_train_lml > 0 .and. subrank == 0) &
         open(newunit=finp, status='scratch', action='readwrite', form='formatted')
      if (dim_stress_train_lml > 0 .and. subrank == 0) &
         open(newunit=sinp, status='scratch', action='readwrite', form='formatted')

      t0 = MY_MPI_WTIME()

      do i = i_start_conf, i_final_conf
         if (.not. config_real(i)%train) cycle

         call test_if_config_is_small(i)
         call calc_neighbours(i)
         call md_allocate_mld_desc(fp_snap)
         post_desc = .true.
         call compute_descriptors(i, post_desc)

         ! --- Energy ---
         if (config_real(i)%has_energy) then
            call md_mld_compute_energy(i, pack_opt=.true.)
            local_ref_energy = config_real(i)%ref_energy
            yh = ene_snap
            yb = config_real(i)%energy(iread_energy) + local_ref_energy
            call accumulate_stats(sums_e, yb, yh)
            call accumulate_stats(sums_epa, yb / dble(config_real(i)%nat), &
                                             yh / dble(config_real(i)%nat))
            if (subrank == 0) write(einp, '(2es20.10," ",(a)," ",(a))') yh, yb, &
               config_real(i)%class, config_real(i)%filename
         end if

         ! --- Forces ---
         if ((config_real(i)%has_force) .and. desc_forces) then
            call md_mld_compute_force(i, pack_opt=.true.)
            do ik = 1, config_real(i)%nat
               do ix = 1, 3
                  call accumulate_stats(sums_f, config_real(i)%force(ix, ik), fp_snap(ix, ik))
               end do
            end do
            if (subrank == 0) then
               do ik = 1, config_real(i)%nat
                  do ix = 1, 3
                     write(finp, '(2es20.10," ",(a)," ",(a))') fp_snap(ix, ik), &
                        config_real(i)%force(ix, ik), config_real(i)%class, config_real(i)%filename
                  end do
               end do
            end if
         end if

         ! --- Stress ---
         if ((config_real(i)%has_stress) .and. desc_forces) then
            call md_mld_compute_stress(i, pack_opt=.true.)
            do ix = 1, 6
               call accumulate_stats(sums_s, config_real(i)%stress(ix), stress_snap(ix))
            end do
            if (subrank == 0) then
               do ix = 1, 6
                  write(sinp, '(2es20.10," ",(a)," ",(a))') stress_snap(ix), &
                     config_real(i)%stress(ix), config_real(i)%class, config_real(i)%filename
               end do
            end if
         end if

         call train_deallocate_desc(i)
      end do

      t1 = MY_MPI_WTIME()

      ! Pack accumulators and reduce across subworlds (32 doubles total)
      reduce_buf(1:NSTAT)             = sums_e
      reduce_buf(NSTAT+1:2*NSTAT)     = sums_epa
      reduce_buf(2*NSTAT+1:3*NSTAT)   = sums_f
      reduce_buf(3*NSTAT+1:4*NSTAT)   = sums_s
#if(PARA)
      call subworlds_allreduce_vect_double(reduce_buf)
#endif
      sums_e   = reduce_buf(1:NSTAT)
      sums_epa = reduce_buf(NSTAT+1:2*NSTAT)
      sums_f   = reduce_buf(2*NSTAT+1:3*NSTAT)
      sums_s   = reduce_buf(3*NSTAT+1:4*NSTAT)

      ! Merge each subworld master's scratch buffer into the single
      ! final output file (ordered, low-memory), then close scratch.
      if (dim_ene_train_lml > 0) then
         call merge_master_scratch(einp, "train_energy.out")
         if (subrank == 0) close(einp)
      end if
      if (dim_force_train_lml > 0) then
         call merge_master_scratch(finp, "train_force.out")
         if (subrank == 0) close(finp)
      end if
      if (dim_stress_train_lml > 0) then
         call merge_master_scratch(sinp, "train_stress.out")
         if (subrank == 0) close(sinp)
      end if

      ! Compute and print statistics from streaming sums
      if (rangml == 0) then
         write(6, '("ML:        ERROR TYPE                             R^2             D             RMSE          MAE    ")')
      end if

      train_rmse_energy = 0.0_kind_double; train_mae_energy = 0.0_kind_double
      train_rmse_force  = 0.0_kind_double; train_mae_force  = 0.0_kind_double
      train_rmse_stress = 0.0_kind_double; train_mae_stress = 0.0_kind_double

      if (dim_ene_train_lml > 0) then
         call stats_from_sums(sums_e, corr, detr, rmse, mae)
         if (rangml == 0) then
            write(6, '("ML: stat  Total_Energy            ",a5,":", 4f15.6)') 'train', corr, detr, rmse, mae
         end if
         call stats_from_sums(sums_epa, corr, detr, rmse, mae)
         if (rangml == 0) then
            write(6, '("ML: stat  Total_Energy_per_atom   ",a5,":", 4f15.6)') 'train', corr, detr, rmse, mae
         end if
         train_rmse_energy = rmse
         train_mae_energy = mae
      end if

      if (dim_force_train_lml > 0) then
         call stats_from_sums(sums_f, corr, detr, rmse, mae)
         if (rangml == 0) then
            write(6, '("ML: stat  Force_each_component    ",a5,":", 4f15.6)') 'train', corr, detr, rmse, mae
         end if
         train_rmse_force = rmse
         train_mae_force = mae
      end if

      if (dim_stress_train_lml > 0) then
         call stats_from_sums(sums_s, corr, detr, rmse, mae)
         if (rangml == 0) then
            write(6, '("ML: stat  Stress_virial           ",a5,":", 4f15.6)') 'train', corr, detr, rmse, mae
         end if
         train_rmse_stress = rmse
         train_mae_stress = mae
      end if

      call log_info("online_fit: train errors computed in " // vtoa(t1-t0) // " seconds")

      _MLD_END_
   end subroutine online_predict_train_errors

end module module_online_predict
