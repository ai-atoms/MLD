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


module module_compute_r_matrix
contains 
subroutine compute_r_matrix(i_start_at, i_final_at, d_n_neigh, d_kind_neigh, iconf)
  ! subroutine compute_r_matrix(i_start_at,i_final_at, iconf)

  USE module_kind_variables, ONLY: double
#ifdef MLD_NDM
  use gen_com_m, ONLY: A2cm, lperiod
  use gen_com_m_ml, ONLY: imm, bg, at
  use tab_imm_m_ml, ONLY: xp
#else
  use ondm_gen_com_m, ONLY: imm, A2cm, lperiod, bg, at
  use ondm_tab_imm_m, ONLY: iwmax2, xp
#endif
  use ml_in_ndm_module, ONLY: imm_neigh, desc_forces, rmat_dim, rangml, dmilady_dim, &
                              power_line, power_coeff_renorm, debug
  use module_neigh_local, only: r_cut 
  use derived_types, only: config_real, config_desc
  use module_neigh_local, only: preallocate_neigh_ja, build_local_neighbours_ja
  use module_chemical_species, only: img_weighted, img_num_ch, tnn_rdist, tcc01_rdist, tcc02_rdist, tcc03_rdist,  tii_rdist
  use time_check_general, only: debug_time, MY_MPI_WTIME
#ifdef MLD_NDM
  use notperiod_mod
#else
  use ondm_transform_coord, only: ondm_notperiod
#endif

  implicit none

  integer, intent(in)  :: i_start_at, i_final_at
  integer, dimension(imm), intent(out)   :: d_n_neigh
  integer, dimension(imm, imm_neigh), intent(out)    :: d_kind_neigh
  real(double), dimension(:, :), allocatable   :: r_matrix, tmp_dxp, tmp_xp               ! , factor_matrix
  real(double), dimension(:, :, :), allocatable   :: r_matrix_ch
  real(double), dimension(:), allocatable      :: inv_r_central, r_central, tmp_central, factor_central, factor_now
  real(double), dimension(:,:), allocatable :: fact_cent_ch, fact_now_ch 
  integer, dimension(:), allocatable     :: i_type, i_type_db, i_central, indx_central, indx_line
  integer, optional    :: iconf

  logical  :: small
  real(double), dimension(:, :), allocatable   :: xpnp
  real(double), dimension(3) :: dxp_ji, ds
  real(double) :: t00, t11, t22, t33, t44, t55, t66, t77, t88, t99 

  integer  :: ja, ia_n, iw2, indx_iline, iline, icol
  integer  :: iia, jja, ia_nn, ich
  double precision     :: r2_ji, r_ji
  double precision     :: factor_ia, factor_ja,  factor_line
  integer  :: max_neigh, max_neigh_local, max_lines_local, nn_g
  logical  :: desc_forces_local
  


  max_neigh = rmat_dim
  if ((i_start_at == 0) .and. (i_final_at == 0)) then
    d_n_neigh(:) = 0
    d_kind_neigh(:, :) = 0
    config_desc(iconf)%energy(:, :) = 0.d0
    ! config_desc(iconf)%force(:,:,:,:)=0.d0
    return
  end if

  desc_forces_local = desc_forces .and. (config_real(iconf)%has_force .or. config_real(iconf)%has_stress)
  
  small = .false.
  if (present(iconf)) then
    small = config_real(iconf)%small
  end if

  ALLOCATE (xpnp(3, imm))
  if (lperiod) then
    xpnp(:, :) = xp(:, :)
  else
#ifdef MLD_NDM
           call notperiod(imm,xp, xpnp,at,bg,.false.)
#else
    call ondm_notperiod(xp, xpnp)
#endif
  end if
  ! call cryst_to_cart (imm, xpnp, bg, -1)


  d_n_neigh(:) = 0
  d_kind_neigh(:, :) = 0
  config_desc(iconf)%energy(:, :) = 0.d0
  if (desc_forces_local) config_desc(iconf)%force(:, :, :, :) = 0.d0



  if (allocated(r_matrix)) deallocate (r_matrix); allocate (r_matrix(max_neigh, max_neigh))
  if (allocated(r_matrix_ch)) deallocate (r_matrix_ch); allocate (r_matrix_ch(max_neigh, max_neigh, img_num_ch))
  r_matrix(:, :) = 0.d0
  r_matrix_ch(:,:,:) = 0.d0


  call preallocate_neigh_ja(r_central, i_type, i_central, tmp_dxp, tmp_xp, imm_neigh)
  ! Loop on all atoms in the system, one image per atom
#ifdef MLD_NDM
  iw2=0
!!$ if (i_start_at == 1) iw2 = 0
!!$  if (i_start_at > 1) iw2 = iwmax2(i_start_at - 1)
  
#else
  if (i_start_at == 1) iw2 = 0
  if (i_start_at > 1) iw2 = iwmax2(i_start_at - 1)
#endif
  do ja = i_start_at, i_final_at
    ! begin small box or not 1/

    if (debug) then
      if (rangml==0) then
        if (mod(ja,20000)==1) then 
        write(*,*) 'ja ==', ja  
        end if 
      end if    
    end if    
    if (debug_time) then 
      t00 =MY_MPI_WTIME()
    end if 
 
    call build_local_neighbours_ja(iconf, ja, imm, xpnp, r_cut, d_n_neigh, d_kind_neigh, &
                                    r_central, i_type, i_type_db, i_central, tmp_dxp, tmp_xp, max_neigh_local, iw2)


    if (img_weighted) then 
      !config_real(iconf)%ntypes
      if (allocated(fact_cent_ch)) deallocate(fact_cent_ch) ; allocate(fact_cent_ch(0:max_neigh_local,img_num_ch))
      if (allocated(fact_now_ch)) deallocate(fact_now_ch) ; allocate(fact_now_ch(1:max_neigh_local,img_num_ch))
      do ich = 1, img_num_ch
        fact_cent_ch(0, ich) = config_real(iconf)%wspecies_per_type_ch(i_type(0), ich)
         do ia_n = 1, max_neigh_local
            fact_cent_ch(ia_n, ich) = config_real(iconf)%wspecies_per_type_ch(i_type(ia_n), ich)
         end do  
      end do   
    else    
      if (allocated(factor_central)) deallocate (factor_central); allocate (factor_central(0:max_neigh_local))
      if (allocated(factor_now)) deallocate (factor_now); allocate (factor_now(1:max_neigh_local))
      factor_ja=1.d0 
      factor_ia=1.d0
      factor_central(0) = 1.d0 
      do ia_n=1, max_neigh_local
         factor_central(ia_n) = factor_ia 
      end do   
  
    end if 
    if (allocated(tmp_central)) deallocate (tmp_central); allocate (tmp_central(max_neigh_local))
       

    nn_g = max_neigh

    !max_neigh_local  gives the order of the graph |G_j| = 1 + max_neigh_local
    !--------------------------------------------------------------- 
    !nn_g = gives the selection for the matrix = n_G
    !           1               2            ...        n_G 
    ! 1       ja:0-1          ja:0-2         ...     ja:0-n_G
    ! 2       ja:1-1          ja:1-2         ...     ja:1-n_G 
    ! :         :               :            ...          :  
    ! :         :               :            ...          :
    ! n_G  ja:(n_G-1)-1    ja:(n_G-1)-2      ...  ja:(n_G-1)-n_G  
    !---------------------------------------------------------------
  
    if (debug_time) then 
      t11 =MY_MPI_WTIME()
      tnn_rdist= tnn_rdist + t11 - t00    
    end if 

    ! here we draw the first line ... with respect the central atom. 
    if (allocated(indx_central)) deallocate (indx_central)
    allocate (indx_central(max_neigh_local))
    if (allocated(indx_line)) deallocate (indx_line)
    allocate (indx_line(max_neigh_local))
    if (allocated(inv_r_central)) deallocate(inv_r_central)
    allocate(inv_r_central(max_neigh_local))
    call indexx(max_neigh_local, r_central(1:max_neigh_local), indx_central(1:max_neigh_local))
    inv_r_central(1:max_neigh_local) = 1.d0/r_central(indx_central(1:max_neigh_local))

    if (debug_time) then 
      t22 =MY_MPI_WTIME()
      tii_rdist= tii_rdist + t22 - t11    
    end if 


    if (max_neigh_local < nn_g) then 
      if (img_weighted) then
        do ich=1, img_num_ch
          r_matrix_ch(1, max_neigh_local + 1:nn_g, ich) = 0.d0
          r_matrix_ch(1,1:max_neigh_local, ich) = fact_cent_ch(0, ich)*fact_cent_ch(indx_central(1:max_neigh_local), ich)*inv_r_central(1:max_neigh_local)**power_line
        end do 
      else   
        r_matrix(1, max_neigh_local + 1:nn_g) = 0.d0
        r_matrix(1,1:max_neigh_local) = factor_central(0)*factor_central(indx_central(1:max_neigh_local))*inv_r_central(1:max_neigh_local)**power_line
      end if 
    else
      if (img_weighted) then
        do ich=1, img_num_ch 
        r_matrix_ch(1,1:nn_G, ich) = fact_cent_ch(0,ich)*fact_cent_ch(indx_central(1:nn_G), ich)*inv_r_central(1:nn_G)**power_line
        end do 
      else    
        r_matrix(1,1:nn_G) = factor_central(0)*factor_central(indx_central(1:nn_G))*inv_r_central(1:nn_G)**power_line
      end if 
    end if  

    if (debug_time) then 
      t33 =MY_MPI_WTIME()
      tcc01_rdist= tcc01_rdist + t33 - t22    
    end if 


    max_lines_local = min(max_neigh_local, nn_g - 1)
    !!            31                 32             32
    ! write(*,*) 'max_lines_local, max_neigh_local, max_neigh, nn_g', max_lines_local, max_neigh_local, max_neigh, nn_g

    ! fill line by line from the second up to min(max_neigh_local, max_neigh - 1)
    do iline = 1, max_lines_local

      if (debug_time) then 
        t44 =MY_MPI_WTIME()
      end if 
  
      indx_iline = indx_central(iline)
      iia = i_central(indx_iline)
      factor_line = inv_r_central(indx_iline)**power_coeff_renorm

      ia_nn = 0
      tmp_central(:) = 0.d0
      do icol = 0, max_neigh_local
        jja = i_central(icol)
        if (indx_iline == icol) cycle
        dxp_ji(1:3) = tmp_xp(1:3, indx_iline) - tmp_xp(1:3, icol)
#ifdef MLD_NDM
!!$        if (small) then
#else
        if (small) then
#endif
          r2_ji = Sum(dxp_ji(1:3)**2)
#ifdef MLD_NDM
!!$        else
!!$          ds = MatMul(dxp_ji, bg)
!!$          WHERE ((ds .GT. 0.5d0) .OR. (ds .LT. -0.5d0))
!!$            ds(1:3) = ds(1:3) - Dble(Nint(ds(1:3)))
!!$          END WHERE
!!$          dxp_ji = MatMul(at, ds)/A2cm
!!$          r2_ji = Sum(dxp_ji(1:3)**2)
!!$        end if
#else
        else
          ds = MatMul(dxp_ji, bg)
          WHERE ((ds .GT. 0.5d0) .OR. (ds .LT. -0.5d0))
            ds(1:3) = ds(1:3) - Dble(Nint(ds(1:3)))
          END WHERE
          dxp_ji = MatMul(at, ds)/A2cm
          r2_ji = Sum(dxp_ji(1:3)**2)
        end if
#endif
        r_ji = dsqrt(r2_ji)
        ia_nn = ia_nn + 1
        tmp_central(ia_nn) = r_ji
        if (img_weighted) then 
          do ich =1, img_num_ch
            fact_now_ch(ia_nn, ich) = fact_cent_ch(icol,ich)
          end do   
        else 
          factor_now(ia_nn) = factor_central(icol) !wtf 
        end if     


        if (tmp_central(ia_nn) == 0.d0) then
          write (*, *) iline, iia, icol, jja
          write (*, *) tmp_xp(1:3, iline)*A2cm
          write (*, *) tmp_xp(1:3, icol)*A2cm
          write (*, *) dxp_ji(:)
          write (*, *) dsqrt(Sum(dxp_ji(1:3)**2)), small
          stop
        end if

      end do                  ! icol= 0, max_neigh_local
      if (ia_nn /= max_neigh_local) then
        stop 'serious trouble in compute_rdist_matrix'
      end if
      if (debug_time) then 
        t55 =MY_MPI_WTIME()
        tcc02_rdist= tcc02_rdist + t55 - t44      
      end if 

      call indexx(max_neigh_local, tmp_central(1:max_neigh_local), indx_line(1:max_neigh_local))
      !write(*,*) 'index1', indx_line (1:max_neigh_local)
      !write(*,*) 'tmp1', tmp_central(indx_line(1:5))
      !now do ii =1,max_neigh_local
      !now   indx_line(ii) = ii
      !now end do   
      !now call dlasrt2('I', max_neigh_local,tmp_central, indx_line, info)
      !write(*,*) 'index2', indx_line(1: max_neigh_local)
      !write(*,*) 'tmp2', tmp_central(1:5)
      !write(*,*) '<Full STOP>'
      !stop 
      if (debug_time) then 
        t66 =MY_MPI_WTIME()
        tii_rdist= tii_rdist + t66 - t55      
      end if 

      if (img_weighted) then 
        do ich=1, img_num_ch
          if (max_neigh_local < nn_g) then
            r_matrix_ch(iline + 1, max_neigh_local + 1:nn_G, ich) = 0.d0
            r_matrix_ch(iline + 1, 1:max_neigh_local,ich) = factor_line*fact_cent_ch(indx_iline,ich)*fact_now_ch(indx_line(1:max_neigh_local),ich)/tmp_central(indx_line(1:max_neigh_local))**power_line
          else 
            r_matrix_ch(iline + 1, 1:nn_g, ich) = factor_line*fact_cent_ch(indx_iline, ich)*fact_now_ch(indx_line(1:nn_g),ich)/tmp_central(indx_line(1:nn_g))**power_line
          end if 
        end do 
      else   
        if (max_neigh_local < nn_g) then
          r_matrix(iline + 1, max_neigh_local + 1:nn_G) = 0.d0
          r_matrix(iline + 1, 1:max_neigh_local) = factor_line*factor_central(indx_iline)*factor_now(indx_line(1:max_neigh_local))/tmp_central(indx_line(1:max_neigh_local))**power_line
          !now r_matrix(iline + 1, 1:max_neigh_local) = factor_line*factor_central(indx_iline)*factor_now(indx_line(1:max_neigh_local))/tmp_central(1:max_neigh_local)**power_line
        else 
          r_matrix(iline + 1, 1:nn_g) = factor_line*factor_central(indx_iline)*factor_now(indx_line(1:nn_g))/tmp_central(indx_line(1:nn_g))**power_line
          !now r_matrix(iline + 1, 1:nn_g) = factor_line*factor_central(indx_iline)*factor_now(indx_line(1:nn_g))/tmp_central(1:nn_g)**power_line
        end if 
      end if 
      if (debug_time) then 
        t77 =MY_MPI_WTIME()
        tcc02_rdist= tcc02_rdist + t77 - t66 
      end if 

    end do                  ! iline=1,max_neigh_local

    if (debug_time) then 
      t88 =MY_MPI_WTIME()     
    end if 

    if (img_weighted) then
      do ich=1, img_num_ch
        r_matrix(:,:) = r_matrix_ch(:,:,ich)
        config_desc(iconf)%channel((ich-1)*dmilady_dim+1:ich*dmilady_dim, ja) = reshape(transpose(r_matrix(:,:)), (/dmilady_dim/)) 
      end do   
    else  
      config_desc(iconf)%energy(:, ja) = reshape(transpose(r_matrix), (/dmilady_dim/))
    end if 
    if (debug_time) then 
      t99 =MY_MPI_WTIME()
      tcc03_rdist= tcc03_rdist + t99 - t88      
    end if 


  end do                  ! ja main loop

  deallocate (xpnp)

end subroutine compute_r_matrix
end module module_compute_r_matrix
