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

#include "../MLD_MACROS.INC"

module module_write_xml
  use iso_fortran_env, only: dp => real64
  implicit none
  private

  public :: write_radialace_to_xml, write_basecnlm_to_xml
   
  contains 

    subroutine write_radialace_to_xml(filename)
        use module_ace_radial, only: radialace, radial_spline, ACE_CHEM_RADIAL_RALF, ACE_CHEM_RADIAL_HSVD, ACE_CHEM_RADIAL_BLOCK_HSVD
        use module_ace_desc, only:  base_params, zetaace_order, ace_numax, ace_dim, delta_zetaace, delta_ace, ace_gencg, & 
                                    ace_radial_poly, ace_chem, ace_chem_low_rank
        use module_chemical_species, only: fix_no_of_elements
        use module_chem_compress, only: chem_compressor
        character(len=*), intent(in) :: filename
        integer :: unit, i, j, qq, ich, ss 
    
        open(newunit=unit, file=filename, status='replace', action='write')
    
        ! Write XML header
        write(unit, '(A)') '<?xml version="1.0"?>'


        !$! integer :: zetaace_order 
        !$! integer, parameter :: MAX_ZETAACE_ORDER = 6 
        !$! real(kind_double), dimension(:,:), allocatable   :: delta_zetaace 
        !$! ! the way in which the CG is generated: 1 for draft way 2 Dusson-Ortner method. Default 2. 
        !$! integer :: ace_gencg
        !$! integer, parameter :: GENCG_DRAFT = 1, GENCG_DUSORT = 2
        !$! real(kind_double), dimension(:,:), allocatable :: delta_ace 

        !-----------------writting the ace_general_params-----------------
        write(unit, '(A)') '<ace_general_params>'
        ! Write scalar parameters
        write(unit, '(A,I0,A)') '  <ace_numax>', ace_numax, '</ace_numax>'
        write(unit, '(A,I0,A)') '  <ace_dim>', ace_dim, '</ace_dim>'
        write(unit, '(A,I0,A)') '  <zetaace_order>', zetaace_order, '</zetaace_order>'
        write(unit, '(A,I0,A)') '  <ace_gencg>', ace_gencg, '</ace_gencg>'
        write(unit, '(A,I0,A)') '  <ace_radial_poly>', ace_radial_poly, '</ace_radial_poly>'
        write(unit, '(A,I0,A)') '  <ace_chem>', ace_chem, '</ace_chem>'
        write(unit, '(A,I0,A)') '  <fix_no_of_elements>', fix_no_of_elements, '</fix_no_of_elements>'
        ! Write delta_zetaace matrix

        if (zetaace_order == 1 ) then 
          if (allocated(delta_zetaace)) deallocate(delta_zetaace)
           allocate(delta_zetaace(1, ace_numax))
           delta_zetaace(1, :) = 1.0_dp
        end if 
        write(unit, '(A)') '  <delta_zetaace>'
        do i = 1, zetaace_order
            write(unit, '(A)', advance="no") '    <row>'
            do j = 1, ace_numax
                write(unit, '(F16.8, A)', advance="no") delta_zetaace(i, j), ' '
            end do
            write(unit, '(A)') '</row>'
        end do
        write(unit, '(A)') '  </delta_zetaace>'
        ! Write delta_ace matrix
        write(unit, '(A)') '  <delta_ace>'
        do i = 1, fix_no_of_elements
            write(unit, '(A)', advance="no") '    <row>'
            do j = 1, fix_no_of_elements
                write(unit, '(F16.8, A)', advance="no") delta_ace(i, j), ' '
            end do
            write(unit, '(A)') '</row>'
        end do
        write(unit, '(A)') '  </delta_ace>'
        ! Close XML file
        write(unit, '(A)') '</ace_general_params>'



        write(unit, '(A)') '<basis_ace_list>'
        do i = 1, size(base_params)
            write(unit, '(A)') '  <basis_ace>'
            write(unit, '(A,I0,A)') '    <nmax>', base_params(i)%nmax, '</nmax>'
            write(unit, '(A,I0,A)') '    <lmax>', base_params(i)%lmax, '</lmax>'
            write(unit, '(A,I0,A)') '    <mumax>', base_params(i)%mumax, '</mumax>'
            write(unit, '(A,I0,A)') '    <kmax>', base_params(i)%kmax, '</kmax>'
            write(unit, '(A,I0,A)') '    <active>', merge(1, 0, base_params(i)%active), '</active>'
            write(unit, '(A,1PE15.8,A)') '    <r_cut_in>', base_params(i)%r_cut_in, '</r_cut_in>'
            write(unit, '(A,1PE15.8,A)') '    <r_cut_out>', base_params(i)%r_cut_out, '</r_cut_out>'
            write(unit, '(A,1PE15.8,A)') '    <r_cut_width_in>', base_params(i)%r_cut_width_in, '</r_cut_width_in>'
            write(unit, '(A,1PE15.8,A)') '    <r_cut_width_out>', base_params(i)%r_cut_width_out, '</r_cut_width_out>'
            write(unit, '(A,1PE15.8,A)') '    <lambda>', base_params(i)%lambda, '</lambda>'
            write(unit, '(A)') '  </basis_ace>'
        end do
        write(unit, '(A)') '</basis_ace_list>'


        write(unit, '(A)') '<ace_radial>'
        
        ! Write radialace properties
        write(unit, '(A,I0,A)') '  <nmax>', radialace%nmax, '</nmax>'
        write(unit, '(A,I0,A)') '  <lmax>', radialace%lmax, '</lmax>'
        write(unit, '(A,I0,A)') '  <dim_mu>', radialace%dim_mu, '</dim_mu>'
        write(unit, '(A,I0,A)') '  <dim_rad>', radialace%dim_rad, '</dim_rad>'
        write(unit, '(A,I0,A)') '  <type_f_radial>', radialace%type_f_radial, '</type_f_radial>'
        write(unit, '(A,I0,A)') '  <type_chem_radial>', radialace%type_chem_radial, '</type_chem_radial>'
        write(unit, '(A,I0,A)') '  <type_fcut_in>', radialace%type_fcut_in, '</type_fcut_in>'
        write(unit, '(A,I0,A)') '  <type_fcut_out>', radialace%type_fcut_out, '</type_fcut_out>'
        write(unit, '(A,I0,A)') '  <npoints>', radialace%npoints, '</npoints>'
        write(unit, '(A,F12.5,A)') '  <r_cut_out>', radialace%r_cut_out, '</r_cut_out>'
        write(unit, '(A,F12.5,A)') '  <r_cut_width_out>', radialace%r_cut_width_out, '</r_cut_width_out>'
        write(unit, '(A,F12.5,A)') '  <r_cut_in>', radialace%r_cut_in, '</r_cut_in>'
        write(unit, '(A,F12.5,A)') '  <r_cut_width_in>', radialace%r_cut_width_in, '</r_cut_width_in>'    
        ! Write general parameters
        !write(unit, '(A,I0,A)') '  <ACE_CHEM_RADIAL_RALF>', ACE_CHEM_RADIAL_RALF, '</ACE_CHEM_RADIAL_RALF>'
        !write(unit, '(A,I0,A)') '  <ACE_CHEM_RADIAL_HSVD>', ACE_CHEM_RADIAL_HSVD, '</ACE_CHEM_RADIAL_HSVD>'
        !write(unit, '(A,I0,A)') '  <ACE_CHEM_RADIAL_BLOCK_HSVD>', ACE_CHEM_RADIAL_BLOCK_HSVD, '</ACE_CHEM_RADIAL_BLOCK_HSVD>'


        write(unit, '(A)') '</ace_radial>'



        ! The compressed representation is used whenever the chemical compressor
        ! holds the data, i.e. HSVD low-rank (ace_chem_low_rank==1). In that case
        ! radial_spline(:) is allocated but its per-spline coefficient arrays were
        ! never filled, so writing them would dereference unallocated components.
        ! The random-projection contraction (ace_radial_chem==5) instead fills the
        ! full radial_spline(:) array like standard HSVD and leaves the compressor
        ! uninitialized, so it correctly falls through to the full-spline writer.
        ! Gate on the compressor state instead of the mode.
        if (chem_compressor%initialized) then
          ! ---- Write compressed prototype splines u^{kl}_q and mixing matrices A ----
          write(unit, '(A)') '<compressed_radial_splines>'
          write(unit, '(A,I0,A)') '  <S>', chem_compressor%S, '</S>'
          write(unit, '(A,I0,A)') '  <Q>', chem_compressor%Q, '</Q>'
          write(unit, '(A,I0,A)') '  <kmax>', chem_compressor%kmax, '</kmax>'
          write(unit, '(A,I0,A)') '  <lmax_p1>', chem_compressor%lmax_p1, '</lmax_p1>'
          write(unit, '(A,I0,A)') '  <n_channels>', chem_compressor%n_channels, '</n_channels>'
          write(unit, '(A,I0,A)') '  <Nr>', chem_compressor%Nr, '</Nr>'

          ! Write mixing matrices A(S, Q, ich) for each channel
          write(unit, '(A)') '  <mixing_matrices>'
          do ich = 1, chem_compressor%n_channels
            write(unit, '(A,I0,A)') '    <channel ich="', ich, '">'
            do ss = 1, chem_compressor%S
              write(unit, '(A)', advance='no') '      <row>'
              do qq = 1, chem_compressor%Q
                write(unit, '(ES15.7E3,1X)', advance='no') chem_compressor%A(ss, qq, ich)
              end do
              write(unit, '(A)') '</row>'
            end do
            write(unit, '(A)') '    </channel>'
          end do
          write(unit, '(A)') '  </mixing_matrices>'

          ! Write prototype splines u_spline(Q, n_channels)
          write(unit, '(A)') '  <prototype_splines>'
          do ich = 1, chem_compressor%n_channels
            do qq = 1, chem_compressor%Q
              write(unit, '(A,I0,A,I0,A)') '    <u_spline q="', qq, '" ich="', ich, '">'
              write(unit, '(A,I0,A)') '      <n>', chem_compressor%u_spline(qq, ich)%n, '</n>'
              write(unit, '(A,1PE15.8,A)') '      <delta>', chem_compressor%u_spline(qq, ich)%delta, '</delta>'
              write(unit, '(A)') '      <coefficients>'
              ! Coefficients a
              write(unit, '(A)', advance='no') '        <a>'
              do j = 1, size(chem_compressor%u_spline(qq, ich)%a)
                write(unit, '(ES15.7E3,1X)', advance='no') chem_compressor%u_spline(qq, ich)%a(j)
              end do
              write(unit, '(A)') '</a>'
              ! Coefficients b
              write(unit, '(A)', advance='no') '        <b>'
              do j = 1, size(chem_compressor%u_spline(qq, ich)%b)
                write(unit, '(ES15.7E3,1X)', advance='no') chem_compressor%u_spline(qq, ich)%b(j)
              end do
              write(unit, '(A)') '</b>'
              ! Coefficients c
              write(unit, '(A)', advance='no') '        <c>'
              do j = 1, size(chem_compressor%u_spline(qq, ich)%c)
                write(unit, '(ES15.7E3,1X)', advance='no') chem_compressor%u_spline(qq, ich)%c(j)
              end do
              write(unit, '(A)') '</c>'
              ! Coefficients d
              write(unit, '(A)', advance='no') '        <d>'
              do j = 1, size(chem_compressor%u_spline(qq, ich)%d)
                write(unit, '(ES15.7E3,1X)', advance='no') chem_compressor%u_spline(qq, ich)%d(j)
              end do
              write(unit, '(A)') '</d>'
              write(unit, '(A)') '      </coefficients>'
              write(unit, '(A)') '    </u_spline>'
            end do
          end do
          write(unit, '(A)') '  </prototype_splines>'
          write(unit, '(A)') '</compressed_radial_splines>'

        else
          ! ---- Write full radial_spline array (S² splines) ----
          write(unit, '(A)') '<radial_splines>'
          do i = 1, size(radial_spline)
              write(unit, '(A)') '  <radial_spline>'
              write(unit, '(A,I0,A)') '    <n>', radial_spline(i)%n, '</n>'
              write(unit, '(A,1PE15.8,A)') '    <delta>', radial_spline(i)%delta, '</delta>'
              write(unit, '(A)') '    <coefficients>'
              write(unit, '(A)', advance='no') '      <a>'
              do j = 1, size(radial_spline(i)%a)
                  write(unit, '(ES15.7E3,1X)', advance='no') radial_spline(i)%a(j)
              end do
              write(unit, '(A)') '</a>'
              write(unit, '(A)', advance='no') '      <b>'
              do j = 1, size(radial_spline(i)%b)
                  write(unit, '(ES15.7E3,1X)', advance='no') radial_spline(i)%b(j)
              end do
              write(unit, '(A)') '</b>'
              write(unit, '(A)', advance='no') '      <c>'
              do j = 1, size(radial_spline(i)%c)
                  write(unit, '(ES15.7E3,1X)', advance='no') radial_spline(i)%c(j)
              end do
              write(unit, '(A)') '</c>'
              write(unit, '(A)', advance='no') '      <d>'
              do j = 1, size(radial_spline(i)%d)
                  write(unit, '(ES15.7E3,1X)', advance='no') radial_spline(i)%d(j)
              end do
              write(unit, '(A)') '</d>'
              write(unit, '(A)') '    </coefficients>'
              write(unit, '(A)') '  </radial_spline>'
          end do
          write(unit, '(A)') '</radial_splines>'
        end if


        write(unit, '(A)') '<ace_parameters>'
        write(unit, '(A)') '</ace_parameters>'
        close(unit)
    end subroutine write_radialace_to_xml


    subroutine write_basecnlm_to_xml(filename)
      use module_base_cnlm, only: base_cnlm
      ! use module_ace_radial, only: radialace, radial_spline, ACE_CHEM_RADIAL_RALF, ACE_CHEM_RADIAL_HOME
      use module_ace_desc, only:  base_params, l_ace_order
      ! use module_ace_desc, only:  base_params, zetaace_order, ace_numax, ace_dim, delta_zetaace, delta_ace, ace_gencg
      ! use module_chemical_species, only: fix_no_of_elements
      character(len=*), intent(in) :: filename
      integer :: unit, i, j, ii, jj
  
      open(newunit=unit, file=filename, status='replace', action='write')
  
      ! Write XML header
      write(unit, '(A)') '<?xml version="1.0"?>'

      !!!! INFO !!!!
      !? type(type_base_dico_mu),   dimension(:),  allocatable:: dico_mu
      !? type(type_base_dico_nn),   dimension(:),  allocatable:: dico_nb  
      ! type(type_base_dico_ll),   dimension(:),  allocatable:: dico_lb0 -> M0mat -> dimr_M0mat
      ! type(type_base_dico_mu_nnll), dimension(:),  allocatable :: dico_mub_nb_lb0
      ! type(type_base_dico_mu_nnll_LLi), dimension(:), allocatable :: dicoB
      ! type(type_dico_cg), dimension(:), allocatable  :: dicoCG 
      ! integer, dimension(:,:), allocatable :: uniqueA_tuples
      ! integer, dimension(:),  allocatable  :: map_large_to_uniqueA 


      ! NECESSARY ??
      ! map_large_to_uniqueA, uniqueA_tuples, dicoCG, dicoB

      !-----------------writting the base_cnlm params-----------------
      write(unit, '(A)') '<base_cnlm_list>'
      do i = 2, size(base_params) !allocate(base_params(ace_numax))
       if (l_ace_order(i)) then
        write(unit, '(A)') '  <base_cnlm>'
            ! write(unit, '(A,I0,A)') '    <nu>', i, '</nu>'
          ! write(unit, '(A,I0,A)') '    <active>', merge(1, 0, base_params(i)%active), '</active>'
          ! if (base_params(i)%active == .True.) then

            !!!!!!! type(type_base_dico_ll),   dimension(:),  allocatable:: dico_lb0
            write(unit, '(A,I0,A)') '      <num_classes_l>', size(base_cnlm(i)%dico_lb0), '</num_classes_l>'
            write(unit, '(A)') '    <dimr_M0mat>'
            do j = 1, size(base_cnlm(i)%dico_lb0)
              write(unit, '(I0, 1X)', advance="no") size(base_cnlm(i)%dico_lb0(j)%M0mat, 1)
            end do
            write(unit, '(A)') '    </dimr_M0mat>'

            !!!!!!! type(type_base_dico_mu_nnll), dimension(:),  allocatable :: dico_mub_nb_lb0
            !   integer :: idx_mub
            !   integer :: idx_nb
            !   integer :: idx_lb0
            !   real(dp), dimension(:,:), allocatable :: rpicg
            !   integer :: idx_cg 
            ! write(unit, '(A)') '    <dico_mub_nb_lb0>'
              write(unit, '(A,I0,A)') '      <dim_munl>', size(base_cnlm(i)%dico_mub_nb_lb0), '</dim_munl>'
              write(unit, '(A)') '      <idx_mub>'
              do j = 1, size(base_cnlm(i)%dico_mub_nb_lb0) !dim_munl
                write(unit, '(I0, 1X)', advance="no") base_cnlm(i)%dico_mub_nb_lb0(j)%idx_mub
              end do
              write(unit, '(A)') '      </idx_mub>'

              write(unit, '(A)') '      <idx_nb>'
              do j = 1, size(base_cnlm(i)%dico_mub_nb_lb0)
                write(unit, '(I0, 1X)', advance="no") base_cnlm(i)%dico_mub_nb_lb0(j)%idx_nb
              end do
              write(unit, '(A)') '      </idx_nb>'

              write(unit, '(A)') '      <idx_lb0>'
              do j = 1, size(base_cnlm(i)%dico_mub_nb_lb0)
                write(unit, '(I0, 1X)', advance="no") base_cnlm(i)%dico_mub_nb_lb0(j)%idx_lb0
              end do
              write(unit, '(A)') '      </idx_lb0>'

              write(unit, '(A)') '      <idx_cg>'
              do j = 1, size(base_cnlm(i)%dico_mub_nb_lb0)
                write(unit, '(I0, 1X)', advance="no") base_cnlm(i)%dico_mub_nb_lb0(j)%idx_cg
              end do
              write(unit, '(A)') '      </idx_cg>'
            ! write(unit, '(A)') '    </dico_mub_nb_lb0>'

            !!!!!!! type(type_base_dico_mu_nnll_LLi), dimension(:), allocatable :: dicoB
            ! integer :: idx_nb       ! this map into dico_nb 
            ! integer :: idx_mub      ! this map into dico_mu
            ! integer :: idx_lb0      ! this map into dico_lb0
            ! integer :: idx_munblb0  ! this map into dico_mu_nnll 
            ! integer :: idx_Li       ! this map into dico_mub_nb_lb0%Lmat(i,:) or dico_lb0%cg(:,i)  or dico_mu_nnll%rpicg(:,i)
            ! write(unit, '(A)') '    <dicoB>'
              write(unit, '(A,I0,A)') '      <dim_baseB>', size(base_cnlm(i)%dicoB), '</dim_baseB>' !dim_baseB = size(base_cnlm(nu)%dicoB)
              write(unit, '(A)') '      <idx_Li>'
              do j = 1, size(base_cnlm(i)%dicoB)
                write(unit, '(I0, 1X)', advance="no") base_cnlm(i)%dicoB(j)%idx_Li
              end do
              write(unit, '(A)') '      </idx_Li>'
            ! write(unit, '(A)') '    </dicoB>'

            !!!!!!! ttype(type_dico_cg), dimension(:), allocatable  :: dicoCG
            ! real(dp), dimension(:,:), allocatable :: cg  
            ! integer :: dimr_cg, dimc_cg
            ! write(unit, '(A)') '    <dicoCG>'
              write(unit, '(A,I0,A)') '      <dim_baseCG>', size(base_cnlm(i)%dicoCG), '</dim_baseCG>'
              write(unit, '(A)') '      <dimr_cg>'
              do j = 1, size(base_cnlm(i)%dicoCG)
                write(unit, '(I0, 1X)', advance="no") size(base_cnlm(i)%dicoCG(j)%cg, 1)
              end do
              write(unit, '(A)') '      </dimr_cg>' 

              write(unit, '(A)') '      <dimc_cg>'
              do j = 1, size(base_cnlm(i)%dicoCG)
                write(unit, '(I0, 1X)', advance="no") size(base_cnlm(i)%dicoCG(j)%cg, 2)
              end do
              write(unit, '(A)') '      </dimc_cg>' 

              write(unit, '(A)') '      <cg>'
              do j = 1, size(base_cnlm(i)%dicoCG)
                do ii = 1, size(base_cnlm(i)%dicoCG(j)%cg, 1) 
                  do jj = 1, size(base_cnlm(i)%dicoCG(j)%cg, 2)
                    write(unit, '(1PE15.8, 1X)', advance='no')  base_cnlm(i)%dicoCG(j)%cg(ii, jj)
                  end do
                end do
              end do
              write(unit, '(A)') '      </cg>'

              ! allocate(this%dicoCG(ii)%cg(size(this%dico_lb0(ii)%cg,1), size(this%dico_lb0(ii)%cg,2)))
            ! write(unit, '(A)') '    </dicoCG>'

            !!!!!!! integer, dimension(:,:), allocatable :: uniqueA_tuples
            ! allocate(this%uniqueA_tuples(4, M * this%nu))
            ! write(unit, '(A,I0,A)') '    <dim_unique_tuples>', int(size(base_cnlm(i)%uniqueA_tuples, 2)/i), '</dim_unique_tuples>' !M
              write(unit, '(A,I0,A)') '    <dim_unique_tuples>', size(base_cnlm(i)%uniqueA_tuples, 2), '</dim_unique_tuples>'
              write(unit, '(A)') '    <uniqueA_tuples>'
            do j = 1, 4
              do ii = 1, size(base_cnlm(i)%uniqueA_tuples, 2)
                write(unit, '(I0, 1X)', advance="no") base_cnlm(i)%uniqueA_tuples(j, ii)
              end do
            end do
            write(unit, '(A)') '    </uniqueA_tuples>'

            !!!!!!! integer, dimension(:),  allocatable  :: map_large_to_uniqueA
            ! allocate(this%map_large_to_uniqueA(M*this%nu))
            write(unit, '(A)') '    <map_large_to_uniqueA>'
            do j = 1, size(base_cnlm(i)%map_large_to_uniqueA)
              write(unit, '(I0, 1X)', advance="no") base_cnlm(i)%map_large_to_uniqueA(j)
            end do
            write(unit, '(A)') '    </map_large_to_uniqueA>'

          ! end if
        write(unit, '(A)') '  </base_cnlm>'
       end if
      end do
      write(unit, '(A)') '</base_cnlm_list>'
      close(unit)
  end subroutine write_basecnlm_to_xml

end module module_write_xml 



module module_write_design_matrix
  use temporary_data_cov, ONLY: dim_xdesc, dim_xdesc_patch
  implicit none 
  contains

  subroutine init_dump_desgin_matrix_test (dunit, iunit)
      use ml_in_ndm_module, only: rangml
      use module_kernel_2b, only: activate_k2b, dim_kernel_2b
      use snap, only : Amat
      use module_ml_scalapack, only :  scalapack_driver, dimc_sca_Amat, dimr_sca_Amat
      use temporary_data_cov, only: dim_data_test
      implicit none 
      integer, intent(inout) :: dunit, iunit
      integer :: nline, ncol, ndata, nk2b
      character(len=80)    :: CHFMT, CHFMT2

      if (activate_k2b) then 
        nk2b = dim_kernel_2b
      else 
        nk2b = 0 
      end if     

      if (scalapack_driver) then
        ! we dump  the transpose of Amat 
        !nline = dimc_sca_Amat
        nline = dim_data_test
        ncol = dimr_sca_Amat
        ndata = dimc_sca_Amat
      else 
        !nline = size(Amat, 2)
        nline = dim_data_test
        ncol = size(Amat, 1)
        !ndata = size(Amat, 2)
        ndata = dim_data_test
      end if

    

    if (rangml == 0) then 
      open (file='design_matrix_test.dat', newunit=dunit, action='write', status='unknown')
      write (CHFMT, *) '( ', ncol, 'e30.18 )'
      write (dunit, '("# design matrix: n_config x dim_descriptor ")')
      write (dunit, '("# design matrix: 1st colum number of atoms for energy_desc, 0 for others")')
      write (dunit, '("# nline ncol nk2b dimD", 4i9)') nline , ncol, nk2b, dim_xdesc
    end if 


    if (rangml == 0 ) then
      open (file='design_matrix_test.info', newunit=iunit, action='write', status='unknown')
      write (CHFMT2, *) '( i6, e30.18, e30.18, a15 )'
      write (iunit, '("# info design matrix: 1st column : 1 if that line comes from energy, 2 from force and 3 from stress ")')
      write (iunit, '("# info design matrix: 2nd column : the ymat target value ")')
      write (iunit, '("# info design matrix: 3nd column : the diagonal part of W weigth matrix ")')
      write (iunit, '("# info design matrix: 4rd column : the original file from where the line comes ")')
      write (iunit, '("# nline ncol", 2i9)') nline, 4
    end if 

  end subroutine init_dump_desgin_matrix_test

  subroutine dump_design_matrix_test (ic, dunit, iunit, dim_wmat, wmat, wtags, wyy)
    use mpi 
    use mld_mpi
    use ml_in_ndm_module, only: rangml 
    use module_kind_variables, only: kind_double
    use derived_types, only: config_real
    use mld_subworld, only: list_of_masters_in_mld, nb_subworlds
    use mld_logger
    implicit none
    integer, intent(in) :: ic, dunit, iunit, dim_wmat 
    real(kind_double), dimension(:,:), intent(in) :: wmat
    integer, dimension(:), intent(in) :: wtags
    real(kind_double), dimension(:), intent(in) :: wyy
    integer  :: ii, info, iproc, iw 
    character(len=80)    :: CHFMT, CHFMT2
    real(kind_double), dimension(:,:), allocatable :: r_wmat
    integer, dimension(:), allocatable :: r_wtags
    real(kind_double), dimension(:), allocatable :: r_wyy
    integer :: r_dim_wmat, icount, r_ic 
  
    ! The design matrix is the transpose of Amat:  Amat (dim_desc, number_of_config).
    _NAMECURRENT_("dump_design_matrix_test")
    _MLD_BEGIN_
    !$! write(*,*) 'ini' , rangml 
    write (CHFMT, *) '( ', size(wmat,1), 'e30.18 )'
    write (CHFMT2, *) '( i6, e30.18, e30.18, a15 )'

    !write(*,*) 'ddebug1', rangml, subrank, dim_wmat, size(wtags)
    ! do ii = 0, nb_subworlds-1
    !   write(*,*) ii, rangml, size(wmat,1), size(wmat,2)
    ! end do  
    icount = 0 

      if (rangml /= 0 ) then 
        do ii = 0, nb_subworlds-1
            iproc = list_of_masters_in_mld(ii)
            if (iproc == 0 ) cycle  
            if  (rangml == iproc ) then 
              icount = icount + 1 
              call MPI_SEND(dim_wmat, 1, MPI_INTEGER, 0, 100000+iproc, mpi_comm_mld, info)
              call MPI_SEND(ic, 1, MPI_INTEGER, 0, 500000+iproc, mpi_comm_mld, info)
              call MPI_SEND(wmat, size(wmat,1)*dim_wmat, MPI_DOUBLE_PRECISION, 0, 200000+iproc, mpi_comm_mld, info)
              call MPI_SEND(wtags,             dim_wmat,          MPI_INTEGER, 0, 300000+iproc, mpi_comm_mld, info)
              call MPI_SEND(wyy,               dim_wmat, MPI_DOUBLE_PRECISION, 0, 400000+iproc, mpi_comm_mld, info)
              
            end if 
        end do 
      end if 

      if (rangml == 0 ) then
        icount = 0 
        do ii = 0, nb_subworlds-1
          iproc = list_of_masters_in_mld(ii)
          if (iproc == 0 ) then 
            do iw = 1, dim_wmat
              write (dunit, CHFMT) wmat(:, iw)
            end do
            do iw = 1, dim_wmat
              write (iunit, CHFMT2) wtags(iw),  wyy(iw), 1.d0, &
                "  "//config_real(ic)%class//"_"//config_real(ic)%klm//"_"//config_real(ic)%cnumber
            end do 
          else 
            icount = icount + 1
            call MPI_RECV(r_dim_wmat, 1, MPI_INTEGER, iproc, 100000+iproc, mpi_comm_mld, MPI_STATUS_IGNORE, info)
            call MPI_RECV(r_ic,       1, MPI_INTEGER, iproc, 500000+iproc, mpi_comm_mld, MPI_STATUS_IGNORE, info)
            
            if (allocated(r_wmat)) deallocate(r_wmat) ; allocate(r_wmat(size(wmat,1), r_dim_wmat))
            if (allocated(r_wtags)) deallocate(r_wtags) ; allocate(r_wtags(r_dim_wmat))
            if (allocated(r_wyy)) deallocate(r_wyy) ; allocate(r_wyy(r_dim_wmat))
            call MPI_RECV(r_wmat,   size(wmat,1)*r_dim_wmat, MPI_DOUBLE_PRECISION, iproc, 200000+iproc, mpi_comm_mld, MPI_STATUS_IGNORE, info)
            call MPI_RECV(r_wtags,               r_dim_wmat,          MPI_INTEGER, iproc, 300000+iproc, mpi_comm_mld, MPI_STATUS_IGNORE, info)
            call MPI_RECV(r_wyy,                 r_dim_wmat, MPI_DOUBLE_PRECISION, iproc, 400000+iproc, mpi_comm_mld, MPI_STATUS_IGNORE, info)
              
            do iw = 1, r_dim_wmat
              write (dunit, CHFMT) r_wmat(:, iw)
            end do  
            do iw = 1, r_dim_wmat
              write (iunit, CHFMT2) r_wtags(iw),  r_wyy(iw), 1.d0, &
                "  "//config_real(r_ic)%class//"_"//config_real(r_ic)%klm//"_"//config_real(r_ic)%cnumber
            end do
          end if 
        end do 
      end if
     call MPI_BARRIER(mpi_comm_mld, info)




    _MLD_END_ 
  end subroutine dump_design_matrix_test


  subroutine dump_design_matrix
    !use mpi 
    use mld_mpi
    use ml_in_ndm_module, only: rangml 
    use module_kind_variables, only: kind_double
    use snap, only: Amat, ymat, fit_snap
    use derived_types, only: config_real
    use module_ml_scalapack, only: sca_Amat, dimr_sca_Amat, dimc_sca_Amat, &
                                   desc_sca_Amat, &
                                   context,  &
                                   scalapack_driver
    use module_scalapack_tools, only: allocation_scalapack_matrix
    use module_kernel_2b, only: activate_k2b, dim_kernel_2b                               
    use mld_logger
    implicit none
    integer  :: dunit, dunit2, ii, jj, tag, ndata, nline, ncol, info 
    character(len=80)    :: CHFMT, CHFMT2
    ! The design matrix is the transpose of Amat:  Amat (dim_desc, number_of_config).
    ! this is dumped in the file design_matrix.dat
    ! there exists also design_matrix.info that contains the information about the
    ! origin of each line of the design matrix: 
    ! tag (1,2,3) for energy, force and stress ;
    ! the y target value ;
    ! the diagonal part of the weight matrix ;
    ! the original file from where the line comes.  
    integer :: size_pack , lld , npack, size_pack_last, icount, rsize_pack, nk2b 
    integer, dimension(9) :: desc_local_matrix
    real(kind_double), dimension(:,:), allocatable  :: local_matrix
    _NAMECURRENT_("dump_design_matrix")
    _MLD_BEGIN_

    if (activate_k2b) then 
      nk2b = dim_kernel_2b
    else 
      nk2b = 0 
    end if     

    if (scalapack_driver) then
      ! we dump  the transpose of Amat 
      nline = dimc_sca_Amat
      ncol = dimr_sca_Amat
      ndata = dimc_sca_Amat
    else 
      nline = size(Amat, 2)
      ncol = size(Amat, 1)
      ndata = size(Amat, 2)
    end if

    

    if (rangml == 0) then 
      open (file='design_matrix.dat', newunit=dunit, action='write', status='unknown')
      write (CHFMT, *) '( ', ncol, 'es30.18e3 )'
      write (dunit, '("# design matrix: n_config x dim_descriptor ")')
      write (dunit, '("# design matrix: 1st colum number of atoms for energy_desc, 0 for others")')
      write (dunit, '("# nline ncol nk2b dimD", 4i9)') nline , ncol, nk2b, dim_xdesc
    end if 

    if (scalapack_driver) then

      ! local matrix local_matrix(dimr_sca_AmatT, size_pack))  
      size_pack=100 
      if (nline < size_pack) size_pack = ncol 
      ! Process 0 prepares the local matrix

      ! Allocate the local matrix to hold the first size_pack columns
      !if (rangml==0) then 
          allocate(local_matrix(dimr_sca_Amat, size_pack))
      !end if 

      ! Local leading dimension
      lld = max(1, dimr_sca_Amat)
      ! Create a scalapack descriptor for the local_matrix
      call descinit(desc_local_matrix, dimr_sca_Amat, size_pack, dimr_sca_Amat, size_pack, 0, 0, context, lld, info)
      if (info /= 0) then
          call log_critical("ML critical: error in descinit in "//NAMECURRENT//" with info = "//vtoa(info))  
          stop "with error in descinit"
      endif

      npack = dimc_sca_Amat/size_pack
      size_pack_last = dimc_sca_Amat - npack*size_pack
      icount = 0
      do ii = 1, dimc_sca_Amat, size_pack 
        if (ii/size_pack < npack) rsize_pack = size_pack
        if (ii/size_pack ==  npack)  then 
          rsize_pack = size_pack_last
          if (size_pack_last == 0) cycle
        end if   
  
        call pdgemr2d(dimr_sca_Amat, rsize_pack, sca_Amat, 1, ii, desc_sca_Amat, local_matrix, 1, 1, desc_local_matrix, context, info)
        if ((info /= 0).and.(rangml==0)) then
            call log_critical("ML critical: error in pdgemr2d in "//NAMECURRENT//" with info = "//vtoa(info))
            stop "with error in pdgemr2d"
        endif

        if (rangml==0) then
          do jj = 1, rsize_pack
            write (dunit, CHFMT) local_matrix(:, jj)
          end do
        end if  
          
      end do             ! ii on dimc_sca_Amat, size_pack
      
      deallocate(local_matrix)
      if (rangml==0) then 
        close (dunit)
      end if

    else 
    
      if (rangml == 0) then
        do ii = 1, size(Amat, 2)
          ! write (*,*) fit_snap(ii)%energy, fit_snap(ii)%force, fit_snap(ii)%stress
          write (dunit, CHFMT) Amat(:, ii)
          !write(*,*) 'aa', Amat(:, 1)
          !stop 
        end do
        close (dunit)
      end if
    end if

    if (rangml == 0 ) then
        open (file='design_matrix.info', newunit=dunit2, action='write', status='unknown')
        write (CHFMT2, *) '( i6, e30.18, e30.18, a15 )'
        write (dunit2, '("# info design matrix: 1st column : 1 if that line comes from energy, 2 from force and 3 from stress ")')
        write (dunit2, '("# info design matrix: 2nd column : the ymat target value ")')
        write (dunit2, '("# info design matrix: 3nd column : the diagonal part of W weigth matrix ")')
        write (dunit2, '("# info design matrix: 4rd column : the original file from where the line comes ")')
        write (dunit2, '("# nline ncol", 2i9)') size(ymat, 1), 4
        do ii = 1,ndata 
          if (fit_snap(ii)%energy) tag = 1
          if (fit_snap(ii)%force) tag = 2
          if (fit_snap(ii)%stress) tag = 3
          write (dunit2, CHFMT2) tag,  ymat(ii, 1), sqrt(fit_snap(ii)%weight), &
            "  "//config_real(fit_snap(ii)%iconf)%class//"_"//config_real(fit_snap(ii)%iconf)%klm//"_"//config_real(fit_snap(ii)%iconf)%cnumber
        end do 
        close (dunit2)
    end if  
    _MLD_END_ 
  end subroutine dump_design_matrix

end module module_write_design_matrix



module module_write_parameters

  use ml_in_ndm_module, only: rangml, debug, char_desc, &
                              mld_order, mld_linear, mld_quadratic, &
                              mld_kernel, mld_polyc, &
                              regularization_name, &
                              mld_regularization_type, mld_regularization_type_home, &
                              descriptor_type, descriptor_afs, &
                              descriptor_g2, descriptor_bispectrum_so4, descriptor_pow_so4, descriptor_pow_so3, &
                              descriptor_pow_so3_3body, &
                              lsoap_diag, lsoap_norm, lsoap_lnorm, alpha_soap, n_soap, l_max, alpha_soap, &
                              eta_max_g2, eta_min_g2, n_g2_eta, n_g2_rs, &
                              j_max, weighted,  &
                              weighted, weighted_3ch, weighted_auto, fix_weighted_for_element, fix_weighted_for_element_3ch, &
                              renorm_mass, renorm_cov, &
                              l_max, descriptor_ftnbody, descriptor_ace, descriptor_tbind
  use module_chemical_species, only: periodic_table_element, fix_type_to_periodic, fix_no_of_elements

  use module_kernel, only: dim_kernel, global_kernel, kernel_type, &
                           kernel_random, kernel_phase_random, &
                           kernel_po, kernel_po_scaled, sigma_kernel, length_kernel, kernel_power, &
                           kernel_random_po, basis_random_po
  use module_mld_quadratic, only: dim_xdesc_quadratic 
  use module_bispectrum_so4, only: lbso4_diag                        
  use module_so3, only: n_rbf_so3, radial_pow_so3
  use snap, only: w_params, dim_xdesc_linear
  use module_afs, only: n_rbf_afs, n_cheb, afs_type
  use temporary_data_cov, ONLY: dim_xdesc, dim_xdesc_patch
  use module_neigh_local, only: r_cut, r_cut_width, r_cut_in, r_cut_width_in, type_fcut
  use module_db_poscar, ONLY: fix_ref_energy_per_element
  use mld_logger
  use module_optimization, only: optimize_weights_db, optimize_weights_chem
  use module_ftnbody, only: delta_rff, length_rff, dim_qbody, ftnbody_omega_2b, ftnbody_omega_3b, ftnbody_omega_4b, ftnbody_omega_5b, &
                            ftnbody_phase_random_2b, ftnbody_phase_random_3b, ftnbody_phase_random_4b, ftnbody_phase_random_5b, &
                            r_cut_ft2b, r_cut_width_ft2b, r_cut_ft3b, r_cut_width_ft3b, &
                            r_cut_ft4b, r_cut_width_ft4b, r_cut_ft5b, r_cut_width_ft5b
  use module_body_desc, only: l_body_order, dim_desc_body, bond_beta, bond_dist_ann, bond_dist_transform 
  use module_zbl, only: zbl_potential, zbl_type, zbl_mode_alone, zbl_mode_default_k2b, r1_zbl, r2_zbl, type_rac_zbl, params_k2b_to_zbl
  use module_kernel_2b, only: activate_k2b, dim_kernel_2b, zpoints_2b, zr_2b, zr_fcut_2b, sigma_2b, delta_2b, & 
                              np_radial_2b, sparse_points_2b_build, r_cut_2b, r_cut_width_2b
  use module_write_xml, only: write_radialace_to_xml, write_basecnlm_to_xml                   

  integer  :: no_of_channels, size_params_to_write, size_supplementary_info_desc, icol, irow, &
              size_zbl_activate_k2b, icount_all

  

contains

  subroutine build_name_of_parameters_file (mld_name, lammps_name)
    implicit none 
    character(len=*), intent(out)  :: mld_name, lammps_name
    character(len=5)     :: tmp_name
    _NAMECURRENT_("build_name_of_parameters_file")
    if (mld_order == mld_linear) tmp_name = "snap1"
    if (mld_order == mld_kernel) tmp_name = "kernl"
    if (mld_order == mld_quadratic) tmp_name = "qnml1"
    if (mld_order == mld_polyc) tmp_name = "polyc"

    if (mld_regularization_type == mld_regularization_type_home) then
      mld_name = trim(char_desc//'_'//tmp_name//'_params.pot_'//regularization_name)
      lammps_name=trim('lammps_'//char_desc//'_'//tmp_name//'_params.pot_'//regularization_name)
    else
      mld_name = trim(char_desc//'_'//tmp_name//'_params.pot')
      lammps_name=trim('lammps_'//char_desc//'_'//tmp_name//'_params.pot')
    end if


  end subroutine build_name_of_parameters_file

  subroutine write_snap_parameters
    use module_ftnbody_potio, only: write_ftnbody_random_xml
    use module_tbind_potio, only: write_tbind_xml
    implicit none
    integer  :: n1_lsoap, n2_lsoap, n3_lsoap
    integer  :: i, ii, jj, ik, itmp
    character(len=80)    :: CHFMT
    integer :: unitpot, unitlammps
    integer :: n_active_patch 
    character(len=120) :: mld_name, lammps_name

  _NAMECURRENT_("write_snap_parameters")



    _MLD_BEGIN_
    !writing only on proc 0
    if (.not. (optimize_weights_db.or.optimize_weights_chem)) then
      if (rangml == 0) then
        write (6, *) 'ML: Writing parameters in write_snap_parameters'
      end if
    end if


    if (rangml == 0) then

      if (mld_order == mld_kernel) then
        if (kernel_type == kernel_random) then
          size_params_to_write = size(w_params, 1) + dim_kernel*dim_xdesc + dim_kernel
        else if (kernel_type == kernel_po) then
          size_params_to_write = size(w_params, 1) + dim_kernel*dim_xdesc + 3
          !write(*,*) size_params_to_write, kernel_type
          !stop 
        else if (kernel_type == kernel_po_scaled) then
          size_params_to_write = size(w_params, 1) + dim_kernel*dim_xdesc + 3  
          !write(*,*) size_params_to_write, kernel_type
          !stop 
        else if (kernel_type == kernel_random_po) then
          itmp = 0
          do ii = 1, dim_kernel
            itmp = itmp + basis_random_po(ii)%dim_omega
          end do
          size_params_to_write = size(w_params, 1) + itmp*dim_xdesc + dim_kernel + 3
        else
          size_params_to_write = size(w_params, 1) + dim_kernel*dim_xdesc
        end if
      else
        size_params_to_write = size(w_params, 1)
      end if

      ! ftnbody random features (omega/phase) are now written to random.xml, so the
      ! .pot stays weights-only (dim_pot == dim_xdesc_full, read back compactly).
      size_supplementary_info_desc  = 0
      size_params_to_write = size_params_to_write + size_supplementary_info_desc

      size_zbl_activate_k2b = 0 
      !TODO sparse_points_build ... 
      sparse_points_2b_build=1
      if (zbl_potential .and. zbl_type == zbl_mode_alone) then
        ! zbl_mode_alone: zbl_flag + zbl_type + r1_zbl + r2_zbl + r_cut_in + r_cut_width_in + type_fcut + r_cut_width
        size_zbl_activate_k2b = 8
      else if ((zbl_potential .and. zbl_type == zbl_mode_default_k2b) .or. activate_k2b) then 
        ! a: 12+5*dim_kernel_2b + np_radial_2b   z:17+iijj+5*dim_kernel_2b + np_radial_2b  
        size_zbl_activate_k2b = 13+5*dim_kernel_2b + np_radial_2b  &
                                 - dim_kernel_2b ! double counting from weights ...  
        if (zbl_potential)  then 
          size_zbl_activate_k2b = 18+5*dim_kernel_2b + np_radial_2b   +  &
                                 size(params_k2b_to_zbl,1)*size(params_k2b_to_zbl,2) &
                                - dim_kernel_2b ! double counting fromweigths 
        end if   
      end if    
      size_params_to_write = size_params_to_write + size_zbl_activate_k2b


      call build_name_of_parameters_file(mld_name, lammps_name)

      !$! tmp_name = "snap1"
      !$! if (mld_order == mld_linear) tmp_name = "snap1"
      !$! if (mld_order == mld_kernel) tmp_name = "kernl"
      !$! if (mld_order == mld_quadratic) tmp_name = "qnml1"
      !$! if (mld_order == mld_polyc) tmp_name = "polyc"
      !$! if (mld_regularization_type == mld_regularization_type_home) then
      !$!   open (file=char_desc//'_'//tmp_name//'_params.pot_'//regularization_name, newunit=unitpot, action='write')
      !$!   open (file='lammps_'//char_desc//'_'//tmp_name//'_params.pot_'//regularization_name, newunit=unitlammps, action='write')
      !$! else
      !$!   open (file=char_desc//'_'//tmp_name//'_params.pot', newunit=unitpot, action='write')
      !$!   open (file='lammps_'//char_desc//'_'//tmp_name//'_params.pot', newunit=unitlammps, action='write')
      !$! end if

      open (file=trim(mld_name), newunit=unitpot, action='write')
      open (file=trim(lammps_name), newunit=unitlammps, action='write')

      write (unitlammps, '("# nchannels nelements element_1 weight_1_for_el_1 weight_2_for_el_1 element_2 weight_1_for_el_2 weight_2_for_el_2 ... etc ")')
      write (unitlammps, '("# descriptor type, nparameters_descriptor, ndetails")')
      write (unitlammps, '("# mass_1 mass_2 ... mass_nelements")')
      write (unitlammps, '("# detail_1 is the mld_order, ..., detail_ndetails (e.g. mld_order no_of_radial_channels no_of_angular_channels)")')
      write (unitlammps, '("# r_cut only one for all descriptors")')
      write (unitlammps, '("# w_i  with i=1, dim_xdesc, dim_xdesc +1 ...  (nchannel-1)*dim_xdesc ... nchannel*dim_xdesc")')

      ! write (CHFMT,*)'(i3, 1x, ',int(fix_no_of_elements),'(a2, " ")')'

      ! Writing the line
      ! no_of_channels type_of_desc element_1 W_1 W_2 element2 W_1 W_2
      ! IFORT_STYLE
      if ((weighted) .and. (weighted_3ch)) then
        no_of_channels = 3
      else if ((weighted) .and. (.not. (weighted_3ch))) then
        no_of_channels = 2
      else
        no_of_channels = 1
      end if
      write (CHFMT, *) '( i7,i7,', fix_no_of_elements, '(  ', ' a4, e20.10, e20.10, e20.10', '))'
      if ((weighted) .and. (.not. weighted_3ch)) then
        !if (weighted_auto) then
        !   write (unitlammps,CHFMT) no_of_channels, fix_no_of_elements, ((periodic_table_element(fix_type_to_periodic(ii))%symbol, fix_weighted_for_element(ii), 1.0 ), ii=1,fix_no_of_elements)
        !else
        !   write (unitlammps,CHFMT) no_of_channels, fix_no_of_elements, ((periodic_table_element(fix_type_to_periodic(ii))%symbol, fix_weighted_for_element(ii), 1.0), ii=1,fix_no_of_elements)
        !end if
        write (unitlammps, '(i7,i7)', advance='no') no_of_channels, fix_no_of_elements
        if (fix_no_of_elements == 1) then
          write (unitlammps, '(a4, e20.10, e20.10, e20.10)') periodic_table_element(fix_type_to_periodic(1))%symbol, fix_weighted_for_element(1), 1.d0, &
            fix_ref_energy_per_element(1)
        else
          do ii = 1, fix_no_of_elements - 1
            write (unitlammps, '(a4, e20.10, e20.10, e20.10)', advance='no') periodic_table_element(fix_type_to_periodic(ii))%symbol, fix_weighted_for_element(ii), 1.d0, &
              fix_ref_energy_per_element(ii)
          end do
        end if
        if (fix_no_of_elements > 1) then
          ii = fix_no_of_elements
          write (unitlammps, '(a4, e20.10, e20.10, e20.10)') periodic_table_element(fix_type_to_periodic(ii))%symbol, fix_weighted_for_element(ii), 1.d0, &
            fix_ref_energy_per_element(ii)
        end if
      end if

      if (weighted_3ch) then
        !write (CHFMT,*) '( i7,',fix_no_of_elements,'(  ',' a4, e20.10, e20.10','))'
        !if (weighted_auto) then
        !write (unitlammps,CHFMT) no_of_channels, fix_no_of_elements, ((periodic_table_element(fix_type_to_periodic(ii))%symbol, fix_weighted_for_element(ii), fix_weighted_for_element_3ch(ii)  ), ii=1,fix_no_of_elements)
        write (unitlammps, '(i7,i7)', advance='no') no_of_channels, fix_no_of_elements
        if (fix_no_of_elements == 1) then
          write (unitlammps, '(a4, e20.10, e20.10, e20.10 )') periodic_table_element(fix_type_to_periodic(1))%symbol, fix_weighted_for_element(1), fix_weighted_for_element_3ch(1), &
            fix_ref_energy_per_element(1)
        else
          do ii = 1, fix_no_of_elements - 1
            write (unitlammps, '(a4, e20.10, e20.10, e20.10)', advance='no') periodic_table_element(fix_type_to_periodic(ii))%symbol, fix_weighted_for_element(ii), &
              fix_weighted_for_element_3ch(ii), fix_ref_energy_per_element(ii)
          end do
        end if
        if (fix_no_of_elements > 1) then
          ii = fix_no_of_elements
          write (unitlammps, '(a4, e20.10, e20.10, e20.10)') periodic_table_element(fix_type_to_periodic(ii))%symbol, fix_weighted_for_element(ii), fix_weighted_for_element_3ch(ii), &
            fix_ref_energy_per_element(ii)
        end if
        !else
        !write (unitlammps,CHFMT) no_of_channels, fix_no_of_elements, ((periodic_table_element(fix_type_to_periodic(ii))%symbol, fix_weighted_for_element(ii), fix_weighted_for_element_3ch(ii)), ii=1,fix_no_of_elements)
        !end if
      end if
      if ((.not. (weighted)) .and. (.not. weighted_3ch)) then
      !!!write (unitlammps,CHFMT) no_of_channels, fix_no_of_elements, ((periodic_table_element(fix_type_to_periodic(ii))%symbol, 1.0,1.0), ii=1,fix_no_of_elements)
        write (unitlammps, '(i7,i7)', advance='no') no_of_channels, fix_no_of_elements
        if (fix_no_of_elements == 1) then
          write (unitlammps, '(a4, e20.10, e20.10, e20.10)') periodic_table_element(fix_type_to_periodic(1))%symbol, 1.d0, 1.d0, fix_ref_energy_per_element(1)
        else
          do ii = 1, fix_no_of_elements - 1
            write (unitlammps, '(a4, e20.10, e20.10, e20.10)', advance='no') periodic_table_element(fix_type_to_periodic(ii))%symbol, 1.d0, 1.d0, fix_ref_energy_per_element(ii)
          end do
        end if

        if (fix_no_of_elements > 1) then
          ii = fix_no_of_elements
          write (unitlammps, '(a4, e20.10, e20.10, e20.10)') periodic_table_element(fix_type_to_periodic(ii))%symbol, 1.d0, 1.d0, fix_ref_energy_per_element(ii)
        end if
      end if

      ! Writing the line
      select case (descriptor_type)

      case (descriptor_g2)
        if (fix_no_of_elements /= 1) then
          write (6, *) 'WARNING G2 NOT YET IMPLEMENTED FOR fix_no_of_elements > 1', fix_no_of_elements
        end if
        if (weighted) then
          !write (unitlammps,'(i4, i6, i4, f20.10)') descriptor_type, size_params_to_write, 4, factor_weight_mass
          write (unitlammps, '(i4, i12, i4)') descriptor_type, size_params_to_write, 5
        else
          write (unitlammps, '(i4, i12, i4)') descriptor_type, size_params_to_write, 5
        end if
        write (unitlammps, '(f20.4)') periodic_table_element(fix_type_to_periodic(fix_no_of_elements))%mass
        write (unitlammps, '(i5, 2f20.6, 2i5)') mld_order, eta_max_g2, eta_min_g2, n_g2_eta, n_g2_rs

      case (descriptor_afs)
        if (mld_order == mld_kernel) then
          write (unitlammps, '(i4, i12, i4 )') descriptor_type, size_params_to_write, 5 + 3
        else
          write (unitlammps, '(i4, i12, i4 )') descriptor_type, size_params_to_write, 5
        end if
        ! line 9: mass_1, mass_2, ...
        write (CHFMT, *) '(', fix_no_of_elements, 'f20.4)'
        write (unitlammps, CHFMT) (periodic_table_element(fix_type_to_periodic(ii))%mass, ii=1, fix_no_of_elements)
        !write (unitlammps,'(f20.4)') periodic_table_element(fix_type_to_periodic(fix_no_of_elements))%mass
        n_active_patch=0 
        if (activate_k2b.or.zbl_potential) then
          n_active_patch=1 
        end if   
        if (mld_order == mld_kernel) then
          write (unitlammps, '(4i5, 4i8)') mld_order, afs_type, n_rbf_afs, n_cheb, n_active_patch, kernel_type, dim_xdesc, dim_kernel
        else
          write (unitlammps, '(5i5)') mld_order, afs_type, n_rbf_afs, n_cheb, n_active_patch
        end if

      case (descriptor_pow_so3)
        ! line 8: descriptor type, nparameters, ndetails
        if (mld_order == mld_kernel) then
          write (unitlammps, '(i4, i12, i4)') descriptor_type, size_params_to_write, 5 + 3
        else
          write (unitlammps, '(i4, i12, i4)') descriptor_type, size_params_to_write, 5
        end if
        ! line 9: mass_1, mass_2, ...
        write (CHFMT, *) '(', fix_no_of_elements, 'f20.4)'
        write (unitlammps, CHFMT) (periodic_table_element(fix_type_to_periodic(ii))%mass, ii=1, fix_no_of_elements)
        !write (unitlammps,'(f20.4)') periodic_table_element(fix_type_to_periodic(fix_no_of_elements))%mass

        n_active_patch=0 
        if (activate_k2b.or.zbl_potential) then
          n_active_patch=1 
        end if   

        if (mld_order == mld_kernel) then
          write (unitlammps, '(i5, 7i8)') mld_order, l_max, n_rbf_so3, radial_pow_so3, n_active_patch, kernel_type, dim_xdesc, dim_kernel
        else
          write (unitlammps, '(i5, 7i8)') mld_order, l_max, n_rbf_so3, radial_pow_so3, n_active_patch
        end if

      case (descriptor_pow_so3_3body)
        ! line 8: descriptor type, nparameters, ndetails
        if (mld_order == mld_kernel) then
          write (unitlammps, '(i4, i12, i4)') descriptor_type, size_params_to_write, 5 + 3
        else
          write (unitlammps, '(i4, i12, i4)') descriptor_type, size_params_to_write, 5
        end if
        ! line 9: mass_1, mass_2, ...
        write (CHFMT, *) '(', fix_no_of_elements, 'f20.4)'
        write (unitlammps, CHFMT) (periodic_table_element(fix_type_to_periodic(ii))%mass, ii=1, fix_no_of_elements)
        !write (unitlammps,'(f20.4)') periodic_table_element(fix_type_to_periodic(fix_no_of_elements))%mass
        n_active_patch=0 
        if (activate_k2b.or.zbl_potential) then
          n_active_patch=1 
        end if   

        if (mld_order == mld_kernel) then
          write (unitlammps, '(i5, 7i8)') mld_order, l_max, n_rbf_so3, radial_pow_so3, n_active_patch, kernel_type, dim_xdesc, dim_kernel
        else
          write (unitlammps, '(i5, 7i8)') mld_order, l_max, n_rbf_so3, radial_pow_so3, n_active_patch
        end if



      case (descriptor_pow_so4)
        ! line 8: descriptor type, nparameters, ndetails
        if (mld_order == mld_kernel) then
          write (unitlammps, '(i4, i12, i4)') descriptor_type, size_params_to_write, 3 + 3
        else
          write (unitlammps, '(i4, i12, i4)') descriptor_type, size_params_to_write, 3
        end if
        ! line 9: mass_1, mass_2, ...
        write (CHFMT, *) '(', fix_no_of_elements, 'f20.4)'
        write (unitlammps, CHFMT) (periodic_table_element(fix_type_to_periodic(ii))%mass, ii=1, fix_no_of_elements)
        !write (unitlammps,'(f20.4)') periodic_table_element(fix_type_to_periodic(fix_no_of_elements))%mass
        n_active_patch=0 
        if (activate_k2b.or.zbl_potential) then
          n_active_patch=1 
        end if   

        if (mld_order == mld_kernel) then
          write (unitlammps, '(i5, f5.1, 3i8)') mld_order, j_max, n_active_patch, kernel_type, dim_xdesc, dim_kernel
        else
          write (unitlammps, '(i5, f5.1)') mld_order, j_max, n_active_patch
        end if

      case (descriptor_bispectrum_so4)
        ! line 8: descriptor type, nparameters, ndetails
        if (mld_order == mld_kernel) then
          write (unitlammps, '(i4, i12, i4)') descriptor_type, size_params_to_write, 4 + 3
        else
          write (unitlammps, '(i4, i12, i4)') descriptor_type, size_params_to_write, 4
        end if
        ! line 9: mass_1, mass_2, ...
        write (CHFMT, *) '(', fix_no_of_elements, 'f20.4)'
        write (unitlammps, CHFMT) (periodic_table_element(fix_type_to_periodic(ii))%mass, ii=1, fix_no_of_elements)
        !line 10: details for descriptor ... for bso4 mld_order -linear, non_linear / j_max / diag or not
        n1_lsoap = 0; if (lbso4_diag) n1_lsoap = 1
        n_active_patch=0 
        if (activate_k2b.or.zbl_potential) then
          n_active_patch=1 
        end if   

        if (mld_order == mld_kernel) then
          write (unitlammps, '(i5, f5.1,i4,4i6)') mld_order, j_max, n1_lsoap, n_active_patch, kernel_type, dim_xdesc, dim_kernel
        else
          write (unitlammps, '(i5, f5.1,2i4)') mld_order, j_max, n1_lsoap, n_active_patch
        end if

      case (descriptor_ace)  
        ! line 8: descriptor type, nparameters, ndetails
        if (mld_order == mld_kernel) then
          write (unitlammps, '(i4, i12, i4)') descriptor_type, size_params_to_write, 2 + 3
        else
          write (unitlammps, '(i4, i12, i4)') descriptor_type, size_params_to_write, 2
        end if
        ! line 9: mass_1, mass_2, ...
        write (CHFMT, *) '(', fix_no_of_elements, 'f20.4)'
        write (unitlammps, CHFMT) (periodic_table_element(fix_type_to_periodic(ii))%mass, ii=1, fix_no_of_elements)

        !line 10: details for descriptor ... for  mld_order -linear, non_linear 
        n_active_patch=0 
        if (activate_k2b.or.zbl_potential) then
          n_active_patch=1 
        end if   

        if (mld_order == mld_kernel) then
          write (unitlammps, '(i5, 4i6)') mld_order, n_active_patch, kernel_type, dim_xdesc, dim_kernel
        else
          write (unitlammps, '(i5, i6)') mld_order,  n_active_patch
        end if

        call write_radialace_to_xml('radialace.xml')
        call write_basecnlm_to_xml('basecnlmace.xml')

      case (descriptor_tbind)
        ! Same compact layout as ACE: the .pot carries only the general model
        ! information and the weights; everything needed to rebuild the
        ! descriptor (parameters + HSVD radial splines) goes to tbind.xml.
        ! line 8: descriptor type, nparameters, ndetails
        if (mld_order == mld_kernel) then
          write (unitlammps, '(i4, i12, i4)') descriptor_type, size_params_to_write, 2 + 3
        else
          write (unitlammps, '(i4, i12, i4)') descriptor_type, size_params_to_write, 2
        end if
        ! line 9: mass_1, mass_2, ...
        write (CHFMT, *) '(', fix_no_of_elements, 'f20.4)'
        write (unitlammps, CHFMT) (periodic_table_element(fix_type_to_periodic(ii))%mass, ii=1, fix_no_of_elements)

        ! line 10: details for descriptor ... for mld_order -linear, non_linear
        n_active_patch=0
        if (activate_k2b.or.zbl_potential) then
          n_active_patch=1
        end if

        if (mld_order == mld_kernel) then
          write (unitlammps, '(i5, 4i6)') mld_order, n_active_patch, kernel_type, dim_xdesc, dim_kernel
        else
          write (unitlammps, '(i5, i6)') mld_order,  n_active_patch
        end if

        call write_tbind_xml('tbind.xml')

      case (descriptor_ftnbody)
        ! line 8: descriptor type, nparameters, ndetails
        ! if (mld_order == mld_kernel) then
        !   write (unitlammps, '(i4, i12, i4)') descriptor_type, size_params_to_write, 28 + 3
        ! else
        !   write (unitlammps, '(i4, i12, i4)') descriptor_type, size_params_to_write, 28
        ! end if
        if (mld_order == mld_kernel) then
          write (unitlammps, '(i4, i12, i4)') descriptor_type, size_params_to_write, 36 + 3
        else
          write (unitlammps, '(i4, i12, i4)') descriptor_type, size_params_to_write, 36
        end if
        ! line 9: mass_1, mass_2, ...
        write (CHFMT, *) '(', fix_no_of_elements, 'f20.4)'
        write (unitlammps, CHFMT) (periodic_table_element(fix_type_to_periodic(ii))%mass, ii=1, fix_no_of_elements)
        
        do ii = 1, size(l_body_order,1) 
          if (.not.(l_body_order(ii))) dim_desc_body(ii) = 0 
        end do   

        n_active_patch=0 
        if (activate_k2b.or.zbl_potential) then
          n_active_patch=1 
        end if   

        if (mld_order == mld_kernel) then
          write (unitlammps, '(i5, 2f10.4, i5, 10i7, 5f10.4, 5f20.10, 3i6, 11f10.4, i5)') mld_order, bond_beta, bond_dist_ann, bond_dist_transform, dim_qbody(1:5), &
                                                                            dim_desc_body(1:5), delta_rff(1:5), length_rff(1:5), &
                                                                            kernel_type, dim_xdesc, dim_kernel , &
                                                                            r_cut_in, r_cut_width , r_cut_width_in, &
                                                                            r_cut_ft2b, r_cut_width_ft2b, r_cut_ft3b, r_cut_width_ft3b, &
                                                                            r_cut_ft4b, r_cut_width_ft4b, r_cut_ft5b, &
                                                                            r_cut_width_ft5b, n_active_patch
        else
          write (unitlammps, '(i5, 2f10.4, i5, 10i7, 5f10.4, 5f20.10, 11f10.4, i5)') mld_order, bond_beta, bond_dist_ann, bond_dist_transform, dim_qbody(1:5), &
                                                                            dim_desc_body(1:5), delta_rff(1:5), length_rff(1:5), &
                                                                            r_cut_in, r_cut_width , r_cut_width_in, &
                                                                            r_cut_ft2b, r_cut_width_ft2b, r_cut_ft3b, r_cut_width_ft3b, &
                                                                            r_cut_ft4b, r_cut_width_ft4b, r_cut_ft5b, r_cut_width_ft5b, n_active_patch
        end if
        
        

      end select
      ! old_version write (unitlammps,'(f20.4)') r_cut*2.d0
      ! line 10: rcut ...
      write (unitlammps, '(f20.4)') r_cut


      if  (mld_order == mld_kernel) then
        write (unitpot, *) size(w_params, 1)
      else
        write (unitpot, *) size_params_to_write
      end if

      icount_all = 0 
      do i = 1, dim_xdesc_linear
        icount_all = icount_all + 1 
        write (unitpot, '(ES25.15E3)') w_params(i, 1)
        write (unitlammps, '(ES25.15E3)') w_params(i, 1)
      end do


      ! init   supplementary info  descriptor ..........................
      if (size_supplementary_info_desc > 0) then
        select case(descriptor_type)

          case(descriptor_ftnbody)

            if (l_body_order(2)) then
              do icol = 1, size(ftnbody_omega_2b,2)
                do irow = 1, size(ftnbody_omega_2b,1)
                  icount_all = icount_all + 1  
                  write(unitlammps, '(ES25.15E3)')   ftnbody_omega_2b(irow, icol)
                  write(unitpot, '(ES25.15E3)')   ftnbody_omega_2b(irow, icol)
                end do 
              end do 
              do icol = 1, size(ftnbody_phase_random_2b,1)
                  icount_all = icount_all + 1 
                  write(unitlammps, '(ES25.15E3)')   ftnbody_phase_random_2b(icol)
                  write(unitpot, '(ES25.15E3)')   ftnbody_phase_random_2b(icol)
              end do       
            end if   
              
            if (l_body_order(3)) then
              do icol = 1, size(ftnbody_omega_3b,2)
                do irow = 1, size(ftnbody_omega_3b,1)
                  icount_all = icount_all + 1 
                  write(unitlammps, '(ES25.15E3)')   ftnbody_omega_3b(irow, icol)
                  write(unitpot, '(ES25.15E3)')   ftnbody_omega_3b(irow, icol)
                end do 
              end do 
              do icol = 1, size(ftnbody_phase_random_3b,1)
                  icount_all = icount_all + 1 
                  write(unitlammps, '(ES25.15E3)')   ftnbody_phase_random_3b(icol)
                  write(unitpot, '(ES25.15E3)')   ftnbody_phase_random_3b(icol)
              end do       
            end if   

            if (l_body_order(4)) then
              do icol = 1, size(ftnbody_omega_4b,2)
                do irow = 1, size(ftnbody_omega_4b,1)
                  icount_all = icount_all + 1 
                  write(unitlammps, '(ES25.15E3)')   ftnbody_omega_4b(irow, icol)
                  write(unitpot, '(ES25.15E3)')   ftnbody_omega_4b(irow, icol)
                end do 
              end do 
              do icol = 1, size(ftnbody_phase_random_4b,1)
                  icount_all = icount_all + 1 
                  write(unitlammps, '(ES25.15E3)')   ftnbody_phase_random_4b(icol)
                    write(unitpot, '(ES25.15E3)')   ftnbody_phase_random_4b(icol)
              end do       
            end if   

            if (l_body_order(5)) then
              do icol = 1, size(ftnbody_omega_5b,2)
                do irow = 1, size(ftnbody_omega_5b,1)
                  icount_all = icount_all + 1 
                  write(unitlammps, '(ES25.15E3)')   ftnbody_omega_5b(irow, icol)
                  write(unitpot, '(ES25.15E3)')   ftnbody_omega_5b(irow, icol)
                end do 
              end do 
              do icol = 1, size(ftnbody_phase_random_5b,1)
                  icount_all = icount_all + 1 
                  write(unitlammps, '(ES25.15E3)')   ftnbody_phase_random_5b(icol)
                  write(unitpot, '(ES25.15E3)')   ftnbody_phase_random_5b(icol)
              end do       
            end if   
              
        end select

      end if
      ! end    supplementary info  descriptor ..........................

      ! random-sampling descriptors: dump the random features to random.xml so a
      ! saved potential can reproduce the descriptor at MD / prediction time.
      if (descriptor_type == descriptor_ftnbody) then
        call write_ftnbody_random_xml('random.xml')
      end if

      ! init zbl and 2b active .........................................
      if (zbl_potential .and. zbl_type == zbl_mode_alone) then
        ! --- Simplified ZBL section for zbl_mode_alone (additive ZBL, no k2b/bridge) ---
        ! 1: zbl flag
        icount_all = icount_all + 1
        write (unitpot, '(ES25.15E3)') 1.d0
        write (unitlammps, '(ES25.15E3)') 1.d0
        ! 2: zbl_type
        icount_all = icount_all + 1
        write (unitpot, '(ES25.15E3)') dble(zbl_type)
        write (unitlammps, '(ES25.15E3)') dble(zbl_type)
        ! 3: r1_zbl
        icount_all = icount_all + 1
        write (unitpot, '(ES25.15E3)') r1_zbl
        write (unitlammps, '(ES25.15E3)') r1_zbl
        ! 4: r2_zbl
        icount_all = icount_all + 1
        write (unitpot, '(ES25.15E3)') r2_zbl
        write (unitlammps, '(ES25.15E3)') r2_zbl
        ! 5: r_cut_in
        icount_all = icount_all + 1
        write (unitpot, '(ES25.15E3)') r_cut_in
        write (unitlammps, '(ES25.15E3)') r_cut_in
        ! 6: r_cut_width_in
        icount_all = icount_all + 1
        write (unitpot, '(ES25.15E3)') r_cut_width_in
        write (unitlammps, '(ES25.15E3)') r_cut_width_in
        ! 7: type_fcut
        icount_all = icount_all + 1
        write (unitpot, '(ES25.15E3)') dble(type_fcut)
        write (unitlammps, '(ES25.15E3)') dble(type_fcut)
        ! 8: r_cut_width
        icount_all = icount_all + 1
        write (unitpot, '(ES25.15E3)') r_cut_width
        write (unitlammps, '(ES25.15E3)') r_cut_width

      else if ((zbl_potential .and. zbl_type == zbl_mode_default_k2b) .or. activate_k2b) then 
        !zbl or not ... 
        ! a: 1   z:1
        if (zbl_potential) then 
          icount_all = icount_all + 1 
          write (unitpot, '(ES25.15E3)') 1.d0
          write (unitlammps, '(ES25.15E3)') 1.d0
        else 
          icount_all = icount_all + 1 
          write (unitpot, '(ES25.15E3)') 0.d0
          write (unitlammps, '(ES25.15E3)') 0.d0
        end if 

        ! a: 2   z:2
        !two dimensions ... 
        icount_all = icount_all + 1 
        write(unitpot, '(ES25.15E3)') dble(dim_xdesc_patch)
        write(unitlammps, '(ES25.15E3)') dble(dim_xdesc_patch)

        ! a: 3   z:3
        icount_all = icount_all + 1 
        write(unitpot, '(ES25.15E3)') dble(dim_kernel_2b)
        write(unitlammps, '(ES25.15E3)') dble(dim_kernel_2b)

        ! a: 3+dim_kernel_2b   z:3+dim_kernel_2b
        !reading parameters ... 
        do i = dim_xdesc_linear + 1, dim_xdesc_linear + dim_kernel_2b 
          icount_all = icount_all + 1 
          write (unitpot, '(ES25.15E3)') w_params(i, 1)
          write (unitlammps, '(ES25.15E3)') w_params(i, 1)
        end do

        !WRITTING SPARSE POINTS BUILD OR NOT ... 
        ! a: 4+dim_kernel_2b   z:4+dim_kernel_2b
        icount_all = icount_all + 1 
        write (unitpot, '(ES25.15E3)') dble(sparse_points_2b_build)
        write (unitlammps, '(ES25.15E3)') dble(sparse_points_2b_build)

        ! a: 5+dim_kernel_2b   z:5+dim_kernel_2b
        icount_all = icount_all + 1 
        write (unitpot, '(ES25.15E3)') dble(np_radial_2b)
        write (unitlammps, '(ES25.15E3)') dble(np_radial_2b)

        ! a: 5+dim_kernel_2b+np_radial_2p   z:5+dim_kernel_2b+np_radial_2p
        do ii = 1, np_radial_2b
          icount_all = icount_all + 1 
          write (unitpot, '(ES25.15E3)') zr_2b(ii)
          write (unitlammps, '(ES25.15E3)') zr_2b(ii)
        end do   

        ! a: 5+5*dim_kernel_2b+np_radial_2p   z:5+5*dim_kernel_2b+np_radial_2p
        do ik = 1, dim_kernel_2b
          icount_all = icount_all + 4 
          !1
          write (unitpot, '(ES25.15E3)') dble(zpoints_2b(ik)%type1_z)
          write (unitlammps, '(ES25.15E3)') dble(zpoints_2b(ik)%type1_z)
          !2
          write (unitpot, '(ES25.15E3)') dble(zpoints_2b(ik)%type2_z)
          write (unitlammps, '(ES25.15E3)') dble(zpoints_2b(ik)%type2_z)
          !3
          write (unitpot, '(ES25.15E3)') dble(zpoints_2b(ik)%irad)
          write (unitlammps, '(ES25.15E3)') dble(zpoints_2b(ik)%irad)
          !4
          write (unitpot, '(ES25.15E3)') zpoints_2b(ik)%rr
          write (unitlammps, '(ES25.15E3)') zpoints_2b(ik)%rr
        end do

        
        ! a: 6+5*dim_kernel_2b + np_radial_2b   z:6+5*dim_kernel_2b + np_radial_2b
        icount_all = icount_all + 1 
        write (unitpot, '(ES25.15E3)') sigma_2b  
        write (unitlammps, '(ES25.15E3)') sigma_2b  

        ! a: 7+5*dim_kernel_2b + np_radial_2b   z:7+5*dim_kernel_2b + np_radial_2b
        icount_all = icount_all + 1 
        write (unitpot, '(ES25.15E3)') delta_2b 
        write (unitlammps, '(ES25.15E3)') delta_2b   

        if (zbl_potential) then 

          ! a: 7+5*dim_kernel_2b + np_radial_2b   z:8+5*dim_kernel_2b + np_radial_2b       
           icount_all = icount_all + 1 
           write (unitpot, '(ES25.15E3)') r1_zbl 
           write (unitlammps, '(ES25.15E3)') r1_zbl 
      
    
          ! a: 7+5*dim_kernel_2b + np_radial_2b   z:9+5*dim_kernel_2b + np_radial_2b  
           icount_all = icount_all + 1 
           write (unitpot, '(ES25.15E3)') r2_zbl 
           write (unitlammps, '(ES25.15E3)') r2_zbl 

          ! a: 7+5*dim_kernel_2b + np_radial_2b   z:10+5*dim_kernel_2b + np_radial_2b      
           icount_all = icount_all + 1 
           write (unitpot, '(ES25.15E3)') dble(type_rac_zbl) 
           write (unitlammps, '(ES25.15E3)') dble(type_rac_zbl)

    
          ! a: 7+5*dim_kernel_2b + np_radial_2b   z:11+5*dim_kernel_2b + np_radial_2b  
           icount_all = icount_all + 1 
           write (unitpot, '(ES25.15E3)') dble(size(params_k2b_to_zbl,1) )
           write (unitlammps, '(ES25.15E3)') dble(size(params_k2b_to_zbl,1))

          ! a: 7+5*dim_kernel_2b + np_radial_2b   z:12+5*dim_kernel_2b + np_radial_2b       
           icount_all = icount_all + 1 
           write (unitpot, '(ES25.15E3)') dble(size(params_k2b_to_zbl,2)) 
           write (unitlammps, '(ES25.15E3)') dble(size(params_k2b_to_zbl,2))
           ! iijj = size(params_k2b_to_zbl,1)*size(params_k2b_to_zbl,2)
           ! IF THIS UPDATED PLEASE UPDATE BEFORE !!!!!!!
          ! a: 7+5*dim_kernel_2b + np_radial_2b   z:12+iijj+5*dim_kernel_2b + np_radial_2b     
           do ii = 1, size(params_k2b_to_zbl,1)
             do jj = 1, size(params_k2b_to_zbl,2)
              icount_all = icount_all + 1 
              write (unitpot, '(ES25.15E3)') params_k2b_to_zbl(ii,jj)
              write (unitlammps, '(ES25.15E3)') params_k2b_to_zbl(ii,jj)
             end do 
           end do     
        end if  ! zbl_potential 


        ! a: 8+5*dim_kernel_2b + np_radial_2b   z:13+iijj+5*dim_kernel_2b + np_radial_2b         
        icount_all = icount_all + 1 
        write (unitpot, '(ES25.15E3)') r_cut_2b
        write (unitlammps, '(ES25.15E3)') r_cut_2b

        ! a: 9+5*dim_kernel_2b + np_radial_2b   z:14+iijj+5*dim_kernel_2b + np_radial_2b         
        icount_all = icount_all + 1 
        write (unitpot, '(ES25.15E3)') r_cut_width_2b
        write (unitlammps, '(ES25.15E3)') r_cut_width_2b

        ! a: 10+5*dim_kernel_2b + np_radial_2b   z:15+iijj+5*dim_kernel_2b + np_radial_2b 
        icount_all = icount_all + 1 
        write (unitpot, '(ES25.15E3)') r_cut_in
        write (unitlammps, '(ES25.15E3)') r_cut_in

        ! a: 11+5*dim_kernel_2b + np_radial_2b   z:16+iijj+5*dim_kernel_2b + np_radial_2b 
        icount_all = icount_all + 1 
        write (unitpot, '(ES25.15E3)') r_cut_width_in
        write (unitlammps, '(ES25.15E3)') r_cut_width_in
        
        ! a: 12+5*dim_kernel_2b + np_radial_2b   z:17+iijj+5*dim_kernel_2b + np_radial_2b 
        icount_all = icount_all + 1 
        write (unitpot, '(ES25.15E3)') dble(type_fcut) 
        write (unitlammps, '(ES25.15E3)') dble(type_fcut) 

        ! a: 13+5*dim_kernel_2b + np_radial_2b   z:18+iijj+5*dim_kernel_2b + np_radial_2b 
        icount_all = icount_all + 1 
        write (unitpot, '(ES25.15E3)') r_cut_width 
        write (unitlammps, '(ES25.15E3)') r_cut_width 

      end if   
      ! end  zbl or  2b active .....................................................

      ! quandratic part ............................................................

      if (mld_order == mld_quadratic) then 
        do i = dim_xdesc_linear + dim_xdesc_patch + 1, dim_xdesc_linear + dim_xdesc_patch + dim_xdesc_quadratic - dim_xdesc_linear
          icount_all = icount_all + 1 
          write (unitpot, '(ES25.15E3)') w_params(i, 1)
          write (unitlammps, '(ES25.15E3)') w_params(i, 1)
        end do   
      end if   
      !end quadratic part ..........................................................


      ! kernel part ................................................................
      if (mld_order == mld_kernel) then
        do i = dim_xdesc_linear + dim_xdesc_patch + 1 , dim_xdesc_linear + dim_xdesc_patch + dim_kernel
          icount_all = icount_all + 1 
          write (unitpot, '(ES25.15E3)') w_params(i, 1)
          write (unitlammps, '(ES25.15E3)') w_params(i, 1)
        end do

        if (kernel_type /= kernel_random_po) then
          do ik = 1, dim_kernel
            do i = 1, dim_xdesc
              icount_all = icount_all + 1 
              !$! write (unitpot, '(ES25.15E3)') global_kernel(i, ik)
              write (unitlammps, '(ES25.15E3)') global_kernel(i, ik)
            end do
          end do
        end if

        if (kernel_type == kernel_random) then
          do ik = 1, dim_kernel
            icount_all = icount_all + 1 
            !$! write (unitpot, '(ES25.15E3)') kernel_phase_random(ik)
            write (unitlammps, '(ES25.15E3)') kernel_phase_random(ik)
          end do
        end if

        if (kernel_type == kernel_po) then
          icount_all = icount_all + 3 
          write (unitlammps, '(e25.15)') sigma_kernel
          write (unitlammps, '(e25.15)') length_kernel
          write (unitlammps, '(e25.15)') kernel_power
        end if

        if (kernel_type == kernel_po_scaled) then
          icount_all = icount_all + 3 
          write (unitlammps, '(e25.15)') sigma_kernel
          write (unitlammps, '(e25.15)') length_kernel
          write (unitlammps, '(e25.15)') kernel_power
        end if

        if (kernel_type == kernel_random_po) then
          do ik = 1, dim_kernel
            !$! write (unitpot, *) basis_random_po(ik)%dim_omega
            !$! do ii = 1, basis_random_po(ik)%dim_omega
            !$!   do i = 1, dim_xdesc
            !$!     icount_all = icount_all + 1 
            !$!     write (unitpot, *) basis_random_po(ik)%omega(i, ii)
            !$!   end do
            !$! end do
            write (unitlammps, *) basis_random_po(ik)%dim_omega
            do ii = 1, basis_random_po(ik)%dim_omega
              do i = 1, dim_xdesc
                icount_all = icount_all + 1 
                write (unitlammps, *) basis_random_po(ik)%omega(i, ii)
              end do
            end do
          end do
          icount_all = icount_all + 3 
          write (unitlammps, '(e25.15)') sigma_kernel
          write (unitlammps, '(e25.15)') length_kernel
          write (unitlammps, '(e25.15)') kernel_power

        end if
      end if
      close (unitpot)
      close (unitlammps)
    end if                  ! rangml==0
    ! end kernel part ...................................................................
    if (icount_all /= size_params_to_write) then
      call log_warning("ML: writting lammps parameters have a  problem. Can have problems with LAMMPS. Ask Master!") 
      call log_warning("ML: size_params_to_write and icount_all "//vtoa(size_params_to_write) // vtoa(icount_all)) 
    end if   
    if (.not. (optimize_weights_db.or.optimize_weights_chem)) then
      !!! if (rangml == 0) then
      !!!   write (6, '("ML: The parameters were wrote ...")')
      !!! end if
      call log_info("ML: The parameters were wrote ...")
    end if

    _MLD_END_

  end subroutine write_snap_parameters


  subroutine dump_kernel_matrix
    use ml_in_ndm_module, only: rangml
    use derived_types, only: config_real, config_desc
    use module_kernel, only: kernel_dump, &
                             kernel_dump_by_mahalanobis_norm, kernel_dump_by_mahalanobis, &
                             kernel_dump_by_cur, kernel_dump_by_cur_maha, info_kernel
    use temporary_data_cov, only: dim_xdesc
    use module_Sigma_matrix, only: Sigma_sample_mcd, Sigma_sample_mcd_inv
    use module_ml_scalapack, only: context, scalapack_driver
    use module_end_ml, only : end_ml
    implicit none
    integer  ::  nunit, ii, ia, icnt, ic
    character(len=80)    :: CHFMT,text 
    ! The design matrix is the transpose of Amat.
    ! Amat (dim_desc+1, number_of_config). The desing matrix is the transpose of Amat.
    !    = size(Amat,1) ! dim_xdesc+1
    !    = size(Amat,2) ! config


    if (rangml == 0) then

      open (file='kernel_matrix.dat', newunit=nunit, action='write', status='unknown')
      write (CHFMT, *) '( i6,', dim_xdesc, 'e20.10, i6, i6, a15 )'

      select case (kernel_dump)

      case (kernel_dump_by_mahalanobis_norm)

        write (nunit, '(2i9)') size(info_kernel), dim_xdesc + 4
        icnt = 0
        do ii = 1, size(info_kernel)
                icnt = icnt + 1
          ic = info_kernel(ii)%iconf
          ia = info_kernel(ii)%ia
                write (nunit, CHFMT) icnt, config_desc(ic)%energy(:, ia), ic, ia, &
                  "  "//config_real(ic)%class//"_"//config_real(ic)%klm//"_"//config_real(ic)%cnumber
              end do

      case (kernel_dump_by_mahalanobis)

        write (nunit, '(2i9)') size(info_kernel), dim_xdesc + 4
        icnt = 0
        do ii = 1, size(info_kernel)
          icnt = icnt + 1
          ic = info_kernel(ii)%iconf
          ia = info_kernel(ii)%ia
          write (nunit, CHFMT) icnt, config_desc(ic)%energy(:, ia), ic, ia, &
            "  "//config_real(ic)%class//"_"//config_real(ic)%klm//"_"//config_real(ic)%cnumber
        end do

      case (kernel_dump_by_cur)

        write (nunit, '(2i9)') size(info_kernel), dim_xdesc + 4
        icnt = 0
        do ii = 1, size(info_kernel)
          icnt = icnt + 1
          ic = info_kernel(ii)%iconf
          ia = info_kernel(ii)%ia
          write (nunit, CHFMT) icnt, config_desc(ic)%energy(:, ia), ic, ia, &
            "  "//config_real(ic)%class//"_"//config_real(ic)%klm//"_"//config_real(ic)%cnumber
        end do

      case (kernel_dump_by_cur_maha)

        write (nunit, '(2i9)') size(info_kernel), dim_xdesc + 4
        icnt = 0
        do ii = 1, size(info_kernel)
          icnt = icnt + 1
          ic = info_kernel(ii)%iconf
          ia = info_kernel(ii)%ia
          write (nunit, CHFMT) icnt, config_desc(ic)%energy(:, ia), ic, ia, &
            "  "//config_real(ic)%class//"_"//config_real(ic)%klm//"_"//config_real(ic)%cnumber
        end do

      case default

        write (6, '("ML: this kernel dump type is not implemented, kernel_dump is set to...",i6)') kernel_dump
        text = ' --- ERROR on kernel dump in dump_kernel_matrix ERROR ---'
        call end_ml(text, scalapack_driver, context, rangml)
      end select
      close (nunit)
    end if                  ! rangml

    if (rangml == 0) then
      open (file='inverse_Sigma_mcd_matrix.mat', newunit=nunit, action='write', status='unknown')
      write (nunit, '(2i9)') dim_xdesc, dim_xdesc
      write (CHFMT, *) '(', dim_xdesc, 'e25.15)'
      do ii = 1, size(Sigma_sample_mcd_inv, 2)
        write (nunit, CHFMT) Sigma_sample_mcd_inv(:, ii)
      end do
      close (nunit)

      open (file='Sigma_mcd_matrix.mat', newunit=nunit, action='write', status='unknown')
      write (nunit, '(2i9)') dim_xdesc, dim_xdesc
      write (CHFMT, *) '(', dim_xdesc, 'e25.15)'
      do ii = 1, size(Sigma_sample_mcd, 2)
        write (nunit, CHFMT) Sigma_sample_mcd(:, ii)
      end do
      close (nunit)
    end if                  ! rangml

  end subroutine dump_kernel_matrix


  subroutine write_zbl_k2b
    use ml_in_ndm_module, only: rangml 
    use module_kind_variables, only: kind_double
    use derived_types, only : typ_species_half
    use module_neigh_local, only: type_fcut
    use module_kernel_2b, only: r_cut_2b, r_cut_width_2b
    use module_zbl, only : rr_k2b
    use module_chemical_species, only: size_species_half
    use module_form_kernel_nbody, only: kernel_2b_end
    use module_potential_zbl, only: full_potential_zbl, potential_zbl, potential_zbl_second
      !cont_test 
    use mesh_grid, only: linear_grid

    implicit none 
        !cont_test
    integer :: ii, is, type_db_ja, type_db_ia  
    integer :: np
    real(kind_double), dimension(:), allocatable :: xp_grid
    real(kind_double) :: rr, ee, ff, ss, Z_ja, Z_ia
    real(kind_double), dimension(:,:), allocatable :: kmat, emat, zimat 
    character(len=100)   :: zbl_filename, k2b_filename 
    character(len=100)   :: CHFMT
    integer :: zunit, kunit, ziunit 

    np=1000
    call linear_grid(np, 0.8d0,3.0d0,xp_grid)

    if (allocated(kmat)) deallocate(kmat) ; allocate(kmat(size_species_half, np))
    if (allocated(emat)) deallocate(emat) ; allocate(emat(size_species_half, np))
    if (allocated(zimat)) deallocate(zimat) ; allocate(zimat(size_species_half, np))
    do is = 1, size_species_half
      type_db_ja= typ_species_half(is)%type1
      type_db_ia= typ_species_half(is)%type2
      Z_ja = dble(typ_species_half(is)%Z1)
      Z_ia = dble(typ_species_half(is)%Z2)
      do ii = 1, np 
         rr = xp_grid(ii)
         if (rr <= rr_k2b) then 
           call full_potential_zbl(.true., is, rr, ee, ff)
           emat(is,ii) = ee 
         end if 
         !if (is==1) write(324, '(3e25.10)') rr, ee, ff
         call kernel_2b_end(type_db_ja, type_db_ia, .true., type_fcut, rr, r_cut_2b, r_cut_width_2b, ee, ff, ss )
         kmat(is,ii) = ee 
         if (rr >= rr_k2b) then
           emat(is,ii) = ee  
         end if  
         !if (is==1) write(325,'(4e25.10)') rr, ee, ff, ss 
         call potential_zbl(.true.,Z_ja, Z_ia, rr, ee, ff)
         zimat(is,ii) = ee 
         !if (is==1) write(327, '(4e25.10)') rr, ee, ff, ss
         call potential_zbl_second(Z_ja, Z_ia, rr, ee, ff, ss)
         !if (is==1) write(326, '(4e25.10)') rr, ee, ff, ss
      end do 
    end do 
    if (rangml == 0) then 
      write (CHFMT, *) '(e20.10, 1x, ', int(size_species_half), 'e20.10)'
      zbl_filename="zbl_potential.dat" 
      k2b_filename="k2b_potential.dat" 
      open (newunit=ziunit, file="zbl_draft.dat", status='unknown')
      open (newunit=zunit, file=zbl_filename, status='unknown')
      open (newunit=kunit, file=k2b_filename, status='unknown')

      ! write header in order
      write(zunit,'("# ", 9x,"r",8x)', advance='no')
      write(ziunit,'("# ", 9x,"r",8x)', advance='no')
      write(kunit,'("# ", 9x,"r",8x)', advance='no')
      do is = 1, size_species_half-1
        write(zunit, '(8x,i3, 9x)', advance='no')  is
        write(ziunit, '(8x,i3, 9x)', advance='no')  is
        write(kunit, '(8x,i3, 9x)', advance='no')  is
      end do
      write(zunit, '(8x,i3, 9x)')  size_species_half   
      write(ziunit, '(8x,i3, 9x)') size_species_half
      write(kunit, '(8x,i3, 9x)')  size_species_half

      !write header in chemical symbols ...
      write(zunit,'("# ", 9x,"r",8x)', advance='no')
      write(ziunit,'("# ", 9x,"r",8x)', advance='no')
      write(kunit,'("# ", 9x,"r",8x)', advance='no')
      do is = 1, size_species_half-1
        write(zunit, '(7x,a2," -", a2,7x)', advance='no')  typ_species_half(is)%ch1,  typ_species_half(is)%ch2
        write(ziunit, '(7x,a2," -", a2,7x)', advance='no')  typ_species_half(is)%ch1,  typ_species_half(is)%ch2
        write(kunit, '(7x,a2," -", a2,7x)', advance='no')  typ_species_half(is)%ch1,  typ_species_half(is)%ch2
      end do
      write(zunit, '(7x,a2,"-", a2,7x)')  typ_species_half(size_species_half)%ch1,  typ_species_half(size_species_half)%ch2   
      write(ziunit, '(7x,a2,"-", a2,7x)')  typ_species_half(size_species_half)%ch1,  typ_species_half(size_species_half)%ch2   
      write(kunit, '(7x,a2,"-", a2,7x)')  typ_species_half(size_species_half)%ch1,  typ_species_half(size_species_half)%ch2   


      ! write header pairs in types ....
      write(zunit,'("# ", 9x,"r",8x)', advance='no')
      write(ziunit,'("# ", 9x,"r",8x)', advance='no')
      write(kunit,'("# ", 9x,"r",8x)', advance='no')
      do is = 1, size_species_half-1
        write(zunit, '(7x,i2," -", i2,7x)', advance='no')  typ_species_half(is)%type1,  typ_species_half(is)%type2
        write(ziunit, '(7x,i2," -", i2,7x)', advance='no')  typ_species_half(is)%type1,  typ_species_half(is)%type2
        write(kunit, '(7x,i2," -", i2,7x)', advance='no')  typ_species_half(is)%type1,  typ_species_half(is)%type2
      end do
      write(zunit, '(7x,i2," -", i2,7x)')  typ_species_half(size_species_half)%type1,  typ_species_half(size_species_half)%type2   
      write(ziunit, '(7x,i2," -", i2,7x)')  typ_species_half(size_species_half)%type1,  typ_species_half(size_species_half)%type2   
      write(kunit, '(7x,i2," -", i2,7x)')  typ_species_half(size_species_half)%type1,  typ_species_half(size_species_half)%type2   

      do ii =1, np
        write(zunit, CHFMT) xp_grid(ii), emat(:,ii)
        write(ziunit, CHFMT) xp_grid(ii), zimat(:,ii)
        write(kunit, CHFMT) xp_grid(ii), kmat(:,ii)
      end do   
    end if   

  end subroutine write_zbl_k2b


end module module_write_parameters



module module_create_descdb

  use ml_in_ndm_module, only: rangml
  implicit none
  integer  :: CSTAT, ESTAT
  character(len=100)    :: namefile, directory_name, CMSG
  logical  :: ok, analysis_get_potential_on_the_fly = .false.

  contains

  subroutine create_descdb(prefix, db_path) 

    use mld_subworld, only: subrank

    character(len=*), intent(in)  :: prefix 
    character(len=*), intent(in)  :: db_path 

    ! Every subworld master creates the directory: each one writes its
    ! own configs, and rangml==0 belongs only to subworld 0. mkdir -p
    ! is idempotent so concurrent calls are safe.
    if (subrank == 0) then
      !test and create the directory ....
      directory_name = trim(adjustl(prefix))//trim(adjustl(db_path))
      !this is better but work only for ifort: inquire(directory=directory_name, exist=ok)
      !this one works correctly only for gfortran. For ifort it executes always mkdir -p
      inquire (file=directory_name, exist=ok)

      if (.not. (ok)) then
        call execute_command_line('mkdir -p '//directory_name, EXITSTAT=ESTAT, CMDSTAT=CSTAT, CMDMSG=CMSG)
        if (CSTAT > 0) then
          write (6, '("Creation of directory", a,   " failed  with error", a)') directory_name, TRIM(CMSG)
        else if (CSTAT < 0) then
          write (6, '("Command execution not supported", i5)') CSTAT
        else
          !write (6,'("Command completed with status", i5)') ESTAT
        end if
      end if
    end if                  ! rangml == 0
  end subroutine create_descdb

  subroutine clean_descdb_energy (iconf, prefix, db_path)
    use ml_in_ndm_module, only:  desc_file_format, eml_type, csv_type, npz_type
    use derived_types, only: config_real
    use mld_subworld, only: subrank
    use mld_logger
    implicit none 
    character(len=*), intent(in)  :: prefix 
    character(len=*), intent(in)  :: db_path 
    integer, intent(in) :: iconf 
    character(len=100)   :: efilename
    logical  :: desc_energy_local
    character(len=80) :: czozo 


    _NAMECURRENT_("clean_descdb_energy")



    _MLD_BEGIN_
    czozo = prefix 
    desc_energy_local = config_real(iconf)%has_energy
    if (.not. desc_energy_local) return
    !if ((i_start_at == 0) .and. (i_final_at == 0)) return
  
    if (subrank == 0) then
      if (desc_file_format == eml_type) then 
        efilename = 'desc'//trim(adjustl(db_path))//config_real(iconf)%class//'_'//config_real(iconf)%klm//'_'//config_real(iconf)%cnumber//'.eml'
      end if 
      
      if (desc_file_format == csv_type) then 
        efilename = 'desc'//trim(adjustl(db_path))//config_real(iconf)%class//'_'//config_real(iconf)%klm//'_'//config_real(iconf)%cnumber//'.csv'
      end if 
      
      if (desc_file_format == npz_type) then 
        efilename = 'desc'//trim(adjustl(db_path))//config_real(iconf)%class//'_'//config_real(iconf)%klm//'_'//config_real(iconf)%cnumber//'.npz'
      end if   

      call execute_command_line('rm -f '//trim(efilename), EXITSTAT=ESTAT, CMDSTAT=CSTAT, CMDMSG=CMSG)
      if (CSTAT > 0) then
        write (6, '("removing the file ", a,   " failed  with error", a)') efilename, TRIM(CMSG)
      else if (CSTAT < 0) then
        write (6, '("Command execution rm not supported", i5)') CSTAT
      end if   
    end if

    _MLD_END_   
     
  end subroutine   

end module module_create_descdb
