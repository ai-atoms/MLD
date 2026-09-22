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

subroutine write_descriptors(iconf)

#if(PARA)
  !use mpi
#endif
  use module_kind_variables, only: kind_double
  use ml_in_ndm_module, only: write_desc, write_desc_dump
  use module_db_setup, only: db_path                            
  use mld_logger
  use module_create_descdb, only: create_descdb, clean_descdb_energy
  use derived_types, only: config_desc, config_real
  use module_chemical_species, only: img_weighted
  use module_kernel_2b, only: dim_kernel_2b, activate_k2b
  use m_npy, only : natoms_npz
  use mld_mpi

  !k2b
  use module_kernel_2b, only : activate_k2b, dim_kernel_2b

  implicit none

  integer, intent(in)  :: iconf
  integer :: ia, ia_start, ia_end, i_pack 
  real(kind_double), dimension(:,:), allocatable :: local_energy 


  _NAMECURRENT_("write_descriptors")

  _MLD_BEGIN_

  if (write_desc .or. write_desc_dump) then
    natoms_npz = config_real(iconf)%nat
    call create_descdb("desc", db_path)
    call clean_descdb_energy(iconf, "desc", db_path)
  end if


  !bricollage for channels - TODO final version 
  if (img_weighted) then
    call write_channels(iconf) 
  else if (activate_k2b) then
    ia_start=1
    ia_end=config_real(iconf)%nat
    i_pack = ia_end - ia_start + 1 
    if (allocated(local_energy)) deallocate(local_energy)
    allocate(local_energy(config_desc(iconf)%dim_desc + dim_kernel_2b + 1, config_real(iconf)%nat))
    do ia = ia_start, ia_end 
        local_energy(1,ia-ia_start + 1) = dble(ia)
        local_energy(2:config_desc(iconf)%dim_desc+1,ia-ia_start + 1) = config_desc(iconf)%energy(:,ia)
        local_energy(config_desc(iconf)%dim_desc+2:config_desc(iconf)%dim_desc+1+dim_kernel_2b,ia-ia_start + 1) = config_desc(iconf)%energy_k2b(:,ia)
    end do   
    call write_energy_descriptors(iconf, config_desc(iconf)%dim_desc + dim_kernel_2b, i_pack, ia_start, ia_end, local_energy)
  else 
    ia_start=1
    ia_end=config_real(iconf)%nat
    i_pack = ia_end - ia_start + 1 

    if (activate_k2b) then
      if (allocated(local_energy)) deallocate(local_energy)
      allocate(local_energy(config_desc(iconf)%dim_desc + 1 + dim_kernel_2b, config_real(iconf)%nat))
      do ia = ia_start, ia_end 
          local_energy(1,ia-ia_start + 1) = dble(ia)
          local_energy(2:config_desc(iconf)%dim_desc+1,ia-ia_start + 1) = config_desc(iconf)%energy(:,ia)
          local_energy(config_desc(iconf)%dim_desc+1:config_desc(iconf)%dim_desc + dim_kernel_2b +1 ,ia-ia_start + 1) = config_desc(iconf)%energy_k2b(:,ia)
      end do   
      call write_energy_descriptors(iconf, config_desc(iconf)%dim_desc + dim_kernel_2b, i_pack, ia_start, ia_end, local_energy)

    else 
      if (allocated(local_energy)) deallocate(local_energy)
      allocate(local_energy(config_desc(iconf)%dim_desc + 1, config_real(iconf)%nat))
      do ia = ia_start, ia_end 
          local_energy(1,ia-ia_start + 1) = dble(ia)
          local_energy(2:config_desc(iconf)%dim_desc+1,ia-ia_start + 1) = config_desc(iconf)%energy(:,ia)
      end do   
      call write_energy_descriptors(iconf, config_desc(iconf)%dim_desc, i_pack, ia_start, ia_end, local_energy)
    end if 

  end if 
  !end bricollage for channels - TODO 

  call write_force_descriptors(iconf)
  if (write_desc_dump) then
    call create_descdb("dump", db_path)
    call write_force_descriptors_dump(iconf)
  end if
  _MLD_END_
end subroutine write_descriptors




subroutine write_energy_descriptors(iconf, ndim, i_pack, ia_start, ia_end, local_shape)
  
  use module_kind_variables, only: kind_double
  use ml_in_ndm_module, only: rangml, write_desc, write_desc_dump, &
                              desc_file_format, eml_type, csv_type, & 
                              npz_type, hdf_type, ml_type_descriptors, ml_type, lmask
  use module_db_setup, only: db_path                             
  use derived_types, only: config_desc, config_real
  use mld_logger
  use m_npy
  use mld_mpi, only: rangml, comm_mld
  use mld_subworld, only: subrank

  implicit none

  integer, intent(in)  :: iconf, i_pack, ia_start, ia_end, ndim
  real(kind_double), dimension(ndim+1,i_pack) :: local_shape 
  integer  :: eunit, ia
  character(len=100)   :: efilename
  character(len=100)   :: CHFMT
  logical  :: desc_energy_local, exist 
  character(len=1), parameter :: comma = char(44), dquote = char(34)

  _NAMECURRENT_("write_energy_descriptors")



  _MLD_BEGIN_
  desc_energy_local = config_real(iconf)%has_energy
  if (.not. desc_energy_local) return

  if (desc_file_format == hdf_type) then
    ! TestingPurposeArnaud hdf5 
    ! Each subworld master writes its own configs (rangml==0 only
    ! belongs to subworld 0, so it would miss every other subworld).
    if (subrank == 0) then
      efilename = config_real(iconf)%class//'_'//config_real(iconf)%klm//'_'//config_real(iconf)%cnumber
      call write_hdf5_energy(iconf, trim(efilename)) 
    end if
    !call comm_mld%barrier

  end if

  if (desc_file_format == eml_type) then 
    efilename = 'desc'//trim(adjustl(db_path))//config_real(iconf)%class//'_'//config_real(iconf)%klm//'_'//config_real(iconf)%cnumber//'.eml'


    if (write_desc .or. write_desc_dump) then
      if (subrank == 0) then
        eunit=42
        inquire(file=efilename, exist=exist)
        if (exist) then
          open(eunit, file=efilename, status="old", position="append", action="write")
        else
          open(eunit, file=efilename, status="new", action="write")
        end if
        write (CHFMT, *) '(es30.18e3, 1x, ', int(ndim), 'es30.18e3)'
        do ia = ia_start, ia_end  
          if (ml_type == ml_type_descriptors) then 
            if (lmask) then 
              if (.not.config_desc(iconf)%amask(ia)) cycle 
            end if
          end if 
          write (eunit, FMT=CHFMT) local_shape(:,ia - ia_start + 1)

        end do
        close (eunit)
      end if                  ! rangml==0
    end if
  end if
  
  if (desc_file_format == csv_type) then 
    ! TestingPurposeArnaud csv 
    efilename = 'desc'//trim(adjustl(db_path))//config_real(iconf)%class//'_'//config_real(iconf)%klm//'_'//config_real(iconf)%cnumber//'.csv'
    !cos ndim = config_desc(iconf)%dim_desc
    if (write_desc .or. write_desc_dump) then
      if (subrank == 0) then
        eunit=42
        inquire(file=efilename, exist=exist)
        if (exist) then
          open(eunit, file=efilename, status="old", position="append", action="write")
        else
          open(eunit, file=efilename, status="new", action="write")
        end if
        !!! open (unit=eunit, file=efilename, status='unknown')
        write (CHFMT, *) '(es30.18e3, 1x, ', int(ndim), '(',dquote, comma, dquote, comma,' es30.18e3))'
        !cos do ia = 1, config_real(iconf)%nat
        do ia = ia_start, ia_end  
          if (ml_type == ml_type_descriptors) then 
            if (lmask) then 
              if (.not.config_desc(iconf)%amask(ia)) cycle
            end if   
          end if    
          !cos write (eunit, FMT=CHFMT) ia, config_desc(iconf)%energy(:, ia)
          write (eunit, FMT=CHFMT) local_shape(:,ia - ia_start + 1)

        end do
        close (eunit)
      end if                  ! rangml==0
    end if
  end if 
  
  if (desc_file_format == npz_type) then 
    ! TestingPurposeArnaud npz 
    efilename = 'desc'//trim(adjustl(db_path))//config_real(iconf)%class//'_'//config_real(iconf)%klm//'_'//config_real(iconf)%cnumber//'.npz'
    !cos ndim = config_desc(iconf)%dim_desc
    if (write_desc .or. write_desc_dump) then
      if (subrank == 0) then
        !cos if (allocated(local_shape)) deallocate(local_shape)
        !cos allocate(local_shape(config_real(iconf)%nat, dim_xdesc+1))
        !cos do ia = 1, config_real(iconf)%nat
        do ia = ia_start, ia_end  
          if (ml_type == ml_type_descriptors) then 
            if (lmask) then 
              if (.not.config_desc(iconf)%amask(ia)) cycle
            end if 
          end if    
          !cos local_shape(ia, 2 : dim_xdesc+1) = config_desc(iconf)%energy(1:dim_xdesc, ia)
        end do
        !write(*,*) 'local_shape npz', size(local_shape,1), size(local_shape,2)
        call save_npy(efilename, transpose(local_shape))  
      end if                  ! rangml==0
    end if
  end if 

  _MLD_END_
end subroutine write_energy_descriptors

subroutine write_channels(iconf)

    use mpi
    use mld_mpi 
    use set_limits, only: i_start_on_proc, i_end_on_proc, i_size_on_proc
    use ml_in_ndm_module, only: rangml, i_start_at, i_final_at
    use derived_types, only: config_desc, config_real
    use module_kind_variables, only: kind_double
    use mld_logger
  
    implicit none
  
    integer, intent(in)  :: iconf
    integer  :: ndim, iproc
    integer  :: ia, ia_start, ia_end, i_pack
    real(kind_double), dimension(:, :), allocatable :: s_channel, r_channel 
    logical  :: desc_energy_local
  
    _NAMECURRENT_("write_channels")
  
    
    _MLD_BEGIN_

    desc_energy_local = config_real(iconf)%has_energy
    if (.not. desc_energy_local) return
    if ((i_start_at == 0) .and. (i_final_at == 0)) return
    ndim = config_desc(iconf)%dim_desc

    !there is only one procs who writing
    if (rangml /= 0) then

      do iproc = 1, nb_procsml - 1
        if (rangml == iproc) then
          ia_start = i_start_on_proc(iproc)
          ia_end = i_end_on_proc(iproc)
          i_pack = i_size_on_proc(iproc)
          if ((ia_start == 0) .and. (ia_end == 0)) cycle
          if (allocated(s_channel)) deallocate (s_channel); allocate (s_channel(ndim+1, i_pack))
          do ia = ia_start, ia_end
            s_channel(1, ia - ia_start + 1) = dble(ia)
            s_channel(2:ndim+1, ia - ia_start + 1) = config_desc(iconf)%channel(1:ndim, ia)
          end do
          ! deallocate(config_desc(iconf)%channel)
          call MPI_SEND(s_channel, size(s_channel, 1)*size(s_channel, 2), MPI_DOUBLE_PRECISION, 0, 1000000 + iproc, mpi_comm_mld, codeml)
          deallocate (s_channel)
        end if
      end do

    end if
  

    if (rangml == 0) then

      do iproc = 0, nb_procsml - 1
        ia_start = i_start_on_proc(iproc)
        ia_end = i_end_on_proc(iproc)
        i_pack = i_size_on_proc(iproc)

        if ((ia_start == 0) .and. (ia_end == 0)) cycle
        
        if (iproc == 0) then
        
          if (allocated(r_channel)) deallocate (r_channel); allocate (r_channel(ndim+1, i_pack))
          do ia = ia_start, ia_end 
            r_channel(1,ia - ia_start + 1) = dble(ia)
            r_channel(2:ndim + 1,ia - ia_start +1 ) = config_desc(iconf)%channel(1:ndim,ia)
          end do 
          call write_energy_descriptors(iconf, ndim, i_pack, ia_start, ia_end, r_channel)
          deallocate(r_channel)
          
        else
        
          if (allocated(r_channel)) deallocate (r_channel); allocate (r_channel(ndim+1, i_pack))
          call MPI_RECV(r_channel, size(r_channel, 1)*size(r_channel, 2), MPI_DOUBLE_PRECISION, iproc, 1000000 + iproc, mpi_comm_mld, statut_ml, codeml)
          call write_energy_descriptors(iconf, ndim, i_pack, ia_start, ia_end, r_channel)
          deallocate(r_channel)
        
        end if
      end do
      
   end if

    _MLD_END_
  
  end subroutine write_channels 
  


subroutine read_energy_descriptors_dump(iconf)

#if(PARA)
  use mpi
  use mld_mpi
#endif

  use set_limits, only: i_start_on_proc, i_end_on_proc, i_size_on_proc
  use ml_in_ndm_module, only: rangml, i_start_at, i_final_at
  use module_db_setup, only: db_path
  use derived_types, only: config_desc, config_real
  use module_kind_variables, only: kind_double
  use mld_logger

  implicit none

  integer, intent(in)  :: iconf
  integer  :: eunit, ip, ia, ndim,  ia_begin, ia_end, ii_count
  character(len=100)   :: efilename
  character(len=100)   :: CHFMT
  logical  :: ok, desc_energy_local
  real(kind_double), dimension(:), allocatable :: etemp
  real(kind_double) :: real_ia_temp
  ! real(kind_double), dimension(:), allocatable :: e_send, e_recieve
  real(kind_double), dimension(:, :), allocatable    :: e_send, e_recieve


  _NAMECURRENT_("read_energy_descriptors_dump")



  _MLD_BEGIN_
#if(PARA)

  if ((i_start_at == 0) .and. (i_final_at == 0)) then
    config_desc(iconf)%energy(:, :) = 0.d0
    return
  end if
  desc_energy_local = config_real(iconf)%has_energy
  if (.not. desc_energy_local) return
  config_desc(iconf)%energy(:, :) = 0.d0
  efilename = 'desc'//trim(adjustl(db_path))//config_real(iconf)%class//'_'//config_real(iconf)%klm//'_'//config_real(iconf)%cnumber//'.eml'
  ndim = config_desc(iconf)%dim_desc
  allocate (etemp(ndim))

  if (rangml == 0) then
    inquire (file=efilename, exist=ok)
    if (.not. (ok)) then
      write (6, '("The file for reading ",a,"is missing. What a bad day!")') efilename
      !call end_ml("stop with reading in read_write_energ")
    end if
  end if                  ! rangml == 0

  if (rangml == 0) then
    write (CHFMT, *) '(es30.18e3, 1x, ', int(ndim), 'es30.18e3)'
    open (newunit=eunit, file=efilename, status='unknown')
    do ip = 0, nb_procsml - 1
      ia_begin = i_start_on_proc(ip)
      ia_end = i_end_on_proc(ip)
      if ((ia_end == 0) .and. (ia_begin == 0)) cycle
      !if (allocated(e_send)) deallocate(e_send) ; allocate(e_send(ndim*i_size_on_proc(ip)))
      if (allocated(e_send)) deallocate (e_send); allocate (e_send(ndim, i_size_on_proc(ip)))
      ii_count = 0
      do ia = ia_begin, ia_end
        read (eunit, FMT=CHFMT) real_ia_temp, etemp(:)
        ii_count = ii_count + 1
        !e_send((ii_count-1)*ndim+1:ii_count*ndim) = etemp(1:ndim)
        e_send(:, ii_count) = etemp(:)
      end do
      if (ip == 0) then
        !do ii_count=ia_begin,ia_end
        !    config_desc(iconf)%energy(:,ii_count) = e_send( (ii_count-ia_begin)*ndim+1 : (ii_count-ia_begin)*ndim + ndim)
        !end do
        do ii_count = ia_begin, ia_end
          config_desc(iconf)%energy(:, ii_count) = e_send(:, ii_count)
        end do

      else
        !call MPI_SEND(e_send,size(e_send,1),MPI_DOUBLE_PRECISION,ip,10000+ip,mpi_comm_mld,codeml)
        call MPI_SEND(e_send, size(e_send, 1)*size(e_send, 2), MPI_DOUBLE_PRECISION, ip, 10000 + ip, mpi_comm_mld, codeml)
      end if
      deallocate (e_send)
    end do                  ! ip - on procs index but on rangml = 0
    close (eunit)
  end if


  do ip = 1, nb_procsml - 1
    ia_begin = i_start_on_proc(ip)
    ia_end = i_end_on_proc(ip)
    if ((ia_end == 0) .and. (ia_begin == 0)) cycle
    if (rangml == ip) then
      !if (allocated(e_recieve)) deallocate(e_recieve) ; allocate(e_recieve(ndim*i_size_on_proc(ip)))
      if (allocated(e_recieve)) deallocate (e_recieve); allocate (e_recieve(ndim, i_size_on_proc(ip)))
      call MPI_RECV(e_recieve, ndim*i_size_on_proc(ip), MPI_DOUBLE_PRECISION, 0, 10000 + ip, mpi_comm_mld, statut_ml, codeml)
      !do ii_count=i_start_on_proc(ip),i_end_on_proc(ip)
      !    config_desc(iconf)%energy(:,ii_count) = e_recieve( (ii_count-i_start_on_proc(ip))*ndim+1 : (ii_count-i_start_on_proc(ip))*ndim + ndim)
      !end do
      do ii_count = i_start_on_proc(ip), i_end_on_proc(ip)
        config_desc(iconf)%energy(:, ii_count) = e_recieve(:, ii_count - i_start_on_proc(ip) + 1)
      end do
      deallocate (e_recieve)
    end if
  end do

  deallocate (etemp)
  _MLD_END_
  return

#else
  write (6, '("ML:.............TODO Houston not MPI in   :",a)') NAMECURRENT
  return
#endif

end subroutine read_energy_descriptors_dump




subroutine write_force_descriptors(iconf)

  use ml_in_ndm_module, only: desc_forces, rangml
  use module_db_setup, only: db_path
  use derived_types, only: config_desc, config_real
  use mld_force_mod, only: pack_force_descriptor
  use mld_logger
  use mld_subworld, only: subrank
  use module_kernel_2b, only : activate_k2b, dim_kernel_2b

  implicit none

  integer, intent(in)  :: iconf
  integer  :: ndim
  character(len=100)   ::  ffilename
  character(len=60)    :: CHFMTf
  integer  :: funit, ia
  logical  :: desc_forces_local

  _NAMECURRENT_("write_force_descriptors")



  _MLD_BEGIN_
  desc_forces_local = desc_forces .and. (config_real(iconf)%has_force .or. config_real(iconf)%has_stress)

  ndim = config_desc(iconf)%dim_desc
  ! efilename='desc'//trim(adjustl(db_path))//config_real(iconf)%class//'_'//config_real(iconf)%klm//'_'//config_real(iconf)%cnumber//'.eml'
  if (desc_forces_local) then
    ffilename = 'desc'//trim(adjustl(db_path))//config_real(iconf)%class//'_'//config_real(iconf)%klm//'_'//config_real(iconf)%cnumber//'.fml'
    !TODOwrite double passage in snap_pack_force_descriptor for forces, this will be call 
    !TODOwrite double passage in snap_pack_force_descriptor for forces, this will be call again when Amat is fill.
    call pack_force_descriptor(iconf)
  end if

  ! Each subworld master writes its own configs (rangml==0 only
  ! belongs to subworld 0, so it would miss every other subworld).
  if (subrank == 0) then
    funit = 42
    if (desc_forces_local) open (funit, file=ffilename, status='unknown')
    write (CHFMTf, *) '(e25.15, 1x, ', int(3*ndim), 'e25.15)'
    do ia = 1, config_real(iconf)%nat
      if (desc_forces_local) then
        if (activate_k2b) then 
          write (funit, FMT=CHFMTf) dble(ia), real(config_desc(iconf)%pack_force_linear(:, 1, ia)), &
                real(config_desc(iconf)%pack_force_k2b(:, 1, ia)) ,&
                real(config_desc(iconf)%pack_force_linear(:, 2, ia)) ,&
                real(config_desc(iconf)%pack_force_k2b(:, 2, ia)) ,&
                real(config_desc(iconf)%pack_force_linear(:, 3, ia)) ,&
                real(config_desc(iconf)%pack_force_k2b(:, 3, ia))
        else
          write (funit, FMT=CHFMTf) dble(ia), real(config_desc(iconf)%pack_force_linear(:, 1, ia)), &
                real(config_desc(iconf)%pack_force_linear(:, 2, ia)), &
                real(config_desc(iconf)%pack_force_linear(:, 3, ia)) 
        end if 
      end if
    end do
    if (desc_forces_local) close (funit)
  end if                  ! rangml==0
  _MLD_END_
end subroutine write_force_descriptors




subroutine write_force_descriptors_dump(iconf)

#if(PARA)
  use mpi
  use mld_mpi 
#endif
  use set_limits, only: i_start_on_proc, i_end_on_proc, i_size_on_proc
  use ml_in_ndm_module, only: desc_forces, rangml, &
                              imm_neigh, i_start_at, i_final_at
  use module_db_setup, only: db_path                            
  use derived_types, only: config_desc, config_real
  use module_kind_variables, only: kind_double
  use mld_logger

  implicit none

  integer, intent(in)  :: iconf
  integer  :: ndim, iproc, inn
  character(len=100)   :: ffilename
  character(len=100)   :: CHFMTf
  integer  :: funit, ia, ia_start, ia_end, i_pack
  real(kind_double), dimension(:, :, :), allocatable :: s_forcex, s_forcey, s_forcez
  real(kind_double), dimension(:, :, :), allocatable :: r_forcex, r_forcey, r_forcez
  logical  :: desc_forces_local

  _NAMECURRENT_("write_force_descriptors_dump")



  _MLD_BEGIN_

#if(PARA)

  if ((i_start_at == 0) .and. (i_final_at == 0)) return
  desc_forces_local = desc_forces .and. (config_real(iconf)%has_force .or. config_real(iconf)%has_stress)
  if (.not. (desc_forces_local)) return
  ndim = config_desc(iconf)%dim_desc
  ffilename = 'dump'//trim(adjustl(db_path))//config_real(iconf)%class//'_'//config_real(iconf)%klm//'_'//config_real(iconf)%cnumber//'.dfml'

  !there is only one procs who writing
  if (rangml /= 0) then
  do iproc = 1, nb_procsml - 1
    if (rangml == iproc) then
      ia_start = i_start_on_proc(iproc)
      ia_end = i_end_on_proc(iproc)
      i_pack = i_size_on_proc(iproc)
      if ((ia_start == 0) .and. (ia_end == 0)) cycle
      if (allocated(s_forcex)) deallocate (s_forcex); allocate (s_forcex(ndim, i_pack, imm_neigh + 1))
      if (allocated(s_forcey)) deallocate (s_forcey); allocate (s_forcey(ndim, i_pack, imm_neigh + 1))
      if (allocated(s_forcez)) deallocate (s_forcez); allocate (s_forcez(ndim, i_pack, imm_neigh + 1))
      do ia = ia_start, ia_end
        do inn = 0, config_desc(iconf)%n_neigh(ia)
          s_forcex(:, ia - ia_start + 1, inn + 1) = config_desc(iconf)%force(:, ia, inn, 1)
          s_forcey(:, ia - ia_start + 1, inn + 1) = config_desc(iconf)%force(:, ia, inn, 2)
          s_forcez(:, ia - ia_start + 1, inn + 1) = config_desc(iconf)%force(:, ia, inn, 3)
        end do
      end do

      call MPI_SEND(s_forcex, size(s_forcex, 1)*size(s_forcex, 2)*size(s_forcex, 3), MPI_DOUBLE_PRECISION, 0, 100 + iproc, mpi_comm_mld, codeml)
      call MPI_SEND(s_forcey, size(s_forcey, 1)*size(s_forcey, 2)*size(s_forcey, 3), MPI_DOUBLE_PRECISION, 0, 200 + iproc, mpi_comm_mld, codeml)
      call MPI_SEND(s_forcez, size(s_forcez, 1)*size(s_forcez, 2)*size(s_forcez, 3), MPI_DOUBLE_PRECISION, 0, 300 + iproc, mpi_comm_mld, codeml)
      deallocate (s_forcex, s_forcey, s_forcez)
    end if
  end do
  end if

  if (rangml == 0) then
    funit = 42
    write (CHFMTf, *) '(i9, 1x, ', int(3*ndim), 'e25.15)'
    open (funit, file=ffilename, status='unknown')   ! , position='write')
    do iproc = 0, nb_procsml - 1
      ia_start = i_start_on_proc(iproc)
      ia_end = i_end_on_proc(iproc)
      i_pack = i_size_on_proc(iproc)
      if ((ia_start == 0) .and. (ia_end == 0)) cycle
      if (iproc == 0) then
        do ia = ia_start, ia_end
          write (funit, '(i9,i6)') ia, config_desc(iconf)%n_neigh(ia)
          do inn = 0, config_desc(iconf)%n_neigh(ia)
            write (funit, FMT=CHFMTf) inn, config_desc(iconf)%force(:, ia, inn, 1), &
              config_desc(iconf)%force(:, ia, inn, 2), &
              config_desc(iconf)%force(:, ia, inn, 3)
          end do                  ! inn
        end do                  ! ia
      else
        if (allocated(r_forcex)) deallocate (r_forcex); allocate (r_forcex(ndim, i_pack, imm_neigh + 1))
        if (allocated(r_forcey)) deallocate (r_forcey); allocate (r_forcey(ndim, i_pack, imm_neigh + 1))
        if (allocated(r_forcez)) deallocate (r_forcez); allocate (r_forcez(ndim, i_pack, imm_neigh + 1))
        call MPI_RECV(r_forcex, size(r_forcex, 1)*size(r_forcex, 2)*size(r_forcex, 3), MPI_DOUBLE_PRECISION, iproc, 100 + iproc, mpi_comm_mld, statut_ml, codeml)
        call MPI_RECV(r_forcey, size(r_forcey, 1)*size(r_forcey, 2)*size(r_forcey, 3), MPI_DOUBLE_PRECISION, iproc, 200 + iproc, mpi_comm_mld, statut_ml, codeml)
        call MPI_RECV(r_forcez, size(r_forcez, 1)*size(r_forcez, 2)*size(r_forcez, 3), MPI_DOUBLE_PRECISION, iproc, 300 + iproc, mpi_comm_mld, statut_ml, codeml)

        do ia = ia_start, ia_end
          write (funit, '(i9,i6)') ia, config_desc(iconf)%n_neigh(ia)
          do inn = 1, config_desc(iconf)%n_neigh(ia) + 1
            write (funit, FMT=CHFMTf) inn - 1, r_forcex(:, ia - ia_start + 1, inn), & ! config_desc(iconf)%force(:,ia,inn,1), &
              r_forcey(:, ia - ia_start + 1, inn), &           ! config_desc(iconf)%force(:,ia,inn,2), &
              r_forcez(:, ia - ia_start + 1, inn)              ! config_desc(iconf)%force(:,ia,inn,3)
          end do                  ! inn
        end do
        deallocate (r_forcex, r_forcey, r_forcez)
      end if
    end do
    close (funit)
 end if
  _MLD_END_
  return

#else
  write (6, '("ML:.............TODO Houston not MPI in   :",a)') _NAMECURRENT_
  return
#endif

end subroutine write_force_descriptors_dump



subroutine read_force_descriptors_dump(iconf)

#if(PARA)
  use mpi
  use mld_mpi 
#endif

  use ml_in_ndm_module, only: desc_forces, rangml, i_start_at, i_final_at
  use module_db_setup, only: db_path
  use derived_types, only: config_desc, config_real
  use module_kind_variables, only: kind_double
  use mld_logger

  implicit none

  integer, intent(in)  :: iconf
  integer  :: ndim, iproc, inn
  character(len=100)   :: ffilename
  character(len=100)   :: CHFMTf
  integer  :: funit, ia
  logical  :: desc_forces_local, ok
  integer  :: ia_temp, inn_temp, nn_temp
  real(kind_double), allocatable, dimension(:) :: tempx, tempy, tempz


  _NAMECURRENT_("read_force_descriptors_dump")



  _MLD_BEGIN_
#if(PARA)

  if ((i_start_at == 0) .and. (i_final_at == 0)) return
  desc_forces_local = desc_forces .and. (config_real(iconf)%has_force .or. config_real(iconf)%has_stress)
  if (.not. (desc_forces_local)) return
  config_desc(iconf)%force(:, :, :, :) = 0.d0
  ndim = config_desc(iconf)%dim_desc
  allocate (tempx(ndim), tempy(ndim), tempz(ndim))
  ffilename = 'dump'//trim(adjustl(db_path))//config_real(iconf)%class//'_'//config_real(iconf)%klm//'_'//config_real(iconf)%cnumber//'.dfml'
  if (rangml == 0) then
    inquire (file=ffilename, exist=ok)
    if (.not. (ok)) then
      write (6, '("The file for reading ",a,"is missing. What a bad day!")') ffilename
      !call end_ml("stop with reading in read_force_descriptors_dump")
    end if
  end if                  ! rangml == 0

  ! there is only one procs who writting ...
  funit = 42
  write (CHFMTf, *) '(i9, 1x, ', int(3*ndim), 'e25.15)'

  do iproc = 0, nb_procsml - 1
    if (rangml == iproc) then
      !do ia=1,config_real(iconf)%nat
      open (funit, file=ffilename, status='unknown')

      do ia = 1, config_real(iconf)%nat
        read (funit, '(i9,i6)') ia_temp, nn_temp         ! config_desc(iconf)%n_neigh(ia)
        if ((ia >= i_start_at) .and. (ia <= i_final_at)) then
          if (nn_temp /= config_desc(iconf)%n_neigh(ia)) then
            write (6, *) 'Problems in the list of neighbours ... possible not good forces', nn_temp, config_desc(iconf)%n_neigh(ia)
          end if
        end if
        do inn = 0, nn_temp
          read (funit, FMT=CHFMTf) inn_temp, tempx(:), tempy(:), tempz(:)
          if ((ia >= i_start_at) .and. (ia <= i_final_at)) then
            config_desc(iconf)%force(:, ia, inn, 1) = tempx(:)
            config_desc(iconf)%force(:, ia, inn, 2) = tempy(:)
            config_desc(iconf)%force(:, ia, inn, 3) = tempz(:)
          end if
        end do                  ! inn
      end do                  ! ia

      close (funit)
    end if
  end do
  deallocate (tempx, tempy, tempz)
  _MLD_END_
  return

#else
  write (6, '("ML:.............TODO Houston not MPI in   :",a)') _NAMECURRENT_
  return
#endif

end subroutine read_force_descriptors_dump


#if(MLD_HDF5)
subroutine write_hdf5_energy(iconf, grpname)
   use mld_hdf5
   use derived_types,   only : config_real, config_desc
   implicit none

   integer, intent(in) :: iconf
   character(len=*), intent(in) :: grpname
   !-----------------------------------------------------------------
   integer               :: nat, ndesc, ii 
   integer, allocatable  :: atom_id(:)
   real(kind=kind(1.d0)), allocatable, dimension(:,:) :: desc
   !-----------------------------------------------------------------


   !$! ! 2. create the configuration group
   !$! write(grpname,'(A,"_",A,"_",A)') trim(config_real(iconf)%class),  &
   !$!                                  trim(config_real(iconf)%klm),    &
   !$!                                  trim(config_real(iconf)%cnumber)
   call hFile%hdf5_create_group(trim(grpname))

   ! 3. gather data from config_real(iconf)
   nat   = config_real(iconf)%nat
   ndesc = size(config_desc(iconf)%energy, 1)
   allocate(atom_id(nat))
   allocate(desc(ndesc, nat))
   do ii = 1, nat 
    atom_id(ii)  = ii 
   end do  
   desc(:,:)    = config_desc(iconf)%energy(:,:)

   ! 4. write attribute and datasets
   call hFile%hdf5_write_integer_group_attribute('natoms', nat)
   call hFile%hdf5_write_integer_group_dataset  ('atom_id',  nat,   atom_id)
   call hFile%hdf5_write_real_group_2Ddataset   ('descriptor', ndesc, nat, desc)

   deallocate(atom_id, desc)


end subroutine write_hdf5_energy
#else
subroutine write_hdf5_energy(iconf, grpname)
  use mld_logger
  implicit none

  integer, intent(in) :: iconf
  character(len=*), intent(in) :: grpname
  call log_critical("ML:.............TODO Houston not HDF5 support ")
  stop "not HDF5 support in write_hdf5_energy"
end subroutine write_hdf5_energy
#endif





! OLD STUFF

subroutine serial_write_descriptors(iconf)

  use ml_in_ndm_module, only: desc_forces, rangml
  use module_db_setup, only: db_path
  use derived_types, only: config_desc, config_real
  use mld_force_mod, only: pack_force_descriptor
  use mld_logger

  implicit none

  integer, intent(in)  :: iconf
  integer  :: ndim
  character(len=100)   :: ffilename
  character(len=60)    :: CHFMTf
  integer  :: funit, ia
  logical  :: desc_forces_local

  _NAMECURRENT_("serial_write_descriptors")



  _MLD_BEGIN_
  desc_forces_local = desc_forces .and. (config_real(iconf)%has_force .or. config_real(iconf)%has_stress)
  ndim = config_desc(iconf)%dim_desc
  !e efilename='desc'//trim(adjustl(db_path))//config_real(iconf)%class//'_'//config_real(iconf)%klm//'_'//config_real(iconf)%cnumber//'.eml'
  if (desc_forces_local) then
    ffilename = 'desc'//trim(adjustl(db_path))//config_real(iconf)%class//'_'//config_real(iconf)%klm//'_'//config_real(iconf)%cnumber//'.fml'
    call pack_force_descriptor(iconf)
  end if

  if (rangml == 0) then   ! there is only one procs who writting
    ! eunit=41
    funit = 42
    ! open(eunit,file=efilename,status='unknown')
    if (desc_forces_local) open (funit, file=ffilename, status='unknown')
    ! write (CHFMT,*)'(i9, 1x, ',int(ndim),'e20.10)'
    !old_not_erase write (CHFMTf,*)'(i6, 1x,i5,1x,i6, ',int(ndim),'e20.10)'
    write (CHFMTf, *) '(e25.15, 1x, ', int(3*ndim), 'e25.15)'
    do ia = 1, config_real(iconf)%nat
      ! write (eunit,FMT=CHFMT) ia, real(config_desc(iconf)%energy(:,ia))
      if (desc_forces_local) then
        !old_not_erase      write (funit,'(i6)') ia
        !old_not_erase      do ix=1,3
        !old_not_erase       write (funit,FMT=CHFMTf) ia, ix, 0,  real(config_desc(iconf)%force(:,ia,0,ix))
        !old_not_erase      end do
        !old_not_erase      write (funit,'(i6)') ia
        write (funit, FMT=CHFMTf) dble(ia), real(config_desc(iconf)%pack_force_linear(:, 1, ia)), &
          real(config_desc(iconf)%pack_force_linear(:, 2, ia)), &
          real(config_desc(iconf)%pack_force_linear(:, 3, ia))
      end if
      !old_not_erase if (desc_forces) then
      !old_not_erase do ja=1,config_desc(iconf)%n_neigh(ia)
      !old_not_erase   do ix=1,3
      !write (funit,FMT=CHFMTf) config_desc(iconf)%kind_neigh(ia,ja), config_real(iconf)%u_ij(ia,ja,ix), config_desc(iconf)%force(:,ia,ja,ix)
      !old_not_erase      write (funit,FMT=CHFMTf) config_desc(iconf)%kind_neigh(ia,ja), ix, ja, config_desc(iconf)%force(:,ia,ja,ix)
      !old_not_erase   end do
      !old_not_erase end do
      !old_not_erase end if
    end do

    ! close(eunit)
    if (desc_forces_local) close (funit)
  end if                  ! rangml==0

  _MLD_END_
end subroutine serial_write_descriptors
