
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

#include "../../MLD_MACROS.INC"


#if(MLD_HDF5)
module mld_hdf5
   !------------------------------------------------------------
   ! Every hdf5 methods needed to read a file
   !------------------------------------------------------------
   use HDF5
   use mpi
   use mld_logger
   implicit none

   public :: HDF5File
   type :: HDF5File
      !private
      character(len=80) :: fname
      integer(hid_t)     :: fid = -1   ! File ID
      integer(hid_t)     :: gid = -1   ! Group ID
      integer(hid_t)     :: did = -1   ! Dataset ID
      integer(hid_t)     :: aid = -1   ! Attribute ID
      integer            :: file_open  ! 0 = closed, 1 = open
      !------------------------------------------------------------
   contains
      ! Initialization and finalization
      procedure :: hdf5_open_file
      procedure :: hdf5_close_file
      procedure :: hdf5_open_group
      ! Read operations
      procedure :: hdf5_get_integer_group_attribute
      procedure :: hdf5_get_real_group_attribute
      procedure :: hdf5_get_integer_group_dataset
      procedure :: hdf5_get_real_group_dataset
      procedure :: hdf5_get_real_group_2Ddataset
      procedure :: hdf5_get_string_group_dataset
      procedure :: hdf5_flush
      procedure :: hdf5_get_filename
      ! Write operations
      procedure :: hdf5_create_file
      procedure :: hdf5_create_group
      procedure :: hdf5_write_integer_group_attribute
      procedure :: hdf5_write_real_group_attribute
      procedure :: hdf5_write_integer_group_dataset
      procedure :: hdf5_write_real_group_dataset
      procedure :: hdf5_write_real_group_2Ddataset
      procedure :: hdf5_write_string_group_dataset
   end type HDF5File

   type(HDF5File) :: hFile
   logical ::  db_hdf5 

contains

   subroutine hdf5_init()
      implicit none
      integer :: hdferr
      ! initialize fortran hdf5
      call h5open_f(hdferr)
   end subroutine hdf5_init

   subroutine hdf5_finalize()
      implicit none
      integer :: hdferr
      ! close FORTRAN interface
      call h5close_f(hdferr)
   end subroutine hdf5_finalize

   subroutine hdf5_get_filename(this, fname)
      implicit none
      class(HDF5File), intent(in) :: this
      character(len=80), intent(inout) :: fname
      if(this%fid == -1) then
         call log_critical('ML HDF5: Fatal Error Trying to get name of a file that has not been opened yet')
         stop "fatal hdf5"
      end if
      fname = this%fname
   end subroutine hdf5_get_filename

   subroutine hdf5_open_file(this, fname)
      implicit none
      class(HDF5File), intent(inout) :: this
      character(len=80), intent(in) :: fname
      integer :: hdferr
      if (this%fid /= -1) then
         call log_critical('ML HDF5: Fatal Error Opening new file when previous one has not been closed')
         stop "fatal hdf5"
      end if
      this%fname = fname
      call h5fopen_f(fname, H5F_ACC_RDONLY_F, this%fid, hdferr)
   end subroutine hdf5_open_file

   subroutine hdf5_open_group(this, gname)
      implicit none
      class(HDF5File), intent(inout) :: this
      character(len=80), intent(in) :: gname
      integer :: hdferr
      call h5gopen_f(this%fid, gname, this%gid, hdferr)
   end subroutine hdf5_open_group

   subroutine hdf5_get_integer_group_attribute(this, aname, attribute)
      implicit none
      class(HDF5File), intent(inout) :: this
      character(len=*), intent(in) :: aname
      integer, intent(out)    :: attribute
      integer(hsize_t), dimension(1) :: dims
      integer :: hdferr
      if(this%gid == -1) then
         call log_critical('ML HDF5: Fatal Error Trying to read attribute on group without opening group first')
         stop "fatal hdf5"
      end if
      dims = 1
      call h5aopen_name_f(this%gid, aname, this%aid, hdferr)
      call h5aread_f(this%aid, H5T_NATIVE_INTEGER, attribute, dims, hdferr)
   end subroutine hdf5_get_integer_group_attribute

   subroutine hdf5_get_real_group_attribute(this, aname, attribute)
      implicit none
      class(HDF5File), intent(inout) :: this
      character(len=*), intent(in) :: aname
      real(kind=kind(1.d0)), intent(inout)    :: attribute
      integer(hsize_t), dimension(1) :: dims
      integer :: hdferr
      if(this%gid == -1) then
         call log_critical('ML HDF5: Fatal Error Trying to read attribute on group without opening group first')
         stop "fatal hdf5"
      end if
      dims = 1
      call h5aopen_name_f(this%gid, aname, this%aid, hdferr)
      call h5aread_f(this%aid, H5T_NATIVE_DOUBLE, attribute, dims, hdferr)
   end subroutine hdf5_get_real_group_attribute

   subroutine hdf5_get_integer_group_dataset(this, dname, dset_sz, dataset)
      implicit none
      class(HDF5File), intent(inout) :: this
      character(len=*), intent(in) :: dname
      integer, intent(in) :: dset_sz
      integer, dimension(:), intent(inout) :: dataset
      integer(hsize_t), dimension(1) :: dims
      integer :: hdferr
      if(this%gid == -1) then
         call log_critical('ML HDF5: Fatal Error Trying to read dataset on group without opening group first')
         stop "fatal hdf5"
      end if
      dims = dset_sz
      call h5dopen_f(this%gid, dname, this%did, hdferr)
      call h5dread_f(this%did, H5T_NATIVE_INTEGER, dataset, dims, hdferr)
   end subroutine hdf5_get_integer_group_dataset

   subroutine hdf5_get_real_group_dataset(this, dname, dset_sz, dataset)
      implicit none
      class(HDF5File), intent(inout) :: this
      character(len=*), intent(in) :: dname
      integer, intent(in) :: dset_sz
      real(kind=kind(1.d0)), dimension(:), intent(inout) :: dataset
      integer(hsize_t), dimension(1) :: dims
      integer :: hdferr
      if(this%gid == -1) then
         call log_critical('ML HDF5: Fatal Error Trying to read dataset on group without opening group first')
         stop "fatal hdf5"
      end if
      dims = dset_sz
      call h5dopen_f(this%gid, dname, this%did, hdferr)
      call h5dread_f(this%did, H5T_NATIVE_DOUBLE, dataset, dims, hdferr)
   end subroutine hdf5_get_real_group_dataset

   subroutine hdf5_get_real_group_2Ddataset(this, dname, dset_nl, dset_nc, dataset)
      implicit none
      class(HDF5File), intent(inout) :: this
      character(len=*), intent(in) :: dname
      integer, intent(in) :: dset_nl, dset_nc
      real(kind=kind(1.d0)), dimension(:,:), intent(inout) :: dataset
      integer(hsize_t), dimension(2) :: dims
      integer :: hdferr
      if(this%gid == -1) then
         call log_critical('ML HDF5: Fatal Error Trying to read dataset on group without opening group first')
         stop "fatal hdf5"
      end if
      dims(1) = dset_nl; dims(2) = dset_nc
      call h5dopen_f(this%gid, dname, this%did, hdferr)
      call h5dread_f(this%did, H5T_NATIVE_DOUBLE, dataset, dims, hdferr)
   end subroutine hdf5_get_real_group_2Ddataset

   subroutine hdf5_get_string_group_dataset(this, dname, dset_sz, string_sz, dataset)
      implicit none
      class(HDF5File), intent(inout) :: this
      character(len=*), intent(in) :: dname
      integer, intent(in) :: dset_sz
      integer(hsize_t), intent(in) :: string_sz
      character(len=string_sz), dimension(:), intent(inout) :: dataset
      integer(hsize_t), dimension(1) :: dims
      integer(hid_t) :: memtype
      integer :: hdferr
      if(this%gid == -1) then
         call log_critical('ML HDF5: Fatal Error Trying to read dataset on group without opening group first')
         stop "fatal hdf5"
      end if
      ! create hdf5 datatype (ie string of size string_sz)
      call h5tcopy_f(H5T_FORTRAN_S1, memtype, hdferr)
      call h5tset_size_f(memtype, string_sz, hdferr)
      ! open and read dataset content
      dims = dset_sz
      call h5dopen_f(this%gid, dname, this%did, hdferr)
      call h5dread_f(this%did, memtype, dataset, dims, hdferr)
   end subroutine hdf5_get_string_group_dataset

   subroutine hdf5_close_file(this)
      implicit none
      class(HDF5File), intent(inout) :: this
      integer :: hdferr
      if(this%fid == -1) then
         call log_critical('ML HDF5: Fatal Error Trying to close file that has not be opened')
         stop "fatal hdf5"
      end if
      if(this%aid /= -1) then 
         call h5aclose_f(this%aid, hdferr)
         this%aid = -1 
      end if    
      if(this%did /= -1) call h5dclose_f(this%did, hdferr)
      if(this%gid /= -1) call h5gclose_f(this%gid,hdferr)
      call h5fclose_f(this%fid,hdferr)
      call this%hdf5_flush()
   end subroutine hdf5_close_file

   subroutine hdf5_flush(this)
      implicit none
      class(HDF5File), intent(inout) :: this
      this%fid = -1
      this%gid = -1
      this%aid = -1
      this%did = -1
   end subroutine hdf5_flush

   !----------------------------------------------------------------
   ! Write interface implementations (updated)
   !----------------------------------------------------------------

  subroutine hdf5_create_file(this, fname_in)
      class(HDF5File), intent(inout) :: this
      character(len=*), intent(in)   :: fname_in
      integer :: hdferr
      this%fname = adjustl(fname_in)
      call h5fcreate_f(this%fname, H5F_ACC_TRUNC_F, this%fid, hdferr)
   end subroutine

   subroutine hdf5_create_group(this, gname)
      class(HDF5File), intent(inout) :: this
      character(len=*), intent(in)   :: gname
      integer :: hdferr
      call h5gcreate_f(this%fid, trim(gname), this%gid, hdferr)
   end subroutine

   subroutine hdf5_write_integer_group_attribute(this, aname, value)
     class(HDF5File), intent(inout) :: this
     character(len=*), intent(in)   :: aname
     integer,           intent(in)  :: value
     integer(hid_t)  :: space_id
     integer(hsize_t) :: dims(1)
     integer          :: hdferr
     !--- 1. build a scalar dataspace ---------------------------------
     dims(1) = 1_hsize_t
     call h5screate_simple_f(1, dims, space_id, hdferr)
     !--- 2. create the attribute using that dataspace ----------------
     call h5acreate_f(this%gid, trim(aname), H5T_NATIVE_INTEGER, space_id, &
                      this%aid, hdferr)
     !--- 3. write the value ------------------------------------------
     !call h5awrite_f(this%aid, H5T_NATIVE_INTEGER, value, hdferr)
     !– write value (legacy 5-arg form) –
     call h5awrite_f(this%aid, H5T_NATIVE_INTEGER, value, dims, hdferr)
     !--- 4. close everything -----------------------------------------
     call h5aclose_f(this%aid, hdferr)
     this%aid = -1
     call h5sclose_f(space_id, hdferr)
   end subroutine


   subroutine hdf5_write_real_group_attribute(this, aname, value)
     class(HDF5File), intent(inout) :: this
     character(len=*), intent(in)   :: aname
     real(kind=kind(1.d0)), intent(in) :: value
     integer(hsize_t) :: dims(1)
     integer(hid_t)   :: space_id          ! <-- new
     integer          :: hdferr
     !--- 1. create a scalar dataspace ---------------------------------
     dims(1) = 1_hsize_t
     call h5screate_simple_f(1, dims, space_id, hdferr)
     !--- 2. create the attribute --------------------------------------
     call h5acreate_f(this%gid, trim(aname), H5T_NATIVE_DOUBLE, space_id, &
                      this%aid, hdferr)
     !--- 3. write its value -------------------------------------------
     !call h5awrite_f(this%aid, H5T_NATIVE_DOUBLE, value, hdferr)
     !– write value (legacy 5-arg form) –
     call h5awrite_f(this%aid, H5T_NATIVE_DOUBLE, value, dims, hdferr)
     ! (use the 5-argument form with dims if your HDF5 headers require it)
     !--- 4. close things ----------------------------------------------
     call h5aclose_f(this%aid, hdferr)
     call h5sclose_f(space_id, hdferr)
   end subroutine

   subroutine hdf5_write_integer_group_dataset(this, dname, count, data)
     class(HDF5File), intent(inout) :: this
     character(len=*), intent(in)   :: dname
     integer,           intent(in)  :: count
     integer, dimension(:), intent(in) :: data
     integer(hsize_t)   :: dims(1)        ! default INTEGER
     integer(hid_t)     :: space_id, dset_id
     integer            :: hdferr
     ! 1. build a 1-D dataspace
     dims(1) = int(count, hsize_t)
     call h5screate_simple_f(1, dims, space_id, hdferr)
     ! 2. create the dataset in the current group
     call h5dcreate_f(this%gid, trim(dname), H5T_NATIVE_INTEGER, space_id,  &
                      dset_id, hdferr)
     ! 3. write the data  (legacy 5-arg form)
     call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data, dims, hdferr)
     ! 4. close handles
     call h5dclose_f(dset_id, hdferr)
     call h5sclose_f(space_id, hdferr)
   end subroutine hdf5_write_integer_group_dataset

   subroutine hdf5_write_real_group_dataset(this, dname, count, data)
     class(HDF5File), intent(inout) :: this
     character(len=*), intent(in)   :: dname
     integer,           intent(in)  :: count
     real(kind=kind(1.d0)), dimension(:), intent(in) :: data
     integer(hsize_t) :: dims(1)          ! HSIZE_T for the dataspace
     integer(hid_t)   :: space_id, dset_id
     integer          :: hdferr
     ! 1. build a 1-D dataspace
     dims(1) = int(count, hsize_t)
     call h5screate_simple_f(1, dims, space_id, hdferr)
     ! 2. create the dataset
     call h5dcreate_f(this%gid, trim(dname), H5T_NATIVE_DOUBLE, space_id,  &
                      dset_id, hdferr)
     ! 3. write the data  (legacy 5-arg form)
     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, hdferr)
     ! 4. close handles
     call h5dclose_f(dset_id, hdferr)
     call h5sclose_f(space_id, hdferr)
   end subroutine hdf5_write_real_group_dataset


   subroutine hdf5_write_real_group_2Ddataset(this, dname, nl, nc, data)
      class(HDF5File), intent(inout) :: this
      character(len=*), intent(in)   :: dname
      integer,           intent(in)  :: nl, nc
      real(kind=kind(1.d0)), dimension(:,:), intent(in) :: data
      integer(hsize_t) :: dims(2)          ! HSIZE_T for the dataspace
      integer(hid_t)   :: space_id, dset_id
      integer          :: hdferr
      ! 1. build a 2-D dataspace
      dims(1) = int(nl, hsize_t)
      dims(2) = int(nc, hsize_t)
      call h5screate_simple_f(2, dims, space_id, hdferr)
      ! 2. create the dataset in the current group
      call h5dcreate_f(this%gid, trim(dname), H5T_NATIVE_DOUBLE, space_id,  &
                       dset_id, hdferr)
      ! 3. write the data  (legacy 5-arg form)
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, hdferr)
      ! 4. close handles
      call h5dclose_f(dset_id, hdferr)
      call h5sclose_f(space_id, hdferr)
   end subroutine hdf5_write_real_group_2Ddataset

   subroutine hdf5_write_string_group_dataset(this, dname, count, str_len, data)
      class(HDF5File), intent(inout) :: this
      character(len=*), intent(in)           :: dname
      integer,           intent(in)          :: count
      integer,           intent(in)          :: str_len
      character(len=str_len), dimension(:), intent(in) :: data
      integer(hsize_t) :: dims(1)          ! HSIZE_T for dataspace
      integer(hid_t)   :: memtype, space_id, dset_id
      integer          :: hdferr
      ! 1. build the fixed-length string datatype
      call h5tcopy_f(H5T_FORTRAN_S1, memtype, hdferr)
      call h5tset_size_f(memtype, int(str_len, hsize_t), hdferr)
      ! 2. build a 1-D dataspace
      dims(1) = int(count, hsize_t)
      call h5screate_simple_f(1, dims, space_id, hdferr)
      ! 3. create the dataset
      call h5dcreate_f(this%gid, trim(dname), memtype, space_id, dset_id, hdferr)
      ! 4. write the data  (legacy 5-arg form)
      call h5dwrite_f(dset_id, memtype, data, dims, hdferr)
      ! 5. close handles
      call h5dclose_f(dset_id, hdferr)
      call h5sclose_f(space_id, hdferr)
      call h5tclose_f(memtype, hdferr)
   end subroutine hdf5_write_string_group_dataset



end module mld_hdf5
#else 

module mld_hdf5
    use, intrinsic :: iso_fortran_env, dp=>real64
    logical ::  db_hdf5 

end module mld_hdf5
#endif 
