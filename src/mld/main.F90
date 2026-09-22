! see READMES and docs. For more details, please contact: mihai-cosmin.marinica@cea.fr

! cmake3 compilation
! cd ${MLD_BUIDIR} ; f_compile_milady ; make -j10
! first small test launch
! cd ${MLD_TESDIR}/NDM_2020/md_cg_ml ; mpirun -n 2 ${MLD_BUIDIR}/bin/milady_main.exe

program milady_main

   use iso_fortran_env
   use mld_mpi
   use mld_string
   use mld_banner
   use mld_logger
   use mld_unit
   use module_end_ml, only: end_ml
   use module_ml_scalapack, only: scalapack_driver, context 
   use module_gin_in_hard, only: gin_in_soft_get_imm
#if(MLD_HDF5)   
   use mld_hdf5
#endif    

   use ondm_gen_com_m, only: fnam
   !use ondm_tab_imm_m, only: xpp

   implicit none

   integer  :: file_size
   logical  :: ok

   call mld_mpi_init()

#if(MLD_HDF5)   
   call hdf5_init()
#endif    
   call log_info(get_banner(compiler_options()))

   ! here comes milady stuff
   ! TODO use MPI io read_all ? etc.
   inquire(file='name.in', exist=ok, size=file_size)
   if (.not. ok) then
      call log_critical("File 'name.in' not found.")
      call log_critical("Please provide 'name.in'.")
      call mld_mpi_finalize("Please provide 'name.in'.")
      stop
   else if (file_size == 0) then
      call log_critical("File 'name.in' is empty.")
      call log_critical("Please provide a valid 'name.in'.")
      call mld_mpi_finalize("Please provide a valid 'name.in'.")
      stop
   endif
   fnam = read_first_line('name.in')

   call ondm_readdm()

   call read_mld_md 
   call gin_in_soft_get_imm

   call prog()

   ! TODO seems stop program in prog
   call log_debug('End regular of program milady_main')
#if(MLD_HDF5)   
   call hdf5_finalize()
#endif    
   !call mld_mpi_finalize()
   call end_ml("----Full end in milady main----", scalapack_driver, context, mld_rank)
end program


subroutine prog()
   use mld_logger
   use ondm_gen_com_m, only: imm, dmtype
   use ondm_tab_imm_m, only: alloc_all_tab_imm
   use mld_mpi, only: mld_mpi_abort, mld_rank
   use ml_in_ndm_module, only: mld_dmtype, MD_MLD_DMTYPE, ML_MLD_DMTYPE
   !use module_ml_scalapack, only: scalapack_driver, context
   use module_end_ml, only: end_ml

   implicit none

   ! Allocation des tableaux dimensionnes sur le nombre d'atomes
   call alloc_all_tab_imm(imm)
   dmtype = mld_dmtype

   call init

   
   select case (dmtype)
    case (ML_MLD_DMTYPE) ! this 18 
      call ml
    case (MD_MLD_DMTYPE) ! this is 181 
       call log_info("ML: MD_mode let's go!")
       !call mld_md !and do your bidule
       call rrrroule_ma_poule
       !call end_ml("----Full MD end----", scalapack_driver, context, mld_rank)
    case default
       call log_critical("ML: critical dmtype. THis dmtype not available: "//vtoa(dmtype))
       call log_critical("We will stop. Meanwhile: Take a cofee. Ask Master. Pray.")
       call mld_mpi_abort("Bye bye in prog!")
   end select

end subroutine prog
