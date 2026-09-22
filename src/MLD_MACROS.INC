

!!!! pay attention that macros are easily tricky, fortran is fortran, not C++


! _NAMECURRENT_() _MLD_BEGIN_ _MLD_END_ are first example of usage preprocessor MACROS
! they are designed for first-intervention-on-problem quick-user-profiling
! they implies very long verbose output log (particularly in debug verbosity)
! if it is not useful, simply modify 'mld_verbose' in mld_logger.F90 (as execution time)
! or modify macros (as compilation time)

#ifndef _NAMECURRENT_
#define _NAMECURRENT_(XX) character(len=*), parameter :: NAMECURRENT=XX
#define _MLD_BEGIN_ if (mld_verbose()) call log_debug("---> Begin "//NAMECURRENT)
#define _MLD_END_ if (mld_verbose()) call log_debug("<--- End   "//NAMECURRENT)
#endif


!!!! pay attention that macros are easily tricky, fortran is fortran, not C++
