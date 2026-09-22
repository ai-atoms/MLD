enable_language( Fortran )

# order no matter.
target_sources( MILADY  PRIVATE
  class_functions.F90
  module_minimizer.F90
  module_cgext.F90
  module_lbfgs.F90
  module_condition_number.F90
)
