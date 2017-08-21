# set_unless_def(ENV_VAL, DEFAULT_VAL)
#   If ENV_VAL is not defined, this command sets ENV_VAL with DEFAULT_VAL.
macro (set_unless_def ENV_VAL DEFAULT_VAL)
  if (NOT DEFINED ${ENV_VAL})
    set(${ENV_VAL} ${DEFAULT_VAL})
  endif ()
endmacro (set_unless_def)

# add_definitions_if(ENV_VAL, DEF_VAL)
#   If ENV_VAL is defined, this command calls add_definitions command.
macro (add_definitions_if ENV_VAL DEF_VAL)
  if (${ENV_VAL})
    add_definitions(${DEF_VAL} ${ARGN})
  endif ()
endmacro (add_definitions_if)

# add_definitions_unless(ENV_VAL, DEF_VAL)
#   If ENV_VAL `is not` defined, this command calls add_definictions command.
macro (add_definitions_unless ENV_VAL DEF_VAL)
  if (NOT ${ENV_VAL})
    add_definitions(${DEF_VAL} ${ARGN})
  endif ()
endmacro (add_definitions_unless)

# check_mpi_compiler(COMPILER_NAME, RESULT)
#   This command checks the prefix of ${COMPILER_NAME} is `mpi`.
#   In almost all of MPI compiler, they takes `mpi` prefix. 
#
#   in:  COMPILER_NAME
#   out: ${RESULT} (TRUE or FALSE)
function (check_mpi_compiler COMPILER_NAME RESULT)
  get_filename_component(COMPILER_NAME ${COMPILER_NAME} NAME)
  string(LENGTH ${COMPILER_NAME} NAME_LEN)
  if (${NAME_LEN} LESS 3)
    set(RET FALSE)
  else ()
    string(SUBSTRING ${COMPILER_NAME} 0 3 COMPILER_HEADER)
    string(TOLOWER ${COMPILER_HEADER} COMPILER_HEADER)
    string(COMPARE EQUAL ${COMPILER_HEADER} "mpi" RET)
  endif ()
  set(${RESULT} ${RET} PARENT_SCOPE)
endfunction (check_mpi_compiler)

# option_set(ENV_VAL, HELP_STR, DEFAULT_VAL)
#   If option is not declare, this command initializes the option with ${DEFAULT_VAL}
macro (option_set ENV_VAL HELP_STR DEFAULT_VAL)
  if (${${ENV_VAL}})
    option(${ENV_VAL} ${HELP_STR} ${${ENV_VAL}})
  else ()
    option(${ENV_VAL} ${HELP_STR} ${DEFAULT_VAL})
  endif()
endmacro (option_set)
