PROGRAM dbcsr_rectangular_benchmark
   USE dbcsr_error_handling,            ONLY: dbcsr_error_handling_setup
   USE dbcsr_api,                       ONLY: dbcsr_finalize_lib,&
                                              dbcsr_init_lib
   USE dbcsr_types, ONLY: dbcsr_mp_obj
   USE dbcsr_methods,                   ONLY: dbcsr_mp_release
   USE dbcsr_rectangular_base,          ONLY: dbcsr_r_destroy
   USE dbcsr_rectangular_base_types,    ONLY: dbcsr_r_type
   USE dbcsr_rectangular_split,         ONLY: dbcsr_r_mp_comm
   USE dbcsr_rectangular_test,          ONLY: dbcsr_r_benchmark_mm,&
                                              dbcsr_r_random_bsizes,&
                                              dbcsr_r_setup_test_matrix
   USE dbcsr_kinds,                           ONLY: int_8,&
                                              real_8,&
                                              default_string_length
   USE dbcsr_machine,                         ONLY: default_output_unit
   USE dbcsr_mpiwrap,                 ONLY: add_mp_perf_env,&
                                              mp_comm_free,&
                                              mp_environ,&
                                              mp_world_finalize,&
                                              mp_world_init,&
                                              rm_mp_perf_env
   !USE string_table,                    ONLY: string_table_allocate,&
   !                                           string_table_deallocate
   USE dbcsr_timings,                         ONLY: add_timer_env,&
                                              rm_timer_env,&
                                              timings_register_hooks
   USE dbcsr_timings_report,                  ONLY: cost_type_time,&
                                              timings_report_print
   USE dbcsr_rectangular_util,          ONLY: dbcsr_mp_environ
   USE dbcsr_base_hooks,                ONLY: timeset, timestop
#include "base/dbcsr_base_uses.f90"

   IMPLICIT NONE

   INTEGER(KIND=int_8) :: dmin, dmax, d1, d2
   CHARACTER(len=default_string_length) :: arg
   INTEGER :: bsize
   REAL(KIND=real_8) :: sparsity
   LOGICAL :: compare_dbcsr = .TRUE.

   INTEGER, DIMENSION(:), ALLOCATABLE :: bsize_max
   INTEGER, DIMENSION(:), ALLOCATABLE :: bsize_min
   TYPE(dbcsr_r_type) :: A, B, C
   INTEGER :: io_unit, mp_comm, mynode, mp_comm_cart, mp_comm_cart_, numnodes
   INTEGER, DIMENSION(2) :: pcoord, pdims
   TYPE(dbcsr_mp_obj)            :: mp_env
   INTEGER :: handle
   REAL(KIND=real_8), PARAMETER :: filter_eps=1.0E-08

   !mp_collect_timings = .TRUE.
   CALL mp_world_init(mp_comm)

   CALL get_command_argument(1, arg)
   READ(arg, *) d1
   CALL get_command_argument(2, arg)
   READ(arg, *) d2
   dmax = MAX(d1, d2)
   dmin = MIN(d1, d2)
   CALL get_command_argument(3, arg)
   READ(arg, *) bsize
   CALL get_command_argument(4, arg)
   READ(arg, *) sparsity
   DBCSR_ASSERT(sparsity .LE. 1.0_real_8)
   CALL get_command_argument(5, arg)
   READ(arg, *) compare_dbcsr

   mp_comm_cart = dbcsr_r_mp_comm(mp_comm, dmax, dmin)
   mp_env = dbcsr_mp_environ(mp_comm_cart)

   !CALL dbcsr_error_handling_setup()
   CALL timings_register_hooks()
   !CALL string_table_allocate()
   CALL add_mp_perf_env()
   CALL add_timer_env()

   CALL mp_environ(numnodes, mynode, mp_comm)

   CALL dbcsr_init_lib()

   io_unit = 0
   IF (mynode .EQ. 0) io_unit = default_output_unit

   IF(io_unit > 0) THEN
      WRITE(io_unit, *) "dimensions", dmax, dmin
      WRITE(io_unit, *) "block sizes", bsize
      WRITE(io_unit, *) "sparsity", sparsity
      WRITE(io_unit, *) "compare with dbcsr?", compare_dbcsr
   ENDIF


   CALL timeset("dbcsr_rectangular_benchmark", handle)


   CALL mp_environ(numnodes, pdims, pcoord, mp_comm_cart)

   ALLOCATE(bsize_max(dmax))
   ALLOCATE(bsize_min(dmin))
   bsize_max(:) = bsize
   bsize_min(:) = bsize
   !CALL dbcsr_r_random_bsizes([bsize], 1, bsize_max)
   !CALL dbcsr_r_random_bsizes([bsize], 1, bsize_min)

   CALL dbcsr_r_setup_test_matrix(A, mp_comm_cart_, mp_comm_cart, dmax, dmin, bsize_max, bsize_min, [5, 1], "A", &
                                  sparsity, reuse_comm=.TRUE.)
   CALL dbcsr_r_setup_test_matrix(B, mp_comm_cart_, mp_comm_cart, dmax, dmin, bsize_max, bsize_min, [5, 1], "B", &
                                  sparsity, reuse_comm=.TRUE.)
   CALL dbcsr_r_setup_test_matrix(C, mp_comm_cart_, mp_comm_cart, dmin, dmin, bsize_min, bsize_min, [2, 3], "C", &
                                  sparsity, reuse_comm=.TRUE.)

   CALL dbcsr_r_benchmark_mm('T', 'N', 'N', A, B, C, compare_dbcsr, filter_eps=filter_eps)

   CALL dbcsr_r_destroy(A)
   CALL dbcsr_r_destroy(B)
   CALL dbcsr_r_destroy(C)

   CALL timestop(handle)

   CALL dbcsr_finalize_lib(mp_comm, io_unit)

   CALL timings_report_print(io_unit, 0.0_real_8, .FALSE., cost_type_time, .TRUE., mp_env)

   CALL dbcsr_mp_release(mp_env)

   CALL mp_comm_free(mp_comm_cart)
   ! Remove timers
   CALL rm_mp_perf_env()
   CALL rm_timer_env()

   CALL mp_world_finalize()

END PROGRAM
