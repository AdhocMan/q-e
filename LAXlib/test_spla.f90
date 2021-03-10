#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
#define __MSGSIZ_MAX 100000

PROGRAM lax_spla
  use laxlib_descriptor
  USE laxlib_parallel_include
  use dspev_module
  use spla
  use iso_c_binding
  IMPLICIT NONE
  include 'laxlib_kinds.fh'
  include 'laxlib_param.fh'
  include 'laxlib_hi.fh'
  include 'laxlib_low.fh'
#if defined(__MPI)
  INTEGER    STATUS(MPI_STATUS_SIZE)
#else
#define MPI_MAX_PROCESSOR_NAME 64
#endif
  INTEGER :: mype, npes, comm, ntgs
  INTEGER :: ierr
  INTEGER :: i
  INTEGER :: nvec
  INTEGER :: ndiag
  INTEGER :: repeats
  INTEGER :: ii
  INTEGER :: kdim, kdmx
  COMPLEX(DP) :: alpha, beta
  type(c_ptr) :: matDis, ctx
  INTEGER, PARAMETER :: blocksize = 256
  INTEGER :: nx ! maximum local block dimension
  LOGICAL :: la_proc ! flag to distinguish procs involved in linear algebra
  INTEGER :: idesc(LAX_DESC_SIZE), idesc_old(LAX_DESC_SIZE)
  INTEGER, ALLOCATABLE :: irc_ip( : )
  INTEGER, ALLOCATABLE :: nrc_ip( : )
  INTEGER, ALLOCATABLE :: rank_ip( :, : )
  !
  COMPLEX(DP), ALLOCATABLE :: A(:,:)
  COMPLEX(DP), ALLOCATABLE :: B(:,:)
  COMPLEX(DP), ALLOCATABLE :: C(:,:)
  integer :: nargs
  CHARACTER(LEN=100) :: arg
  CHARACTER(LEN=300) :: output_file_name
  CHARACTER(LEN=40) :: timer_name

#if defined(__MPI)

#if defined(_OPENMP)
  CALL MPI_Init_thread(MPI_THREAD_FUNNELED, PROVIDED, ierr)
#else
  CALL MPI_Init(ierr)
#endif

  CALL mpi_comm_rank(MPI_COMM_WORLD,mype,ierr)
  CALL mpi_comm_size(MPI_COMM_WORLD,npes,ierr)
  comm = MPI_COMM_WORLD
  ntgs = 1

#else

  mype = 0
  npes = 1
  comm = 0

#endif

  alpha = ONE
  beta = ZERO

  repeats = 20
  nvec = 2000
  kdim = 10
  output_file_name = "timer.json"

  nargs = command_argument_count()
  do i = 1, nargs - 1
    CALL get_command_argument(i, arg)
    IF (TRIM(arg) == '-r') THEN
      CALL get_command_argument(i + 1, arg)
      READ (arg, *) repeats
    END IF
    IF (TRIM(arg) == '-k') THEN
      CALL get_command_argument(i + 1, arg)
      READ (arg, *) kdim
    END IF
    IF (TRIM(arg) == '-n') THEN
      CALL get_command_argument(i + 1, arg)
      READ (arg, *) nvec
    END IF
    IF (TRIM(arg) == '-o') THEN
      CALL get_command_argument(i + 1, arg)
      output_file_name = trim(arg)
    END IF
  end do


  kdmx = kdim

  ALLOCATE(  A( kdmx, nvec ), STAT=ierr )
  ALLOCATE(  B( kdmx, nvec ), STAT=ierr )
  ALLOCATE(  C( nvec, nvec ), STAT=ierr )

  A = ONE
  B = 2*ONE
  C = 3*ONE

  ndiag = 0
  CALL laxlib_start_drv(ndiag, MPI_COMM_WORLD, MPI_COMM_WORLD, .false.)
  CALL desc_init( nvec, nx, la_proc, idesc, rank_ip, irc_ip, nrc_ip )

  ierr = spla_ctx_create(ctx, SPLA_PU_HOST)
  ierr = spla_mat_dis_create_block_cyclic(matDis, comm, 'C', idesc(LAX_DESC_NPR), &
            idesc(LAX_DESC_NPC), idesc(LAX_DESC_NRCX), idesc(LAX_DESC_NRCX))

  if( mype == 0 ) then
    write(6,*) 'n  = ', nvec
    write(6,*) 'k  = ', kdim
    write(6,*) 'repeats  = ', repeats
    write(6,*) '=================================='
    write(6,*) 'num. procs  = ', npes
    write(6,*) 'nx = ', nx
    write(6,*) 'la_proc = ', la_proc
    write(6,*) 'idesc = ', idesc
    write(6,*) 'rank_ip = ', rank_ip
    write(6,*) 'irc_ip = ', irc_ip
    write(6,*) 'nrc_ip = ', nrc_ip
    write(6,*) 'idesc(LAX_DESC_NPC)= ', idesc(LAX_DESC_NPC)
    write(6,*) 'idesc(LAX_DESC_IR)= ', idesc(LAX_DESC_IR)
    write(6,*) 'idesc(LAX_DESC_NR)= ', idesc(LAX_DESC_NR)
    write(6,*) 'idesc(LAX_DESC_IC)= ', idesc(LAX_DESC_IC)
    write(6,*) 'idesc(LAX_DESC_NC)= ', idesc(LAX_DESC_NC)
    write(6,*) 'idesc(LAX_DESC_NX)= ', idesc(LAX_DESC_NC)
    write(6,*) 'idesc(LAX_DESC_N)= ', idesc(LAX_DESC_N)
    write(6,*) 'idesc(LAX_DESC_NRCX)= ', idesc(LAX_DESC_NRCX)
    write(6,*) 'idesc(LAX_DESC_NPR)= ', idesc(LAX_DESC_NPR)
    write(6,*) 'idesc(LAX_DESC_NPC)= ', idesc(LAX_DESC_NPC)
  endif

  ! allocate( proc_name( npes ) )
  ! allocate( node_name( npes ) )
  ! allocate( proc2node( npes ) )


  CALL init_clocks( .true. )


  CALL compute_distmat(C, A, B) ! warm up
  CALL init_clocks( .true. )
  timer_name = "qe"
  ierr = spla_timer_start(len(trim(timer_name), timer_name)
  DO ii = 1, repeats
    timer_name = "iter"
    ierr = spla_timer_start(len(trim(timer_name), timer_name)
    CALL start_clock( 'compute_distmat' )
    CALL compute_distmat(C, A, B)
    CALL stop_clock( 'compute_distmat' )
    timer_name = "iter"
    ierr = spla_timer_stop(len(trim(timer_name), timer_name)
  END DO

  ! warm up
  ierr = spla_pzgemm_ssbtr(nvec, nvec, kdim, SPLA_OP_CONJ_TRANSPOSE, alpha, A, &
                         kdim, B, kdim, beta, C, nvec, 0, 0, SPLA_FILL_MODE_UPPER,&
                         matDis, ctx)
  timer_name = "spla"
  ierr = spla_timer_start(len(trim(timer_name), timer_name)
  DO ii = 1, repeats
    timer_name = "iter"
    ierr = spla_timer_start(len(trim(timer_name), timer_name)
    CALL start_clock( 'spla_pzgemm' )
    ierr = spla_pzgemm_ssbtr(nvec, nvec, kdim, SPLA_OP_CONJ_TRANSPOSE, alpha, A, &
                           kdim, B, kdim, beta, C, nvec, 0, 0, SPLA_FILL_MODE_UPPER,&
                           matDis, ctx)
    CALL stop_clock( 'spla_pzgemm' )
    timer_name = "iter"
    ierr = spla_timer_stop(len(trim(timer_name), timer_name)
  END DO
  timer_name = "spla"
  ierr = spla_timer_stop(len(trim(timer_name), timer_name)


  if( mype == 0 ) then
    ierr = spla_timer_print(len()
    ierr = spla_timer_export_json(len(trim(output_file_name), output_file_name)
    write(6,*)  'compute_distmat' 
    CALL print_clock( 'compute_distmat' )
    write(6,*)  'compute matrix blocks' 
    CALL print_clock( 'compute matrix blocks' )
    write(6,*)  'sym' 
    CALL print_clock( 'sym' )
    write(6,*)  'spla_pzgemm_ssb' 
    CALL print_clock( 'spla_pzgemm' )
  endif

#if defined(__MPI)
  CALL mpi_finalize(ierr)
#endif

 DEALLOCATE( A )
 DEALLOCATE( B )
 DEALLOCATE( C )


CONTAINS
!----------------------------------------------------------------------------
  SUBROUTINE reduce_base_real_to( dim, ps, psout, comm, root )
  !----------------------------------------------------------------------------
  !
  ! ... sums a distributed variable ps(dim) over the processors,
  ! ... and store the results in variable psout.
  ! ... This version uses a fixed-length buffer of appropriate (?) length
  !
  USE util_param, ONLY : DP
  USE parallel_include  
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)  :: dim
  COMPLEX(DP), INTENT(IN)  :: ps(dim)
  COMPLEX(DP)              :: psout(dim)
  INTEGER,  INTENT(IN)  :: comm    ! communecator
  INTEGER,  INTENT(IN)  :: root    ! if root <  0 perform a reduction to all procs
                                   ! if root >= 0 perform a reduce only to root proc.
  !
#if defined (__MPI)  
  !
  INTEGER            :: info, n, nbuf, nproc, myid
  INTEGER, PARAMETER :: maxb = __MSGSIZ_MAX
  !
#if defined __TRACE
  write(*,*) 'reduce_base_real_to IN'
#endif

  CALL mpi_comm_size( comm, nproc, info )
  IF( info /= 0 ) CALL errore( 'reduce_base_real_to', 'error in mpi_comm_size', info )

  CALL mpi_comm_rank( comm, myid, info )
  IF( info /= 0 ) CALL errore( 'reduce_base_real_to', 'error in mpi_comm_rank', info )
  !
  IF ( dim > 0 .AND. nproc <= 1 ) THEN
     psout = ps
  END IF
  IF( dim <= 0 .OR. nproc <= 1 ) GO TO 1 ! go to the end of the subroutine
  !
  ! ... synchronize processes
  !
#if defined __USE_BARRIER
  CALL mp_synchronize( comm )
#endif
  !
  nbuf = dim / maxb
  !
  DO n = 1, nbuf
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps(1+(n-1)*maxb), psout(1+(n-1)*maxb), 2*maxb, MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_real_to', 'error in mpi_reduce 1', info )
     ELSE
        CALL MPI_ALLREDUCE( ps(1+(n-1)*maxb), psout(1+(n-1)*maxb), 2*maxb, MPI_DOUBLE_PRECISION, MPI_SUM, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_real_to', 'error in mpi_allreduce 1', info )
     END IF
     !                    
  END DO
  !
  ! ... possible remaining elements < maxb
  !
  IF ( ( dim - nbuf * maxb ) > 0 ) THEN
     !
     IF( root >= 0 ) THEN
        CALL MPI_REDUCE( ps(1+nbuf*maxb), psout(1+nbuf*maxb), 2*(dim-nbuf*maxb), MPI_DOUBLE_PRECISION, MPI_SUM, root, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_real_to', 'error in mpi_reduce 2', info )
     ELSE
        CALL MPI_ALLREDUCE( ps(1+nbuf*maxb), psout(1+nbuf*maxb), 2*(dim-nbuf*maxb), MPI_DOUBLE_PRECISION, MPI_SUM, comm, info )
        IF( info /= 0 ) CALL errore( 'reduce_base_real_to', 'error in mpi_allreduce 2', info )
     END IF
     !
  END IF
  !
1 CONTINUE
  !
#if defined __TRACE
  write(*,*) 'reduce_base_real_to OUT'
#endif
  !
#endif
  !
  RETURN
  !
END SUBROUTINE reduce_base_real_to

SUBROUTINE mp_root_sum_cm( msg, res, root, gid )
  IMPLICIT NONE
  COMPLEX (DP), INTENT (IN)  :: msg(:,:)
  COMPLEX (DP), INTENT (OUT) :: res(:,:)
  INTEGER,   INTENT (IN)  :: root
  INTEGER,  INTENT (IN) :: gid
#if defined(__MPI)
  INTEGER :: msglen, ierr, taskid

  msglen = size(msg)

  CALL mpi_comm_rank( gid, taskid, ierr)
  IF( ierr /= 0 ) ERROR STOP

  IF( taskid == root ) THEN
     IF( msglen > size(res) ) ERROR STOP
  END IF

  CALL reduce_base_real_to( msglen, msg, res, gid, root )

#else

        res = msg

#endif

END SUBROUTINE mp_root_sum_cm



SUBROUTINE compute_distmat( dm, v, w )
   !
   !  This subroutine compute <vi|wj> and store the
   !  result in distributed matrix dm 
   !
   INTEGER :: ipc, ipr
   INTEGER :: nr, nc, ir, ic, root
   COMPLEX(DP), INTENT(OUT) :: dm( :, : )
   COMPLEX(DP) :: v(:,:), w(:,:)
   COMPLEX(DP), ALLOCATABLE :: work( :, : )
   !
   CALL start_clock( 'compute matrix blocks' )
   ALLOCATE( work( nx, nx ) )
   !
   work = ZERO
   !
   !  Only upper triangle is computed, then the matrix is hermitianized
   !
   DO ipc = 1, idesc(LAX_DESC_NPC) !  loop on column procs 
      !
      nc = nrc_ip( ipc )
      ic = irc_ip( ipc )
      !
      DO ipr = 1, ipc ! idesc(LAX_DESC_NPR) ! ipc ! use symmetry for the loop on row procs
         !
         nr = nrc_ip( ipr )
         ir = irc_ip( ipr )
         !
         !  rank of the processor for which this block (ipr,ipc) is destinated
         !
         root = rank_ip( ipr, ipc )

         ! use blas subs. on the matrix block

         CALL ZGEMM( 'C', 'N', nr, nc, kdim, ONE , &
                     v(1,ir), kdmx, w(1,ic), kdmx, ZERO, work, nx )

         ! accumulate result on dm of root proc.

         CALL mp_root_sum_cm( work, dm, root, comm )

      END DO
      !
   END DO
   CALL stop_clock( 'compute matrix blocks' )
   !
   !  The matrix is hermitianized using upper triangle
   !
   ! CALL start_clock( 'sym' )
   ! CALL laxlib_zsqmher( nvec, dm, nx, idesc )
   ! CALL stop_clock( 'sym' )
   !
   DEALLOCATE( work )
   !
   RETURN
END SUBROUTINE compute_distmat




end program lax_spla


