module m_solver
  use iso_fortran_env, only: stderr => error_unit
  use mpi

  use m_allocator, only: allocator_t, field_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, get_argument, &
                      RDR_X2Y, RDR_X2Z, RDR_Y2X, RDR_Y2Z, RDR_Z2X, RDR_Z2Y, &
                      RDR_Z2C, RDR_C2Z, &
                      DIR_X, DIR_Y, DIR_Z, DIR_C, VERT, CELL
  use m_config, only: solver_config_t
  use m_tdsops, only: tdsops_t, dirps_t
  use m_time_integrator, only: time_intg_t
  use m_vector_calculus, only: vector_calculus_t
  use m_mesh, only: mesh_t

  implicit none

  type :: solver_t
      !! solver class defines the Incompact3D algorithm at a very high level.
      !!
      !! Procedures defined here that are part of the Incompact3D algorithm
      !! are: transeq, divergence, poisson, and gradient.
      !!
      !! The operations these high level procedures require are provided by
      !! the relavant backend implementations.
      !!
      !! transeq procedure obtains the derivations in x, y, and z directions
      !! using the transeq_x, transeq_y, and transeq_z operations provided by
      !! the backend.
      !! There are two different algorithms available for this operation, a
      !! distributed algorithm and the Thomas algorithm. At the solver class
      !! level it isn't known which algorithm will be executed, that is decided
      !! at run time and therefore backend implementations are responsible for
      !! executing the right subroutines.
      !!
      !! Allocator is responsible from giving us a field sized array when
      !! requested. For example, when the derivations in x direction are
      !! completed and we are ready for the y directional derivatives, we need
      !! three fields to reorder and store the velocities in y direction. Also,
      !! we need three more fields for storing the results, and the get_block
      !! method of the allocator is used to arrange all these memory
      !! assignments. Later, when a field is no more required, release_block
      !! method of the allocator can be used to make this field available
      !! for later use.

    real(dp) :: dt, nu
    integer :: n_iters, n_output
    integer :: ngrid

    class(field_t), pointer :: u, v, w

    class(base_backend_t), pointer :: backend
    type(mesh_t), pointer :: mesh
    type(time_intg_t) :: time_integrator
    type(allocator_t), pointer :: host_allocator
    type(dirps_t), pointer :: xdirps, ydirps, zdirps
    type(vector_calculus_t) :: vector_calculus
    procedure(poisson_solver), pointer :: poisson => null()
  contains
    procedure :: transeq
    procedure :: pressure_correction
    procedure :: divergence_v2p
    procedure :: gradient_p2v
    procedure :: curl
  end type solver_t

  abstract interface
    subroutine poisson_solver(self, pressure, div_u)
      import :: solver_t
      import :: field_t
      implicit none

      class(solver_t) :: self
      class(field_t), intent(inout) :: pressure
      class(field_t), intent(in) :: div_u
    end subroutine poisson_solver
  end interface

  interface solver_t
    module procedure init
  end interface solver_t

contains

  function init(backend, mesh, host_allocator) result(solver)
    implicit none

    class(base_backend_t), target, intent(inout) :: backend
    type(mesh_t), target, intent(inout) :: mesh
    type(allocator_t), target, intent(inout) :: host_allocator
    type(solver_t) :: solver

    type(solver_config_t) :: solver_cfg

    solver%backend => backend
    solver%mesh => mesh
    solver%host_allocator => host_allocator

    allocate (solver%xdirps, solver%ydirps, solver%zdirps)
    solver%xdirps%dir = DIR_X
    solver%ydirps%dir = DIR_Y
    solver%zdirps%dir = DIR_Z

    solver%vector_calculus = vector_calculus_t(solver%backend)

    solver%u => solver%backend%allocator%get_block(DIR_X)
    solver%v => solver%backend%allocator%get_block(DIR_X)
    solver%w => solver%backend%allocator%get_block(DIR_X)

    call solver_cfg%read(nml_file=get_argument(1))

    solver%time_integrator = time_intg_t(solver%backend, &
                                         solver%backend%allocator, &
                                         solver_cfg%time_intg)
    if (solver%mesh%par%is_root()) then
      print *, solver_cfg%time_intg//' time integrator instantiated'
    end if

    solver%dt = solver_cfg%dt
    solver%backend%nu = 1._dp/solver_cfg%Re
    solver%n_iters = solver_cfg%n_iters
    solver%n_output = solver_cfg%n_output
    solver%ngrid = product(solver%mesh%get_global_dims(VERT))

    ! Allocate and set the tdsops
    call allocate_tdsops( &
      solver%xdirps, solver%backend, solver%mesh, solver_cfg%der1st_scheme, &
      solver_cfg%der2nd_scheme, solver_cfg%interpl_scheme, &
      solver_cfg%stagder_scheme &
      )
    call allocate_tdsops( &
      solver%ydirps, solver%backend, solver%mesh, solver_cfg%der1st_scheme, &
      solver_cfg%der2nd_scheme, solver_cfg%interpl_scheme, &
      solver_cfg%stagder_scheme &
      )
    call allocate_tdsops( &
      solver%zdirps, solver%backend, solver%mesh, solver_cfg%der1st_scheme, &
      solver_cfg%der2nd_scheme, solver_cfg%interpl_scheme, &
      solver_cfg%stagder_scheme &
      )

    select case (trim(solver_cfg%poisson_solver_type))
    case ('FFT')
      if (solver%mesh%par%is_root()) print *, 'Poisson solver: FFT'
      call solver%backend%init_poisson_fft(solver%mesh, solver%xdirps, &
                                           solver%ydirps, solver%zdirps)
      solver%poisson => poisson_fft
    case ('CG')
      if (solver%mesh%par%is_root()) &
        print *, 'Poisson solver: CG, not yet implemented'
      solver%poisson => poisson_cg
    case default
      error stop 'poisson_solver_type is not valid. Use "FFT" or "CG".'
    end select

  end function init

  subroutine allocate_tdsops(dirps, backend, mesh, der1st_scheme, &
                             der2nd_scheme, interpl_scheme, stagder_scheme)
    type(dirps_t), intent(inout) :: dirps
    class(base_backend_t), intent(in) :: backend
    type(mesh_t), intent(in) :: mesh
    character(*), intent(in) :: der1st_scheme, der2nd_scheme, &
                                interpl_scheme, stagder_scheme

    integer :: dir, bc_start, bc_end, n_vert, n_cell
    real(dp) :: d

    dir = dirps%dir
    bc_start = mesh%grid%BCs(dir, 1)
    bc_end = mesh%grid%BCs(dir, 2)
    d = mesh%geo%d(dir)

    n_vert = mesh%get_n(dir, VERT)
    n_cell = mesh%get_n(dir, CELL)

    call backend%alloc_tdsops(dirps%der1st, n_vert, d, 'first-deriv', &
                              der1st_scheme, bc_start, bc_end)
    call backend%alloc_tdsops(dirps%der1st_sym, n_vert, d, 'first-deriv', &
                              der1st_scheme, bc_start, bc_end)
    call backend%alloc_tdsops(dirps%der2nd, n_vert, d, 'second-deriv', &
                              der2nd_scheme, bc_start, bc_end)
    call backend%alloc_tdsops(dirps%der2nd_sym, n_vert, d, 'second-deriv', &
                              der2nd_scheme, bc_start, bc_end)
    call backend%alloc_tdsops(dirps%interpl_v2p, n_cell, d, 'interpolate', &
                              interpl_scheme, bc_start, bc_end, from_to='v2p')
    call backend%alloc_tdsops(dirps%interpl_p2v, n_vert, d, 'interpolate', &
                              interpl_scheme, bc_start, bc_end, from_to='p2v')
    call backend%alloc_tdsops(dirps%stagder_v2p, n_cell, d, 'stag-deriv', &
                              stagder_scheme, bc_start, bc_end, from_to='v2p')
    call backend%alloc_tdsops(dirps%stagder_p2v, n_vert, d, 'stag-deriv', &
                              stagder_scheme, bc_start, bc_end, from_to='p2v')

    if (dir == 2) &
    print*, 'dirps%stagder_v2p', dirps%stagder_v2p%n_tds, dirps%stagder_v2p%n_rhs
  end subroutine

  subroutine transeq(self, du, dv, dw, u, v, w)
    !! Skew-symmetric form of convection-diffusion terms in the
    !! incompressible Navier-Stokes momemtum equations, excluding
    !! pressure terms.
    !! Inputs from velocity grid and outputs to velocity grid.
    implicit none

    class(solver_t) :: self
    class(field_t), intent(inout) :: du, dv, dw
    class(field_t), intent(in) :: u, v, w

    class(field_t), pointer :: u_y, v_y, w_y, u_z, v_z, w_z, &
      du_y, dv_y, dw_y, du_z, dv_z, dw_z

    ! -1/2(nabla u curl u + u nabla u) + nu nablasq u

    ! call derivatives in x direction. Based on the run time arguments this
    ! executes a distributed algorithm or the Thomas algorithm.
    call self%backend%transeq_x(du, dv, dw, u, v, w, self%xdirps)

    ! request fields from the allocator
    u_y => self%backend%allocator%get_block(DIR_Y)
    v_y => self%backend%allocator%get_block(DIR_Y)
    w_y => self%backend%allocator%get_block(DIR_Y)
    du_y => self%backend%allocator%get_block(DIR_Y)
    dv_y => self%backend%allocator%get_block(DIR_Y)
    dw_y => self%backend%allocator%get_block(DIR_Y)

    ! reorder data from x orientation to y orientation
    call self%backend%reorder(u_y, u, RDR_X2Y)
    call self%backend%reorder(v_y, v, RDR_X2Y)
    call self%backend%reorder(w_y, w, RDR_X2Y)

    ! similar to the x direction, obtain derivatives in y.
    call self%backend%transeq_y(du_y, dv_y, dw_y, u_y, v_y, w_y, self%ydirps)

    ! we don't need the velocities in y orientation any more, so release
    ! them to open up space.
    ! It is important that this doesn't actually deallocate any memory,
    ! it just makes the corresponding memory space available for use.
    call self%backend%allocator%release_block(u_y)
    call self%backend%allocator%release_block(v_y)
    call self%backend%allocator%release_block(w_y)

    call self%backend%sum_yintox(du, du_y)
    call self%backend%sum_yintox(dv, dv_y)
    call self%backend%sum_yintox(dw, dw_y)

    call self%backend%allocator%release_block(du_y)
    call self%backend%allocator%release_block(dv_y)
    call self%backend%allocator%release_block(dw_y)

    ! just like in y direction, get some fields for the z derivatives.
    u_z => self%backend%allocator%get_block(DIR_Z)
    v_z => self%backend%allocator%get_block(DIR_Z)
    w_z => self%backend%allocator%get_block(DIR_Z)
    du_z => self%backend%allocator%get_block(DIR_Z)
    dv_z => self%backend%allocator%get_block(DIR_Z)
    dw_z => self%backend%allocator%get_block(DIR_Z)

    ! reorder from x to z
    call self%backend%reorder(u_z, u, RDR_X2Z)
    call self%backend%reorder(v_z, v, RDR_X2Z)
    call self%backend%reorder(w_z, w, RDR_X2Z)

    ! get the derivatives in z
    call self%backend%transeq_z(du_z, dv_z, dw_z, u_z, v_z, w_z, self%zdirps)

    ! there is no need to keep velocities in z orientation around, so release
    call self%backend%allocator%release_block(u_z)
    call self%backend%allocator%release_block(v_z)
    call self%backend%allocator%release_block(w_z)

    ! gather all the contributions into the x result array
    call self%backend%sum_zintox(du, du_z)
    call self%backend%sum_zintox(dv, dv_z)
    call self%backend%sum_zintox(dw, dw_z)

    ! release all the unnecessary blocks.
    call self%backend%allocator%release_block(du_z)
    call self%backend%allocator%release_block(dv_z)
    call self%backend%allocator%release_block(dw_z)

  end subroutine transeq

  subroutine divergence_v2p(self, div_u, u, v, w)
    !! Wrapper for divergence_v2p
    implicit none

    class(solver_t) :: self
    class(field_t), intent(inout) :: div_u
    class(field_t), intent(in) :: u, v, w


    integer :: dims(3), ierr
    real(dp) :: u_max, u_min, u_mean
    class(field_t), pointer :: u_out
    class(field_t), pointer :: du_x, dv_x, dw_x, &
      u_y, v_y, w_y, du_y, dv_y, dw_y, &
      u_z, w_z, dw_z


    call self%vector_calculus%divergence_v2c( &
      div_u, u, v, w, &
      self%xdirps%stagder_v2p, self%xdirps%interpl_v2p, &
      self%ydirps%stagder_v2p, self%ydirps%interpl_v2p, &
      self%zdirps%stagder_v2p, self%zdirps%interpl_v2p &
      )

if (.false.) then
    if (div_u%dir /= DIR_Z .or. u%dir /= DIR_X .or. v%dir /= DIR_X &
        .or. w%dir /= DIR_X) then
      error stop 'Error in divergence_v2c input/output field dirs: &
                  &output must be in DIR_Z, inputs must be in DIR_X layout.'
    end if

    du_x => self%backend%allocator%get_block(DIR_X)
    dv_x => self%backend%allocator%get_block(DIR_X)
    dw_x => self%backend%allocator%get_block(DIR_X)

    ! Staggared der for u field in x
    ! Interpolation for v field in x
    ! Interpolation for w field in x
    call self%backend%tds_solve(du_x, u, self%xdirps%stagder_v2p)
    call self%backend%tds_solve(dv_x, v, self%xdirps%interpl_v2p)
    call self%backend%tds_solve(dw_x, w, self%xdirps%interpl_v2p)


    u_out => self%host_allocator%get_block(DIR_C)
    call self%backend%get_field_data(u_out%data, v)

    dims = self%mesh%get_dims(v%data_loc)
    print*, 'v data loc', v%data_loc, dims
    u_max = maxval(abs(u_out%data(1:dims(1), 1:dims(2), 1:dims(3))))
    u_min = minval(u_out%data(1:dims(1), 1:dims(2), 1:dims(3)))
    u_mean = sum(abs(u_out%data(1:dims(1), 1:dims(2), 1:dims(3)))) &
                 /self%ngrid

    call self%host_allocator%release_block(u_out)

    call MPI_Allreduce(MPI_IN_PLACE, u_max, 1, MPI_DOUBLE_PRECISION, &
                       MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, u_mean, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, MPI_COMM_WORLD, ierr)
    if (self%mesh%par%is_root()) &
      print *, 'v max min mean:', u_max, u_min, u_mean


    u_out => self%host_allocator%get_block(DIR_C)
    call self%backend%get_field_data(u_out%data, du_x)

    dims = self%mesh%get_dims(du_x%data_loc)
    print*, 'du_x data loc', du_x%data_loc, dims
    u_max = maxval(abs(u_out%data(1:dims(1), 1:dims(2), 1:dims(3))))
    u_min = minval(u_out%data(1:dims(1), 1:dims(2), 1:dims(3)))
    u_mean = sum(abs(u_out%data(1:dims(1), 1:dims(2), 1:dims(3)))) &
                 /self%ngrid

    call self%host_allocator%release_block(u_out)

    call MPI_Allreduce(MPI_IN_PLACE, u_max, 1, MPI_DOUBLE_PRECISION, &
                       MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, u_mean, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, MPI_COMM_WORLD, ierr)
    if (self%mesh%par%is_root()) &
      print *, 'du_x max min mean:', u_max, u_min, u_mean


    u_out => self%host_allocator%get_block(DIR_C)
    call self%backend%get_field_data(u_out%data, dv_x)

    dims = self%mesh%get_dims(dv_x%data_loc)
    print*, 'dv_x data loc', dv_x%data_loc, dims
    u_max = maxval(abs(u_out%data(1:dims(1), 1:dims(2), 1:dims(3))))
    u_min = minval(u_out%data(1:dims(1), 1:dims(2), 1:dims(3)))
    u_mean = sum(abs(u_out%data(1:dims(1), 1:dims(2), 1:dims(3)))) &
                 /self%ngrid

    call self%host_allocator%release_block(u_out)

    call MPI_Allreduce(MPI_IN_PLACE, u_max, 1, MPI_DOUBLE_PRECISION, &
                       MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, u_mean, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, MPI_COMM_WORLD, ierr)
    if (self%mesh%par%is_root()) &
      print *, 'dv_x max min mean:', u_max, u_min, u_mean



    ! request fields from the allocator
    u_y => self%backend%allocator%get_block(DIR_Y)
    v_y => self%backend%allocator%get_block(DIR_Y)
    w_y => self%backend%allocator%get_block(DIR_Y)

    ! reorder data from x orientation to y orientation
    call self%backend%reorder(u_y, du_x, RDR_X2Y)
    call self%backend%reorder(v_y, dv_x, RDR_X2Y)
    call self%backend%reorder(w_y, dw_x, RDR_X2Y)

    call self%backend%allocator%release_block(du_x)
    call self%backend%allocator%release_block(dv_x)
    call self%backend%allocator%release_block(dw_x)

    du_y => self%backend%allocator%get_block(DIR_Y)
    dv_y => self%backend%allocator%get_block(DIR_Y)
    dw_y => self%backend%allocator%get_block(DIR_Y)




    u_out => self%host_allocator%get_block(DIR_C)
    call self%backend%get_field_data(u_out%data, v_y)

    dims = self%mesh%get_dims(v_y%data_loc)
    print*, 'v_y data loc', v_y%data_loc, dims
    u_max = maxval(abs(u_out%data(1:dims(1), 1:dims(2), 1:dims(3))))
    u_min = minval(u_out%data(1:dims(1), 1:dims(2), 1:dims(3)))
    u_mean = sum(abs(u_out%data(1:dims(1), 1:dims(2), 1:dims(3)))) &
                 /self%ngrid
!print*, 'v_y(1, :, 1)', u_out%data(1, 1:dims(2), 1)
    call self%host_allocator%release_block(u_out)

    call MPI_Allreduce(MPI_IN_PLACE, u_max, 1, MPI_DOUBLE_PRECISION, &
                       MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, u_mean, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, MPI_COMM_WORLD, ierr)
    if (self%mesh%par%is_root()) &
      print *, 'v_y before max min mean:', u_max, u_min, u_mean




    ! similar to the x direction, obtain derivatives in y.
    call self%backend%tds_solve(du_y, u_y, self%ydirps%interpl_v2p)
    call self%backend%tds_solve(dv_y, v_y, self%ydirps%stagder_v2p)
    call self%backend%tds_solve(dw_y, w_y, self%ydirps%interpl_v2p)





    u_out => self%host_allocator%get_block(DIR_C)
    call self%backend%get_field_data(u_out%data, dv_y)

    dims = self%mesh%get_dims(dv_y%data_loc)
    print*, 'dv_y data loc', dv_y%data_loc, dims
    u_max = maxval(abs(u_out%data(1:dims(1), 1:dims(2), 1:dims(3))))
    u_min = minval(u_out%data(1:dims(1), 1:dims(2), 1:dims(3)))
    u_mean = sum(abs(u_out%data(1:dims(1), 1:dims(2), 1:dims(3)))) &
                 /self%ngrid
!print*, 'dv_y(1, :, 1)', u_out%data(1, 1:dims(2), 1)
    call self%host_allocator%release_block(u_out)

    call MPI_Allreduce(MPI_IN_PLACE, u_max, 1, MPI_DOUBLE_PRECISION, &
                       MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, u_mean, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, MPI_COMM_WORLD, ierr)
    if (self%mesh%par%is_root()) &
      print *, 'dv_y before max min mean:', u_max, u_min, u_mean




    ! we don't need the velocities in y orientation any more, so release
    ! them to open up space.
    ! It is important that this doesn't actually deallocate any memory,
    ! it just makes the corresponding memory space available for use.
    call self%backend%allocator%release_block(u_y)
    call self%backend%allocator%release_block(v_y)
    call self%backend%allocator%release_block(w_y)

    ! just like in y direction, get some fields for the z derivatives.
    u_z => self%backend%allocator%get_block(DIR_Z)
    w_z => self%backend%allocator%get_block(DIR_Z)

    ! du_y = dv_y + du_y
    call self%backend%vecadd(1._dp, dv_y, 1._dp, du_y)

    ! reorder from y to z
    call self%backend%reorder(u_z, du_y, RDR_Y2Z)
    call self%backend%reorder(w_z, dw_y, RDR_Y2Z)

    ! release all the unnecessary blocks.
    call self%backend%allocator%release_block(du_y)
    call self%backend%allocator%release_block(dv_y)
    call self%backend%allocator%release_block(dw_y)

    dw_z => self%backend%allocator%get_block(DIR_Z)

    ! get the derivatives in z
    call self%backend%tds_solve(div_u, u_z, self%zdirps%interpl_v2p)
    call self%backend%tds_solve(dw_z, w_z, self%zdirps%stagder_v2p)

    ! div_u = div_u + dw_z
    call self%backend%vecadd(1._dp, dw_z, 1._dp, div_u)

    ! div_u array is in z orientation

    ! there is no need to keep velocities in z orientation around, so release
    call self%backend%allocator%release_block(u_z)
    call self%backend%allocator%release_block(w_z)
    call self%backend%allocator%release_block(dw_z)
end if
  end subroutine divergence_v2p

  subroutine gradient_p2v(self, dpdx, dpdy, dpdz, pressure)
    !! Wrapper for gradient_p2v
    implicit none

    class(solver_t) :: self
    class(field_t), intent(inout) :: dpdx, dpdy, dpdz
    class(field_t), intent(in) :: pressure

    call self%vector_calculus%gradient_c2v( &
      dpdx, dpdy, dpdz, pressure, &
      self%xdirps%stagder_p2v, self%xdirps%interpl_p2v, &
      self%ydirps%stagder_p2v, self%ydirps%interpl_p2v, &
      self%zdirps%stagder_p2v, self%zdirps%interpl_p2v &
      )

  end subroutine gradient_p2v

  subroutine curl(self, o_i_hat, o_j_hat, o_k_hat, u, v, w)
    !! Wrapper for curl
    implicit none

    class(solver_t) :: self
    !> Vector components of the output vector field Omega
    class(field_t), intent(inout) :: o_i_hat, o_j_hat, o_k_hat
    class(field_t), intent(in) :: u, v, w

    call self%vector_calculus%curl( &
      o_i_hat, o_j_hat, o_k_hat, u, v, w, &
      self%xdirps%der1st, self%ydirps%der1st, self%zdirps%der1st &
      )

  end subroutine curl

  subroutine poisson_fft(self, pressure, div_u)
    implicit none

    class(solver_t) :: self
    class(field_t), intent(inout) :: pressure
    class(field_t), intent(in) :: div_u

    class(field_t), pointer :: p_temp, temp

    ! reorder into 3D Cartesian data structure
    p_temp => self%backend%allocator%get_block(DIR_C)
    call self%backend%reorder(p_temp, div_u, RDR_Z2C)

    temp => self%backend%allocator%get_block(DIR_C)
    ! solve poisson equation with FFT based approach
    call self%backend%poisson_fft%solve_poisson(p_temp, temp)
    call self%backend%allocator%release_block(temp)

    ! reorder back to our specialist data structure from 3D Cartesian
    call self%backend%reorder(pressure, p_temp, RDR_C2Z)

    call self%backend%allocator%release_block(p_temp)

  end subroutine poisson_fft

  subroutine poisson_cg(self, pressure, div_u)
    implicit none

    class(solver_t) :: self
    class(field_t), intent(inout) :: pressure
    class(field_t), intent(in) :: div_u

    call pressure%fill(0._dp)
  end subroutine poisson_cg

  subroutine pressure_correction(self, u, v, w)
    implicit none

    class(solver_t) :: self
    class(field_t), intent(inout) :: u, v, w

    class(field_t), pointer :: div_u, pressure, dpdx, dpdy, dpdz
    class(field_t), pointer :: u_out
    real(dp) :: div_u_max, div_u_mean, u_min
    integer :: ierr, dims(3), i, j, k
    class(field_t), pointer :: u_host

    u_host => self%host_allocator%get_block(DIR_C)
    call self%backend%get_field_data(u_host%data, u)
    dims = self%mesh%get_dims(VERT)
    !print*, 'dims of divu', dims
    !print*, 'data loc', u%data_loc
    u_host%data(:, 1, :) = 0
    u_host%data(:, dims(2), :) = 0
    call self%backend%set_field_data(u, u_host%data)
    call u%set_data_loc(VERT)

    call self%backend%get_field_data(u_host%data, v)
    dims = self%mesh%get_dims(VERT)
    !print*, 'dims of divu', dims
    u_host%data(:, 1, :) = 0
    u_host%data(:, dims(2), :) = 0
    call self%backend%set_field_data(v, u_host%data)
    call v%set_data_loc(VERT)

    call self%backend%get_field_data(u_host%data, w)
    dims = self%mesh%get_dims(VERT)
    !print*, 'dims of divu', dims
    u_host%data(:, 1, :) = 0
    u_host%data(:, dims(2), :) = 0
    call self%backend%set_field_data(w, u_host%data)
    call w%set_data_loc(VERT)
!    print*, 'BC set' 



    call self%host_allocator%release_block(u_host)




    div_u => self%backend%allocator%get_block(DIR_Z)

    call self%divergence_v2p(div_u, u, v, w)




!    u_out => self%host_allocator%get_block(DIR_C)
!    call self%backend%get_field_data(u_out%data, div_u)

!    dims = self%mesh%get_dims(div_u%data_loc)
!    div_u_max = maxval(abs(u_out%data(1:dims(1), 1:dims(2), 1:dims(3))))
!    u_min = minval(u_out%data(1:dims(1), 1:dims(2), 1:dims(3)))
!    div_u_mean = sum(abs(u_out%data(1:dims(1), 1:dims(2), 1:dims(3)))) &
!                 /self%ngrid

!print*, 'div_u', u_out%data(1, 1:dims(2), 1)
!    call self%host_allocator%release_block(u_out)
!
!    call MPI_Allreduce(MPI_IN_PLACE, div_u_max, 1, MPI_DOUBLE_PRECISION, &
!                       MPI_MAX, MPI_COMM_WORLD, ierr)
!    call MPI_Allreduce(MPI_IN_PLACE, div_u_mean, 1, MPI_DOUBLE_PRECISION, &
!                       MPI_SUM, MPI_COMM_WORLD, ierr)
!    if (self%mesh%par%is_root()) &
!      print *, 'divu max min mean:', div_u_max, u_min, div_u_mean




    pressure => self%backend%allocator%get_block(DIR_Z)

    call self%poisson(pressure, div_u)






!    u_out => self%host_allocator%get_block(DIR_C)
!    call self%backend%get_field_data(u_out%data, pressure)

!    dims = self%mesh%get_dims(pressure%data_loc)
!    div_u_max = maxval(abs(u_out%data(1:dims(1), 1:dims(2), 1:dims(3))))
!    u_min = minval(u_out%data(1:dims(1), 1:dims(2), 1:dims(3)))
!    div_u_mean = sum(abs(u_out%data(1:dims(1), 1:dims(2), 1:dims(3)))) &
!                 /self%ngrid

!print*, 'pressure', u_out%data(1, 1:dims(2), 1)
!    call self%host_allocator%release_block(u_out)

!    call MPI_Allreduce(MPI_IN_PLACE, div_u_max, 1, MPI_DOUBLE_PRECISION, &
!                       MPI_MAX, MPI_COMM_WORLD, ierr)
!    call MPI_Allreduce(MPI_IN_PLACE, div_u_mean, 1, MPI_DOUBLE_PRECISION, &
!                       MPI_SUM, MPI_COMM_WORLD, ierr)
!    if (self%mesh%par%is_root()) &
!      print *, 'pressure max min mean:', div_u_max, u_min, div_u_mean
!





    call self%backend%allocator%release_block(div_u)

    dpdx => self%backend%allocator%get_block(DIR_X)
    dpdy => self%backend%allocator%get_block(DIR_X)
    dpdz => self%backend%allocator%get_block(DIR_X)

    call self%gradient_p2v(dpdx, dpdy, dpdz, pressure)

    call self%backend%allocator%release_block(pressure)


!    u_out => self%host_allocator%get_block(DIR_C)
!    call self%backend%get_field_data(u_out%data, dpdx)
!
!    dims = self%mesh%get_dims(dpdx%data_loc)
!    print*, 'dpdx data loc', dpdx%data_loc, dims
!    div_u_max = maxval(abs(u_out%data(1:dims(1), 1:dims(2), 1:dims(3))))
!    u_min = minval(u_out%data(1:dims(1), 1:dims(2), 1:dims(3)))
!    div_u_mean = sum(abs(u_out%data(1:dims(1), 1:dims(2), 1:dims(3)))) &
!                 /self%ngrid

!    call self%host_allocator%release_block(u_out)

!    call MPI_Allreduce(MPI_IN_PLACE, div_u_max, 1, MPI_DOUBLE_PRECISION, &
!                       MPI_MAX, MPI_COMM_WORLD, ierr)
!    call MPI_Allreduce(MPI_IN_PLACE, div_u_mean, 1, MPI_DOUBLE_PRECISION, &
!                       MPI_SUM, MPI_COMM_WORLD, ierr)
!    if (self%mesh%par%is_root()) &
!      print *, 'dpdx max min mean:', div_u_max, u_min, div_u_mean




!    u_out => self%host_allocator%get_block(DIR_C)
!    call self%backend%get_field_data(u_out%data, dpdy)
!
!    dims = self%mesh%get_dims(dpdy%data_loc)
!    print*, 'dpdy data loc', dpdy%data_loc, dims
!    div_u_max = maxval(abs(u_out%data(1:dims(1), 1:dims(2), 1:dims(3))))
!    u_min = minval(u_out%data(1:dims(1), 1:dims(2), 1:dims(3)))
!    div_u_mean = sum(abs(u_out%data(1:dims(1), 1:dims(2), 1:dims(3)))) &
!                 /self%ngrid

!print*, 'dpdy', u_out%data(1, 1:dims(2), 1)
!    call self%host_allocator%release_block(u_out)
!
!    call MPI_Allreduce(MPI_IN_PLACE, div_u_max, 1, MPI_DOUBLE_PRECISION, &
!                       MPI_MAX, MPI_COMM_WORLD, ierr)
!    call MPI_Allreduce(MPI_IN_PLACE, div_u_mean, 1, MPI_DOUBLE_PRECISION, &
!                       MPI_SUM, MPI_COMM_WORLD, ierr)
!    if (self%mesh%par%is_root()) &
!      print *, 'dpdy max min mean:', div_u_max, u_min, div_u_mean


!    u_out => self%host_allocator%get_block(DIR_C)
!    call self%backend%get_field_data(u_out%data, dpdz)
!
!    dims = self%mesh%get_dims(dpdz%data_loc)
!    print*, 'dpdz data loc', dpdz%data_loc, dims
!    div_u_max = maxval(abs(u_out%data(1:dims(1), 1:dims(2), 1:dims(3))))
!    u_min = minval(u_out%data(1:dims(1), 1:dims(2), 1:dims(3)))
!    div_u_mean = sum(abs(u_out%data(1:dims(1), 1:dims(2), 1:dims(3)))) &
!                 /self%ngrid
!
!    call self%host_allocator%release_block(u_out)
!
!    call MPI_Allreduce(MPI_IN_PLACE, div_u_max, 1, MPI_DOUBLE_PRECISION, &
!                       MPI_MAX, MPI_COMM_WORLD, ierr)
!    call MPI_Allreduce(MPI_IN_PLACE, div_u_mean, 1, MPI_DOUBLE_PRECISION, &
!                       MPI_SUM, MPI_COMM_WORLD, ierr)
!    if (self%mesh%par%is_root()) &
!      print *, 'dpdz max min mean:', div_u_max, u_min, div_u_mean




    ! velocity correction
    call self%backend%vecadd(-1._dp, dpdx, 1._dp, u)
    call self%backend%vecadd(-1._dp, dpdy, 1._dp, v)
    call self%backend%vecadd(-1._dp, dpdz, 1._dp, w)

    call self%backend%allocator%release_block(dpdx)
    call self%backend%allocator%release_block(dpdy)
    call self%backend%allocator%release_block(dpdz)

  end subroutine pressure_correction

end module m_solver
