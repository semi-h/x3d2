module m_case_channel
  use iso_fortran_env, only: stderr => error_unit
  use mpi

  use m_allocator, only: allocator_t, field_t
  use m_base_backend, only: base_backend_t
  use m_base_case, only: base_case_t
  use m_common, only: dp, DIR_C, VERT, CELL
  use m_mesh, only: mesh_t
  use m_solver, only: init

  implicit none

  type, extends(base_case_t) :: case_channel_t
  contains
    procedure :: boundary_conditions => boundary_conditions_channel
    procedure :: initial_conditions => initial_conditions_channel
    procedure :: forcings => forcings_channel
    procedure :: postprocess => postprocess_channel
  end type case_channel_t

  interface case_channel_t
    module procedure case_channel_init
  end interface case_channel_t

contains

  function case_channel_init(backend, mesh, host_allocator) result(flow_case)
    implicit none

    class(base_backend_t), target, intent(inout) :: backend
    type(mesh_t), target, intent(inout) :: mesh
    type(allocator_t), target, intent(inout) :: host_allocator
    type(case_channel_t) :: flow_case

    call flow_case%case_init(backend, mesh, host_allocator)

  end function case_channel_init

  subroutine boundary_conditions_channel(self)
    implicit none

    class(case_channel_t) :: self

    class(field_t), pointer :: u_host
    real(dp) :: can, ub, coeff, dy, L_y
    integer :: dims(3), i, j, k, ierr

    u_host => self%solver%host_allocator%get_block(DIR_C)
    call self%solver%backend%get_field_data(u_host%data, self%solver%u)

    dims = self%solver%mesh%get_dims(VERT)
    ub = 0._dp
    do k = 1, dims(3)
      do j = 1, dims(2)
       do i = 1, dims(1)
          ub = ub + u_host%data(i, j, k)
        end do
      end do
    end do

    call self%solver%host_allocator%release_block(u_host)

    dy = self%solver%mesh%geo%d(2)
    L_y = self%solver%mesh%geo%L(2)
    dims = self%solver%mesh%get_global_dims(CELL)
    coeff = dy/(L_y*dims(1)*dims(3))
    !coeff = 1._dp/(product(self%solver%mesh%get_global_dims(VERT)))

    ub = ub*coeff
    call MPI_Allreduce(MPI_IN_PLACE, ub, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, MPI_COMM_WORLD, ierr)

    can = 2._dp/3._dp - ub
!print*, 'can', can, coeff, dy, L_y, dims
    call self%solver%backend%field_shift(self%solver%u, can)

  end subroutine boundary_conditions_channel

  subroutine initial_conditions_channel(self)
    implicit none

    class(case_channel_t) :: self

    class(field_t), pointer :: u_init, v_init, w_init

    integer :: i, j, k, dims(3), ii, code
    real(dp) :: xloc(3), y, noise, um

    dims = self%solver%mesh%get_dims(VERT)
    u_init => self%solver%host_allocator%get_block(DIR_C)
    v_init => self%solver%host_allocator%get_block(DIR_C)
    w_init => self%solver%host_allocator%get_block(DIR_C)

!    call system_clock(count=code)
!    call random_seed(size=ii)
!    print*, 'seed', ii, code
!    call random_seed(put=code+63946*[(i - 1, i = 1, ii)])

    call random_number(u_init%data(1:dims(1), 1:dims(2), 1:dims(3)))
    call random_number(v_init%data(1:dims(1), 1:dims(2), 1:dims(3)))
    call random_number(w_init%data(1:dims(1), 1:dims(2), 1:dims(3)))

    noise = 0.125_dp
    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          xloc = self%solver%mesh%get_coordinates(i, j, k)
          y = xloc(2) - self%solver%mesh%geo%L(2)/2._dp
          um = exp(-0.2_dp*y*y)

          u_init%data(i, j, k) = 1._dp - y*y &
                                 + noise*um*(2*u_init%data(i, j, k) - 1._dp)
          v_init%data(i, j, k) = noise*um*(2*v_init%data(i, j, k) - 1._dp)
          w_init%data(i, j, k) = noise*um*(2*w_init%data(i, j, k) - 1._dp)
        end do
      end do
    end do

    u_init%data(:, 1, :) = 0
    v_init%data(:, 1, :) = 0
    w_init%data(:, 1, :) = 0
    u_init%data(:, dims(2), :) = 0
    v_init%data(:, dims(2), :) = 0
    w_init%data(:, dims(2), :) = 0

    call self%solver%backend%set_field_data(self%solver%u, u_init%data)
    call self%solver%backend%set_field_data(self%solver%v, v_init%data)
    call self%solver%backend%set_field_data(self%solver%w, w_init%data)

    call self%solver%host_allocator%release_block(u_init)
    call self%solver%host_allocator%release_block(v_init)
    call self%solver%host_allocator%release_block(w_init)

    call self%solver%u%set_data_loc(VERT)
    call self%solver%v%set_data_loc(VERT)
    call self%solver%w%set_data_loc(VERT)

  end subroutine initial_conditions_channel

  subroutine postprocess_channel(self, iter, t)
    implicit none

    class(case_channel_t) :: self
    integer, intent(in) :: iter
    real(dp), intent(in) :: t

    integer :: iounit, dims(3), i, j, k
    character(len=20) :: name, iterchar
    class(field_t), pointer :: u_host

    if (self%solver%mesh%par%is_root()) print *, 'time =', t, 'iteration =', iter
    call self%print_enstrophy(self%solver%u, self%solver%v, self%solver%w)
    call self%print_div_max_mean(self%solver%u, self%solver%v, self%solver%w)

    u_host => self%solver%host_allocator%get_block(DIR_C)
    call self%solver%backend%get_field_data(u_host%data, self%solver%u)
    write(iterchar, '(i0)') iter
    name = 'u'//trim(iterchar)//'.vtr'
    open(newunit=iounit, file=trim(name), status='replace')

    dims = self%solver%mesh%get_dims(VERT)
    write(iounit, '(a)') '# vtk DataFile Version 2.0'
    write(iounit, '(a)') 'title'
    write(iounit, '(a)') 'ASCII'
    write(iounit, '(a)') 'DATASET STRUCTURED_POINTS'
    write(iounit, '("DIMENSIONS ", i0, " ", i0, " ", i0)') dims(1), dims(2), dims(3)
    write(iounit, '(a)') 'ORIGIN 0 0 0'
    write(iounit, '(a)', advance='no') 'SPACING '
    write(iounit, '(3f0.15)') self%solver%mesh%geo%d(1), self%solver%mesh%geo%d(2), self%solver%mesh%geo%d(3)
    write(iounit, '(a)', advance='no') 'POINT_DATA '
    write(iounit, '(i0)') product(dims)
    write(iounit, '(a)') 'SCALARS vol double 1'
    write(iounit, '(a)') 'LOOKUP_TABLE default'
    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          write(iounit, '(f0.15)') u_host%data(i, j, k)
        end do
      end do
    end do
    close(iounit)
    call self%solver%host_allocator%release_block(u_host)

  end subroutine postprocess_channel

  subroutine forcings_channel(self, du, dv, dw, i)
    implicit none

    class(case_channel_t) :: self
    class(field_t), intent(inout) :: du, dv, dw
    integer, intent(in) :: i
    integer :: dims(3)

    if (i < 5000) then
      call self%solver%backend%vecadd(-0.12_dp, self%solver%v, 1._dp, du)
      call self%solver%backend%vecadd(0.12_dp, self%solver%u, 1._dp, dv)
      !print*, 'u v shape', self%solver%u%get_shape(), self%solver%v%get_shape()
      !print*, 'du dv shape', du%get_shape(), dv%get_shape()
    end if

  end subroutine forcings_channel

end module m_case_channel
