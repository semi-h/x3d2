module m_cuda_poisson_fft
  use iso_c_binding, only: c_loc, c_ptr, c_f_pointer
  use iso_fortran_env, only: stderr => error_unit
  use cudafor
  use cufftXt
  use cufft
  use mpi

  use m_allocator, only: field_t
  use m_common, only: dp, DIR_C, CELL
  use m_mesh, only: mesh_t
  use m_poisson_fft, only: poisson_fft_t
  use m_tdsops, only: dirps_t

  use m_cuda_allocator, only: cuda_field_t
  use m_cuda_spectral, only: memcpy3D, &
                             process_spectral_000, process_spectral_010, &
                             enforce_periodicity_y, undo_periodicity_y

  implicit none

  type, extends(poisson_fft_t) :: cuda_poisson_fft_t
    !! FFT based Poisson solver

    !> Local domain sized array storing the spectral equivalence constants
    complex(dp), device, allocatable, dimension(:, :, :) :: waves_dev
    !> Wave numbers in x, y, and z
    real(dp), device, allocatable, dimension(:) :: ax_dev, bx_dev, &
                                                   ay_dev, by_dev, &
                                                   az_dev, bz_dev
    !> Forward and backward FFT transform plans
    integer :: plan3D_fw, plan3D_bw

    !> cuFFTMp object manages decomposition and data storage
    type(cudaLibXtDesc), pointer :: xtdesc
    complex(dp), device, allocatable, dimension(:, :, :) :: temp_dev

    real(dp), allocatable, dimension(:, :, :) :: host_data, host_data_2
    complex(dp), allocatable, dimension(:, :, :) :: host_cdata, host_cdata_2
    real(dp), device, allocatable, dimension(:, :, :) :: device_data_in, device_data_out
  contains
    procedure :: fft_forward => fft_forward_cuda
    procedure :: fft_backward => fft_backward_cuda
    procedure :: fft_postprocess_000 => fft_postprocess_000_cuda
    procedure :: fft_postprocess_010 => fft_postprocess_010_cuda
    procedure :: enforce_periodicity_y => enforce_periodicity_y_cuda
    procedure :: undo_periodicity_y => undo_periodicity_y_cuda
  end type cuda_poisson_fft_t

  interface cuda_poisson_fft_t
    module procedure init
  end interface cuda_poisson_fft_t

  private :: init

contains

  function init(mesh, xdirps, ydirps, zdirps) result(poisson_fft)
    implicit none

    type(mesh_t), intent(in) :: mesh
    type(dirps_t), intent(in) :: xdirps, ydirps, zdirps

    type(cuda_poisson_fft_t) :: poisson_fft

    integer :: nx, ny, nz

    integer :: ierr, i
    type(dim3) :: blocks, threads
    integer(int_ptr_kind()) :: worksize

    integer :: dims_glob(3), dims_loc(3), n_spec(3), n_sp_st(3)

    ! 1D decomposition along Z in real domain, and along Y in spectral space
    if (mesh%par%nproc_dir(2) /= 1) print *, 'nproc_dir in y-dir must be 1'

    ! Work out the spectral dimensions in the permuted state
    dims_glob = mesh%get_global_dims(CELL)
    dims_loc = mesh%get_dims(CELL)

    n_spec(1) = dims_loc(1)/2 + 1
    n_spec(2) = dims_loc(2)/mesh%par%nproc_dir(3)
    n_spec(3) = dims_glob(3)

    n_sp_st(1) = 0
    n_sp_st(2) = dims_loc(2)/mesh%par%nproc_dir(3)*mesh%par%nrank_dir(3)
    n_sp_st(3) = 0

    call poisson_fft%base_init(mesh, xdirps, ydirps, zdirps, n_spec, n_sp_st)

    dims_loc = mesh%get_padded_dims(DIR_C)
    allocate(poisson_fft%host_data(dims_loc(1), dims_loc(2), dims_loc(3)))
    allocate(poisson_fft%host_data_2(dims_loc(1), dims_loc(2), dims_loc(3)))
    allocate(poisson_fft%device_data_in(dims_loc(1), dims_loc(2), dims_loc(3)))
    allocate(poisson_fft%device_data_out(dims_loc(1), dims_loc(2), dims_loc(3)))
    print*, 'host data allocated', dims_loc

    if (.false.) then
    do i = 1, poisson_fft%ny_loc
      poisson_fft%host_data(:, i, :) = i
    end do
    poisson_fft%device_data_in = poisson_fft%host_data

    blocks = dim3(poisson_fft%nz_spec, 1, 1)
    threads = dim3(poisson_fft%nx_spec, 1, 1)
    call enforce_periodicity_y<<<blocks, threads>>>( & !&
      poisson_fft%device_data_out, poisson_fft%device_data_in, poisson_fft%ny_spec &
      )

    poisson_fft%host_data = poisson_fft%device_data_out
    poisson_fft%device_data_in = poisson_fft%host_data
    print*, 'host data after enforcing', poisson_fft%host_data(1, :, 1)

    blocks = dim3(poisson_fft%nz_spec, 1, 1)
    threads = dim3(poisson_fft%nx_spec, 1, 1)
    call undo_periodicity_y<<<blocks, threads>>>( & !&
      poisson_fft%device_data_out, poisson_fft%device_data_in, poisson_fft%ny_spec &
      )

    poisson_fft%host_data = poisson_fft%device_data_out

    print*, 'host data after undo', poisson_fft%host_data(1, :, 1)
    stop
    end if


    nx = poisson_fft%nx_glob
    ny = poisson_fft%ny_glob
    nz = poisson_fft%nz_glob

    allocate (poisson_fft%waves_dev(poisson_fft%nx_spec, &
                                    poisson_fft%ny_spec, &
                                    poisson_fft%nz_spec))
    poisson_fft%waves_dev = poisson_fft%waves
    allocate (poisson_fft%temp_dev(poisson_fft%nx_spec, &
                                    poisson_fft%ny_spec, &
                                    poisson_fft%nz_spec))
    allocate (poisson_fft%host_cdata(poisson_fft%nx_spec, &
                                    poisson_fft%ny_spec, &
                                    poisson_fft%nz_spec))
    allocate (poisson_fft%host_cdata_2(poisson_fft%nx_spec, &
                                    poisson_fft%ny_spec, &
                                    poisson_fft%nz_spec))

    allocate (poisson_fft%ax_dev(nx), poisson_fft%bx_dev(nx))
    allocate (poisson_fft%ay_dev(ny), poisson_fft%by_dev(ny))
    allocate (poisson_fft%az_dev(nz), poisson_fft%bz_dev(nz))
    poisson_fft%ax_dev = poisson_fft%ax; poisson_fft%bx_dev = poisson_fft%bx
    poisson_fft%ay_dev = poisson_fft%ay; poisson_fft%by_dev = poisson_fft%by
    poisson_fft%az_dev = poisson_fft%az; poisson_fft%bz_dev = poisson_fft%bz

    ! 3D plans
    ierr = cufftCreate(poisson_fft%plan3D_fw)
    ierr = cufftMpAttachComm(poisson_fft%plan3D_fw, CUFFT_COMM_MPI, &
                             MPI_COMM_WORLD)
    ierr = cufftMakePlan3D(poisson_fft%plan3D_fw, nz, ny, nx, CUFFT_D2Z, &
                           worksize)
    if (ierr /= 0) then
      write (stderr, *), 'cuFFT Error Code: ', ierr
      error stop 'Forward 3D FFT plan generation failed'
    end if

    ierr = cufftCreate(poisson_fft%plan3D_bw)
    ierr = cufftMpAttachComm(poisson_fft%plan3D_bw, CUFFT_COMM_MPI, &
                             MPI_COMM_WORLD)
    ierr = cufftMakePlan3D(poisson_fft%plan3D_bw, nz, ny, nx, CUFFT_Z2D, &
                           worksize)
    if (ierr /= 0) then
      write (stderr, *), 'cuFFT Error Code: ', ierr
      error stop 'Backward 3D FFT plan generation failed'
    end if

    ! allocate storage for cuFFTMp
    ierr = cufftXtMalloc(poisson_fft%plan3D_fw, poisson_fft%xtdesc, &
                         CUFFT_XT_FORMAT_INPLACE)
    if (ierr /= 0) then
      write (stderr, *), 'cuFFT Error Code: ', ierr
      error stop 'cufftXtMalloc failed'
    end if

  end function init

  subroutine fft_forward_cuda(self, f)
    implicit none

    class(cuda_poisson_fft_t) :: self
    class(field_t), intent(in) :: f

    real(dp), device, pointer :: padded_dev(:, :, :), d_dev(:, :, :)

    type(cudaXtDesc), pointer :: descriptor

    integer :: tsize, ierr
    type(dim3) :: blocks, threads

    select type (f)
    type is (cuda_field_t)
      padded_dev => f%data_d
    end select

    !self%host_data = padded_dev
    !print*, 'host data before forward fft', poisson_fft%host_data(1, :, 1)

    call c_f_pointer(self%xtdesc%descriptor, descriptor)
    call c_f_pointer(descriptor%data(1), d_dev, &
                     [self%nx_loc + 2, self%ny_loc, self%nz_loc])

    ! tsize is different than SZ, because here we work on a 3D Cartesian
    ! data structure, and free to specify any suitable thread/block size.
    tsize = 16
    blocks = dim3((self%ny_loc - 1)/tsize + 1, self%nz_loc, 1)
    threads = dim3(tsize, 1, 1)

    call memcpy3D<<<blocks, threads>>>(d_dev, padded_dev, & !&
                                       self%nx_loc, self%ny_loc, self%nz_loc)

    ierr = cufftXtExecDescriptor(self%plan3D_fw, self%xtdesc, self%xtdesc, &
                                 CUFFT_FORWARD)

    if (ierr /= 0) then
      write (stderr, *), 'cuFFT Error Code: ', ierr
      error stop 'Forward 3D FFT execution failed'
    end if

  end subroutine fft_forward_cuda

  subroutine fft_backward_cuda(self, f)
    implicit none

    class(cuda_poisson_fft_t) :: self
    class(field_t), intent(inout) :: f

    real(dp), device, pointer :: padded_dev(:, :, :), d_dev(:, :, :)

    type(cudaXtDesc), pointer :: descriptor

    integer :: tsize, ierr
    type(dim3) :: blocks, threads

    ierr = cufftXtExecDescriptor(self%plan3D_bw, self%xtdesc, self%xtdesc, &
                                 CUFFT_INVERSE)
    if (ierr /= 0) then
      write (stderr, *), 'cuFFT Error Code: ', ierr
      error stop 'Backward 3D FFT execution failed'
    end if

    select type (f)
    type is (cuda_field_t)
      padded_dev => f%data_d
    end select

    call c_f_pointer(self%xtdesc%descriptor, descriptor)
    call c_f_pointer(descriptor%data(1), d_dev, &
                     [self%nx_loc + 2, self%ny_loc, self%nz_loc])

    tsize = 16
    blocks = dim3((self%ny_loc - 1)/tsize + 1, self%nz_loc, 1)
    threads = dim3(tsize, 1, 1)
    call memcpy3D<<<blocks, threads>>>(padded_dev, d_dev, & !&
                                       self%nx_loc, self%ny_loc, self%nz_loc)

  end subroutine fft_backward_cuda

  subroutine fft_postprocess_000_cuda(self)
    implicit none

    class(cuda_poisson_fft_t) :: self

    type(cudaXtDesc), pointer :: descriptor

    complex(dp), device, dimension(:, :, :), pointer :: c_dev
    type(dim3) :: blocks, threads
    integer :: tsize

    ! obtain a pointer to descriptor so that we can carry out postprocessing
    call c_f_pointer(self%xtdesc%descriptor, descriptor)
    call c_f_pointer(descriptor%data(1), c_dev, &
                     [self%nx_spec, self%ny_spec, self%nz_spec])

    ! tsize is different than SZ, because here we work on a 3D Cartesian
    ! data structure, and free to specify any suitable thread/block size.
    tsize = 16
    blocks = dim3((self%ny_spec - 1)/tsize + 1, self%nz_spec, 1)
    threads = dim3(tsize, 1, 1)

    ! Postprocess div_u in spectral space
    call process_spectral_000<<<blocks, threads>>>( & !&
      c_dev, self%waves_dev, self%nx_spec, self%ny_spec, self%y_sp_st, &
      self%nx_glob, self%ny_glob, self%nz_glob, &
      self%ax_dev, self%bx_dev, self%ay_dev, self%by_dev, &
      self%az_dev, self%bz_dev &
      )

  end subroutine fft_postprocess_000_cuda

  subroutine fft_postprocess_010_cuda(self)
    implicit none

    class(cuda_poisson_fft_t) :: self

    type(cudaXtDesc), pointer :: descriptor

    complex(dp), device, dimension(:, :, :), pointer :: c_dev
    type(dim3) :: blocks, threads
    integer :: tsize

    integer :: i, j, k, ix, iy, iz
    real(dp) :: tmp_r, tmp_c, div_r, div_c, l_r, l_c, r_r, r_c, &
                l1, l2, l3, l4, r1, r2, r3, r4

    ! obtain a pointer to descriptor so that we can carry out postprocessing
    call c_f_pointer(self%xtdesc%descriptor, descriptor)
    call c_f_pointer(descriptor%data(1), c_dev, &
                     [self%nx_spec, self%ny_spec, self%nz_spec])

    ! tsize is different than SZ, because here we work on a 3D Cartesian
    ! data structure, and free to specify any suitable thread/block size.
    tsize = 16
!    blocks = dim3((self%nx_spec - 1)/tsize + 1, self%nz_spec, 1)
!    threads = dim3(tsize, 1, 1)

    blocks = dim3(self%nz_spec, 1, 1)
    !threads = dim3(self%nx_spec, 1, 1)
    threads = dim3(256, 1, 1)
!print*, 'postprocess', blocks, threads
    ! Postprocess div_u in spectral space
!    call process_spectral_010<<<blocks, threads>>>( & !&
!      c_dev, self%waves_dev, self%nx_spec, self%ny_spec, self%y_sp_st, &
!      self%nx_glob, self%ny_glob, self%nz_glob, &
!      self%ax_dev, self%bx_dev, self%ay_dev, self%by_dev, &
!      self%az_dev, self%bz_dev &
!      )

if (.false.) then
    call process_spectral_010_xz_fw<<<blocks, threads>>>( & !&
      c_dev, self%waves_dev, self%nx_spec, self%ny_spec, self%y_sp_st, &
      self%nx_glob, self%ny_glob, self%nz_glob, &
      self%ax_dev, self%bx_dev, self%ay_dev, self%by_dev, &
      self%az_dev, self%bz_dev &
      )
    call process_spectral_010_ysy<<<blocks, threads>>>( & !&
      c_dev, self%temp_dev, self%waves_dev, self%nx_spec, self%ny_spec, self%y_sp_st, &
      self%nx_glob, self%ny_glob, self%nz_glob, &
      self%ax_dev, self%bx_dev, self%ay_dev, self%by_dev, &
      self%az_dev, self%bz_dev &
      )
    call process_spectral_010_xz_bw<<<blocks, threads>>>( & !&
      c_dev, self%waves_dev, self%nx_spec, self%ny_spec, self%y_sp_st, &
      self%nx_glob, self%ny_glob, self%nz_glob, &
      self%ax_dev, self%bx_dev, self%ay_dev, self%by_dev, &
      self%az_dev, self%bz_dev &
      )
end if

if (.true.) then
    !sort things at host side
    self%host_cdata = c_dev
    !normalise
    self%host_cdata = self%host_cdata/(self%nx_glob*self%ny_glob*self%nz_glob)
    ! hack
    !self%host_cdata = 0
    ! postproces in z
    do k = 1, self%nz_spec
      do j = 1, self%ny_spec
        do i = 1, self%nx_spec
          ix = i; iy = j + self%y_sp_st; iz = k
          div_r = real(self%host_cdata(i, j, k), kind=dp)
          div_c = aimag(self%host_cdata(i, j, k))

          self%host_cdata(i, j, k) = cmplx(div_r*self%bz(iz) + div_c*self%az(iz), &
                                           div_c*self%bz(iz) - div_r*self%az(iz), kind=dp)
          if (iz > self%nz_glob/2 + 1) self%host_cdata(i, j, k) = -self%host_cdata(i, j, k)
!if (abs(self%host_cdata(i, j, k)) > 1.0e-4) print*, 'bw z, >e-4 at', i, j, k, self%host_cdata(i, j, k)
        end do
      end do
    end do
    ! postproces in x
    do k = 1, self%nz_spec
      do j = 1, self%ny_spec
        do i = 1, self%nx_spec
          ix = i; iy = j + self%y_sp_st; iz = k
          div_r = real(self%host_cdata(i, j, k), kind=dp)
          div_c = aimag(self%host_cdata(i, j, k))

          self%host_cdata(i, j, k) = cmplx(div_r*self%bx(ix) + div_c*self%ax(ix), &
                                           div_c*self%bx(ix) - div_r*self%ax(ix), kind=dp)
!if (abs(self%host_cdata(i, j, k)) > 1.0e-4) print*, 'bw x, >e-4 at', i, j, k, self%host_cdata(i, j, k)
        end do
      end do
    end do
    ! postprocess in y
    do k = 1, self%nz_spec
      do i = 1, self%nx_spec
        self%host_cdata_2(i, 1, k) = self%host_cdata(i, 1, k)
        do j = 2, self%ny_spec
          ix = i; iy = j + self%y_sp_st; iz = k

          l_r = real(self%host_cdata(i, j, k), kind=dp)
          l_c = aimag(self%host_cdata(i, j, k))
          r_r = real(self%host_cdata(i, self%ny_spec - j + 2, k), kind=dp)
          r_c = aimag(self%host_cdata(i, self%ny_spec - j + 2, k))
          l1 = l_r*self%by(iy)
          l2 = l_r*self%ay(iy)
          l3 = l_c*self%by(iy)
          l4 = l_c*self%ay(iy)
          r1 = r_r*self%by(iy)
          r2 = r_r*self%ay(iy)
          r3 = r_c*self%by(iy)
          r4 = r_c*self%ay(iy)

          ! update the entry
          self%host_cdata_2(i, j, k) = 0.5_dp*cmplx(l1 - l4 + r1 - r4, &
                                                    -l2 + l3 + r2 + r3, kind=dp)
!if (abs(self%host_cdata_2(i, j, k)) > 1.0e-4) print*, 'fw y, >e-4 at', i, j, k, self%host_cdata_2(i, j, k)
        end do
      end do
    end do
    !solve poisson
    do k = 1, self%nz_spec
      do j = 1, self%ny_spec
        do i = 1, self%nx_spec
          div_r = real(self%host_cdata_2(i, j, k), kind=dp)
          div_c = aimag(self%host_cdata_2(i, j, k))

          tmp_r = real(self%waves(i, j, k), kind=dp)
          tmp_c = aimag(self%waves(i, j, k))
          if (abs(tmp_r) < 1.e-14_dp) then
            div_r = 0._dp
          else
            div_r = -div_r/tmp_r
          end if
          if (abs(tmp_c) < 1.e-14_dp) then
            div_c = 0._dp
          else
            div_c = -div_c/tmp_c
          end if

          ! update the entry
          self%host_cdata_2(i, j, k) = cmplx(div_r, div_c, kind=dp)
          if (ix == self%nx_glob/2 + 1 .and. iz == self%nz_glob/2 + 1) self%host_cdata_2(i, j, k) = 0._dp
!if (abs(self%host_cdata_2(i, j, k)) > 1.0e-4) print*, 'solve, >e-4 at', i, j, k, self%host_cdata_2(i, j, k)
        end do
      end do
    end do

    !post-process backward
    ! post process in y
    do k = 1, self%nz_spec
      do i = 1, self%nx_spec
        self%host_cdata(i, 1, k) = self%host_cdata_2(i, 1, k)
        do j = 2, self%ny_spec
          ix = i; iy = j + self%y_sp_st; iz = k

          l_r = real(self%host_cdata_2(i, j, k), kind=dp)
          l_c = aimag(self%host_cdata_2(i, j, k))
          r_r = real(self%host_cdata_2(i, self%ny_spec - j + 2, k), kind=dp)
          r_c = aimag(self%host_cdata_2(i, self%ny_spec - j + 2, k))
          l1 = l_r*self%by(iy)
          l2 = l_r*self%ay(iy)
          l3 = l_c*self%by(iy)
          l4 = l_c*self%ay(iy)
          r1 = r_r*self%by(iy)
          r2 = r_r*self%ay(iy)
          r3 = r_c*self%by(iy)
          r4 = r_c*self%ay(iy)

          ! update the entry
          self%host_cdata(i, j, k) = cmplx(l1 - l4 + r2 + r3, &
                                           l2 + l3 - r1 + r4, kind=dp)
!if (abs(self%host_cdata(i, j, k)) > 1.0e-4) print*, 'bw y, >e-4 at', i, j, k, self%host_cdata(i, j, k)
        end do
      end do
    end do

    ! postproces in x
    do k = 1, self%nz_spec
      do j = 1, self%ny_spec
        do i = 1, self%nx_spec
          ix = i; iy = j + self%y_sp_st; iz = k
          div_r = real(self%host_cdata(i, j, k), kind=dp)
          div_c = aimag(self%host_cdata(i, j, k))

          self%host_cdata(i, j, k) = cmplx(div_r*self%bx(ix) - div_c*self%ax(ix), &
                                           div_c*self%bx(ix) + div_r*self%ax(ix), kind=dp)
!if (abs(self%host_cdata(i, j, k)) > 1.0e-4) print*, 'bw x, >e-4 at', i, j, k, self%host_cdata(i, j, k)
        end do
      end do
    end do
    ! postproces in z
    do k = 1, self%nz_spec
      do j = 1, self%ny_spec
        do i = 1, self%nx_spec
          ix = i; iy = j + self%y_sp_st; iz = k
          div_r = real(self%host_cdata(i, j, k), kind=dp)
          div_c = aimag(self%host_cdata(i, j, k))

          self%host_cdata(i, j, k) = cmplx(div_r*self%bz(iz) - div_c*self%az(iz), &
                                           div_c*self%bz(iz) + div_r*self%az(iz), kind=dp)
          if (iz > self%nz_glob/2 + 1) self%host_cdata(i, j, k) = -self%host_cdata(i, j, k)
!if (abs(self%host_cdata(i, j, k)) > 1.0e-4) print*, 'bw z, >e-4 at', i, j, k, self%host_cdata(i, j, k)
        end do
      end do
    end do

    ! move back to the descriptor so that it can be inverse fft'd
    c_dev = self%host_cdata
end if

  end subroutine fft_postprocess_010_cuda

  subroutine enforce_periodicity_y_cuda(self, f_out, f_in)
    implicit none

    class(cuda_poisson_fft_t) :: self
    class(field_t), intent(inout) :: f_out
    class(field_t), intent(in) :: f_in

    real(dp), device, pointer, dimension(:, :, :) :: f_out_dev, f_in_dev
    type(dim3) :: blocks, threads
    integer :: i, j, k

    select type (f_out)
    type is (cuda_field_t)
      f_out_dev => f_out%data_d
    end select
    select type (f_in)
    type is (cuda_field_t)
      f_in_dev => f_in%data_d
    end select

    !self%host_data = f_in_dev
    !print*, 'line before enforcing periodicity'
    !print*, self%host_data(1,:,1)
!    blocks = dim3(self%nz_loc, 1, 1)
!    threads = dim3(self%nx_loc, 1, 1)
!    call enforce_periodicity_y<<<blocks, threads>>>( & !&
!      f_out_dev, f_in_dev, self%ny_glob &
!      )
    !self%host_data = f_out_dev
    !print*, 'line after enforcing periodicity'
    !print*, self%host_data(1,:,1)

    self%host_data = f_in_dev
    print*, 'before enforce', self%host_data(1, :, 1)
    do k = 1, self%nz_loc
       do i = 1, self%nx_loc
          do j = 1, self%ny_glob/2
             self%host_data_2(i,j,k) = self%host_data(i,2*(j-1)+1,k)
          enddo
          do j = self%ny_glob/2 + 1, self%ny_glob
             self%host_data_2(i,j,k) = self%host_data(i,2*self%ny_glob-2*j+2,k)
          enddo
       enddo
    end do
    f_out_dev = self%host_data_2
    print*, 'after enforce', self%host_data_2(1, :, 1)

  end subroutine enforce_periodicity_y_cuda

  subroutine undo_periodicity_y_cuda(self, f_out, f_in)
    implicit none

    class(cuda_poisson_fft_t) :: self
    class(field_t), intent(inout) :: f_out
    class(field_t), intent(in) :: f_in

    real(dp), device, pointer, dimension(:, :, :) :: f_out_dev, f_in_dev
    type(dim3) :: blocks, threads
    integer :: i, j, k

    select type (f_out)
    type is (cuda_field_t)
      f_out_dev => f_out%data_d
    end select
    select type (f_in)
    type is (cuda_field_t)
      f_in_dev => f_in%data_d
    end select

    !self%host_data = f_in_dev
    !print*, 'line before undo periodicity'
    !print*, self%host_data(1,:,1)
!    blocks = dim3(self%nz_loc, 1, 1)
!    threads = dim3(self%nx_loc, 1, 1)
!    call undo_periodicity_y<<<blocks, threads>>>( & !&
!!      f_out_dev, f_in_dev, self%ny_glob &
!      )
    !self%host_data = f_out_dev
    !print*, 'line after undo periodicity'
    !print*, self%host_data(1,:,1)
!
    self%host_data = f_in_dev
    do k = 1, self%nz_loc
       do i = 1, self%nx_loc
          do j = 1, self%ny_glob/2
             self%host_data_2(i,2*j-1,k) = self%host_data(i,j,k)
          enddo
          do j=1,self%ny_glob/2
             self%host_data_2(i,2*j,k) = self%host_data(i,self%ny_glob-j+1,k)
          enddo
       enddo
    end do
    f_out_dev = self%host_data_2
    print*, 'after undo', self%host_data_2(1, :, 1)

  end subroutine undo_periodicity_y_cuda

end module m_cuda_poisson_fft
