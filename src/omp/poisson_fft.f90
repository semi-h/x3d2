module m_omp_poisson_fft

  use decomp_2d_constants, only: PHYSICAL_IN_X
  use decomp_2d_fft, only: decomp_2d_fft_init, decomp_2d_fft_3d, &
                           decomp_2d_fft_get_size
  use m_allocator, only: field_t
  use m_common, only: dp
  use m_poisson_fft, only: poisson_fft_t
  use m_tdsops, only: dirps_t
  use m_mesh, only: mesh_t
  use m_omp_spectral, only: process_spectral_div_u

  implicit none

  type, extends(poisson_fft_t) :: omp_poisson_fft_t
      !! FFT based Poisson solver
    complex(dp), allocatable, dimension(:, :, :) :: c_x, c_y, c_z
  contains
    procedure :: fft_forward => fft_forward_omp
    procedure :: fft_backward => fft_backward_omp
    procedure :: fft_postprocess_000 => fft_postprocess_000_omp
    procedure :: fft_postprocess_010 => fft_postprocess_010_omp
    procedure :: enforce_periodicity_y => enforce_periodicity_y_omp
    procedure :: undo_periodicity_y => undo_periodicity_y_omp
  end type omp_poisson_fft_t

  interface omp_poisson_fft_t
    module procedure init
  end interface omp_poisson_fft_t

  private :: init

contains

  function init(mesh, xdirps, ydirps, zdirps) result(poisson_fft)
    implicit none

    type(mesh_t), intent(in) :: mesh
    class(dirps_t), intent(in) :: xdirps, ydirps, zdirps
    integer, dimension(3) :: istart, iend, isize

    type(omp_poisson_fft_t) :: poisson_fft

    if (mesh%par%is_root()) then
      print *, "Initialising 2decomp&fft"
    end if

    ! Work out the spectral dimensions in the permuted state
    call decomp_2d_fft_init(PHYSICAL_IN_X)
    call decomp_2d_fft_get_size(istart, iend, isize)
    ! Converts a start position into an offset
    istart(:) = istart(:) - 1

    call poisson_fft%base_init(mesh, xdirps, ydirps, zdirps, isize, istart)

    allocate (poisson_fft%c_x(poisson_fft%nx_spec, poisson_fft%ny_spec, &
                              poisson_fft%nz_spec))

  end function init

  subroutine fft_forward_omp(self, f_in)
    implicit none

    class(omp_poisson_fft_t) :: self
    class(field_t), intent(in) :: f_in

    call decomp_2d_fft_3d(f_in%data, self%c_x)

  end subroutine fft_forward_omp

  subroutine fft_backward_omp(self, f_out)
    implicit none

    class(omp_poisson_fft_t) :: self
    class(field_t), intent(inout) :: f_out

    call decomp_2d_fft_3d(self%c_x, f_out%data)

  end subroutine fft_backward_omp

  subroutine fft_postprocess_000_omp(self)
    implicit none

    class(omp_poisson_fft_t) :: self

    call process_spectral_div_u( &
      self%c_x, self%waves, self%nx_spec, self%ny_spec, self%nz_spec, &
      self%x_sp_st, self%y_sp_st, self%z_sp_st, &
      self%nx_glob, self%ny_glob, self%nz_glob, &
      self%ax, self%bx, self%ay, self%by, self%az, self%bz &
      )

  end subroutine fft_postprocess_000_omp

  subroutine fft_postprocess_010_omp(self)
    implicit none

    class(omp_poisson_fft_t) :: self

    real(dp) :: dir_r, div_c, temp_r, temp_c
    integer :: i, j, k, ix, iy, iz

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
          self%host_cdata_2(i, j, k) = 0.5_dp*cmplx(l1 + l4 + r1 - r4, &
                                                    -l2 + l3 + r2 + r3, kind=dp)
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
          if (abs(tmp_r) < 1.e-16_dp) then
            div_r = 0._dp
          else
            div_r = -div_r/tmp_r
          end if
          if (abs(tmp_c) < 1.e-16_dp) then
            div_c = 0._dp
          else
            div_c = -div_c/tmp_c
          end if

          ! update the entry
          self%host_cdata_2(i, j, k) = cmplx(div_r, div_c, kind=dp)
          if (ix == self%nx_glob/2 + 1 .and. iz == self%nz_glob/2 + 1) then
            self%host_cdata_2(i, j, k) = cmplx(0._dp, 0._dp, kind=dp)
          end if
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
        end do
      end do
    end do

  end subroutine fft_postprocess_010_omp

  subroutine enforce_periodicity_y_omp(self, f_out, f_in)
    implicit none

    class(omp_poisson_fft_t) :: self
    class(field_t), intent(inout) :: f_out
    class(field_t), intent(in) :: f_in

    integer :: i, j, k

    do k = 1, self%nz_loc
       do i = 1, self%nx_loc
          do j = 1, self%ny_glob/2
             f_out%data(i, j, k) = f_in%data(i, 2*(j - 1) + 1, k)
          end do
          do j = self%ny_glob/2 + 1, self%ny_glob
             f_out%data(i, j, k) = f_in%data(i, 2*self%ny_glob - 2*j + 2, k)
          end do
       end do
    end do

  end subroutine enforce_periodicity_y_omp

  subroutine undo_periodicity_y_omp(self, f_out, f_in)
    implicit none

    class(omp_poisson_fft_t) :: self
    class(field_t), intent(inout) :: f_out
    class(field_t), intent(in) :: f_in

    integer :: i, j, k

    do k = 1, self%nz_loc
       do i = 1, self%nx_loc
          do j = 1, self%ny_glob/2
             f_out%data(i, 2*j - 1, k) = f_in%data(i, j, k)
          end do
          do j = 1, self%ny_glob/2
             f_out%data(i, 2*j, k) = f_in%data(i,self%ny_glob - j + 1, k)
          end do
       end do
    end do

  end subroutine undo_periodicity_y_omp

end module m_omp_poisson_fft
