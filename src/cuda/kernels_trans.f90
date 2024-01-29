module m_cuda_kernels_trans
   use cudafor

   use m_common, only: dp
   use m_cuda_common, only: SZ

contains

   attributes(global) subroutine trans_x2y_k(u_y, u_x, nz)
      implicit none

      real(dp), device, intent(out), dimension(:, :, :) :: u_y
      real(dp), device, intent(in), dimension(:, :, :) :: u_x
      integer, value, intent(in) :: nz

      real(dp), shared :: tile(SZ, SZ)
      integer :: i, j, b_i, b_j, b_k

      i = threadIdx%x; j = threadIdx%y
      b_i = blockIdx%x; b_j = blockIdx%y; b_k = blockIdx%z

      ! copy into shared
      tile(i, j) = u_x(i, j + (b_i - 1)*SZ, b_j + (b_k - 1)*nz)

      call syncthreads()

      ! copy into output array from shared
      u_y(i, j + (b_k - 1)*SZ, b_j + (b_i - 1)*nz) = tile(j, i)

   end subroutine trans_x2y_k

   attributes(global) subroutine trans_x2z_k(u_z, u_x, nz)
      implicit none

      real(dp), device, intent(out), dimension(:, :, :) :: u_z
      real(dp), device, intent(in), dimension(:, :, :) :: u_x
      integer, value, intent(in) :: nz

      integer :: i, j, b_i, b_j, nx!, nz

      i = threadIdx%x; b_i = blockIdx%x; b_j = blockIdx%y
      nx = gridDim%x

      do j = 1, nz
         u_z(i, j, b_i + (b_j - 1)*nx) = u_x(i, b_i, j + (b_j - 1)*nz)
      end do

   end subroutine trans_x2z_k

   attributes(global) subroutine trans_y2x_k(u_x, u_y, nz)
      implicit none

      real(dp), device, intent(out), dimension(:, :, :) :: u_x
      real(dp), device, intent(in), dimension(:, :, :) :: u_y
      integer, value, intent(in) :: nz

      real(dp), shared :: tile(SZ, SZ)
      integer :: i, j, b_i, b_j, b_k

      i = threadIdx%x; j = threadIdx%y
      b_i = blockIdx%x; b_j = blockIdx%y; b_k = blockIdx%z

      ! copy into shared
      tile(i, j) = u_y(i, (b_j - 1)*SZ + j, (b_i - 1)*nz + b_k)

      call syncthreads()

      ! copy into output array from shared
      u_x(i, (b_i - 1)*SZ + j, (b_j - 1)*nz + b_k) = tile(j, i)

   end subroutine trans_y2x_k

   attributes(global) subroutine trans_y2z_k(u_z, u_y, nx, nz)
      implicit none

      real(dp), device, intent(out), dimension(:, :, :) :: u_z
      real(dp), device, intent(in), dimension(:, :, :) :: u_y
      integer, value, intent(in) :: nx, nz

      real(dp), shared :: tile(SZ, SZ)
      integer :: i, j, b_i, b_j, b_k

      i = threadIdx%x; j = threadIdx%y
      b_i = blockIdx%x; b_j = blockIdx%y; b_k = blockIdx%z

      ! copy into shared
      tile(i, j) = u_y(i, (b_j - 1)*SZ + j, (b_i - 1)*nz + b_k)

      call syncthreads()

      ! copy into output array from shared
      u_z(i, b_k, (b_i - 1)*SZ + j + (b_j - 1)*nx) = tile(j, i)

   end subroutine trans_y2z_k

   attributes(global) subroutine trans_z2y_k(u_y, u_z, nx, nz)
      implicit none

      real(dp), device, intent(out), dimension(:, :, :) :: u_y
      real(dp), device, intent(in), dimension(:, :, :) :: u_z
      integer, value, intent(in) :: nx, nz

      real(dp), shared :: tile(SZ, SZ)
      integer :: i, j, b_i, b_j, b_k

      i = threadIdx%x; j = threadIdx%y
      b_i = blockIdx%x; b_j = blockIdx%y; b_k = blockIdx%z

      ! copy into shared
      tile(i, j) = u_z(i, b_k, (b_i - 1)*SZ + j + (b_j - 1)*nx)

      call syncthreads()

      ! copy into output array from shared
      u_y(i, (b_j - 1)*SZ + j, (b_i - 1)*nz + b_k) = tile(j, i)

   end subroutine trans_z2y_k

end module m_cuda_kernels_trans
