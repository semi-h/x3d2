module m_cuda_kernels_thom
  use cudafor

  use m_common, only: dp
  use m_cuda_common, only: SZ

  implicit none

contains

  attributes(global) subroutine der_univ_thom( &
    du, u, coeffs_s, coeffs_e, coeffs, n, thom_f, thom_s, thom_w &
    )
    implicit none

    real(dp), device, intent(out), dimension(:, :, :) :: du
    real(dp), device, intent(in), dimension(:, :, :) :: u
    real(dp), device, intent(in), dimension(:, :) :: coeffs_s, coeffs_e
    real(dp), device, intent(in), dimension(:) :: coeffs
    integer, value, intent(in) :: n
    real(dp), device, intent(in), dimension(:) :: thom_f, thom_s, thom_w

    integer :: i, j, b

    real(dp) :: c_m4, c_m3, c_m2, c_m1, c_j, c_p1, c_p2, c_p3, c_p4, temp_du

    i = threadIdx%x
    b = blockIdx%x

    ! store bulk coeffs in the registers
    c_m4 = coeffs(1); c_m3 = coeffs(2); c_m2 = coeffs(3); c_m1 = coeffs(4)
    c_j = coeffs(5)
    c_p1 = coeffs(6); c_p2 = coeffs(7); c_p3 = coeffs(8); c_p4 = coeffs(9)

    du(i, 1, b) = coeffs_s(5, 1)*u(i, 1, b) &
                  + coeffs_s(6, 1)*u(i, 2, b) &
                  + coeffs_s(7, 1)*u(i, 3, b) &
                  + coeffs_s(8, 1)*u(i, 4, b) &
                  + coeffs_s(9, 1)*u(i, 5, b)
    du(i, 1, b) = du(i, 1, b)
    du(i, 2, b) = coeffs_s(4, 2)*u(i, 1, b) &
                  + coeffs_s(5, 2)*u(i, 2, b) &
                  + coeffs_s(6, 2)*u(i, 3, b) &
                  + coeffs_s(7, 2)*u(i, 4, b) &
                  + coeffs_s(8, 2)*u(i, 5, b) &
                  + coeffs_s(9, 2)*u(i, 6, b)
    du(i, 2, b) = du(i, 2, b) - du(i, 1, b)*thom_s(2)
    du(i, 3, b) = coeffs_s(3, 3)*u(i, 1, b) &
                  + coeffs_s(4, 3)*u(i, 2, b) &
                  + coeffs_s(5, 3)*u(i, 3, b) &
                  + coeffs_s(6, 3)*u(i, 4, b) &
                  + coeffs_s(7, 3)*u(i, 5, b) &
                  + coeffs_s(8, 3)*u(i, 6, b) &
                  + coeffs_s(9, 3)*u(i, 7, b)
    du(i, 3, b) = du(i, 3, b) - du(i, 2, b)*thom_s(3)
    du(i, 4, b) = coeffs_s(2, 4)*u(i, 1, b) &
                  + coeffs_s(3, 4)*u(i, 2, b) &
                  + coeffs_s(4, 4)*u(i, 3, b) &
                  + coeffs_s(5, 4)*u(i, 4, b) &
                  + coeffs_s(6, 4)*u(i, 5, b) &
                  + coeffs_s(7, 4)*u(i, 6, b) &
                  + coeffs_s(8, 4)*u(i, 7, b) &
                  + coeffs_s(9, 4)*u(i, 8, b)
    du(i, 4, b) = du(i, 4, b) - du(i, 3, b)*thom_s(4)

    do j = 5, n - 4
      temp_du = c_m4*u(i, j - 4, b) + c_m3*u(i, j - 3, b) &
                + c_m2*u(i, j - 2, b) + c_m1*u(i, j - 1, b) &
                + c_j*u(i, j, b) &
                + c_p1*u(i, j + 1, b) + c_p2*u(i, j + 2, b) &
                + c_p3*u(i, j + 3, b) + c_p4*u(i, j + 4, b)
      du(i, j, b) = temp_du - du(i, j - 1, b)*thom_s(j)
    end do

    j = n - 3
    du(i, j, b) = coeffs_e(1, 1)*u(i, j - 4, b) &
                  + coeffs_e(2, 1)*u(i, j - 3, b) &
                  + coeffs_e(3, 1)*u(i, j - 2, b) &
                  + coeffs_e(4, 1)*u(i, j - 1, b) &
                  + coeffs_e(5, 1)*u(i, j, b) &
                  + coeffs_e(6, 1)*u(i, j + 1, b) &
                  + coeffs_e(7, 1)*u(i, j + 2, b) &
                  + coeffs_e(8, 1)*u(i, j + 3, b)
    du(i, j, b) = du(i, j, b) - du(i, j - 1, b)*thom_s(j)
    j = n - 2
    du(i, j, b) = coeffs_e(1, 2)*u(i, j - 4, b) &
                  + coeffs_e(2, 2)*u(i, j - 3, b) &
                  + coeffs_e(3, 2)*u(i, j - 2, b) &
                  + coeffs_e(4, 2)*u(i, j - 1, b) &
                  + coeffs_e(5, 2)*u(i, j, b) &
                  + coeffs_e(6, 2)*u(i, j + 1, b) &
                  + coeffs_e(7, 2)*u(i, j + 2, b)
    du(i, j, b) = du(i, j, b) - du(i, j - 1, b)*thom_s(j)
    j = n - 1
    du(i, j, b) = coeffs_e(1, 3)*u(i, j - 4, b) &
                  + coeffs_e(2, 3)*u(i, j - 3, b) &
                  + coeffs_e(3, 3)*u(i, j - 2, b) &
                  + coeffs_e(4, 3)*u(i, j - 1, b) &
                  + coeffs_e(5, 3)*u(i, j, b) &
                  + coeffs_e(6, 3)*u(i, j + 1, b)
    du(i, j, b) = du(i, j, b) - du(i, j - 1, b)*thom_s(j)
    j = n
    du(i, j, b) = coeffs_e(1, 4)*u(i, j - 4, b) &
                  + coeffs_e(2, 4)*u(i, j - 3, b) &
                  + coeffs_e(3, 4)*u(i, j - 2, b) &
                  + coeffs_e(4, 4)*u(i, j - 1, b) &
                  + coeffs_e(5, 4)*u(i, j, b)
    du(i, j, b) = du(i, j, b) - du(i, j - 1, b)*thom_s(j)

    ! Backward pass of the Thomas algorithm
    du(i, n, b) = du(i, n, b)*thom_w(n)
    do j = n - 1, 1, -1
      du(i, j, b) = (du(i, j, b) - thom_f(j)*du(i, j + 1, b))*thom_w(j)
    end do

  end subroutine der_univ_thom

  attributes(global) subroutine der_univ_thom_per( &
    du, u, coeffs, n, alpha, thom_f, thom_s, thom_w, thom_p &
    )
    implicit none

    real(dp), device, intent(out), dimension(:, :, :) :: du
    real(dp), device, intent(in), dimension(:, :, :) :: u
    real(dp), device, intent(in), dimension(:) :: coeffs
    integer, value, intent(in) :: n
    real(dp), value, intent(in) :: alpha
    real(dp), device, intent(in), dimension(:) :: thom_f, thom_s, thom_w, &
                                                  thom_p

    integer :: i, j, b
    integer :: jm4, jm3, jm2, jm1, jp1, jp2, jp3, jp4

    real(dp) :: c_m4, c_m3, c_m2, c_m1, c_j, c_p1, c_p2, c_p3, c_p4
    real(dp) :: temp_du, ss

    i = threadIdx%x
    b = blockIdx%x

    ! store bulk coeffs in the registers
    c_m4 = coeffs(1); c_m3 = coeffs(2); c_m2 = coeffs(3); c_m1 = coeffs(4)
    c_j = coeffs(5)
    c_p1 = coeffs(6); c_p2 = coeffs(7); c_p3 = coeffs(8); c_p4 = coeffs(9)

    do j = 1, n
      jm4 = modulo(j - 5, n) + 1
      jm3 = modulo(j - 4, n) + 1
      jm2 = modulo(j - 3, n) + 1
      jm1 = modulo(j - 2, n) + 1
      jp1 = modulo(j - n, n) + 1
      jp2 = modulo(j - n + 1, n) + 1
      jp3 = modulo(j - n + 2, n) + 1
      jp4 = modulo(j - n + 3, n) + 1

      temp_du = c_m4*u(i, jm4, b) + c_m3*u(i, jm3, b) &
                + c_m2*u(i, jm2, b) + c_m1*u(i, jm1, b) &
                + c_j*u(i, j, b) &
                + c_p1*u(i, jp1, b) + c_p2*u(i, jp2, b) &
                + c_p3*u(i, jp3, b) + c_p4*u(i, jp4, b)
      du(i, j, b) = temp_du - du(i, jm1, b)*thom_s(j)
    end do

    ! Backward pass of the Thomas algorithm
    du(i, n, b) = du(i, n, b)*thom_w(n)
    do j = n - 1, 1, -1
      du(i, j, b) = (du(i, j, b) - thom_f(j)*du(i, j + 1, b))*thom_w(j)
    end do

    ! Periodic final pass
    ss = (du(i, 1, b) - alpha*du(i, n, b)) &
         /(1._dp + thom_p(1) - alpha*thom_p(n))
    do j = 1, n
      du(i, j, b) = du(i, j, b) - ss*thom_p(j)
    end do

  end subroutine der_univ_thom_per

  attributes(global) subroutine der_univ_thom_shared( &
    du, u, coeffs_s, coeffs_e, coeffs, n, thom_f, thom_s, thom_w &
    )
    implicit none

    real(dp), device, intent(out), dimension(:, :, :) :: du
    real(dp), device, intent(in), dimension(:, :, :) :: u
    real(dp), device, intent(in), dimension(:, :) :: coeffs_s, coeffs_e
    real(dp), device, intent(in), dimension(:) :: coeffs
    integer, value, intent(in) :: n
    real(dp), device, intent(in), dimension(:) :: thom_f, thom_s, thom_w

    real(8), shared :: s(*)
    integer :: i, j, b

    real(dp) :: c_m4, c_m3, c_m2, c_m1, c_j, c_p1, c_p2, c_p3, c_p4, temp_du, temp

    i = threadIdx%x
    b = blockIdx%x

    ! store bulk coeffs in the registers
    c_m4 = coeffs(1); c_m3 = coeffs(2); c_m2 = coeffs(3); c_m1 = coeffs(4)
    c_j = coeffs(5)
    c_p1 = coeffs(6); c_p2 = coeffs(7); c_p3 = coeffs(8); c_p4 = coeffs(9)

    s(i + 0*SZ) = coeffs_s(5, 1)*u(i, 1, b) &
                  + coeffs_s(6, 1)*u(i, 2, b) &
                  + coeffs_s(7, 1)*u(i, 3, b) &
                  + coeffs_s(8, 1)*u(i, 4, b) &
                  + coeffs_s(9, 1)*u(i, 5, b)
    s(i + 1*SZ) = coeffs_s(4, 2)*u(i, 1, b) &
                  + coeffs_s(5, 2)*u(i, 2, b) &
                  + coeffs_s(6, 2)*u(i, 3, b) &
                  + coeffs_s(7, 2)*u(i, 4, b) &
                  + coeffs_s(8, 2)*u(i, 5, b) &
                  + coeffs_s(9, 2)*u(i, 6, b) &
                  - s(i + 0*SZ)*thom_s(2)
    s(i + 2*SZ) = coeffs_s(3, 3)*u(i, 1, b) &
                  + coeffs_s(4, 3)*u(i, 2, b) &
                  + coeffs_s(5, 3)*u(i, 3, b) &
                  + coeffs_s(6, 3)*u(i, 4, b) &
                  + coeffs_s(7, 3)*u(i, 5, b) &
                  + coeffs_s(8, 3)*u(i, 6, b) &
                  + coeffs_s(9, 3)*u(i, 7, b) &
                  - s(i + 1*SZ)*thom_s(3)
    s(i + 3*SZ) = coeffs_s(2, 4)*u(i, 1, b) &
                  + coeffs_s(3, 4)*u(i, 2, b) &
                  + coeffs_s(4, 4)*u(i, 3, b) &
                  + coeffs_s(5, 4)*u(i, 4, b) &
                  + coeffs_s(6, 4)*u(i, 5, b) &
                  + coeffs_s(7, 4)*u(i, 6, b) &
                  + coeffs_s(8, 4)*u(i, 7, b) &
                  + coeffs_s(9, 4)*u(i, 8, b) &
                  - s(i + 2*SZ)*thom_s(4)

    do j = 5, n - 4
      s(i + (j - 1)*SZ) = c_m4*u(i, j - 4, b) + c_m3*u(i, j - 3, b) &
                          + c_m2*u(i, j - 2, b) + c_m1*u(i, j - 1, b) &
                          + c_j*u(i, j, b) &
                          + c_p1*u(i, j + 1, b) + c_p2*u(i, j + 2, b) &
                          + c_p3*u(i, j + 3, b) + c_p4*u(i, j + 4, b) &
                          - s(i + (j - 2)*SZ)*thom_s(j)
    end do

    j = n - 3
    s(i + (j - 1)*SZ) = coeffs_e(1, 1)*u(i, j - 4, b) &
                        + coeffs_e(2, 1)*u(i, j - 3, b) &
                        + coeffs_e(3, 1)*u(i, j - 2, b) &
                        + coeffs_e(4, 1)*u(i, j - 1, b) &
                        + coeffs_e(5, 1)*u(i, j, b) &
                        + coeffs_e(6, 1)*u(i, j + 1, b) &
                        + coeffs_e(7, 1)*u(i, j + 2, b) &
                        + coeffs_e(8, 1)*u(i, j + 3, b) &
                        - s(i + (j - 2)*SZ)*thom_s(j)
    j = n - 2
    s(i + (j - 1)*SZ) = coeffs_e(1, 2)*u(i, j - 4, b) &
                        + coeffs_e(2, 2)*u(i, j - 3, b) &
                        + coeffs_e(3, 2)*u(i, j - 2, b) &
                        + coeffs_e(4, 2)*u(i, j - 1, b) &
                        + coeffs_e(5, 2)*u(i, j, b) &
                        + coeffs_e(6, 2)*u(i, j + 1, b) &
                        + coeffs_e(7, 2)*u(i, j + 2, b) &
                        - s(i + (j - 2)*SZ)*thom_s(j)
    j = n - 1
    s(i + (j - 1)*SZ) = coeffs_e(1, 3)*u(i, j - 4, b) &
                        + coeffs_e(2, 3)*u(i, j - 3, b) &
                        + coeffs_e(3, 3)*u(i, j - 2, b) &
                        + coeffs_e(4, 3)*u(i, j - 1, b) &
                        + coeffs_e(5, 3)*u(i, j, b) &
                        + coeffs_e(6, 3)*u(i, j + 1, b) &
                        - s(i + (j - 2)*SZ)*thom_s(j)
    j = n
    s(i + (j - 1)*SZ) = coeffs_e(1, 4)*u(i, j - 4, b) &
                        + coeffs_e(2, 4)*u(i, j - 3, b) &
                        + coeffs_e(3, 4)*u(i, j - 2, b) &
                        + coeffs_e(4, 4)*u(i, j - 1, b) &
                        + coeffs_e(5, 4)*u(i, j, b) &
                        - s(i + (j - 2)*SZ)*thom_s(j)

    ! Backward pass of the Thomas algorithm
    !du(i, n, b) = s(i + (n - 1)*SZ)*thom_w(n)
    !do j = n - 1, 1, -1
    !  du(i, j, b) = (s(i + (j - 1)*SZ) - thom_f(j)*s(i + (j + 1 - 1)*SZ))*thom_w(j)
    !end do

    !du(i, n, b) = s(i + (n - 1)*SZ)*thom_w(n)

    !s(i + (n - 1)*SZ) = s(i + (n - 1)*SZ)*thom_w(n)
    !du(i, n, b) = s(i + (n - 1)*SZ)
    temp = s(i + (n - 1)*SZ)*thom_w(n)
    du(i, n, b) = temp
    do j = n - 1, 1, -1
    !  !!du(i, j, b) = (s(i + (j - 1)*SZ) - thom_f(j)*s(i + (j + 1 - 1)*SZ))*thom_w(j)
    !  s(i + (j - 1)*SZ) = (s(i + (j - 1)*SZ) - thom_f(j)*s(i + (j + 1 - 1)*SZ))*thom_w(j)
    !  du(i, j, b) = s(i + (j - 1)*SZ)
      
      temp = (s(i + (j - 1)*SZ) - thom_f(j)*temp)*thom_w(j)
      du(i, j, b) = temp
    end do

    ! Backward pass of the Thomas algorithm
    !du(i, n, b) = du(i, n, b)*thom_w(n)
    !do j = n - 1, 1, -1
    !  du(i, j, b) = (du(i, j, b) - thom_f(j)*du(i, j + 1, b))*thom_w(j)
    !end do

  end subroutine der_univ_thom_shared

  attributes(global) subroutine der_univ_thom_shared_per( &
    du, u, coeffs, n, alpha, thom_f, thom_s, thom_w, thom_p &
    )
    implicit none

    real(dp), device, intent(out), dimension(:, :, :) :: du
    real(dp), device, intent(in), dimension(:, :, :) :: u
    real(dp), device, intent(in), dimension(:) :: coeffs
    integer, value, intent(in) :: n
    real(dp), value, intent(in) :: alpha
    real(dp), device, intent(in), dimension(:) :: thom_f, thom_s, thom_w, &
                                                  thom_p

    real(8), shared :: s(*)

    integer :: i, j, b
    integer :: jm4, jm3, jm2, jm1, jp1, jp2, jp3, jp4

    real(dp) :: c_m4, c_m3, c_m2, c_m1, c_j, c_p1, c_p2, c_p3, c_p4
    real(dp) :: temp_du, ss, temp_j, temp_jp1
    real(dp) :: temp_u_0, temp_u_1, temp_u_2, temp_u_3, temp_u_4

    i = threadIdx%x
    b = blockIdx%x

    ! store bulk coeffs in the registers
    c_m4 = coeffs(1); c_m3 = coeffs(2); c_m2 = coeffs(3); c_m1 = coeffs(4)
    c_j = coeffs(5)
    c_p1 = coeffs(6); c_p2 = coeffs(7); c_p3 = coeffs(8); c_p4 = coeffs(9)

!    do j = 1, n
!      jm4 = modulo(j - 5, n) + 1
!      jm3 = modulo(j - 4, n) + 1
!      jm2 = modulo(j - 3, n) + 1
!      jm1 = modulo(j - 2, n) + 1
!      jp1 = modulo(j - n, n) + 1
!      jp2 = modulo(j - n + 1, n) + 1
!      jp3 = modulo(j - n + 2, n) + 1
!      jp4 = modulo(j - n + 3, n) + 1
!
!      !temp_du = c_m4*u(i, jm4, b) + c_m3*u(i, jm3, b) &
!      !          + c_m2*u(i, jm2, b) + c_m1*u(i, jm1, b) &
!      !          + c_j*u(i, j, b) &
!      !          + c_p1*u(i, jp1, b) + c_p2*u(i, jp2, b) &
!      !          + c_p3*u(i, jp3, b) + c_p4*u(i, jp4, b)
!      !s(i + (j - 1)*SZ) = temp_du - s(i + (jm1 - 1)*SZ)*thom_s(j)
!      s(i + (j - 1)*SZ) = c_m4*u(i, jm4, b) + c_m3*u(i, jm3, b) &
!                          + c_m2*u(i, jm2, b) + c_m1*u(i, jm1, b) &
!                          + c_j*u(i, j, b) &
!                          + c_p1*u(i, jp1, b) + c_p2*u(i, jp2, b) &
!                          + c_p3*u(i, jp3, b) + c_p4*u(i, jp4, b) &
!                          - s(i + (jm1 - 1)*SZ)*thom_s(j)
!    end do

    !! test  shared copy
    !do j = 1, n
    !  s(i + (j - 1)*SZ) = u(i, j, b)
    !end do

    temp_u_3 = u(i, n-3, b)
    temp_u_2 = u(i, n-2, b)
    temp_u_1 = u(i, n-1, b)
    temp_u_0 = u(i, n, b)

    s(i + 0*SZ) = c_m4*temp_u_3 + c_m3*temp_u_2 &
                  + c_m2*temp_u_1 + c_m1*temp_u_0 &
                  + c_j*u(i, 1, b) &
                  + c_p1*u(i, 2, b) + c_p2*u(i, 3, b) &
                  + c_p3*u(i, 4, b) + c_p4*u(i, 5, b)
    s(i + 1*SZ) = c_m4*temp_u_2 + c_m3*temp_u_1 &
                  + c_m2*temp_u_0 + c_m1*u(i, 1, b) &
                  + c_j*u(i, 2, b) &
                  + c_p1*u(i, 3, b) + c_p2*u(i, 4, b) &
                  + c_p3*u(i, 5, b) + c_p4*u(i, 6, b) &
                  - s(i + 0*SZ)*thom_s(2)
    s(i + 2*SZ) = c_m4*temp_u_1 + c_m3*temp_u_0 &
                  + c_m2*u(i, 1, b) + c_m1*u(i, 2, b) &
                  + c_j*u(i, 3, b) &
                  + c_p1*u(i, 4, b) + c_p2*u(i, 5, b) &
                  + c_p3*u(i, 6, b) + c_p4*u(i, 7, b) &
                  - s(i + 1*SZ)*thom_s(3)
    s(i + 3*SZ) = c_m4*temp_u_0 + c_m3*u(i, 1, b) &
                  + c_m2*u(i, 2, b) + c_m1*u(i, 3, b) &
                  + c_j*u(i, 4, b) &
                  + c_p1*u(i, 5, b) + c_p2*u(i, 6, b) &
                  + c_p3*u(i, 7, b) + c_p4*u(i, 8, b) &
                  - s(i + 2*SZ)*thom_s(4)

!    s(i + 0*SZ) = c_m4*u(i, n-3, b) + c_m3*u(i, n-2, b) &
!                  + c_m2*u(i, n-1, b) + c_m1*u(i, n, b) &
!                  + c_j*u(i, 1, b) &
!                  + c_p1*u(i, 2, b) + c_p2*u(i, 3, b) &
!                  + c_p3*u(i, 4, b) + c_p4*u(i, 5, b)
!    s(i + 1*SZ) = c_m4*u(i, n-2, b) + c_m3*u(i, n-1, b) &
!                  + c_m2*u(i, n, b) + c_m1*u(i, 1, b) &
!                  + c_j*u(i, 2, b) &
!                  + c_p1*u(i, 3, b) + c_p2*u(i, 4, b) &
!                  + c_p3*u(i, 5, b) + c_p4*u(i, 6, b) &
!                  - s(i + 0*SZ)*thom_s(2)
!    s(i + 2*SZ) = c_m4*u(i, n-1, b) + c_m3*u(i, n, b) &
!                  + c_m2*u(i, 1, b) + c_m1*u(i, 2, b) &
!                  + c_j*u(i, 3, b) &
!                  + c_p1*u(i, 4, b) + c_p2*u(i, 5, b) &
!                  + c_p3*u(i, 6, b) + c_p4*u(i, 7, b) &
!                  - s(i + 1*SZ)*thom_s(3)
!    s(i + 3*SZ) = c_m4*u(i, n, b) + c_m3*u(i, 1, b) &
!                  + c_m2*u(i, 2, b) + c_m1*u(i, 3, b) &
!                  + c_j*u(i, 4, b) &
!                  + c_p1*u(i, 5, b) + c_p2*u(i, 6, b) &
!                  + c_p3*u(i, 7, b) + c_p4*u(i, 8, b) &
!                  - s(i + 2*SZ)*thom_s(4)
    do j = 5, n - 4
      s(i + (j - 1)*SZ) = c_m4*u(i, j - 4, b) + c_m3*u(i, j - 3, b) &
                          + c_m2*u(i, j - 2, b) + c_m1*u(i, j - 1, b) &
                          + c_j*u(i, j, b) &
                          + c_p1*u(i, j + 1, b) + c_p2*u(i, j + 2, b) &
                          + c_p3*u(i, j + 3, b) + c_p4*u(i, j + 4, b) &
                          - s(i + (j - 2)*SZ)*thom_s(j)
    end do

    !j = n - 3
    s(i + (n - 4)*SZ) = c_m4*u(i, n-7, b) + c_m3*u(i, n-6, b) &
                        + c_m2*u(i, n-5, b) + c_m1*u(i, n-4, b) &
                        + c_j*u(i, n-3, b) &
                        + c_p1*u(i, n-2, b) + c_p2*u(i, n-1, b) &
                        + c_p3*u(i, n, b) + c_p4*u(i, 1, b) &
                        - s(i + (n - 5)*SZ)*thom_s(n-3)
    !j = n - 2
    s(i + (n - 3)*SZ) = c_m4*u(i, n-6, b) + c_m3*u(i, n-5, b) &
                        + c_m2*u(i, n-4, b) + c_m1*u(i, n-3, b) &
                        + c_j*u(i, n-2, b) &
                        + c_p1*u(i, n-1, b) + c_p2*u(i, n, b) &
                        + c_p3*u(i, 1, b) + c_p4*u(i, 2, b) &
                        - s(i + (n - 4)*SZ)*thom_s(n-2)
    !j = n - 1
    s(i + (n - 2)*SZ) = c_m4*u(i, n-5, b) + c_m3*u(i, n-4, b) &
                        + c_m2*u(i, n-3, b) + c_m1*u(i, n-2, b) &
                        + c_j*u(i, n-1, b) &
                        + c_p1*u(i, n, b) + c_p2*u(i, 1, b) &
                        + c_p3*u(i, 2, b) + c_p4*u(i, 3, b) &
                        - s(i + (n - 3)*SZ)*thom_s(n-1)
    !j = n
    s(i + (n - 1)*SZ) = c_m4*u(i, n-4, b) + c_m3*u(i, n-3, b) &
                        + c_m2*u(i, n-2, b) + c_m1*u(i, n-1, b) &
                        + c_j*u(i, n, b) &
                        + c_p1*u(i, 1, b) + c_p2*u(i, 2, b) &
                        + c_p3*u(i, 3, b) + c_p4*u(i, 4, b) &
                        - s(i + (n - 2)*SZ)*thom_s(n)

    ! Backward pass of the Thomas algorithm
    s(i + (n - 1)*SZ) = s(i + (n - 1)*SZ)*thom_w(n)
    !temp_jp1 = s(i + (n - 1)*SZ)*thom_w(n)
    !s(i + (n - 1)*SZ) = temp_jp1
    do j = n - 1, 1, -1
      s(i + (j - 1)*SZ) = (s(i + (j - 1)*SZ) - thom_f(j)*s(i + (j + 1 - 1)*SZ))*thom_w(j)

      !s(i + (j - 1)*SZ) = (s(i + (j - 1)*SZ) - thom_f(j)*temp_jp1)*thom_w(j)
      !temp_jp1 = s(i + (j - 1)*SZ)

      !temp_jp1 = (s(i + (j - 1)*SZ) - thom_f(j)*temp_jp1)*thom_w(j)
      !s(i + (j - 1)*SZ) = temp_jp1
    end do

    ! Periodic final pass
    ss = (s(i) - alpha*s(i + (n - 1)*SZ)) &
         /(1._dp + thom_p(1) - alpha*thom_p(n))
    do j = 1, n
      du(i, j, b) = s(i + (j - 1)*SZ) - ss*thom_p(j)
    end do

     ! test bw
     !do j = 1, n
     !  du(i, j, b) = u(i, j, b)
     !end do

    ! Backward pass of the Thomas algorithm
!    du(i, n, b) = s(i + (n - 1)*SZ)*thom_w(n)
!    do j = n - 1, 1, -1
!      du(i, j, b) = (s(i + (j - 1)*SZ) - thom_f(j)*s(i + (j + 1 - 1)*SZ))*thom_w(j)
!    end do

  end subroutine der_univ_thom_shared_per

  attributes(global) subroutine copycuda(u_, u, n)
    implicit none

    real(dp), device, intent(out), dimension(:, :, :) :: u_
    real(dp), device, intent(in), dimension(:, :, :) :: u
    integer, value :: n

    integer :: i, j, b

    i = threadIdx%x
    b = blockIdx%x

    do j = 1, n
      u_(i, j, b) = u(i, j, b)
    end do

  end subroutine copycuda

end module m_cuda_kernels_thom
