module quadrature
  implicit none
contains

  ! Compute 1D Gauss quadrature points and weights for 1, 2, or 3 points
  subroutine get_1D_integration_points(n, points, weights)
    ! Inputs:
    ! n       : Number of Gauss points (1, 2, or 3)
    ! Outputs:
    ! points  : Gauss quadrature points
    ! weights : Gauss quadrature weights
    integer, intent(in) :: n
    real(8), allocatable, intent(inout) :: points(:), weights(:)

    if (n == 1) then
      ! Single-point rule
      points(1) = 0.0d0
      weights(1) = 2.0d0
    elseif (n == 2) then
      ! Two-point rule
      points(1) = -0.5773502691896258d0
      points(2) =  0.5773502691896258d0
      weights(1) = 1.0d0
      weights(2) = 1.0d0
    elseif (n == 3) then
      ! Three-point rule
      points(1) = -0.7745966692414834d0
      points(2) =  0.0d0
      points(3) =  0.7745966692414834d0
      weights(1) = 0.5555555555555556d0
      weights(2) = 0.8888888888888888d0
      weights(3) = 0.5555555555555556d0
    else
      print *, "Error: Only 1, 2, or 3 points are supported."
      stop
    end if
  end subroutine get_1D_integration_points

  ! Compute 2D Gauss quadrature points and weights using tensor product
  subroutine get_2D_integration_points(n, points_2D, weights_2D)
    ! Inputs:
    ! n         : Number of Gauss points in each dimension (1, 2, or 3)
    ! Outputs:
    ! points_2D : Gauss quadrature points (x, y) in 2D
    ! weights_2D: Gauss quadrature weights in 2D
    integer, intent(in) :: n
    real(8), intent(inout) :: points_2D(:, :), weights_2D(:)
    real(8), allocatable :: points_1D(:), weights_1D(:)
    integer :: i, j, index

    ! Get 1D points and weights
    allocate(points_1D(n), weights_1D(n))
    call get_1D_integration_points(n, points_1D, weights_1D)
    !allocate(points_2D(n*n, 2), weights_2D(n*n))

    ! Tensor product for 2D quadrature
    index = 1
    do i = 1, n
      do j = 1, n
        points_2D(index, 1) = points_1D(i)  ! x-coordinate
        points_2D(index, 2) = points_1D(j)  ! y-coordinate
        weights_2D(index) = weights_1D(i) * weights_1D(j)
        index = index + 1
      end do
    end do

    deallocate(points_1D, weights_1D)

  end subroutine get_2D_integration_points

  
    ! Compute 3D Gauss quadrature points and weights
  subroutine get_3D_integration_points(n, points_3D, weights_3D)
    integer, intent(in) :: n
    real(8), intent(inout) :: points_3D(:, :), weights_3D(:)
    real(8), allocatable :: points_1D(:), weights_1D(:)
    integer :: i, j, k, index

    allocate(points_1D(n), weights_1D(n))
    call get_1D_integration_points(n, points_1D, weights_1D)
    !allocate(points_3D(n**3, 3), weights_3D(n**3))

    index = 1
    do i = 1, n
      do j = 1, n
        do k = 1, n
          points_3D(index, 1) = points_1D(i)  ! x-coordinate
          points_3D(index, 2) = points_1D(j)  ! y-coordinate
          points_3D(index, 3) = points_1D(k)  ! z-coordinate
          weights_3D(index) = weights_1D(i) * weights_1D(j) * weights_1D(k)
          index = index + 1
        end do
      end do
    end do
    
    deallocate(points_1D, weights_1D)
  end subroutine get_3D_integration_points
  
    ! Function to get 3D integration points on a specific face of a hexahedral element
  subroutine get_face_integration_points(num_2D_points, face_id, points_2D, face_points_3D)
    integer, intent(in) :: num_2D_points                       ! Number of Gauss points per direction
    integer, intent(in) :: face_id                          ! Face ID (1 to 6)
    real(8), intent(in) :: points_2D(num_2D_points, 2)         ! input 2D integration points  
    real(8), intent(inout) :: face_points_3D(num_2D_points, 3) ! Output 3D points
    integer :: i

    ! Compute 2D Gauss quadrature points and weights
    ! num_2D_points = num_points * num_points

    ! Map 2D points to 3D natural coordinates based on the face
    select case (face_id)
    case (4) ! Face: xi = -1
      do i = 1, num_2D_points
        face_points_3D(i, :) = [-1.0d0, points_2D(i, 1), points_2D(i, 2)]
      end do
    case (3) ! Face: xi = +1
      do i = 1, num_2D_points
        face_points_3D(i, :) = [1.0d0, points_2D(i, 1), points_2D(i, 2)]
      end do
    case (5) ! Face: eta = -1
      do i = 1, num_2D_points
        face_points_3D(i, :) = [points_2D(i, 1), -1.0d0, points_2D(i, 2)]
      end do
    case (6) ! Face: eta = +1
      do i = 1, num_2D_points
        face_points_3D(i, :) = [points_2D(i, 1), 1.0d0, points_2D(i, 2)]
      end do
    case (1) ! Face: zeta = -1
      do i = 1, num_2D_points
        face_points_3D(i, :) = [points_2D(i, 1), points_2D(i, 2), -1.0d0]
      end do
    case (2) ! Face: zeta = +1
      do i = 1, num_2D_points
        face_points_3D(i, :) = [points_2D(i, 1), points_2D(i, 2), 1.0d0]
      end do
    case default
      print *, "Invalid face ID. Must be between 1 and 6."
      stop
    end select

    ! deallocate(points_2D, weights_2D)
  end subroutine get_face_integration_points
  
  
  subroutine natural_coord(coord, numNodes, xbar, natCoords)
    ! Computes the natural coordinates (r, s, t) of a point within an element.
    ! Arguments:
    !    coord - 2D array (3x8) containing the nodal coordinates of the hexahedral element.
    !    xbar  - 1D array (3) containing the global coordinates of the point.
    ! Outputs:
    !    r, s, t - Natural (local) coordinates of the point xbar.
    
    implicit none

    ! Input arguments
    real(8), intent(in) :: coord(3, numNodes)
    real(8), intent(in) :: xbar(3)
    integer, intent(in) :: numNodes
    
    ! Output arguments
    real(8), dimension(3), intent(out) :: natCoords
    
    ! Local variables
    real(8) :: tol, error, norm_residual
    real(8) :: d_rst(3), residual(3), x(3)
    real(8) :: J(3, 3), invJ(3, 3), detJ
    integer :: iter
    real(8) :: shapeFunc(numNodes), dN_dXi(numNodes, 3)
    
    ! Constants
    tol = 1.0d-10
    error = 1.0d0
    ! nnode = 8
    
    ! Initial guess
    natCoords = 0.0d0
    iter = 0
    do while (error >= tol .and. iter <= 15)
        ! Compute x = shape_val3D * coord
        call solidShapeFunctions(natCoords, numNodes, shapeFunc, dN_dXi)
        
        x = matmul(coord, shapeFunc)
        !x(2) = matmul(coord(2, :), shapeFunc)
        !x(3) = matmul(coord(3, :), shapeFunc)
        
        ! Compute derivatives (Nr, Ns, Nt)
        ! call dN_dr_ds_dt(numNodes, rst(1), rst(2), rst(3), Nr, Ns, Nt)
        
        ! Compute residual
        residual = x - xbar
        
        ! Compute Jacobian
        J = matmul(transpose(dN_dXi), transpose(coord))
        !J(2, :) = matmul(dN_dXi(:, 2), coord)
        !J(3, :) = matmul(dN_dXi(:, 3), coord)
        
        call innve(J, invJ, detJ)
        
        if (detJ <= 0.0d0) then
            print *, "Singular system when evaluating natural coordinates."
            stop
        end if
        
        ! Solve for d_rst
        d_rst = -matmul(invJ, residual)
        
        ! Update rst
        natCoords = natCoords + d_rst
        
        ! Compute error as the norm of the residual
        iter = iter + 1
        error = sqrt(sum(d_rst**2))
    end do
    
    if (iter > 15) then
        print *, "Newton-Raphson did not converge when evaluating natural coordinates."
        ! stop
    end if
    
    ! set close to face to 1 or -1
    !if (abs(natCoords(3)) <= 0.999 ) natCoords(3) = sign(0.9999, natCoords(3))

end subroutine natural_coord

  
  
  
  
end module quadrature
