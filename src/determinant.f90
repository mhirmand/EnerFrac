    subroutine calculate_determinant(mat, determinant)
    
    implicit none
    
    real(8), intent(in) :: mat(3, 3)
    real(8), intent(out) :: determinant
    real :: minor1, minor2, minor3

    ! Calculate the minors and determinant
    minor1 = mat(2, 2) * mat(3, 3) - mat(2, 3) * mat(3, 2)
    minor2 = mat(2, 1) * mat(3, 3) - mat(2, 3) * mat(3, 1)
    minor3 = mat(2, 1) * mat(3, 2) - mat(2, 2) * mat(3, 1)

    determinant = mat(1, 1) * minor1 - mat(1, 2) * minor2 + mat(1, 3) * minor3
    
    return
    end
    ! end subroutine calculate_determinant