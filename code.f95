PROGRAM main
    IMPLICIT NONE
    real, parameter::pi = 3.14159265359
    real*8 a(3,3), expo_mat(3,3)
    integer i
    
    do i=1,3
        read *, a(i, 1:3)
    end do
    call exponent_matrix(a,expo_mat,3)
    do i=1,3
        print *, expo_mat(i, 1:3)
    end do

    
END PROGRAM main

!roots of quadratic equation
subroutine roots(a,b,c)
    integer a,b,c
    real*8 root1, root2, discriminant

    discriminant = b**2 - (4*a*c)
    if(discriminant <0.0) then
        print *, "No real roots! "
        return
    endif
    root1 = (-b + sqrt(discriminant))/(2.0 * a)
    root2 = (-b - sqrt(discriminant))/(2.0 * a)
    print *, "Root 1: ", root1
    print *, "Root 2: ", root2
end subroutine roots

!polar coordinate to cartesian coordinate
subroutine polartocart(r, theta,x,y)
    real*8 x,y,r,theta
    x = r * cos(theta)
    y = r * sin(theta)
    print *, x," ",y
end subroutine polartocart

!cartesian to polar
subroutine carttopolar(x, y)
    real*8 x,y,theta,r,theta_in_degree,z
    r = sqrt(x**2 + y**2)
    z = x/r
    theta = acos(z)
    theta_in_degree = (theta * 180) / pi
    print *, "r = ", r," theta in radian = ", theta
end subroutine carttopolar

! matrix multiplication (standard lib -> MATMUL(matrix_1, matrix_2))
! mat_mul -> matrix_a(n,m) matrix_b(m,l) product(n,l)
subroutine mat_mul(matrix_a,matrix_b, product,n,m,l)
    real*8 matrix_a(n,m), matrix_b(m,l), product(n,l)
    integer i,j,k
    do i=1,n
        do j=1,l
            product(i,j) = 0
            do k=1,m
                product(i,j) = product(i,j) + matrix_a(i,k) * matrix_b(k,j)
            end do
        end do
    end do
end subroutine mat_mul

!exponent (standard lib -> exp())
subroutine exponent(x,ex,term)
    real*8 x,ex,term
    do i = 1,5000
        term = term * x/(i)
        ex = ex + term
        if(term < 10**(-8.0)) then
            exit
        endif
    end do
    
end subroutine exponent

! exponent_matrix
subroutine exponent_matrix(a,expo_mat,n)
    real*8 a(n,n), prod(n,n), expo_mat(n,n), term(n,n)
    do i=1,n
        do j=1,n
            if(i .eq. j) then
                expo_mat(i,j) =1
                term(i,j) = 1
            else
                expo_mat(i,j) = 0
                term(i,j) = 0
            end if
        end do
    end do
    
    do i=1,5000
        call mat_mul(term,a,prod,n,n,n)
        term = prod/i
        expo_mat = expo_mat + term
    end do
end subroutine exponent_matrix