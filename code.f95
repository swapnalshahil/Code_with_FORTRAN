!SWAPNAL SHAHIL 190122051
PROGRAM main
    IMPLICIT NONE
    real, parameter::pi = acos(-1.0)
    real*8 x0,y0,h,xfinal
    x0 = 0.0
    y0 = 0.0
    h = 0.25
    xfinal = 1.0
    !call rk1(x0,y0,h,xfinal)

    !call rk2(x0,y0,h,xfinal)

    call rk4(x0,y0,h,xfinal)
    
    !OPEN(UNIT = outunit, FILE = 'input.txt', FORM = "formatted" )
    !CLOSE(outunit)
    
    
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
subroutine mat_mul(matrix_a,matrix_b, prod,n,m,l)
    real*8 matrix_a(n,m), matrix_b(m,l), prod(n,l)
    integer i,j,k
    do i=1,n
        do j=1,l
            prod(i,j) = 0
            do k=1,m
                prod(i,j) = prod(i,j) + matrix_a(i,k) * matrix_b(k,j)
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

! derivative
! 1.2 : fdx(i) = ( -f(i-1)* 1/2 + f(i+1)*1/2)/dx**2
! 1.4 : fdx(i) = (f(i-2)*1/12 - f(i-1)*2/3 + f(i+1)*2/3 - f(i+2)*1/12)/dx**2
! 1.6 : fdx(i) = (-f(i-3)*1/60 + f(i-2)*3/20 - f(i-1)*3/4 + f(i+1)*3/4 - f(i+2)*3/20 + f(i+3)*1/60)/dx**2
! 1.8 : fdx(i) = (f(i-4)*1/280 -f(i-3)*4/105 + f(i-2)*1/5 - f(i-1)*4/5 + f(i+1)*4/5 - f(i+2)*1/5 + f(i+3)*4/105 - f(i+4)*1/280)/dx**2
! 2.2 : fdx(i) = (f(i-1)*1 - f(i)*2 + f(i+1))/dx**2
! 2.4 : fdx(i) = (-f(i-2)*1/12 + f(i-1)*4/3 - f(i)*5/2 + f(i+1)*4/3 - f(i+2)*1/12)/dx**2
! 2.6 : fdx(i) = (f(i-3)* 1/90 - f(i-2) * 3/20 + f(i-1)*3/2 - f(i) * 49/18 + f(i+1)*3/2 - f(i+2)*3/20 + f(i+3)*1/90)/dx**2
! 2.8 : fdx(i) = (-f(i-4)*1/560 + f(i-3)*8/315 - f(i-2)*1/5 + f(i-1)*8/5 - f(i)*205/72 + f(i+1)*8/5 - f(i+2)*1/5 + f(i+3)*8/315 - f(i+4)*1/560)/dx**2
! 3.2 : fdx(i) = (-f(i-2) * 1/2 + f(i-1) - f(i+1) + f(i+2)*1/2)/dx**2
! 3.4 : fdx(i) = (f(i-3)*1/8 - f(i-2)*1 + f(i-1) * 13/8 - f(i+1) * 13/8 + f(i+2) - f(i+3) *1/8)/dx**2
! 3.6 : fdx(i) = (-f(i-4)* 7/240 + f(i-3)*3/10 - f(i-2)*169/120 + f(i-1)*61/30 - f(i+1)*61/30 + f(i+2)*169/120 - f(i+3)*3/10 + f(i+4)*7/240)/dx**2;

subroutine derivative(f,dx,fdx,n)
    integer n
    real*8 f(200), dx, fdx(200)
    real*8 s
    fdx =0
    do i=4,n-3
        fdx(i) = (f(i-3)* 1/90 - f(i-2) * 3/20 + f(i-1)*3/2 - f(i) * 49/18 + f(i+1)*3/2 - f(i+2)*3/20 + f(i+3)*1/90)/dx**2
    end do
    s =0
    do i=1,200
        s = s+f(i)*fdx(i)*dx
    end do
    print *,s 
end subroutine derivative

!Range Kutta methods
!rk1(x0,y0,h,xfinal) xfinal is the point at which we need to find out.

real(kind = 8) function func(x0,y0)
    real*8 x0,y0
    !change function here
    func = x0 +2*y0
end function func


subroutine rk1(x0,y0,h,xfinal)
    real(kind=8),external :: func
    real*8 x0,y0,h,xfinal,x1,y1
    integer i, n
    n = int((xfinal-x0)/h + 0.5)
    do i=1,n
        x1 = x0 + h
        y1 = y0 + h * func(x0,y0)
        
        x0 = x1
        y0 = y1
        print *, x1," ",y1
    end do
  
end subroutine rk1

subroutine rk2(x0,y0,h,xfinal)
    real(kind=8),external :: func
    real*8 x0,y0,h,xfinal,k1,k2
    integer i, n
    n = int((xfinal-x0)/h + 0.5)
    do i=1,n
        k1 = h * func(x0,y0)
        k2 = h * func(x0+h,y0 +k1)
        y0 = y0 + (k1 + k2)/2.0
        x0 = x0 + h
        print *, x0," ",y0
    end do
    
end subroutine rk2

subroutine rk4(x0,y0,h,xfinal)
    real(kind=8),external :: func
    real*8 x0,y0,h,xfinal,k1,k2,k3,k4
    integer i, n
    n = int((xfinal-x0)/h + 0.5)
    do i=1,n
        k1 = h * func(x0,y0)
        k2 = h * func(x0+h/2.0,y0 + k1/2.0)
        k3 = h * func(x0+h/2.0,y0 + k2/2.0)
        k4 = h * func(x0+h,y0 + k3)
        y0 = y0 + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0
        x0 = x0 + h
        print *, x0," ",y0
    end do
    

end subroutine rk4

!cosine of matrix
subroutine cos_subroutine(a,n)
    real*8 a(2,2),pro(2,2),term(2,2),temp(2,2),prod(2,2)
    integer i,n,x
    n=2
    x =1
    do i=1,n
        do j=1,n
            if(i .eq. j) then
                pro(i,j) =1
                term(i,j) = 1
                prod(i,j)=1
            else
                pro(i,j) = 0
                term(i,j) = 0
                prod(i,j) = 0
            end if
        end do
    end do

    call mat_mul(a,a,temp,2,2,2)
    
    
    do i=1,5000,2
        call mat_mul(term,temp, prod,2,2,2)
        term = prod/(i * (i+1))
        if(x .eq. 1)then
            pro = pro - prod/(i * (i+1))
            x =0
        else
            pro = pro + prod/(i * (i+1))
            x = 1
        end if
    end do
    do i=1,2
            print *, pro(i,1:2)
    end do

end subroutine cos_subroutine
