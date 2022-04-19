!SWAPNAL SHAHIL 190122051
PROGRAM main
    IMPLICIT NONE
    real*8, parameter::pi = acos(-1.d0)
    integer, parameter::outunit = 6
    integer n,i ,nt
    real*8 y(2),dt,t,tmax,Et,k
    !complex*16 f(200), g(200) , func(200)
    OPEN(UNIT = outunit, FILE = 'output.txt')
    n = 2
    t = 0
    dt = 0.0001
    tmax = 60
    nt = int(tmax/dt)
    
    y = 0
    y(1) = 1.d0
    y(2) = 0
    do i =1,nt
        t = t + dt
        call rk4(y,n,dt,t)
        Et = (1.d0/2.d0)*y(1)**2 + (1.d0/2.d0)*y(2)**2
        write(outunit,*) t, Et
    end do
    
    
    CLOSE(outunit)
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



! Range Kutta methods

subroutine fn_for_rk(temp, f, n,dt,tt)
    integer n
    real*8 temp(n), f(n), dt, tt,w,k
    k = 1.d0
    w = 1.d0
    f(1) = temp(2)
    f(2) = -k*temp(1) + 0.01d0*cos(w*tt)
   
    dt = dt
    tt = tt
end subroutine fn_for_rk


subroutine rk4(y, n, dt, t)
    real*8 y(n), k1(n), k2(n), k3(n), k4(n), dt , temp(n),f(n),t,tt
    integer n
    temp = y
    tt = t
    call fn_for_rk(temp, f, n, dt, tt)
    k1 = f*dt
    temp = y + k1/2.0
    tt = tt + dt/2.0              !t = t + dt/2.0
    call fn_for_rk(temp, f, n, dt, tt)
    k2 = f*dt
    temp = y + k2/2.0
    call fn_for_rk(temp, f, n, dt,tt)
    k3 = f*dt
    temp = y + k3
    tt = tt + dt/2.0              !t = t + dt/2.0 (t = t + h)
    call fn_for_rk(temp, f, n, dt,tt)
    k4 = f*dt
    y = y + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0
    
end subroutine rk4



!cosine of matrix
subroutine cos_of_matrix(a,n)
    real*8 a(n,n),pro(n,n),term(n,n),temp(n,n),prod(n,n)
    integer i,n,x
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

    call mat_mul(a,a,temp,n,n,n)
    
    
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
    do i=1,n
            print *, pro(i,1:n)
    end do

end subroutine cos_of_matrix
