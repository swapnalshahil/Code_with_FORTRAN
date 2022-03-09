!SWAPNAL SHAHIL 190122051
PROGRAM main
    IMPLICIT NONE
    real*8, parameter::pi = acos(-1.d0)
    integer, parameter::outunit = 6
    integer n,i 
    real*8 dk,dx,xi,a ,k
    complex*16 f(200), g(200) , func(200)
    !OPEN(UNIT = outunit, FILE = 'output.txt')
    n = 200
    dx = (10.0 + 10.0)/n
    dk = 2.0*pi/(n*dx)
    a = 1.51
    do i=1,n 
        xi = -10 + i*dx 
        f(i) = exp(-xi**2)
    end do
     call dft(f,g,n,dk,dx)
    do i=1,n 
        xi = -10 + i*dx 
        k = -pi/dx + i*dk
        g(i) = g(i)
        !print *,xi, realpart(f(i)), realpart(g(i))
    end do

    
    func = f
    f = g
    g=0
    call dft_inverse(f,g,n,dk,dx)
    do i=1,n 
        xi = -10 + i*dx 
        print *,xi, realpart(func(i)),realpart(f(i)), realpart(g(i))
    end do

    
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

! Range Kutta methods

subroutine fn_for_rk(temp, f, n,dt,t)
    integer n
    real*8 temp(n), f(n), dt, t,g,m 
    g = 6.67408E-11
    m = 5.972E24
    f(1) = temp(3)
    f(2) = temp(4)
    f(3) =  -g*m* 1.0*temp(1) / (sqrt(temp(1)**2 + temp(2)**2) * (temp(1)**2 + temp(2)**2))
    f(4) =  -g*m*1.0*temp(1) / (sqrt(temp(1)**2 + temp(2)**2) * (temp(1)**2 + temp(2)**2))
    dt = dt
    t = t
end subroutine fn_for_rk

subroutine rk1(y, n, dt, t)
    integer n
    real*8 y(n),temp(n), dt, f(n),t
    temp = y
    call fn_for_rk(temp, f, n, dt,t)
    y = y + dt*f
end subroutine rk1

subroutine rk4(y, n, dt, t)
    real*8 y(n), k1(n), k2(n), k3(n), k4(n), dt , temp(n),f(n),t
    integer n
    temp = y
    call fn_for_rk(temp, f, n, dt, t)
    k1 = f*dt
    temp = y + k1/2.0
    t = t + dt/2.0              !t = t + dt/2.0
    call fn_for_rk(temp, f, n, dt, t)
    k2 = f*dt
    temp = y + k2/2.0
    call fn_for_rk(temp, f, n, dt,t)
    k3 = f*dt
    temp = y + k3
    t = t + dt/2.0              !t = t + dt/2.0 (t = t + h)
    call fn_for_rk(temp, f, n, dt,t)
    k4 = f*dt
    y = y + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0
    
end subroutine rk4

subroutine rk2(y, n, dt, t)
    real*8 y(n), k1(n), k2(n), dt, temp(n),f(n),t
    integer n
    temp = y
    call fn_for_rk(temp, f, n, dt, t)
    k1 = f*dt
    temp = y + k1/2.0
    t = t + dt/2.0              !t = t + dt/2.0
    call fn_for_rk(temp, f, n, dt,t)
    k2 = f*dt
    y = y + k2
end subroutine rk2

subroutine dft(f, g, n, dk, dx)
    complex*16 f(n), g(n),iota
    integer k, n
    real*8 dk, dx , xi,ki ,pi,a
    a = 1.51
    pi = acos(-1.d0)
    iota = (0.d0,1.d0)
    do k=1,n
        ki = -pi/dx + k*dk
        g(k) = 0.0
        do i=1,n
            xi = -10.0 + i*dx
            g(k) = g(k) + f(i)*exp(-iota*ki*xi)*dx
        end do
    end do
    g = g * (1.0/sqrt(2.d0*pi))
end subroutine dft

subroutine dft_inverse(f,g,n,dk,dx)
    complex*16 f(n), g(n),iota
    integer k, n
    real*8 dk, dx , xi,ki ,pi
    pi = acos(-1.d0)
    iota = (0.d0,1.d0)
    do i=1,n
        xi = -10.0 + i*dx
        g(i) = 0.0
        do k=1,n
            ki = -pi/dx + k*dk
            g(i) = g(i) + f(k)*exp(iota*ki*xi)*dk
        end do
    end do
    g = g * (1.0/sqrt(2.d0*pi))

end subroutine dft_inverse

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

!sin of matrix
subroutine sin_of_matrix(a,n)
    real*8 a(n,n),pro(n,n),term(n,n),temp(n,n),prod(n,n)
    integer i,n,x
    x=1
    pro = a
    term = 0.0
    prod = 0.0

    call mat_mul(a,a,temp,n,n,n)
    term = a
    do i=2,5000,2
        call mat_mul(term,temp, prod,n,n,n)
        term = prod/(i * (i+1))
         if(x .eq. 1)then
            pro = pro - prod/((i * (i+1)))
            x =0
        else
            pro = pro + prod/(i * (i+1))
            x = 1
        end if
    end do
    
    do i=1,n
            print *, pro(i,1:n)
    end do
end subroutine sin_of_matrix

! subroutine for log(1+x) -> do change x while using this subroutine
subroutine log_of_matrix(a,n)
    real*8 a(n,n),pro(n,n),term(n,n),prod(n,n)
    integer i,n,x
    x = 1
    pro = a
    
    prod = 0.0
    term = a
    do i=2,5000
        call mat_mul(term,a, prod,n,n,n)
        term = prod
         if(x .eq. 1)then
            pro = pro - prod/(i)
            x =0
        else
            pro = pro + prod/(i)
            x = 1
        end if
    end do
    do i=1,n
            print *, pro(i,1:n)
    end do
end subroutine log_of_matrix