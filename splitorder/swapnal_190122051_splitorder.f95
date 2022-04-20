!SWAPNAL SHAHIL 190122051
PROGRAM main
    IMPLICIT NONE
    real*8, parameter::pi = acos(-1.d0)
    !real*8, parameter::g = 0.51d0
    integer, parameter::outunit = 6
    integer, parameter::nx = 200
    integer, parameter::nb = 50
    integer i,k
    real*8 dx,dt,xst,xend,x(nx),l, dk,ki
    complex*16 iota,v(nx),psi(nx),f(nx),g(nx)
    OPEN(UNIT = outunit, FILE = 'output.txt')
    iota = (0.d0,1.d0)
    xst = -10.d0
    xend = 10.d0
    l = (xend-xst)  !length of the box = xend - xstart
    dx = l/nx
    dt = 0.1
    dk = 2.d0*pi/(nx*dx)

    
   
    call calc_xi(xst,dx,x,nx)
    call calc_psi(psi,x,nx)
    call calc_potential(v,x,nx)
    do i=1,nx
       v(i) = exp(-iota*v(i)*dt/2.d0) * psi(i)
    end do
    f = v
     call dft_inverse(f,g,nx,dk,dx)
        
        do k=1,nx
            ki = -pi/dx + k*dk
            g(i) = g(i)*exp((-iota*(iota*ki)**2)/2.d0)
        end do
        f =0
        call dft(g, f, nx, dk, dx)
   
    do i=1,nx
       f(i) = exp(-iota*v(i)*dt/2.d0) * f(i)
    end do

    do i=1,nx
       write(6,*)x(i), real(f(i))
    end do
    

    CLOSE(outunit)
END PROGRAM main




subroutine calc_xi(xst,dx,x,nx)
  integer i,nx
  real*8 xst,dx,x(nx)
  do i=1,nx
     x(i) = xst + (i)*dx
  end do
end subroutine calc_xi


subroutine calc_potential(v,x,nx)
  integer i,nx
  real*8 x(nx)
  complex*16 v(nx)
  v = 0.d0
  do i=1,nx
     v(i) = v(i) + (1/2.d0)*x(i)**2
  end do
    
end subroutine calc_potential


subroutine calc_psi(psi,x,nx)
  integer i,nx
  real*8 x(nx)
  complex*16 psi(nx)
  do i=1,nx
        psi(i) = exp(-(x(i)**2)/2.d0)
  end do
end subroutine calc_psi































































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


subroutine complex_mat_mul(matrix_a,matrix_b, prod,n,m,l)
    complex*16 matrix_a(n,m), matrix_b(m,l), prod(n,l)
    integer i,j,k
    do i=1,n
        do j=1,l
            prod(i,j) = 0
            do k=1,m
                prod(i,j) = prod(i,j) + matrix_a(i,k) * matrix_b(k,j)
            end do
        end do
    end do
end subroutine complex_mat_mul
subroutine complex_exponent_matrix(a,expo_mat,n)
    complex*16 a(n,n), prod(n,n), expo_mat(n,n), term(n,n)
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
        call complex_mat_mul(term,a,prod,n,n,n)
        term = prod/i
        expo_mat = expo_mat + term
    end do
end subroutine complex_exponent_matrix

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
