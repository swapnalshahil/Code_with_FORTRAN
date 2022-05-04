! 
! 
!                             Online Fortran Compiler.
!                 Code, Compile, Run and Debug Fortran program online.
! Write your code in this editor and press "Run" button to execute it.
! 
! 


subroutine e_A(A,n,ea)
    implicit none
    integer i,j,n
    real*8 A(n,n),ea(n,n),t(n,n),Id(n,n),t2(n,n),zi    
    real*8 dt 
    dt=0.1
    zi=(0.d0,1.d0)
    Id=0
    do i=1,n
        Id(i,i)=1
    enddo  

    t=Id     
    ea=Id
    do i=1,100
        t2 = matmul(t,A)/i
        ea=ea+t2
        t=t2
    enddo
end


program main
implicit none
integer, parameter::nb=50,nx=400,nt=200
real*8 dx,l,dt,pi,x,xend,xst,sum
real*8 V(nb,nb),A(nb,nb),C(nb,1),X_matrix(nb,nb),ea(nb,nb),Xt(1,1),vv(nx),vvv(nx),T(nb,1),Temp(nb,nb)
complex*16 zi
integer m,n,k,i

zi=(0.d0,1.d0)
pi = dacos(-1.d0)
xst=-10
xend=+10
sum=0
dx=(xend-xst)/nx
dt=0.1d0
l=(xend-xst)
call pot(vv,xst,dx,nb,nx)

do m=1,nb 
    do n=1,nb 
        V(m,n)=0
        do k=1,nx 
           x=xst+k*dx
           V(m,n)=V(m,n)+(2/l)*vv(k)*sin((m*pi*(x-xst))/l)*sin((n*pi*(x-xst))/l)*dx
        enddo
    enddo
enddo
!write(6,*) "V(nb,nb) , ", V

! Adding T to V i.e. V=V+T
do m=1,nb 
    V(m,m)=V(m,m)+0.5d0*m**2*pi**2/l**2 ! 
enddo 
!write(6,*) "H matrix ,", V(nb,nb)

! calculate A
do m=1,nb 
    do n=1,nb 
        A(m,n)=-dt*V(m,n)
    enddo
enddo
!write(6,*)"A matrix, ",A

! calculate C col matrix at t=0
do m=1,nb 
    C(m,1)=0
    do i=1,nx 
        x=xst+i*dx 
        C(m,1)=C(m,1)+sqrt(2/l)*sin((m*pi*(x-xst))/l)*(x+1)*exp(-1*x**2/2)*dx
     enddo
     sum=sum+C(m,1)**2
  enddo
C=C/sqrt(sum) 
 
 
!write(6,*)"C matrix, ",C
!stop
! calculate x matrix 
X_matrix=0
do m=1,nb
    do n=1,nb
        X_matrix(m,n)=0
        do k=1,nx
            x=xst+k*dx
            X_matrix(m,n)=X_matrix(m,n)+2/l*sin((m*pi*(x-xst))/l)*sin((n*pi*(x-xst))/l)*dx*vv(k)
        enddo
        ! print*,v(i,j)
    enddo
 enddo

do m=1,nb 
    X_matrix(m,m)=X_matrix(m,m)+0.5d0*m**2*pi**2/l**2 ! 
enddo

Temp=0
do m=1,nb
    do n=1,nb
        Temp(m,n)=0
        do k=1,nx
            x=xst+k*dx
            Temp(m,n)=Temp(m,n)+2/l*sin((m*pi*(x-xst))/l)*sin((n*pi*(x-xst))/l)*dx*vv(k)*vv(k)
        enddo
    enddo
 enddo
 
do m=1,nb 
    Temp(m,m)=Temp(m,m)+(0.25d0*m**4*pi**4/l**4) 
enddo

! calculate C(i) at t=ndt 
call e_A(A,nb,ea)

! write(6,*)"ea matrix, ",ea

do i=0,nt
   C=matmul(ea,C)
   sum=0
   do m=1,nb
      sum=sum+C(m,1)**2
     enddo
   C=C/sqrt(sum)
   Xt = matmul(transpose(C),matmul(X_matrix,C))
    if(i.eq.(0)) then
    write(6,*)i*dt, "Xt matrix, ",Xt
 endif
    if(i.eq.(20)) then
    write(6,*)i*dt, "Xt matrix, ",Xt
 endif
    if(i.eq.(50)) then
    write(6,*)i*dt, "Xt matrix, ",Xt
 endif
    if(i.eq.(100)) then
    write(6,*)i*dt, "Xt matrix, ",Xt
 endif
    if(i.eq.(150)) then
    write(6,*)i*dt, "Xt matrix, ",Xt
 endif
    if(i.eq.(200)) then
    write(6,*)i*dt, "Xt matrix, ",Xt
 endif
 Xt = matmul(transpose(C),matmul(Temp,C))
    if(i.eq.(0)) then
    write(6,*)i*dt, "Xt_2 matrix, ",Xt
 endif
    if(i.eq.(20)) then
    write(6,*)i*dt, "Xt_2 matrix, ",Xt
 endif
    if(i.eq.(50)) then
    write(6,*)i*dt, "Xt_2 matrix, ",Xt
 endif
    if(i.eq.(100)) then
    write(6,*)i*dt, "Xt_2 matrix, ",Xt
 endif
    if(i.eq.(150)) then
    write(6,*)i*dt, "Xt_2 matrix, ",Xt
 endif
    if(i.eq.(200)) then
    write(6,*)i*dt, "Xt_2 matrix, ",Xt
 endif
 enddo

 
end program main


subroutine pot(vv,xst,dx,nb,nx)

  integer nb,nx ,i,k
  real*8 vv(nx),zi
  real*8 xst,dx,g
  g=0.29d0
  zi=(0.d0,1.d0)
  do k=1,nx 
     x=xst+k*dx
     vv(k)= ((1+g)/2)*x**2
    enddo
  end subroutine pot



