PROGRAM fibonacci
    IMPLICIT NONE
    integer*8 n,f1,f2,i,temp
    read *, n
    f1=0
    f2=1
    if(n==0)then
        print *,f1
    endif
    if(n==1)then
        !print*,f1
        print *,f2
    endif
    if(n>1)then
        !print*,f1
        !print *,f2
        do i=2,n
            temp = f1 + f2
            f1 = f2
            f2 = temp
            
        enddo
        print *, temp
    endif
    
END PROGRAM fibonacci