PROGRAM multiplyComplex
    IMPLICIT NONE
    COMPLEX:: x,y
    PRINT *, "Enter X (Complex Number)"
    READ *, x
    PRINT *, "Enter Y (Complex Number)"
    READ *, y
    PRINT *, x, " x ", y, "is equal to ", x*y
    STOP
    
END PROGRAM multiplyComplex