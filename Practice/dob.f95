PROGRAM dateofbirth
    IMPLICIT NONE
    CHARACTER(LEN = 20):: name, title
    INTEGER:: d,m,y
    PRINT *, "What is your name? (20 characters only)"
    READ *, name,title
    PRINT *, "Enter your birth date, month and year (separate by space)"
    READ *, d,m,y
    PRINT *, d,".",m,".",y," is the date of birth of ", name," ",title
    PRINT *, name," ",title
    STOP
    
END PROGRAM dateofbirth