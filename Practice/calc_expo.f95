PROGRAM cal_expo
    IMPLICIT NONE
    real value,g,temp,final
    integer total
    READ(5,*)g
    value=1 +g/100
    final=1+value
    temp = value
    DO total= 2,10
        temp = temp*value/total
        final = final + temp
        write(6,*)final
    ENDDO

   
    
END PROGRAM cal_expo