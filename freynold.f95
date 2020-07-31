!c**********************c
        subroutine reno_stress
!c**********************c
        include 'table.prc'
        include 'table.cre'
        include 'table.lei'
        
        open(20,file='..\input\RANS-vw\ww_Re400_DNS.txt')
        read(20,*) wpwp
        close(20)

        open(21,file='..\input\RANS-vw\vw_Re400_DNS.txt')
        read(21,*) vpwp
        close(21)
        
        open(22,file='..\input\RANS-vw\vv_Re400_DNS.txt')
        read(22,*) vpvp
        close(22)

        end
