!c===============================c
        subroutine initiv
!c===============================c
        use varalc
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
        include 'table.lei'
!c-------------------------------c
        m21=m11-1
        n21=n11-1
!c-------------------------------c
        j11=1
        k11=1
        do j=j11,m11
            do k=k11,n11
                aiw1(j,k)=0.0d+00
                aie1(j,k)=0.0d+00
                ajs1(j,k)=0.0d+00
                ajn1(j,k)=0.0d+00
                akb1(j,k)=0.0d+00
                akt1(j,k)=0.0d+00
                apc1(j,k)=0.0d+00
                bpp1(j,k)=0.0d+00
                vpvp1(j,k)=vpvp(j,k)
                vpwp1(j,k)=vpwp(j,k)
                wpwp1(j,k)=wpwp(j,k)
            enddo
        enddo
        return
        end

!c===============================c
        subroutine initiw
!c===============================c
        use varalc
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
        include 'table.lei'
!c-------------------------------c
        m21=m11-1
        n21=n11-1
!c-------------------------------c
        j11=1
        k11=1
        do j=j11,m11
            do k=k11,n11
                aiw2(j,k)=0.0d+00
                aie2(j,k)=0.0d+00
                ajs2(j,k)=0.0d+00
                ajn2(j,k)=0.0d+00
                akb2(j,k)=0.0d+00
                akt2(j,k)=0.0d+00
                apc2(j,k)=0.0d+00
                bpp2(j,k)=0.0d+00
            enddo
        enddo
        return
        end

!c===============================c
        subroutine initip
!c===============================c
        use varalc
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
!c-------------------------------c
        j11=1
        k11=1
        do j=j11,m11
            do k=k11,n11
                aiw3(j,k)=0.0d+00
                aie3(j,k)=0.0d+00
                ajs3(j,k)=0.0d+00
                ajn3(j,k)=0.0d+00
                akb3(j,k)=0.0d+00
                akt3(j,k)=0.0d+00
                apc3(j,k)=0.0d+00
                bpp3(j,k)=0.0d+00
            enddo
        enddo
        return
        end
