!c===============================c
        subroutine equatp
!c===============================c
        use varalc
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
!c-------------------------------c
        j11=1;k11=1
        j21=j11+1;k21=k11+1
        m21=m11-1;n21=n11-1
!c-------------------------------c
        do j=j21,m21
            do k=k21,n21
                if((j.ge.jb1.and.j.lt.je1).and.(k.ge.kb1.and.k.lt.ke1)) cycle
                aiw3(j,k)=0.0d+00
                aie3(j,k)=0.0d+00
                ajs3(j,k)=dzm1(k)*dzm1(k)/apc1(j-1,k)
                ajn3(j,k)=dzm1(k)*dzm1(k)/apc1(j,k)
                akb3(j,k)=dym1(j)*dym1(j)/apc2(j,k-1)
                akt3(j,k)=dym1(j)*dym1(j)/apc2(j,k)
            enddo
         enddo
         
        do j=j21,m21
            do k=k21,n21                
                if(j.eq.m21.or.(j.eq.(jb1-1).and.(k.ge.kb1.and.k.lt.ke1))) ajn1(j,k)=0.0d+00
                if(k.eq.n21.or.(k.eq.(kb1-1).and.(j.ge.jb1.and.j.lt.je1))) akt1(j,k)=0.0d+00
                if(j.eq.j21.or.(j.eq.je1.and.(k.ge.kb1.and.k.lt.ke1))) ajs1(j,k)=0.0d+00
                if(k.eq.k21.or.(k.eq.ke1.and.(j.ge.jb1.and.j.lt.je1))) akb1(j,k)=0.0d+00
                
                apc3(j,k)=aiw3(j,k)+aie3(j,k)+ajs3(j,k)+ajn3(j,k)+akb3(j,k)+akt3(j,k)
        
                bpp3(j,k)= - (vv1(j,k) - vv1(j-1,k)) * dzm1(k) - (ww1(j,k) - ww1(j,k-1)) * dym1(j)
            enddo
        enddo


!        apc3(jr1,kr1)=1.0d+00
!        aiw3(jr1,kr1)=0.0d+00
!        aie3(jr1,kr1)=0.0d+00
!        ajs3(jr1,kr1)=0.0d+00
!        ajn3(jr1,kr1)=0.0d+00
!        akb3(jr1,kr1)=0.0d+00
!        akt3(jr1,kr1)=0.0d+00
!        bpp3(jr1,kr1)=0.0d+00
        
!        mode=2
        mode=1
        if(levgrd.eq.1) then
            call solve1(j11,k11,m11,n11,jb1,je1,kb1,ke1,ni1,nj1,nk1,mode,epsp11,epsp21,epsp31,epsp41,ressum0,ressum, &
                        pre2,res1,apc3,bpp3,aiw3,aie3,ajs3,ajn3,akb3,akt3)
        end if
        if(levgrd.eq.2) then
            call solve2(j11,k11,m11,n11,jb1,je1,kb1,ke1,ni1,nj1,nk1,mode,epsp11,epsp21,epsp31,epsp41,ressum0,ressum, &
                        pre2,res1,apc3,bpp3,aiw3,aie3,ajs3,ajn3,akb3,akt3)
        end if
        if(levgrd.eq.3) then
            call solve3(j11,k11,m11,n11,jb1,je1,kb1,ke1,ni1,nj1,nk1,mode,epsp11,epsp21,epsp31,epsp41,ressum0,ressum, &
                        pre2,res1,apc3,bpp3,aiw3,aie3,ajs3,ajn3,akb3,akt3)
        end if
        if(levgrd.eq.4) then
            call solve4(j11,k11,m11,n11,jb1,je1,kb1,ke1,ni1,nj1,nk1,mode,epsp11,epsp21,epsp31,epsp41,ressum0,ressum, &
                        pre2,res1,apc3,bpp3,aiw3,aie3,ajs3,ajn3,akb3,akt3)
        end if
        if(levgrd.eq.5) then
            call solve5(j11,k11,m11,n11,jb1,je1,kb1,ke1,ni1,nj1,nk1,mode,epsp11,epsp21,epsp31,epsp41,ressum0,ressum, &
                        pre2,res1,apc3,bpp3,aiw3,aie3,ajs3,ajn3,akb3,akt3)
        end if
        if(levgrd.eq.6) then
            call solve6(j11,k11,m11,n11,jb1,je1,kb1,ke1,ni1,nj1,nk1,mode,epsp11,epsp21,epsp31,epsp41,ressum0,ressum, &
                        pre2,res1,apc3,bpp3,aiw3,aie3,ajs3,ajn3,akb3,akt3)
        end if
        if(levgrd.eq.7) then
            call solve7(j11,k11,m11,n11,jb1,je1,kb1,ke1,ni1,nj1,nk1,mode,epsp11,epsp21,epsp31,epsp41,ressum0,ressum,&
                        pre2,res1,apc3,bpp3,aiw3,aie3,ajs3,ajn3,akb3,akt3)
        end if
        if(levgrd.eq.8) then
            call solve8(j11,k11,m11,n11,jb1,je1,kb1,ke1,ni1,nj1,nk1,mode,epsp11,epsp21,epsp31,epsp41,ressum0,ressum,&
                        pre2,res1,apc3,bpp3,aiw3,aie3,ajs3,ajn3,akb3,akt3)
        end if
!c-------pressure boundary value-------c

        do j=j21,m21
            pre2(j,k11)=pre2(j,k21)+dzw1(k21)/dzw1(k21+1)*(pre2(j,k21)-pre2(j,k21+1))
            pre2(j,n11)=pre2(j,n21)+dzw1(n21+1)/dzw1(n21)*(pre2(j,n21)-pre2(j,n21-1))
        enddo

        do k=k21,n21
            pre2(j11,k)=pre2(j21,k)+dyv1(j21)/dyv1(j21+1)*(pre2(j21,k)-pre2(j21+1,k))
            pre2(m11,k)=pre2(m21,k)+dyv1(m21+1)/dyv1(m21)*(pre2(m21,k)-pre2(m21-1,k))
        enddo

        prej=pre2(j21,k11)+dyv1(j21)/dyv1(j21+1)*(pre2(j21,k11)-pre2(j21+1,k11))
        prek=pre2(j11,k21)+dzw1(k21)/dzw1(k21+1)*(pre2(j11,k21)-pre2(j11,k21+1))
        pre2(j11,k11)=5.0d-01*(prej+prek)

        prej=pre2(j21,n11)+dyv1(j21)/dyv1(j21+1)*(pre2(j21,n11)-pre2(j21+1,n11))
        prek=pre2(j11,n21)+dzw1(n21+1)/dzw1(n21)*(pre2(j11,n21)-pre2(j11,n21-1))
        pre2(j11,n11)=5.0d-01*(prej+prek)

        prej=pre2(m21,k11)+dyv1(m21+1)/dyv1(m21)*(pre2(m21,k11)-pre2(m21-1,k11))
        prek=pre2(m11,k21)+dzw1(k21)/dzw1(k21+1)*(pre2(m11,k21)-pre2(m11,k21+1))
        pre2(m11,k11)=5.0d-01*(prej+prek)

        prej=pre2(m21,n11)+dyv1(m21+1)/dyv1(m21)*(pre2(m21,n11)-pre2(m21-1,n11))
        prek=pre2(m11,n21)+dzw1(n21+1)/dzw1(n21)*(pre2(m11,n21)-pre2(m11,n21-1))
        pre2(m11,n11)=5.0d-01*(prej+prek)

        do j=jb1+1,je1-2
            pre2(j,ke1-1)=pre2(j,ke1)+dzw1(ke1)/dzw1(ke1+1)*(pre2(j,ke1)-pre2(j,ke1+1))
            pre2(j,kb1)=pre2(j,kb1-1)+dzw1(kb1)/dzw1(kb1-1)*(pre2(j,kb1-1)-pre2(j,kb1-2))
        enddo

        do k=kb1+1,ke1-2
            pre2(je1-1,k)=pre2(je1,k)+dyv1(je1)/dyv1(je1+1)*(pre2(je1,k)-pre2(je1+1,k))
            pre2(jb1,k)=pre2(jb1-1,k)+dyv1(jb1)/dyv1(jb1-1)*(pre2(jb1-1,k)-pre2(jb1-2,k))
        enddo

        prej=pre2(je1,ke1-1)+dyv1(je1)/dyv1(je1+1)*(pre2(je1,ke1-1)-pre2(je1+1,ke1-1))
        prek=pre2(je1-1,ke1)+dzw1(ke1)/dzw1(ke1+1)*(pre2(je1-1,ke1)-pre2(je1-1,ke1+1))
        pre2(je1-1,ke1-1)=5.0d-01*(prej+prek)

        prej=pre2(je1,kb1)+dyv1(je1)/dyv1(je1+1)*(pre2(je1,kb1)-pre2(je1+1,kb1))
        prek=pre2(je1,kb1-1)+dzw1(kb1)/dzw1(kb1-1)*(pre2(je1,kb1-1)-pre2(je1,kb1-2))
        pre2(je1-1,kb1)=5.0d-01*(prej+prek)

        prej=pre2(jb1-1,ke1)+dyv1(jb1)/dyv1(jb1-1)*(pre2(jb1-1,ke1)-pre2(jb1-2,ke1))
        prek=pre2(jb1,ke1)+dzw1(ke1)/dzw1(ke1+1)*(pre2(jb1,ke1)-pre2(jb1,ke1+1))
        pre2(jb1,ke1-1)=5.0d-01*(prej+prek)

        prej=pre2(jb1-1,kb1)+dyv1(jb1)/dyv1(jb1-1)*(pre2(jb1-1,kb1)-pre2(jb1-2,kb1))
        prek=pre2(jb1,kb1-1)+dzw1(kb1)/dzw1(kb1-1)*(pre2(jb1,kb1-1)-pre2(jb1,kb1-2))
        pre2(jb1,kb1)=5.0d-01*(prej+prek)


!c-------correction of v by p-------c
        j11=2
        k11=1
        j21=j11+1
        k21=k11+1
        do j=j21,m21
            do k=k21,n21
                if((j.ge.jb1.and.j.le.je1).and.(k.ge.kb1.and.k.lt.ke1)) cycle
                vv2(j,k)=vv1(j,k)-(pre2(j+1,k)-pre2(j,k))*dzm1(k)/apc1(j,k)
            enddo
        enddo
!c-------correction of w by p-------c
        j11=1
        k11=2
        j21=j11+1
        k21=k11+1
        do j=j21,m21
            do k=k21,n21
                if((j.ge.jb1.and.j.lt.je1).and.(k.ge.kb1.and.k.le.ke1)) cycle
                ww2(j,k)=ww1(j,k)-(pre2(j,k+1)-pre2(j,k))*dym1(j)/apc2(j,k)
            enddo
        enddo

        j11=1
        k11=1
        j21=j11+1
        k21=k11+1
        do j=j11,m11
            do k=k11,n11
                if((j.ge.jb1.and.j.lt.je1).and.(k.ge.kb1.and.k.le.ke1)) cycle
                pre2(j,k)=pre1(j,k)+pre2(j,k)
            enddo
        enddo
!c-------continuity satisfaction level check-------c
        divsum=0.0d+00
        do j=j21,m21
            do k=k21,n21
                if((j.ge.jb1.and.j.lt.je1).and.(k.ge.kb1.and.k.lt.ke1)) cycle
                dvy=(vv1(j+1,k)-vv1(j,k))*dzm1(k)
                dwz=(ww1(j,k+1)-ww1(j,k))*dym1(j)
                divsum=divsum+dabs(dvy+dwz)
            enddo
        enddo
        
        write(*,601) divsum0,divsum
601     format(' *',8x,'divergence drop',1pe9.2,'/',1pe9.2,8x,'*')
        return
        end
