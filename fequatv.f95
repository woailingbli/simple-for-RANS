!c===============================c
        subroutine equatv
!c===============================c
        use varalc
        include 'table.prc'
        include 'table.cre'
        include 'table.v_gd11'
        include 'table.lei'

        dimension Fn(nj1,nk1),Fs(nj1,nk1),Ft(nj1,nk1),Fb(nj1,nk1),vpwpt(nj1,nk1),vpwpb(nj1,nk1)
        dimension Dn(nj1,nk1),Ds(nj1,nk1),Dt(nj1,nk1),Db(nj1,nk1),deltaF(nj1,nk1),deltaA(nj1,nk1)

!c-------------------------------c
        j11=2;k11=1
        j21=j11+1;k21=k11+1
        m21=m11-1;n21=n11-1
!c-------------------------------c
!c-------doing vyy,ajs,ajn-------c

        do j=j21,m21
            do k=k21,n21
                if((j.ge.jb1.and.j.lt.je1).and.(k.ge.kb1.and.k.lt.ke1)) cycle
                       
                        Fn(j,k) = dzm1(k)*(vi1(j+1,k)+vi1(j,k))/2
                        Fs(j,k) = dzm1(k)*(vi1(j-1,k)+vi1(j,k))/2
                        Ft(j,k) = dyv1(j)*(wi1(j,k+1)+wi1(j-1,k+1))/2
                        Fb(j,k) = dyv1(j)*(wi1(j,k)+wi1(j-1,k))/2
                        
                        Dn(j,k) = dzm1(k)/(reno*dym1(j))
                        Ds(j,k) = dzm1(k)/(reno*dym1(j-1))
                        Dt(j,k) = dyv1(j)/(reno*dzw1(k+1))
                        Db(j,k) = dyv1(j)/(reno*dzw1(k))
                        
                        ajs1(j,k) = Ds(j,k) + max(0.0d-01 , Fs(j,k))
                        ajn1(j,k) = Dn(j,k) + max(0.0d-01 , -Fn(j,k))
                        akb1(j,k) = Db(j,k) + max(0.0d-01 , Fb(j,k))
                        akt1(j,k) = Dt(j,k) + max(0.0d-01 , -Ft(j,k))
                        deltaF(j,k) = Fn(j,k) - Fs(j,k) + Ft(j,k) - Fb(j,k)
                      
                        apc1(j,k) = ajs1(j,k) + ajn1(j,k) + akb1(j,k) + akt1(j,k)
                        
                        deltaA(j,k) = (ajs1(j,k)+ajn1(j,k)+akb1(j,k)+akt1(j,k))-apc1(j,k)
                        vpwpt(j,k) = (vpwp(j-1,k+1)+vpwp(j,k+1))*fz21(k+1)*0.5 + (vpwp(j-1,k)+vpwp(j,k))*fz11(k+1)*0.5
                        vpwpb(j,k) = (vpwp(j-1,k)+vpwp(j,k))*fz21(k)*0.5 + (vpwp(j-1,k-1)+vpwp(j,k-1))*fz11(k)*0.5
                        aiw1(j,k)=0
                        aie1(j,k)=0
                        bpp1(j,k) = - (pre1(j,k)-pre1(j-1,k)) * dzm1(k)  
                        bpp1(j,k) = bpp1(j,k) - (vpvp(j,k) - vpvp(j-1,k)) * dzm1(k) - (vpwpt(j,k) - vpwpb(j,k)) * dyv1(j)
            enddo
        enddo
!!!c-------------------------------c
!!!c********************************c
        call bondv2(j11,k11,m11,n11,jb1,je1,kb1,ke1,ni1,nj1,nk1,vv1,apc1,bpp1,aiw1,ajs1,akb1,aie1,ajn1,akt1)
!!!c********************************c
        mode=1

        if(levgrd.eq.1) then
            call v_solve1(j11,k11,m11,n11,jb1,je1,kb1,ke1,ni1,nj1,nk1,mode,epsp11,epsp21,epsp31,epsp41,ressum0,ressum, &
                        vv1,res1,apc1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1)
        end if
        if(levgrd.eq.2) then
            call v_solve2(j11,k11,m11,n11,jb1,je1,kb1,ke1,ni1,nj1,nk1,mode,epsp11,epsp21,epsp31,epsp41,ressum0,ressum, &
                        vv1,res1,apc1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1)
        end if
        if(levgrd.eq.3) then
            call v_solve3(j11,k11,m11,n11,jb1,je1,kb1,ke1,ni1,nj1,nk1,mode,epsp11,epsp21,epsp31,epsp41,ressum0,ressum, &
                        vv1,res1,apc1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1)
        end if
        if(levgrd.eq.4) then
            call v_solve4(j11,k11,m11,n11,jb1,je1,kb1,ke1,ni1,nj1,nk1,mode,epsp11,epsp21,epsp31,epsp41,ressum0,ressum, &
                        vv1,res1,apc1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1)
        end if
        if(levgrd.eq.5) then
            call v_solve5(j11,k11,m11,n11,jb1,je1,kb1,ke1,ni1,nj1,nk1,mode,epsp11,epsp21,epsp31,epsp41,ressum0,ressum, &
                        vv1,res1,apc1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1)
        end if
        if(levgrd.eq.6) then
            call v_solve6(j11,k11,m11,n11,jb1,je1,kb1,ke1,ni1,nj1,nk1,mode,epsp11,epsp21,epsp31,epsp41,ressum0,ressum, &
                        vv1,res1,apc1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1)
        end if
        if(levgrd.eq.7) then
            call v_solve7(j11,k11,m11,n11,jb1,je1,kb1,ke1,ni1,nj1,nk1,mode,epsp11,epsp21,epsp31,epsp41,ressum0,ressum,&
                        vv1,res1,apc1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1)
        end if
        if(levgrd.eq.8) then
            call v_solve8(j11,k11,m11,n11,jb1,je1,kb1,ke1,ni1,nj1,nk1,mode,epsp11,epsp21,epsp31,epsp41,ressum0,ressum,&
                        vv1,res1,apc1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1)
        end if
        if(levgrd.eq.9) then
            call v_solve9(j11,k11,m11,n11,jb1,je1,kb1,ke1,ni1,nj1,nk1,mode,epsp11,epsp21,epsp31,epsp41,ressum0,ressum,&
                        vv1,res1,apc1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1)
        end if

1000    return
!        call residum(j21,k21,m21,n21,jb1,je1+1,kb1,ke1,ni1,nj1,nk1,ressum0,vv1,apc1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1,res1)
!        ressum1=ressum0
!        write(*,200) ressum0
!        iter=0
!        mode=2
!!c--------------------------------c
!100     iter=iter+1
!!c--------------------------------c
!
!!c--------------------------------c
!        do k=k21,n21
!            if(k.ge.kb1.and.k.lt.ke1) then
!                call trdgmj(j11,jb1,k,nj1,nk1,vv1,apc1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1)
!                call trdgmj(je1,m11,k,nj1,nk1,vv1,apc1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1)
!            else
!                call trdgmj(j11,m11,k,nj1,nk1,vv1,apc1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1)
!            end if
!        enddo
!!c--------------------------------c
!        do j=j21,m21
!            if(j.ge.jb1.and.j.le.je1) then
!                call trdgmk(k11,kb1,j,nj1,nk1,vv1,apc1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1)
!                call trdgmk(ke1-1,n11,j,nj1,nk1,vv1,apc1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1)
!            else
!                call trdgmk(k11,n11,j,nj1,nk1,vv1,apc1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1)
!            end if
!        enddo
!!c--------------------------------c
!154     call residum(j21,k21,m21,n21,jb1,je1+1,kb1,ke1,ni1,nj1,nk1,ressum,vv1,apc1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1,res1)
!        write(*,201) ressum
!
!        if(ressum.ge.ressum1) then
!            write(*,202)
!            goto 1000
!        end if
!
!!c mode=1: control residual precision
!        if(mode.eq.1) then
!            if(ressum.le.epsm21) then
!                write(*,203)
!                goto 1000
!            else
!                ressum1=ressum
!                write(*,204)
!                goto 100
!            end if
!        end if
!
!!c mode=2: control residual decreasing level
!        if(mode.eq.2) then
!            if(ressum.le.epsm21) then
!                write(*,203)
!                goto 1000
!            end if
!            declev=ressum/ressum0
!            if(declev.gt.epsm31) then
!                ressum1=ressum
!                write(*,205)
!                goto 100
!            else
!                write(*,206)
!                goto 1000
!            end if
!        end if
!
!!c mode=3: control residual decreasing rate
!        if(mode.eq.3) then
!            if(ressum.le.epsm21) then
!                write(*,203)
!                goto 1000
!            end if
!            decrat=ressum/ressum1
!            if(decrat.le.epsm41) then
!                ressum1=ressum
!                write(*,207)
!                goto 100
!            else
!                write(*,208)
!                goto 1000
!            end if
!        end if
!
!200     format(' *',8x,'initial residual ressum=',1pe10.3,8x,'*')
!201     format(' *',8x,'total  residual  ressum=',1pe10.3,8x,'*')
!202     format(' *        residual not improving,     return        *')
!203     format(' *        total residual satisfied,   return        *')
!204     format(' *        residual not satisfied, continuing        *')
!205     format(' *        declev not satisfied,   continuing        *')
!206     format(' *        declev satisfied,  go out of L-B-L        *')
!207     format(' *        decrat small, go on next iteration        *')
!208     format(' *        decrat large, stop the LBL process        *')
!
!1000    return
        end
