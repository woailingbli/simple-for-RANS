!c===============================c
        subroutine equatw
!c===============================c
        use varalc
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
        include 'table.lei'
!c-------------------------------c
        j11=1;k11=2
        j21=j11+1;k21=k11+1
        m21=m11-1;n21=n11-1
!c-------------------------------c
!c-------doing wzz,akb,akt-------c
        do j=j21,m21
            do k=k21,n21
                if((j.ge.jb1.and.j.lt.je1).and.(k.ge.kb1.and.k.le.ke1)) cycle
                       beta_w_plus1 = max(0.0d-01 , (5.0d-01 * (wi1(j,k) + wi1(j,k-1))))!j,k,+
                       beta_w_minus1 = min(0.0d-01 , (5.0d-01 * (wi1(j,k) + wi1(j,k-1))))!j,k,-
                       beta_v_plus1 = max(0.0d-01 , (5.0d-01 * (vi1(j+1,k) + vi1(j+1,k-1))))!j+0.5,k+0.5,+
                       beta_v_minus1 = min(0.0d-01 , (5.0d-01 * (vi1(j+1,k) + vi1(j+1,k-1))))!j+0.5,k+0.5,-
                       
                       beta_w_plus2 = max(0.0d-01 , (5.0d-01 * (wi1(j,k+1) + wi1(j,k))))!j,k+1,+
                       beta_w_minus2 = min(0.0d-01 , (5.0d-01 * (wi1(j,k+1) + wi1(j,k))))!j,k+1,-
                       beta_v_plus2 = max(0.0d-01 , (5.0d-01 * (vi1(j,k) + vi1(j,k+1))))!j-0.5,k+0.5,+
                       beta_v_minus2 = min(0.0d-01 , (5.0d-01 * (vi1(j,k) + vi1(j,k+1))))!j-0.5,k+0.5,-
                        
                       ajs2(j,k)=dzw1(k) * (beta_v_plus2 + 1/(reno*dyv1(j)))
                       ajn2(j,k)=dzw1(k)* (-beta_v_minus1 + 1/(reno*dyv1(j+1)))!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                       akb2(j,k)=dym1(j) * (beta_w_plus1 + 1/(reno*dzm1(k-1)))
                       akt2(j,k)=dym1(j) * (-beta_w_minus2 + 1/(reno*dzm1(k)))
                        
                       apc2(j,k)=dzw1(k) * (beta_v_plus1 + 1/(reno*dyv1(j)) -beta_v_minus2 + 1/(reno*dyv1(j+1)))!!!!!!!!!!!!!!!!!
                       apc2(j,k)=apc2(j,k)+dym1(j) * (beta_w_plus2 + 1/(reno*dzm1(k)) -beta_w_minus1 + 1/(reno*dzm1(k-1)))
                       
                       aiw2(j,k)=0
                       aie2(j,k)=0
                       
                       vpwpn = (vpwp(j+1,k)+vpwp(j+1,k-1))*fy21(j+1)*0.5 + (vpwp(j,k)+vpwp(j,k-1))*fy11(j+1)*0.5!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                       vpwps = (vpwp(j,k)+vpwp(j,k-1))*fy21(j)*0.5 + (vpwp(j-1,k)+vpwp(j-1,k-1))*fy11(j)*0.5
                       bpp2(j,k) = - (pre1(j,k)-pre1(j,k-1)) * dym1(j)  
                       bpp2(j,k) = bpp2(j,k) - (vpwpn - vpwps) * dzw1(k) - (wpwp(j,k) - wpwp(j,k-1)) * dym1(j)
            enddo
        enddo

!c********************************c
        call bondw2(j11,k11,m11,n11,jb1,je1,kb1,ke1,ni1,nj1,nk1,ww1,bpp2,apc2,aiw2,ajs2,akb2,aie2,ajn2,akt2)
!c********************************c
        mode=1
        call w_solve8(j11,k11,m11,n11,jb1,je1,kb1,ke1,ni1,nj1,nk1,mode,epsp11,epsp21,epsp31,epsp41,ressum0,ressum,&
                        vv1,res1,apc2,bpp2,aiw2,aie2,ajs2,ajn2,akb2,akt2)

!!        call residum(j21,k21,m21,n21,jb1,je1,kb1,ke1+1,ni1,nj1,nk1,ressum0,ww1,apc2,bpp2,aiw2,ajs2,akb2,aie2,ajn2,akt2,res2)
!        ressum1=ressum0
!        write(*,200) ressum0
!        iter=0
!        mode=2
!!c--------------------------------c
!100     iter=iter+1
!!c--------------------------------c
!
!        do k=k21,n21
!            if(k.ge.kb1.and.k.le.ke1) then
!                call trdgmj(j11,jb1,k,nj1,nk1,ww1,apc2,bpp2,aiw2,ajs2,akb2,aie2,ajn2,akt2)
!                call trdgmj(je1-1,m11,k,nj1,nk1,ww1,apc2,bpp2,aiw2,ajs2,akb2,aie2,ajn2,akt2)
!            else
!                call trdgmj(j11,m11,k,nj1,nk1,ww1,apc2,bpp2,aiw2,ajs2,akb2,aie2,ajn2,akt2)
!            end if
!        enddo
!!c--------------------------------c
!        do j=j21,m21
!            if(j.ge.jb1.and.j.lt.je1) then
!                call trdgmk(k11,kb1,j,nj1,nk1,ww1,apc2,bpp2,aiw2,ajs2,akb2,aie2,ajn2,akt2)
!                call trdgmk(ke1,n11,j,nj1,nk1,ww1,apc2,bpp2,aiw2,ajs2,akb2,aie2,ajn2,akt2)
!            else
!                call trdgmk(k11,n11,j,nj1,nk1,ww1,apc2,bpp2,aiw2,ajs2,akb2,aie2,ajn2,akt2)
!            end if
!        enddo
!!c--------------------------------c
!154     call residum(j21,k21,m21,n21,jb1,je1,kb1,ke1+1,ni1,nj1,nk1,ressum,ww1,apc2,bpp2,aiw2,ajs2,akb2,aie2,ajn2,akt2,res2)
!        write(*,201) ressum
!
!        if(ressum.ge.ressum1) then
!        write(*,202)
!        goto 1000
!        end if
!
!!c mode=1: control residual precision
!        if(mode.eq.1) then
!        if(ressum.le.epsm21) then
!        write(*,203)
!        goto 1000
!        else
!        ressum1=ressum
!        write(*,204)
!        goto 100
!        end if
!        end if
!
!!c mode=2: control residual decreasing level
!        if(mode.eq.2) then
!                if(ressum.le.epsm21) then
!                        write(*,203)
!                        goto 1000
!                        end if
!                        declev=ressum/ressum0
!                        if(declev.gt.epsm31) then
!                                ressum1=ressum
!                                write(*,205)
!                                goto 100
!                        else
!                                write(*,206)
!                                goto 1000
!                end if
!        end if
!
!!c mode=3: control residual decreasing rate
!        if(mode.eq.3) then
!        if(ressum.le.epsm21) then
!        write(*,203)
!        goto 1000
!        end if
!        decrat=ressum/ressum1
!        if(decrat.le.epsm41) then
!        ressum1=ressum
!        write(*,207)
!        goto 100
!        else
!        write(*,208)
!        goto 1000
!        end if
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
1000    return
        end
