!c================================c
        subroutine v_solve0(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode,eps1,eps2,eps3,eps4,ressum, &
                          ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!c================================c
        include  'table.prc'
        include  'table.dmw'
        dimension ff(nj,nk),res(nj,nk)
        dimension app(nj,nk),bpp(nj,nk),aiw(nj,nk),aie(nj,nk),ajs(nj,nk),ajn(nj,nk),akb(nj,nk),akt(nj,nk)
        dimension bj(njw),bk(nkw),fw(njw,nkw)
!c--------------------------------c
        j2=j1+1
        k2=k1+1
        m2=m1-1
        n2=n1-1
        iter=0
!c--------------------------------c
        call residup(j2,k2,m2,n2,jb,je,kb,ke,nj,nk,ressum0,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt,res)
        ressum1=ressum0
!        write(*,10) ressum0
10      format(' *',8x,'initial residual ressum=',1pe10.3,8x,'*')
!c--------------------------------c
100     iter=iter+1
        do j=j2,m2
            do k=k2,n2
                if((j.ge.jb.and.j.lt.je).and.(k.ge.kb.and.k.lt.ke)) cycle
                fw(j,k)=ff(j,k)
            enddo
        enddo

!c--------------------------------c
        do k=k2,n2
            if(k.ge.kb.and.k.lt.ke) then
                call trdgmj(j1,jb,k,nj,nk,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt)
                call trdgmj(je-1,m1,k,nj,nk,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt)
            else
                call trdgmj(j1,m1,k,nj,nk,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt)
            end if
        enddo
         call residup(j2,k2,m2,n2,jb,je,kb,ke,nj,nk,ressum,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt,res)
!        write(*,200) ressum
!c--------------------------------c
        do j=j2,m2
            if(j.ge.jb.and.j.lt.je) then
                call trdgmk(k1,kb,j,nj,nk,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt)
                call trdgmk(ke-1,n1,j,nj,nk,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt)
            else
                call trdgmk(k1,n1,j,nj,nk,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt)
            end if
        enddo
!c--------------------------------c
154     call residup(j2,k2,m2,n2,jb,je,kb,ke,nj,nk,ressum,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt,res)
!        write(*,200) ressum

        if(ressum.ge.ressum1) then
        do 160 j=j2,m2
        do 160 k=k2,n2
        if((j.ge.jb.and.j.lt.je).and.(k.ge.kb.and.k.lt.ke)) goto 160
        ff(j,k)=fw(j,k)
160     continue
        write(*,201)
        goto 1000
        end if

!c mode=1: control residual precision
        if(mode.eq.1) then
        if(ressum.le.eps2) then
!        write(*,202)
        goto 1000
        else
        ressum1=ressum
!        write(*,203)
        goto 100
        end if
        end if

!c mode=2: control residual decreasing level
        if(mode.eq.2) then
        if(ressum.le.eps2) then
!        write(*,202)
        goto 1000
        end if
        declev=ressum/ressum0
        if(declev.gt.eps3) then
        ressum1=ressum
!        write(*,206)
        goto 100
        else
!        write(*,207)
        goto 1000
        end if
        end if

!c mode=3: control residual decreasing rate
        if(mode.eq.3) then
        if(ressum.le.eps2) then
!        write(*,202)
        goto 1000
        end if
        decrat=ressum/ressum1
        if(decrat.le.eps4) then
        ressum1=ressum
!        write(*,204)
        goto 100
        else
!        write(*,205)
        goto 1000
        end if
        end if
200     format(' *',8x,'total  residual  ressum=',1pe10.3,8x,'*')
201     format(' *        residual not improving,     return        *')
202     format(' *        total residual satisfied,   return        *')
203     format(' *        residual not satisfied, continuing        *')
204     format(' *        decrat small, go on next iteration        *')
205     format(' *        decrat large, stop the LBL process        *')
206     format(' *        declev not satisfied,   continuing        *')
207     format(' *        declev satisfied,  go out of L-B-L        *')
1000    return
        end

!c================================c
        subroutine v_solve01(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode,eps1,eps2,eps3,eps4,ressum, &
                          ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!c================================c
        include  'table.prc'
        include  'table.dmw'
        dimension ff(nj,nk),res(nj,nk)
        dimension app(nj,nk),bpp(nj,nk),aiw(nj,nk),aie(nj,nk),ajs(nj,nk),ajn(nj,nk),akb(nj,nk),akt(nj,nk)
        dimension bj(njw),bk(nkw),fw(njw,nkw),p(nj),q(nj)
!c--------------------------------c
        j2=j1+1
        k2=k1+1
        m2=m1-1
        n2=n1-1
        iter=0
!c--------------------------------c
        call residup(j2,k2,m2,n2,jb,je,kb,ke,nj,nk,ressum0,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt,res)
        ressum1=ressum0
!        write(*,10) ressum0
10      format(' *',8x,'initial residual ressum=',1pe10.3,8x,'*')
!c--------------------------------c
100     iter=iter+1
        do j=j2,m2
            do k=k2,n2
                if((j.ge.jb.and.j.lt.je).and.(k.ge.kb.and.k.lt.ke)) cycle
                fw(j,k)=ff(j,k)
            enddo
        enddo
        do k=2,3
                j3=j2+1
                m3=m2-1
                denom=1.0d+00/app(j2,k)
                p(j2)=ajn(j2,k)*denom
                temp=bpp(j2,k)+akb(j2,k)*ff(j2,k-1)+akt(j2,k)*ff(j2,k+1)
                q(j2)=temp*denom
                do j=j3,m2
                denom=1.0d+00/(app(j,k)-p(j-1)*ajs(j,k))
                p(j)=ajn(j,k)*denom
                temp=bpp(j,k)+akb(j,k)*ff(j,k-1)+akt(j,k)*ff(j,k+1)
                q(j)=(temp+ajs(j,k)*q(j-1))*denom
                enddo
                ff(m2,k)=q(m2)
                do j=m3,j2,-1
                ff(j,k)=ff(j+1,k)*p(j)+q(j)
                enddo
        enddo

        call residup(j2,k2,m2,n2,jb,je,kb,ke,nj,nk,ressum,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt,res)
        !        write(*,11) ressum
        do j=j2,m2
                b2=bpp(j,2)+ajs(j,2)*ff(j-1,2)+ajn(j,2)*ff(j+1,2)
                b3=bpp(j,3)+ajs(j,3)*ff(j-1,3)+ajn(j,3)*ff(j+1,3)
                d1=app(j,2)*app(j,3)-akt(j,2)*akb(j,3)
                d2=app(j,3)*b2+akt(j,2)*b3
                d3=app(j,2)*b3+akb(j,3)*b2
                ff(j,2)=d2/d1
                ff(j,3)=d3/d1
        enddo
        call residup(j2,k2,m2,n2,jb,je,kb,ke,nj,nk,ressum,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt,res)
!        write(*,11) ressum
        
!c--------------------------------c

!c--------------------------------c
    ! call residup(j2,k2,m2,n2,jb,je,kb,ke,nj,nk,ressum,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt,res)
154        ressum=0.0d+00
        do j=j2,m2
        do k=2,3
        res(j,k)=bpp(j,k)-app(j,k)*ff(j,k)+ajs(j,k)*ff(j-1,k)+ajn(j,k)*ff(j+1,k)+ &
                 akb(j,k)*ff(j,k-1)+akt(j,k)*ff(j,k+1)
        ressum=ressum+dabs(res(j,k))
        enddo
        enddo
!        write(*,200) ressum

        if(ressum.ge.ressum1) then
        do j=j2,m2
        do k=k2,n2
        ff(j,k)=fw(j,k)
        enddo
        enddo
        goto 1000
        end if

!c mode=1: control residual precision
        if(mode.eq.1) then
        if(ressum.le.eps2) then
!        write(*,202)
        goto 1000
        else
        ressum1=ressum
!        write(*,203)
        goto 100
        end if
        end if

!c mode=2: control residual decreasing level
        if(mode.eq.2) then
        if(ressum.le.eps2) then
!        write(*,202)
        goto 1000
        end if
        declev=ressum/ressum0
        if(declev.gt.eps3) then
        ressum1=ressum
!        write(*,206)
        goto 100
        else
!        write(*,207)
        goto 1000
        end if
        end if

!c mode=3: control residual decreasing rate
        if(mode.eq.3) then
        if(ressum.le.eps2) then
!        write(*,202)
        goto 1000
        end if
        decrat=ressum/ressum1
        if(decrat.le.eps4) then
        ressum1=ressum
!        write(*,204)
        goto 100
        else
!        write(*,205)
        goto 1000
        end if
        end if
200     format(' *',8x,'total  residual  ressum=',1pe10.3,8x,'*')
201     format(' *        residual not improving,     return        *')
202     format(' *        total residual satisfied,   return        *')
203     format(' *        residual not satisfied, continuing        *')
204     format(' *        decrat small, go on next iteration        *')
205     format(' *        decrat large, stop the LBL process        *')
206     format(' *        declev not satisfied,   continuing        *')
207     format(' *        declev satisfied,  go out of L-B-L        *')
1000    return
        end
        
!c================================c
        subroutine v_solve1(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode,eps1,eps2,eps3,eps4,ressum0,ressum,&
                          ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!c================================c
        include  'table.prc'
        include  'table.dmw'
        dimension app(nj,nk),bpp(nj,nk),aiw(nj,nk),aie(nj,nk),ajs(nj,nk),ajn(nj,nk),akb(nj,nk),akt(nj,nk)
        dimension ff(nj,nk),res(nj,nk)
        dimension p(nj),q(nj)
        dimension v_app0(nj,3),v_bpp0(nj,3),v_aiw0(nj,3),v_aie0(nj,3),v_ajs0(nj,3),v_ajn0(nj,3),v_akb0(nj,3),v_akt0(nj,3),res2(nj,3)
        dimension v_ff0(nj,3)
!c--------------------------------c
        j2=j1+1
        k2=k1+1
        m2=m1-1
        n2=n1-1
        mode1=3
        call residup(j2,k2,m2,n2,jb,je,kb,ke,nj,nk,ressum0,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt,res)
        write(*,11) ressum0
11      format(' *','  l-1   ','initial residual ressum=',1pe10.3,8x,'*')
!c--------------------------------c
100     call v_solve0(j1,k1,m1,n1,1000,je,kb,ke,ni,nj,nk,mode1,eps1,eps2,eps3,eps4,ressum,ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
       write(*,12) ressum
12      format(' *','  l-1   ','postLBL residual ressum=',1pe10.3,8x,'*')
        mode=1
!c mode=1: control residual precision
        if(mode.eq.1) then
        if(ressum.le.eps2) then
!        write(*,201)
        goto 1000
        else
!        write(*,202)
        goto 1000
        end if
        end if

!c mode=2: control residual decreasing level
        if(mode.eq.2) then
        if(ressum.le.eps2) then
!c        write(*,201)
        goto 1000
        end if
        declev=ressum/ressum0
        if(declev.gt.eps3) then
!c        write(*,203)
        goto 200
        else
!c        write(*,204)
        goto 1000
        end if
        end if
     
200     v_app0=0.0d+00
        v_ajs0=0.0d+00
        v_ajn0=0.0d+00
        v_akt0=0.0d+00
        v_akb0=0.0d+00
        v_aie0=0.0d+00
        v_aiw0=0.0d+00
        v_bpp0=0.0d+00
        p=0.0d+00
        q=0.0d+00
        v_ff0=0.0d+00
        res2=0.0d+00
        do j=j2,m2
                v_app0(j,2)=app(j,2)+app(j,3)-akt(j,2)-akb(j,3)
                v_ajs0(j,2)=ajs(j,2)+ajs(j,3)
                v_ajn0(j,2)=ajn(j,2)+ajn(j,3)
                v_akt0(j,2)=akt(j,3)
                v_akb0(j,2)=akb(j,2)
                v_aie0(j,2)=0.0d+00
                v_aiw0(j,2)=0.0d+00
                v_bpp0(j,2)=res(j,2)+res(j,3)
        enddo
        
        j3=j2+1
        m3=m2-1
        denom=1.0d+00/v_app0(j2,2)
        p(j2)=v_ajn0(j2,2)*denom
        temp=v_bpp0(j2,2)
        q(j2)=temp*denom
        do j=j3,m2
        denom=1.0d+00/(v_app0(j,2)-p(j-1)*v_ajs0(j,2))
        p(j)=v_ajn0(j,2)*denom
        temp=v_bpp0(j,2)
        q(j)=(temp+v_ajs0(j,2)*q(j-1))*denom
        enddo
        v_ff0(m2,2)=q(m2)
        do j=m3,j2,-1
        v_ff0(j,2)=v_ff0(j+1,2)*p(j)+q(j)
        enddo

        call residup(j2,k2,m2,n2,jb,je,kb,ke,nj,nk,ressum2,v_ff0,v_app0,v_bpp0,v_aiw0,v_aie0,v_ajs0,v_ajn0,v_akb0,v_akt0,res2)
!        write(*,11) ressum2
        do j=j2,m2
        ff(j,2)=ff(j,2)+v_ff0(j,2)
        ff(j,3)=ff(j,3)+v_ff0(j,2)      
        enddo
        
        call residup(j2,k2,m2,n2,jb,je,kb,ke,nj,nk,ressum1,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt,res)
!        write(*,11) ressum1        
        goto 100
201     format(' *        total residual satisfied,   return        *')
202     format(' *        residual not satisfied, continuing        *')
203     format(' *        declev not reached,decrat exceeded        *')
204     format(' *        declev satisfied,  go out of L-B-L        *')
1000    return
        end


!c================================c
        subroutine v_solve2(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode,eps1,eps2,eps3,eps4,ressum0,ressum, &
                          ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!c================================c
        include  'table.prc'
        include  'table.v_gd2'
        dimension app(nj,nk),bpp(nj,nk),aiw(nj,nk),aie(nj,nk),ajs(nj,nk),ajn(nj,nk),akb(nj,nk),akt(nj,nk)
        dimension ff(nj,nk),res(nj,nk)
!c--------------------------------c
        j2=j1+1
        k2=k1+1
        m2=m1-1
        n2=n1-1
        mode1=3
        mode2=2
        call residup(j2,k2,m2,n2,jb,je,kb,ke,nj,nk,ressum0,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt,res)
!        write(*,11) ressum0
11      format(' *','  l-2   ','initial residual ressum=',1pe10.3,8x,'*')
!c--------------------------------c
100     call v_solve0(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode1,eps1,eps2,eps3,eps4,ressum,ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!        write(*,12) ressum
12      format(' *','  l-2   ','postLBL residual ressum=',1pe10.3,8x,'*')

!c mode=1: control residual precision
        if(mode.eq.1) then
        if(ressum.le.eps2) then
!        write(*,201)
        goto 1000
        else
!        write(*,202)
        goto 200
        end if
        end if

!c mode=2: control residual decreasing level
        if(mode.eq.2) then
        if(ressum.le.eps2) then
!        write(*,201)
        goto 1000
        end if
        declev=ressum/ressum0
        if(declev.gt.eps3) then
!        write(*,203)
        goto 200
        else
!        write(*,204)
        goto 1000
        end if
        end if

200     j12=2
        k12=1
        j22=j12+1
        k22=k12+1
        m22=m12-1
        n22=n12-1
        call v_amgcoi(nj,nk,nj2,nk2,j22,k22,m22,n22,jb2,je2,kb2,ke2,app,res,aiw,aie,ajs,ajn,akb,akt, &
                    v_app2,v_bpp2,v_aiw2,v_aie2,v_ajs2,v_ajn2,v_akb2,v_akt2)
        call v_solve1(j12,k12,m12,n12,jb2,je2,kb2,ke2,ni2,nj2,nk2,mode2,eps12,eps22,eps32,eps42,ressum02,ressum2, &
                    v_ff2,v_res2,v_app2,v_bpp2,v_aiw2,v_aie2,v_ajs2,v_ajn2,v_akb2,v_akt2)
        call v_amgcri_2_to_4(nj,nk,nj2,nk2,j22,k22,m22,n22,jb2,je2,kb2,ke2,ff,v_ff2)
        goto 100
201     format(' *        total residual satisfied,   return        *')
202     format(' *        residual not satisfied, continuing        *')
203     format(' *        declev not reached,decrat exceeded        *')
204     format(' *        declev satisfied,  go out of L-B-L        *')
1000    return
        end


!c================================c
        subroutine v_solve3(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode,eps1,eps2,eps3,eps4,ressum0,ressum, &
                          ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!c================================c
        include  'table.prc'
        include  'table.v_gd3'
        dimension app(nj,nk),bpp(nj,nk),aiw(nj,nk),aie(nj,nk),ajs(nj,nk),ajn(nj,nk),akb(nj,nk),akt(nj,nk)
        dimension ff(nj,nk),res(nj,nk)
!c--------------------------------c
        j2=j1+1
        k2=k1+1
        m2=m1-1
        n2=n1-1
        mode1=3
        mode2=2
        call residup(j2,k2,m2,n2,jb,je,kb,ke,nj,nk,ressum0,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt,res)
!        write(*,11) ressum0
11      format(' *','  l-3   ','initial residual ressum=',1pe10.3,8x,'*')
!c--------------------------------c
100     call v_solve0(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode1,eps1,eps2,eps3,eps4,ressum,ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!        write(*,12) ressum
12      format(' *','  l-3   ','postLBL residual ressum=',1pe10.3,8x,'*')

!c mode=1: control residual precision
        if(mode.eq.1) then
        if(ressum.le.eps2) then
        write(*,201)
        goto 1000
        else
        write(*,202)
        goto 200
        end if
        end if

!c mode=2: control residual decreasing level
        if(mode.eq.2) then
        if(ressum.le.eps2) then
!        write(*,201)
        goto 1000
        end if
        declev=ressum/ressum0
        if(declev.gt.eps3) then
!        write(*,203)
        goto 200
        else
!        write(*,204)
        goto 1000
        end if
        end if

200     j13=2
        k13=1
        j23=j13+1
        k23=k13+1
        m23=m13-1
        n23=n13-1
        call v_amgcoi(nj,nk,nj3,nk3,j23,k23,m23,n23,jb3,je3,kb3,ke3,app,res,aiw,aie,ajs,ajn,akb,akt, &
                    v_app3,v_bpp3,v_aiw3,v_aie3,v_ajs3,v_ajn3,v_akb3,v_akt3)
        call v_solve2(j13,k13,m13,n13,jb3,je3,kb3,ke3,ni3,nj3,nk3,mode2,eps13,eps23,eps33,eps43,ressum03,ressum3, &
                    v_ff3,v_res3,v_app3,v_bpp3,v_aiw3,v_aie3,v_ajs3,v_ajn3,v_akb3,v_akt3)
        call v_amgcri(nj,nk,nj3,nk3,j23,k23,m23,n23,jb3,je3,kb3,ke3,ff,v_ff3)
        goto 100
201     format(' *        total residual satisfied,   return        *')
202     format(' *        residual not satisfied, continuing        *')
203     format(' *        declev not reached,decrat exceeded        *')
204     format(' *        declev satisfied,  go out of L-B-L        *')
1000    return
        end


!c================================c
        subroutine v_solve4(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode,eps1,eps2,eps3,eps4,ressum0,ressum, &
                          ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!c================================c
        include  'table.prc'
        include  'table.v_gd4'
        dimension app(nj,nk),bpp(nj,nk),aiw(nj,nk),aie(nj,nk),ajs(nj,nk),ajn(nj,nk),akb(nj,nk),akt(nj,nk)
        dimension ff(nj,nk),res(nj,nk)
!c--------------------------------c
        j2=j1+1
        k2=k1+1
        m2=m1-1
        n2=n1-1
        mode1=3
        mode2=2
        call residup(j2,k2,m2,n2,jb,je,kb,ke,nj,nk,ressum0,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt,res)
     
!        write(*,11) ressum0
11      format(' *','  l-4   ','initial residual ressum=',1pe10.3,8x,'*')
!c--------------------------------c
100     call v_solve0(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode1,eps1,eps2,eps3,eps4,ressum,ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!        write(*,12) ressum
12      format(' *','  l-4   ','postLBL residual ressum=',1pe10.3,8x,'*')

!c mode=1: control residual precision
        if(mode.eq.1) then
        if(ressum.le.eps2) then
!        write(*,201)
        goto 1000
        else
!        write(*,202)
        goto 200
        end if
        end if

!c mode=2: control residual decreasing level
        if(mode.eq.2) then
        if(ressum.le.eps2) then
!        write(*,201)
        goto 1000
        end if
        declev=ressum/ressum0
        if(declev.gt.eps3) then
!        write(*,203)
        goto 200
        else
!        write(*,204)
        goto 1000
        end if
        end if

200     j14=2
        k14=1
        j24=j14+1
        k24=k14+1
        m24=m14-1
        n24=n14-1
        call v_amgcoi(nj,nk,nj4,nk4,j24,k24,m24,n24,jb4,je4,kb4,ke4,app,res,aiw,aie,ajs,ajn,akb,akt, &
                    v_app4,v_bpp4,v_aiw4,v_aie4,v_ajs4,v_ajn4,v_akb4,v_akt4)
        call v_solve3(j14,k14,m14,n14,jb4,je4,kb4,ke4,ni4,nj4,nk4,mode2,eps14,eps24,eps34,eps44,ressum04,ressum4, &
                    v_ff4,v_res4,v_app4,v_bpp4,v_aiw4,v_aie4,v_ajs4,v_ajn4,v_akb4,v_akt4)
        call v_amgcri(nj,nk,nj4,nk4,j24,k24,m24,n24,jb4,je4,kb4,ke4,ff,v_ff4)
        goto 100
201     format(' *        total residual satisfied,   return        *')
202     format(' *        residual not satisfied, continuing        *')
203     format(' *        declev not reached,decrat exceeded        *')
204     format(' *        declev satisfied,  go out of L-B-L        *')
1000    return
        end


!c================================c
        subroutine v_solve5(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode,eps1,eps2,eps3,eps4,ressum0,ressum, &
                          ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!c================================c
        include  'table.prc'
        include  'table.v_gd5'
        dimension app(nj,nk),bpp(nj,nk),aiw(nj,nk),aie(nj,nk),ajs(nj,nk),ajn(nj,nk),akb(nj,nk),akt(nj,nk)
        dimension ff(nj,nk),res(nj,nk)
!c--------------------------------c
        j2=j1+1
        k2=k1+1
        m2=m1-1
        n2=n1-1
        mode1=3
        mode2=2
        call residup(j2,k2,m2,n2,jb,je,kb,ke,nj,nk,ressum0,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt,res)
!        write(*,11) ressum0
11      format(' *','  l-5   ','initial residual ressum=',1pe10.3,8x,'*')
        goto 200
!c--------------------------------c
100     call v_solve0(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode1,eps1,eps2,eps3,eps4,ressum,ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!        write(*,12) ressum
12      format(' *','  l-5   ','postLBL residual ressum=',1pe10.3,8x,'*')

!c mode=1: control residual precision
        if(mode.eq.1) then
        if(ressum.le.eps2) then
!        write(*,201)
        goto 1000
        else
!        write(*,202)
        goto 200
        end if
        end if

!c mode=2: control residual decreasing level
        if(mode.eq.2) then
        if(ressum.le.eps2) then
!        write(*,201)
        goto 1000
        end if
        declev=ressum/ressum0
        if(declev.gt.eps3) then
!        write(*,203)
        goto 200
        else
!        write(*,204)
        goto 1000
        end if
        end if

200     j15=2
        k15=1
        j25=j15+1
        k25=k15+1
        m25=m15-1
        n25=n15-1
        call v_amgcoi(nj,nk,nj5,nk5,j25,k25,m25,n25,jb5,je5,kb5,ke5,app,res,aiw,aie,ajs,ajn,akb,akt, &
                    v_app5,v_bpp5,v_aiw5,v_aie5,v_ajs5,v_ajn5,v_akb5,v_akt5)
        call v_solve4(j15,k15,m15,n15,jb5,je5,kb5,ke5,ni5,nj5,nk5,mode2,eps15,eps25,eps35,eps45,ressum05,ressum5, &
                    v_ff5,v_res5,v_app5,v_bpp5,v_aiw5,v_aie5,v_ajs5,v_ajn5,v_akb5,v_akt5)
        call v_amgcri(nj,nk,nj5,nk5,j25,k25,m25,n25,jb5,je5,kb5,ke5,ff,v_ff5)
        goto 100
201     format(' *        total residual satisfied,   return        *')
202     format(' *        residual not satisfied, continuing        *')
203     format(' *        declev not reached,decrat exceeded        *')
204     format(' *        declev satisfied,  go out of L-B-L        *')
1000    return
        end


!c================================c
        subroutine v_solve6(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode,eps1,eps2,eps3,eps4,ressum0,ressum, &
                          ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!c================================c
        include  'table.prc'
        include  'table.v_gd6'
        dimension app(nj,nk),bpp(nj,nk),aiw(nj,nk),aie(nj,nk),ajs(nj,nk),ajn(nj,nk),akb(nj,nk),akt(nj,nk)
        dimension ff(nj,nk),res(nj,nk)
!c--------------------------------c
        j2=j1+1
        k2=k1+1
        m2=m1-1
        n2=n1-1
        mode1=3
        mode2=2
        call residup(j2,k2,m2,n2,jb,je,kb,ke,nj,nk,ressum0,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt,res)
!        write(*,11) ressum0
11      format(' *','  l-6   ','initial residual ressum=',1pe10.3,8x,'*')
!c--------------------------------c
100     call v_solve0(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode1,eps1,eps2,eps3,eps4,ressum,ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!        write(*,12) ressum
12      format(' *','  l-6   ','postLBL residual ressum=',1pe10.3,8x,'*')

!c mode=1: control residual precision
        if(mode.eq.1) then
        if(ressum.le.eps2) then
!        write(*,201)
        goto 1000
        else
!        write(*,202)
        goto 200
        end if
        end if

!c mode=2: control residual decreasing level
        if(mode.eq.2) then
        if(ressum.le.eps2) then
!        write(*,201)
        goto 1000
        end if
        declev=ressum/ressum0
        if(declev.gt.eps3) then
!        write(*,203)
        goto 200
        else
!        write(*,204)
        goto 1000
        end if
        end if

200     j16=2
        k16=1
        j26=j16+1
        k26=k16+1
        m26=m16-1
        n26=n16-1
        call v_amgcoi(nj,nk,nj6,nk6,j26,k26,m26,n26,jb6,je6,kb6,ke6,app,res,aiw,aie,ajs,ajn,akb,akt, &
                    v_app6,v_bpp6,v_aiw6,v_aie6,v_ajs6,v_ajn6,v_akb6,v_akt6)
        call v_solve5(j16,k16,m16,n16,jb6,je6,kb6,ke6,ni6,nj6,nk6,mode2,eps16,eps26,eps36,eps46,ressum06,ressum6, &
                    v_ff6,v_res6,v_app6,v_bpp6,v_aiw6,v_aie6,v_ajs6,v_ajn6,v_akb6,v_akt6)
        call v_amgcri(nj,nk,nj6,nk6,j26,k26,m26,n26,jb6,je6,kb6,ke6,ff,v_ff6)
!        write(*,*)'*        The end of accelerating process           *'
!        write(*,*)'*                                                  *'
        goto 100
201     format(' *        total residual satisfied,   return        *')
202     format(' *        residual not satisfied, continuing        *')
203     format(' *        declev not satisfied,   continuing        *')
204     format(' *        declev satisfied,  go out of L-B-L        *')
1000    return
        end

!c================================c
        subroutine v_solve7(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode,eps1,eps2,eps3,eps4,ressum0,ressum, &
                          ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!c================================c
        include  'table.prc'
        include  'table.v_gd7'
        dimension app(nj,nk),bpp(nj,nk),aiw(nj,nk),aie(nj,nk),ajs(nj,nk),ajn(nj,nk),akb(nj,nk),akt(nj,nk)
        dimension ff(nj,nk),res(nj,nk)
!c--------------------------------c
        j2=j1+1
        k2=k1+1
        m2=m1-1
        n2=n1-1
        mode1=3
        mode2=2
        call residup(j2,k2,m2,n2,jb,je,kb,ke,nj,nk,ressum0,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt,res)
        write(*,11) ressum0
11      format(' *','  l-7   ','initial residual ressum=',1pe10.3,8x,'*')
!c--------------------------------c
100     call v_solve0(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode1,eps1,eps2,eps3,eps4,ressum,ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
        write(*,12) ressum
12      format(' *','  l-7   ','postLBL residual ressum=',1pe10.3,8x,'*')

!c mode=1: control residual precision
        if(mode.eq.1) then
        if(ressum.le.eps2) then
!        write(*,201)
        goto 1000
        else
!        write(*,202)
        goto 200
        end if
        end if

!c mode=2: control residual decreasing level
        if(mode.eq.2) then
        if(ressum.le.eps2) then
        write(*,201)
        goto 1000
        end if
        declev=ressum/ressum0
        if(declev.gt.eps3) then
        write(*,203)
        goto 200
        else
        write(*,204)
        goto 1000
        end if
        end if

200     j17=2
        k17=1
        j27=j17+1
        k27=k17+1
        m27=m17-1
        n27=n17-1
        call v_amgcoi(nj,nk,nj7,nk7,j27,k27,m27,n27,jb7,je7,kb7,ke7,app,res,aiw,aie,ajs,ajn,akb,akt, &
                    v_app7,v_bpp7,v_aiw7,v_aie7,v_ajs7,v_ajn7,v_akb7,v_akt7)
        call v_solve6(j17,k17,m17,n17,jb7,je7,kb7,ke7,ni7,nj7,nk7,mode2,eps17,eps27,eps37,eps47,ressum07,ressum7, &
                    v_ff7,v_res7,v_app7,v_bpp7,v_aiw7,v_aie7,v_ajs7,v_ajn7,v_akb7,v_akt7)
        call v_amgcri(nj,nk,nj7,nk7,j27,k27,m27,n27,jb7,je7,kb7,ke7,ff,v_ff7)
!        write(*,*)'*        The end of accelerating process           *'
!        write(*,*)'*                                                  *'
        goto 100
201     format(' *        total residual satisfied,   return        *')
202     format(' *        residual not satisfied, continuing        *')
203     format(' *        declev not satisfied,   continuing        *')
204     format(' *        declev satisfied,  go out of L-B-L        *')
1000    return
        end

!c================================c
        subroutine v_solve8(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode,eps1,eps2,eps3,eps4,ressum0,ressum, &
                          ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!c================================c
        include  'table.prc'
        include  'table.v_gd8'
        dimension app(nj,nk),bpp(nj,nk),aiw(nj,nk),aie(nj,nk),ajs(nj,nk),ajn(nj,nk),akb(nj,nk),akt(nj,nk)
        dimension ff(nj,nk),res(nj,nk)
!c--------------------------------c
        j2=j1+1
        k2=k1+1
        m2=m1-1
        n2=n1-1
        mode1=3
        mode2=2
        call residup(j2,k2,m2,n2,jb,je,kb,ke,nj,nk,ressum0,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt,res)
       write(*,11) ressum0
11      format(' *','  l-8   ','initial residual ressum=',1pe10.3,8x,'*')
!c--------------------------------c
100     call v_solve0(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode1,eps1,eps2,eps3,eps4,ressum,ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
        write(*,12) ressum
12      format(' *','  l-8   ','postLBL residual ressum=',1pe10.3,8x,'*')

!c mode=1: control residual precision
        if(mode.eq.1) then
        if(ressum.le.eps2) then
!        write(*,201)
        goto 1000
        else
!        write(*,202)
        goto 200
        end if
        end if

!c mode=2: control residual decreasing level
        if(mode.eq.2) then
        if(ressum.le.eps2) then
        write(*,201)
        goto 1000
        end if
        declev=ressum/ressum0
        if(declev.gt.eps3) then
        write(*,203)
        goto 200
        else
        write(*,204)
        goto 1000
        end if
        end if

200     j18=2
        k18=1
        j28=j18+1
        k28=k18+1
        m28=m18-1
        n28=n18-1
        call v_amgcoi(nj,nk,nj8,nk8,j28,k28,m28,n28,jb8,je8,kb8,ke8,app,res,aiw,aie,ajs,ajn,akb,akt, &
                    v_app8,v_bpp8,v_aiw8,v_aie8,v_ajs8,v_ajn8,v_akb8,v_akt8)
        call v_solve7(j18,k18,m18,n18,jb8,je8,kb8,ke8,ni8,nj8,nk8,mode2,eps18,eps28,eps38,eps48,ressum08,ressum8, &
                    v_ff8,v_res8,v_app8,v_bpp8,v_aiw8,v_aie8,v_ajs8,v_ajn8,v_akb8,v_akt8)
        call v_amgcri(nj,nk,nj8,nk8,j28,k28,m28,n28,jb8,je8,kb8,ke8,ff,v_ff8)
!        write(*,*)'*        The end of accelerating process           *'
!        write(*,*)'*                                                  *'
        goto 100
201     format(' *        total residual satisfied,   return        *')
202     format(' *        residual not satisfied, continuing        *')
203     format(' *        declev not satisfied,   continuing        *')
204     format(' *        declev satisfied,  go out of L-B-L        *')
1000    return
        end

!c================================c
        subroutine v_solve9(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode,eps1,eps2,eps3,eps4,ressum0,ressum, &
                          ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!c================================c
        include  'table.prc'
        include  'table.v_gd9'
        dimension app(nj,nk),bpp(nj,nk),aiw(nj,nk),aie(nj,nk),ajs(nj,nk),ajn(nj,nk),akb(nj,nk),akt(nj,nk)
        dimension ff(nj,nk),res(nj,nk)
!c--------------------------------c
        j2=j1+1
        k2=k1+1
        m2=m1-1
        n2=n1-1
        mode1=3
        mode2=2
        m19=m1
        call residup(j2,k2,m2,n2,jb,je,kb,ke,nj,nk,ressum0,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt,res)
        write(*,11) ressum0
11      format(' *','  l-9   ','initial residual ressum=',1pe10.3,8x,'*')
!c--------------------------------c
100     call v_solve0(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode1,eps1,eps2,eps3,eps4,ressum,ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
        write(*,12) ressum
12      format(' *','  l-9   ','postLBL residual ressum=',1pe10.3,8x,'*')

!c mode=1: control residual precision
        if(mode.eq.1) then
        if(ressum.le.eps2) then
!        write(*,201)
        goto 1000
        else
!        write(*,202)
        goto 200
        end if
        end if

!c mode=2: control residual decreasing level
        if(mode.eq.2) then
        if(ressum.le.eps2) then
        write(*,201)
        goto 1000
        end if
        declev=ressum/ressum0
        if(declev.gt.eps3) then
        write(*,203)
        goto 200
        else
        write(*,204)
        goto 1000
        end if
        end if

200     j19=2
        k19=1
        j29=j19+1
        k29=k19+1
        m29=m19-1
        n29=n19-1
        call v_amgcoi(nj,nk,nj9,nk9,j29,k29,m29,n29,jb9,je9,kb9,ke9,app,res,aiw,aie,ajs,ajn,akb,akt, &
                    v_app9,v_bpp9,v_aiw9,v_aie9,v_ajs9,v_ajn9,v_akb9,v_akt9)
        call v_solve8(j19,k19,m19,n19,jb9,je9,kb9,ke9,ni9,nj9,nk9,mode2,eps19,eps29,eps39,eps49,ressum09,ressum9, &
                    v_ff9,v_res9,v_app9,v_bpp9,v_aiw9,v_aie9,v_ajs9,v_ajn9,v_akb9,v_akt9)
        call v_amgcri(nj,nk,nj9,nk9,j29,k29,m29,n29,jb9,je9,kb9,ke9,ff,v_ff9)
        write(*,*)'*        The end of accelerating process           *'
        write(*,*)'*                                                  *'
        goto 100
201     format(' *        total residual satisfied,   return        *')
202     format(' *        residual not satisfied, continuing        *')
203     format(' *        declev not satisfied,   continuing        *')
204     format(' *        declev satisfied,  go out of L-B-L        *')
1000    return
        end




!c===============================c
        subroutine v_amgcoi(nj1,nk1,nj2,nk2,j22,k22,m22,n22,jb2,je2,kb2,ke2,v_app1,v_res1,v_aiw1,v_aie1,&
                           v_ajs1,v_ajn1,v_akb1,v_akt1, v_app2,v_bpp2,v_aiw2,v_aie2,v_ajs2,v_ajn2,v_akb2,v_akt2)
!c===============================c
        include  'table.prc'
        dimension v_app1(nj1,nk1),v_res1(nj1,nk1),v_aiw1(nj1,nk1),v_aie1(nj1,nk1),v_ajs1(nj1,nk1),v_ajn1(nj1,nk1), &
                  v_akb1(nj1,nk1),v_akt1(nj1,nk1)
        dimension v_app2(nj2,nk2),v_bpp2(nj2,nk2),v_aiw2(nj2,nk2),v_aie2(nj2,nk2),v_ajs2(nj2,nk2),v_ajn2(nj2,nk2), &
                  v_akb2(nj2,nk2),v_akt2(nj2,nk2),delta(nj2,nk2)
!c-------------------------------c
        delta=0.0
        do 100 j=j22,m22
        do 100 k=k22,n22
        if((j.ge.jb2.and.j.lt.je2).and.(k.ge.kb2.and.k.lt.ke2)) goto 100
        v_aiw2(j,k)=0
        v_aie2(j,k)=0
        v_ajs2(j,k)=v_ajs1(j,2*k-2)+v_ajs1(j,2*k-1)
        v_ajn2(j,k)=v_ajn1(j,2*k-2)+v_ajn1(j,2*k-1)
        v_akb2(j,k)=v_akb1(j,2*k-2)
        v_akt2(j,k)=v_akt1(j,2*k-1)
        v_app2(j,k)=v_app1(j,2*k-1)+v_app1(j,2*k-2)-v_akb1(j,2*k-1)-v_akt1(j,2*k-2)
                
        delta(j,k)= v_app2(j,k)-v_ajs2(j,k)-v_ajn2(j,k)-v_akb2(j,k)- v_akt2(j,k)
        v_bpp2(j,k)=v_res1(j,2*k-2)+v_res1(j,2*k-1)
100     continue
        return
        end

!c===============================c
        subroutine v_amgcri(nj1,nk1,nj2,nk2,j22,k22,m22,n22,jb2,je2,kb2,ke2,ff1,ff2)
!c===============================c
        include  'table.prc'
        dimension ff1(nj1,nk1),ff2(nj2,nk2)
!c-------------------------------c
        do 100 j=j22,m22
        do 100 k=k22,n22
        if((j.ge.jb2.and.j.lt.je2).and.(k.ge.kb2.and.k.lt.ke2)) goto 100
        ff1(j,2*k-2)=ff1(j,2*k-2)+ff2(j,k)
        ff1(j,2*k-1)=ff1(j,2*k-1)+ff2(j,k)
100     continue
        return
        end

!c===============================c
        subroutine v_amgcri_2_to_4(nj1,nk1,nj2,nk2,j22,k22,m22,n22,jb2,je2,kb2,ke2,ff1,ff2)
!c===============================c
        include  'table.prc'
        dimension ff1(nj1,nk1),ff2(nj2,nk2)
!c-------------------------------c
        do j=j22,m22
                ff1(j,2)=ff1(j,2)+ff2(j,2)
                ff1(j,5)=ff1(j,5)+ff2(j,3)
        enddo 
        do j=j22,jb2-1
                ff1(j,3)=ff1(j,3)+ff2(j,2)
                ff1(j,4)=ff1(j,4)+ff2(j,3) 
        enddo
        do j=je2,m22
                ff1(j,3)=ff1(j,3)+ff2(j,2)
                ff1(j,4)=ff1(j,4)+ff2(j,3) 
        enddo
        return
        end
