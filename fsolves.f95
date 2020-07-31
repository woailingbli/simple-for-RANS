!c================================c
        subroutine solve0(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode,eps1,eps2,eps3,eps4,ressum, &
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
!c        write(*,10) ressum0
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
                call trdgpj(j1,jb,k,nj,nk,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt)
                call trdgpj(je-1,m1,k,nj,nk,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt)
            else
                call trdgpj(j1,m1,k,nj,nk,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt)
            end if
        enddo
!c--------------------------------c
        do j=j2,m2
            if(j.ge.jb.and.j.lt.je) then
                call trdgpk(k1,kb,j,nj,nk,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt)
                call trdgpk(ke-1,n1,j,nj,nk,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt)
            else
                call trdgpk(k1,n1,j,nj,nk,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt)
            end if
        enddo
!c--------------------------------c
154     call residup(j2,k2,m2,n2,jb,je,kb,ke,nj,nk,ressum,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt,res)
!c        write(*,200) ressum

        if(ressum.ge.ressum1) then
        do 160 j=j2,m2
        do 160 k=k2,n2
        if((j.ge.jb.and.j.lt.je).and.(k.ge.kb.and.k.lt.ke)) goto 160
        ff(j,k)=fw(j,k)
160     continue
!c        write(*,201)
        goto 1000
        end if

!c mode=1: control residual precision
        if(mode.eq.1) then
        if(ressum.le.eps2) then
!c        write(*,202)
        goto 1000
        else
        ressum1=ressum
!c        write(*,203)
        goto 100
        end if
        end if

!c mode=2: control residual decreasing level
        if(mode.eq.2) then
        if(ressum.le.eps2) then
!c        write(*,202)
        goto 1000
        end if
        declev=ressum/ressum0
        if(declev.gt.eps3) then
        ressum1=ressum
!c        write(*,206)
        goto 100
        else
!c        write(*,207)
        goto 1000
        end if
        end if

!c mode=3: control residual decreasing rate
        if(mode.eq.3) then
        if(ressum.le.eps2) then
!c        write(*,202)
        goto 1000
        end if
        decrat=ressum/ressum1
        if(decrat.le.eps4) then
        ressum1=ressum
!c        write(*,204)
        goto 100
        else
!c        write(*,205)
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
        subroutine solve1(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode,eps1,eps2,eps3,eps4,ressum0,ressum,&
                          ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!c================================c
        include  'table.prc'
        include  'table.dmw'
        dimension app(nj,nk),bpp(nj,nk),aiw(nj,nk),aie(nj,nk),ajs(nj,nk),ajn(nj,nk),akb(nj,nk),akt(nj,nk)
        dimension ff(nj,nk),res(nj,nk)
!c--------------------------------c
        j2=j1+1
        k2=k1+1
        m2=m1-1
        n2=n1-1
        mode1=3
        call residup(j2,k2,m2,n2,jb,je,kb,ke,nj,nk,ressum0,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt,res)
!c        write(*,11) ressum0
11      format(' *','  l-1   ','initial residual ressum=',1pe10.3,8x,'*')
!c--------------------------------c
100     call solve0(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode1,eps1,eps2,eps3,eps4,ressum,ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!c        write(*,12) ressum
12      format(' *','  l-1   ','postLBL residual ressum=',1pe10.3,8x,'*')

!c mode=1: control residual precision
        if(mode.eq.1) then
        if(ressum.le.eps2) then
!c        write(*,201)
        goto 1000
        else
!c        write(*,202)
        goto 200
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
        goto 100
        else
!c        write(*,204)
        goto 1000
        end if
        end if
     
200     app0=0.0d+00
        bpp0=0.0d+00
        aiw0=0.0d+00
        aie0=0.0d+00

        do 120 j=j2,m2
        do 120 k=k2,n2
        if((j.ge.jb.and.j.lt.je).and.(k.ge.kb.and.k.lt.ke)) goto 120
        aiw0=aiw0+aiw(j,k)
        aie0=aie0+aie(j,k)
        app0=app0+app(j,k)-(ajs(j,k)+ajn(j,k)+akb(j,k)+akt(j,k))
122     bpp0=bpp0+res(j,k)
120     continue
        do 130 j=j2,m2
        do 130 k=k2,n2
        if((j.ge.jb.and.j.lt.je).and.(k.ge.kb.and.k.lt.ke)) goto 130

131     ff(j,k)=ff(j,k)+bi0
130     continue
        goto 100
201     format(' *        total residual satisfied,   return        *')
202     format(' *        residual not satisfied, continuing        *')
203     format(' *        declev not reached,decrat exceeded        *')
204     format(' *        declev satisfied,  go out of L-B-L        *')
1000    return
        end


!c================================c
        subroutine solve2(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode,eps1,eps2,eps3,eps4,ressum0,ressum,&
                          ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!c================================c
        include  'table.prc'
        include  'table.gd2'
        dimension app(nj,nk),bpp(nj,nk),aiw(nj,nk),aie(nj,nk),ajs(nj,nk),ajn(nj,nk),akb(nj,nk),akt(nj,nk)
        dimension ff(nj,nk) ,res(nj,nk)
!c--------------------------------c
        j2=j1+1
        k2=k1+1
        m2=m1-1
        n2=n1-1
        mode1=3
        mode2=2
        call residup(j2,k2,m2,n2,jb,je,kb,ke,nj,nk,ressum0,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt,res)
!c        write(*,11) ressum0
11      format(' *','  l-2   ','initial residual ressum=',1pe10.3,8x,'*')
!c--------------------------------c
100     call solve0(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode1,eps1,eps2,eps3,eps4,ressum,ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!c        write(*,12) ressum
12      format(' *','  l-2   ','postLBL residual ressum=',1pe10.3,8x,'*')

!c mode=1: control residual precision
        if(mode.eq.1) then
        if(ressum.le.eps2) then
!c        write(*,201)
        goto 1000
        else
!c        write(*,202)
        goto 200
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

200     j12=1
        k12=1
        j22=j12+1
        k22=k12+1
        m22=m12-1
        n22=n12-1
        call amgcoi(nj,nk,nj2,nk2,j22,k22,m22,n22,jb2,je2,kb2,ke2,app,res,aiw,aie,ajs,ajn,akb,akt,app2,bpp2, &
                    aiw2,aie2,ajs2,ajn2,akb2,akt2)
        call solve1(j12,k12,m12,n12,jb2,je2,kb2,ke2,ni2,nj2,nk2,mode2,eps12,eps22,eps32,eps42,ressum20,ressum2, &
                    ff2,res2,app2,bpp2,aiw2,aie2,ajs2,ajn2,akb2,akt2)
        call amgcri(nj,nk,nj2,nk2,j22,k22,m22,n22,jb2,je2,kb2,ke2,ff,ff2)
        goto 100
201     format(' *        total residual satisfied,   return        *')
202     format(' *        residual not satisfied, continuing        *')
203     format(' *        declev not reached,decrat exceeded        *')
204     format(' *        declev satisfied,  go out of L-B-L        *')
1000    return
        end


!c================================c
        subroutine solve3(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode,eps1,eps2,eps3,eps4,ressum0,ressum, &
                          ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!c================================c
        include  'table.prc'
        include  'table.gd3'
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
!c        write(*,11) ressum0
11      format(' *','  l-3   ','initial residual ressum=',1pe10.3,8x,'*')
!c--------------------------------c
100     call solve0(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode1,eps1,eps2,eps3,eps4,ressum,ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!c        write(*,12) ressum
12      format(' *','  l-3   ','postLBL residual ressum=',1pe10.3,8x,'*')

!c mode=1: control residual precision
        if(mode.eq.1) then
        if(ressum.le.eps2) then
!c        write(*,201)
        goto 1000
        else
!c        write(*,202)
        goto 200
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

200     j13=1
        k13=1
        j23=j13+1
        k23=k13+1
        m23=m13-1
        n23=n13-1
        call amgcoi(nj,nk,nj3,nk3,j23,k23,m23,n23,jb3,je3,kb3,ke3,app,res,aiw,aie,ajs,ajn,akb,akt, &
                    app3,bpp3,aiw3,aie3,ajs3,ajn3,akb3,akt3)
        call solve2(j13,k13,m13,n13,jb3,je3,kb3,ke3,ni3,nj3,nk3,mode2,eps13,eps23,eps33,eps43,ressum03,ressum3, &
                    ff3,res3,app3,bpp3,aiw3,aie3,ajs3,ajn3,akb3,akt3)
        call amgcri(nj,nk,nj3,nk3,j23,k23,m23,n23,jb3,je3,kb3,ke3,ff,ff3)
        goto 100
201     format(' *        total residual satisfied,   return        *')
202     format(' *        residual not satisfied, continuing        *')
203     format(' *        declev not reached,decrat exceeded        *')
204     format(' *        declev satisfied,  go out of L-B-L        *')
1000    return
        end


!c================================c
        subroutine solve4(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode,eps1,eps2,eps3,eps4,ressum0,ressum, &
                          ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!c================================c
        include  'table.prc'
        include  'table.gd4'
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
     
!c        write(*,11) ressum0
11      format(' *','  l-4   ','initial residual ressum=',1pe10.3,8x,'*')
!c--------------------------------c
100     call solve0(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode1,eps1,eps2,eps3,eps4,ressum,ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!c        write(*,12) ressum
12      format(' *','  l-4   ','postLBL residual ressum=',1pe10.3,8x,'*')

!c mode=1: control residual precision
        if(mode.eq.1) then
        if(ressum.le.eps2) then
!c        write(*,201)
        goto 1000
        else
!c        write(*,202)
        goto 200
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

200     j14=1
        k14=1
        j24=j14+1
        k24=k14+1
        m24=m14-1
        n24=n14-1
        call amgcoi(nj,nk,nj4,nk4,j24,k24,m24,n24,jb4,je4,kb4,ke4,app,res,aiw,aie,ajs,ajn,akb,akt, &
                    app4,bpp4,aiw4,aie4,ajs4,ajn4,akb4,akt4)
        call solve3(j14,k14,m14,n14,jb4,je4,kb4,ke4,ni4,nj4,nk4,mode2,eps14,eps24,eps34,eps44,ressum04,ressum4, &
                    ff4,res4,app4,bpp4,aiw4,aie4,ajs4,ajn4,akb4,akt4)
        call amgcri(nj,nk,nj4,nk4,j24,k24,m24,n24,jb4,je4,kb4,ke4,ff,ff4)
        goto 100
201     format(' *        total residual satisfied,   return        *')
202     format(' *        residual not satisfied, continuing        *')
203     format(' *        declev not reached,decrat exceeded        *')
204     format(' *        declev satisfied,  go out of L-B-L        *')
1000    return
        end


!c================================c
        subroutine solve5(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode,eps1,eps2,eps3,eps4,ressum0,ressum, &
                          ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!c================================c
        include  'table.prc'
        include  'table.gd5'
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
!c        write(*,11) ressum0
11      format(' *','  l-5   ','initial residual ressum=',1pe10.3,8x,'*')
!c        goto 200
!c--------------------------------c
100     call solve0(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode1,eps1,eps2,eps3,eps4,ressum,ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!c        write(*,12) ressum
12      format(' *','  l-5   ','postLBL residual ressum=',1pe10.3,8x,'*')

!c mode=1: control residual precision
        if(mode.eq.1) then
        if(ressum.le.eps2) then
!c        write(*,201)
        goto 1000
        else
!c        write(*,202)
        goto 200
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

200     j15=1
        k15=1
        j25=j15+1
        k25=k15+1
        m25=m15-1
        n25=n15-1
        call amgcoi(nj,nk,nj5,nk5,j25,k25,m25,n25,jb5,je5,kb5,ke5,app,res,aiw,aie,ajs,ajn,akb,akt, &
                    app5,bpp5,aiw5,aie5,ajs5,ajn5,akb5,akt5)
        call solve4(j15,k15,m15,n15,jb5,je5,kb5,ke5,ni5,nj5,nk5,mode2,eps15,eps25,eps35,eps45,ressum05,ressum5, &
                    ff5,res5,app5,bpp5,aiw5,aie5,ajs5,ajn5,akb5,akt5)
        call amgcri(nj,nk,nj5,nk5,j25,k25,m25,n25,jb5,je5,kb5,ke5,ff,ff5)
        goto 100
201     format(' *        total residual satisfied,   return        *')
202     format(' *        residual not satisfied, continuing        *')
203     format(' *        declev not reached,decrat exceeded        *')
204     format(' *        declev satisfied,  go out of L-B-L        *')
1000    return
        end


!c================================c
        subroutine solve6(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode,eps1,eps2,eps3,eps4,ressum0,ressum, &
                          ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!c================================c
        include  'table.prc'
        include  'table.gd6'
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
100     call solve0(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode1,eps1,eps2,eps3,eps4,ressum,ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!c        write(*,12) ressum
12      format(' *','  l-6   ','postLBL residual ressum=',1pe10.3,8x,'*')

!c mode=1: control residual precision
        if(mode.eq.1) then
        if(ressum.le.eps2) then
!c        write(*,201)
        goto 1000
        else
!c        write(*,202)
        goto 200
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

200     j16=1
        k16=1
        j26=j16+1
        k26=k16+1
        m26=m16-1
        n26=n16-1
        call amgcoi(nj,nk,nj6,nk6,j26,k26,m26,n26,jb6,je6,kb6,ke6,app,res,aiw,aie,ajs,ajn,akb,akt, &
                    app6,bpp6,aiw6,aie6,ajs6,ajn6,akb6,akt6)
        call solve5(j16,k16,m16,n16,jb6,je6,kb6,ke6,ni6,nj6,nk6,mode2,eps16,eps26,eps36,eps46,ressum06,ressum6, &
                    ff6,res6,app6,bpp6,aiw6,aie6,ajs6,ajn6,akb6,akt6)
        call amgcri(nj,nk,nj6,nk6,j26,k26,m26,n26,jb6,je6,kb6,ke6,ff,ff6)
!c        write(*,*)'*        The end of accelerating process           *'
!c        write(*,*)'*                                                  *'
        goto 100
201     format(' *        total residual satisfied,   return        *')
202     format(' *        residual not satisfied, continuing        *')
203     format(' *        declev not satisfied,   continuing        *')
204     format(' *        declev satisfied,  go out of L-B-L        *')
1000    return
        end

!c================================c
        subroutine solve7(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode,eps1,eps2,eps3,eps4,ressum0,ressum, &
                          ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!c================================c
        include  'table.prc'
        include  'table.gd7'
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
!       write(*,11) ressum0
11      format(' *','  l-7   ','initial residual ressum=',1pe10.3,8x,'*')
!c--------------------------------c
100     call solve0(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode1,eps1,eps2,eps3,eps4,ressum,ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!        write(*,12) ressum
12      format(' *','  l-7   ','postLBL residual ressum=',1pe10.3,8x,'*')

!c mode=1: control residual precision
        if(mode.eq.1) then
        if(ressum.le.eps2) then
!c        write(*,201)
        goto 1000
        else
!c        write(*,202)
        goto 200
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

200     j17=1
        k17=1
        j27=j17+1
        k27=k17+1
        m27=m17-1
        n27=n17-1
        call amgcoi(nj,nk,nj7,nk7,j27,k27,m27,n27,jb7,je7,kb7,ke7,app,res,aiw,aie,ajs,ajn,akb,akt, &
                    app7,bpp7,aiw7,aie7,ajs7,ajn7,akb7,akt7)
        call solve6(j17,k17,m17,n17,jb7,je7,kb7,ke7,ni7,nj7,nk7,mode2,eps17,eps27,eps37,eps47,ressum07,ressum7, &
                    ff7,res7,app7,bpp7,aiw7,aie7,ajs7,ajn7,akb7,akt7)
        call amgcri(nj,nk,nj7,nk7,j27,k27,m27,n27,jb7,je7,kb7,ke7,ff,ff7)
!c        write(*,*)'*        The end of accelerating process           *'
!c        write(*,*)'*                                                  *'
        goto 100
201     format(' *        total residual satisfied,   return        *')
202     format(' *        residual not satisfied, continuing        *')
203     format(' *        declev not satisfied,   continuing        *')
204     format(' *        declev satisfied,  go out of L-B-L        *')
1000    return
        end

!c================================c
        subroutine solve8(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode,eps1,eps2,eps3,eps4,ressum0,ressum, &
                          ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!c================================c
        include  'table.prc'
        include  'table.gd8'
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
100     call solve0(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode1,eps1,eps2,eps3,eps4,ressum,ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!        write(*,12) ressum
12      format(' *','  l-8   ','postLBL residual ressum=',1pe10.3,8x,'*')

!c mode=1: control residual precision
        if(mode.eq.1) then
        if(ressum.le.eps2) then
!c        write(*,201)
        goto 1000
        else
!c        write(*,202)
        goto 200
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

200     j18=1
        k18=1
        j28=j18+1
        k28=k18+1
        m28=m18-1
        n28=n18-1
        call amgcoi(nj,nk,nj8,nk8,j28,k28,m28,n28,jb8,je8,kb8,ke8,app,res,aiw,aie,ajs,ajn,akb,akt, &
                    app8,bpp8,aiw8,aie8,ajs8,ajn8,akb8,akt8)
        call solve7(j18,k18,m18,n18,jb8,je8,kb8,ke8,ni8,nj8,nk8,mode2,eps18,eps28,eps38,eps48,ressum08,ressum8, &
                    ff8,res8,app8,bpp8,aiw8,aie8,ajs8,ajn8,akb8,akt8)
        call amgcri(nj,nk,nj8,nk8,j28,k28,m28,n28,jb8,je8,kb8,ke8,ff,ff8)
!c        write(*,*)'*        The end of accelerating process           *'
!c        write(*,*)'*                                                  *'
        goto 100
201     format(' *        total residual satisfied,   return        *')
202     format(' *        residual not satisfied, continuing        *')
203     format(' *        declev not satisfied,   continuing        *')
204     format(' *        declev satisfied,  go out of L-B-L        *')
1000    return
        end


!c===============================c
        subroutine residum(j2,k2,m2,n2,jb,je,kb,ke,ni,nj,nk,ressum,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt,res)
!c===============================c
        include 'table.prc'
        dimension app(nj,nk),bpp(nj,nk),aiw(nj,nk),aie(nj,nk),ajs(nj,nk),ajn(nj,nk),akb(nj,nk),akt(nj,nk)
        dimension ff(nj,nk),res(nj,nk)
!c*******************************c
        ressum=0.0d+00
        do 101 j=j2,m2
        do 101 k=k2,n2
        if((j.ge.jb.and.j.lt.je).and.(k.ge.kb.and.k.lt.ke)) goto 101
        res(j,k)=bpp(j,k)-app(j,k)*ff(j,k)+aiw(j,k)*ff(j,k)+aie(j,k)*ff(j,k)+ajs(j,k)*ff(j-1,k)+ajn(j,k)*ff(j+1,k)+ &
                 akb(j,k)*ff(j,k-1)+akt(j,k)*ff(j,k+1)
        ressum=ressum+dabs(res(j,k))
101     continue
        return
        end
!c===============================c
        subroutine residup(j2,k2,m2,n2,jb,je,kb,ke,nj,nk,ressum,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt,res)
!c===============================c
        include 'table.prc'
        dimension app(nj,nk),bpp(nj,nk),aiw(nj,nk),aie(nj,nk),ajs(nj,nk),ajn(nj,nk),akb(nj,nk),akt(nj,nk)
        dimension ff(nj,nk),res(nj,nk)
!c*******************************c
        ressum=0.0d+00
        do 101 j=j2,m2
        do 101 k=k2,n2
        if((j.ge.jb.and.j.lt.je).and.(k.ge.kb.and.k.lt.ke)) goto 101
        res(j,k)=bpp(j,k)-app(j,k)*ff(j,k)+aiw(j,k)*ff(j,k)+aie(j,k)*ff(j,k)+ajs(j,k)*ff(j-1,k)+ajn(j,k)*ff(j+1,k)+ &
                 akb(j,k)*ff(j,k-1)+akt(j,k)*ff(j,k+1)
103     ressum=ressum+dabs(res(j,k))

101     continue
        return
        end
        
!c===============================c
        subroutine residual_all(j2,k2,m2,n2,jb,je,kb,ke,nj,nk,ressum,ff1,ff2,res)
!c===============================c
        include 'table.prc'
        dimension ff1(nj,nk),ff2(nj,nk),res(nj,nk)
!c*******************************c
        ressum=0.0d+00
        do 101 j=j2,m2
        do 101 k=k2,n2
        if((j.ge.jb.and.j.lt.je).and.(k.ge.kb.and.k.lt.ke)) goto 101
        res(j,k)=ff1(j,k)-ff2(j,k)
        ressum=ressum+dabs(res(j,k))
101     continue
        return
        end

!c********************************c
        subroutine trdgmj(j1,m1,k,nj,nk,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!c********************************c
        include 'table.prc'
        include 'table.dmw'
        dimension p(njw),q(njw)
        dimension ff(nj,nk)
        dimension app(nj,nk),bpp(nj,nk),aiw(nj,nk),aie(nj,nk),ajs(nj,nk),ajn(nj,nk),akb(nj,nk),akt(nj,nk)
!c--------------------------------c
        j2=j1+1
        j3=j2+1
        m2=m1-1
        m3=m2-1
        denom=1.0d+00/app(j2,k)
        p(j2)=ajn(j2,k)*denom
        temp=bpp(j2,k)+akb(j2,k)*ff(j2,k-1)+akt(j2,k)*ff(j2,k+1)
        q(j2)=temp*denom
        do 100 j=j3,m2
        denom=1.0d+00/(app(j,k)-p(j-1)*ajs(j,k))
        p(j)=ajn(j,k)*denom
        temp=bpp(j,k)+akb(j,k)*ff(j,k-1)+akt(j,k)*ff(j,k+1)
        q(j)=(temp+ajs(j,k)*q(j-1))*denom
100     continue
        ff(m2,k)=q(m2)
        do 101 j=m3,j2,-1
101     ff(j,k)=ff(j+1,k)*p(j)+q(j)
        return
        end

!c********************************c
        subroutine trdgmk(k1,n1,j,nj,nk,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!c********************************c
        include 'table.prc'               
        include 'table.dmw'
        dimension p(nkw),q(nkw)
        dimension ff(nj,nk)
        dimension app(nj,nk),bpp(nj,nk),aiw(nj,nk),aie(nj,nk),ajs(nj,nk),ajn(nj,nk),akb(nj,nk),akt(nj,nk)
!c--------------------------------c
        k2=k1+1
        k3=k2+1
        n2=n1-1
        n3=n2-1
        denom=1.0d+00/app(j,k2)
        p(k2)=akt(j,k2)*denom
        temp=bpp(j,k2)+ajs(j,k2)*ff(j-1,k2)+ajn(j,k2)*ff(j+1,k2)
        q(k2)=temp*denom
        do 100 k=k3,n2
        denom=1.0d+00/(app(j,k)-p(k-1)*akb(j,k))
        p(k)=akt(j,k)*denom
        temp=bpp(j,k)+ajs(j,k)*ff(j-1,k)+ajn(j,k)*ff(j+1,k)
        q(k)=(temp+akb(j,k)*q(k-1))*denom
100     continue
        ff(j,n2)=q(n2)
        do 101 k=n3,k2,-1
101     ff(j,k)=ff(j,k+1)*p(k)+q(k)
        return
        end


!c********************************c
        subroutine trdgpj(j1,m1,k,nj,nk,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!c********************************c
        include 'table.prc'
        include 'table.dmw'
        dimension p(njw),q(njw)
        dimension ff(nj,nk)
        dimension app(nj,nk),bpp(nj,nk),aiw(nj,nk),aie(nj,nk),ajs(nj,nk),ajn(nj,nk),akb(nj,nk),akt(nj,nk)
!c--------------------------------c 
        j2=j1+1
        j3=j2+1
        m2=m1-1
        m3=m2-1


        denom=1.0d+00/app(j2,k)
        p(j2)=ajn(j2,k)*denom
        temp=bpp(j2,k)+akb(j2,k)*ff(j2,k-1)+akt(j2,k)*ff(j2,k+1)
        q(j2)=temp*denom
        do 100 j=j3,m2
        denom=1.0d+00/(app(j,k)-p(j-1)*ajs(j,k))
        p(j)=ajn(j,k)*denom
        temp=bpp(j,k)+akb(j,k)*ff(j,k-1)+akt(j,k)*ff(j,k+1)
        q(j)=(temp+ajs(j,k)*q(j-1))*denom
100     continue

102     ff(m2,k)=q(m2)
        do 101 j=m3,j2,-1
101     ff(j,k)=ff(j+1,k)*p(j)+q(j)
        return
        end

!c********************************c
        subroutine trdgpk(k1,n1,j,nj,nk,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!c********************************c
        include 'table.prc'
        include 'table.dmw'
        dimension p(nkw),q(nkw)
        dimension ff(nj,nk)
        dimension app(nj,nk),bpp(nj,nk),aiw(nj,nk),aie(nj,nk),ajs(nj,nk),ajn(nj,nk),akb(nj,nk),akt(nj,nk)
!c--------------------------------c
        k2=k1+1
        k3=k2+1
        n2=n1-1
        n3=n2-1


        denom=1.0d+00/app(j,k2)
        p(k2)=akt(j,k2)*denom
        temp=bpp(j,k2)+ajs(j,k2)*ff(j-1,k2)+ajn(j,k2)*ff(j+1,k2)
        q(k2)=temp*denom
        do 100 k=k3,n2
        denom=1.0d+00/(app(j,k)-p(k-1)*akb(j,k))
        p(k)=akt(j,k)*denom
        temp=bpp(j,k)+ajs(j,k)*ff(j-1,k)+ajn(j,k)*ff(j+1,k)
        q(k)=(temp+akb(j,k)*q(k-1))*denom
100     continue

102     ff(j,n2)=q(n2)
        do 101 k=n3,k2,-1
101     ff(j,k)=ff(j,k+1)*p(k)+q(k)
        return
        end


!c===============================c
        subroutine amgcoi(nj1,nk1,nj2,nk2,j22,k22,m22,n22,jb2,je2,kb2,ke2,app1,res1,aiw1,aie1,ajs1,ajn1,akb1,akt1, &
                          app2,bpp2,aiw2,aie2,ajs2,ajn2,akb2,akt2)
!c===============================c
        include  'table.prc'
        dimension app1(nj1,nk1),res1(nj1,nk1),aiw1(nj1,nk1),aie1(nj1,nk1),ajs1(nj1,nk1),ajn1(nj1,nk1), &
                  akb1(nj1,nk1),akt1(nj1,nk1)
        dimension app2(nj2,nk2),bpp2(nj2,nk2),aiw2(nj2,nk2),aie2(nj2,nk2),ajs2(nj2,nk2),ajn2(nj2,nk2), &
                  akb2(nj2,nk2),akt2(nj2,nk2)
!c-------------------------------c
        do 100 j=j22,m22
        do 100 k=k22,n22
        if((j.ge.jb2.and.j.lt.je2).and.(k.ge.kb2.and.k.lt.ke2)) goto 100
        app2(j,k)=app1(2*j-2,2*k-2)+app1(2*j-1,2*k-2)+app1(2*j-2,2*k-1)+app1(2*j-1,2*k-1)-ajs1(2*j-1,2*k-2)-ajs1(2*j-1,2*k-1) &
                  -ajn1(2*j-2,2*k-2)-ajn1(2*j-2,2*k-1)-akb1(2*j-2,2*k-1)-akb1(2*j-1,2*k-1)-akt1(2*j-2,2*k-2)-akt1(2*j-1,2*k-2)
        aiw2(j,k)=aiw1(2*j-2,2*k-2)+aiw1(2*j-1,2*k-2)+aiw1(2*j-2,2*k-1)+aiw1(2*j-1,2*k-1)
        aie2(j,k)=aie1(2*j-2,2*k-2)+aie1(2*j-1,2*k-2)+aie1(2*j-2,2*k-1)+aie1(2*j-1,2*k-1)
        ajs2(j,k)=ajs1(2*j-2,2*k-2)+ajs1(2*j-2,2*k-1)
        ajn2(j,k)=ajn1(2*j-1,2*k-2)+ajn1(2*j-1,2*k-1)
        akb2(j,k)=akb1(2*j-2,2*k-2)+akb1(2*j-1,2*k-2)
        akt2(j,k)=akt1(2*j-2,2*k-1)+akt1(2*j-1,2*k-1)

        bpp2(j,k)=res1(2*j-2,2*k-2)+res1(2*j-1,2*k-2)+res1(2*j-2,2*k-1)+res1(2*j-1,2*k-1)
100     continue
        return
        end

!c===============================c
        subroutine amgcri(nj1,nk1,nj2,nk2,j22,k22,m22,n22,jb2,je2,kb2,ke2,ff1,ff2)
!c===============================c
        include  'table.prc'
        dimension ff1(nj1,nk1),ff2(nj2,nk2)
!c-------------------------------c
        do 100 j=j22,m22
        do 100 k=k22,n22
        if((j.ge.jb2.and.j.lt.je2).and.(k.ge.kb2.and.k.lt.ke2)) goto 100
        ff1(2*j-2,2*k-2)=ff1(2*j-2,2*k-2)+ff2(j,k)
        ff1(2*j-1,2*k-2)=ff1(2*j-1,2*k-2)+ff2(j,k)
        ff1(2*j-2,2*k-1)=ff1(2*j-2,2*k-1)+ff2(j,k)
        ff1(2*j-1,2*k-1)=ff1(2*j-1,2*k-1)+ff2(j,k)
100     continue
        return
        end
