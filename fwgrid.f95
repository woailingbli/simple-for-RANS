!c================================c
        subroutine w_solve0(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode,eps1,eps2,eps3,eps4,ressum, &
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
!c--------------------------------c
        do j=j2,m2
            if(j.ge.jb.and.j.le.je) then
                call trdgmk(k1,kb,j,nj,nk,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt)
                call trdgmk(ke,n1,j,nj,nk,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt)
            else
                call trdgmk(k1,n1,j,nj,nk,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt)
            end if
        enddo
!c--------------------------------c
154     call residup(j2,k2,m2,n2,jb,je,kb,ke,nj,nk,ressum,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt,res)
        write(*,200) ressum

        if(ressum.ge.ressum1) then
        do 160 j=j2,m2
        do 160 k=k2,n2
        if((j.ge.jb.and.j.lt.je).and.(k.ge.kb.and.k.lt.ke)) goto 160
        ff(j,k)=fw(j,k)
160     continue
!        write(*,201)
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
        subroutine w_solve1(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode,eps1,eps2,eps3,eps4,ressum0,ressum,&
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
!        write(*,11) ressum0
11      format(' *','  l-1   ','initial residual ressum=',1pe10.3,8x,'*')
!c--------------------------------c
100     call w_solve0(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode1,eps1,eps2,eps3,eps4,ressum,ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!        write(*,12) ressum
12      format(' *','  l-1   ','postLBL residual ressum=',1pe10.3,8x,'*')

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
        goto 100
        else
!        write(*,204)
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
        subroutine w_solve2(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode,eps1,eps2,eps3,eps4,ressum0,ressum,&
                          ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!c================================c
        include  'table.prc'
        include  'table.w_gd2'
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
!        write(*,11) ressum0
11      format(' *','  l-2   ','initial residual ressum=',1pe10.3,8x,'*')
!c--------------------------------c
100     call w_solve0(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode1,eps1,eps2,eps3,eps4,ressum,ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
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

200     j12=1
        k12=2
        j22=j12+1
        k22=k12+1
        m22=m12-1
        n22=n12-1
        call w_amgcoi(nj,nk,nj2,nk2,j22,k22,m22,n22,jb2,je2,kb2,ke2,app,res,aiw,aie,ajs,ajn,akb,akt,w_app2,w_bpp2, &
                    w_aiw2,w_aie2,w_ajs2,w_ajn2,w_akb2,w_akt2)
        call w_solve1(j12,k12,m12,n12,jb2,je2,kb2,ke2,ni2,nj2,nk2,mode2,eps12,eps22,eps32,eps42,ressum20,ressum2, &
                    w_ff2,w_res2,w_app2,w_bpp2,w_aiw2,w_aie2,w_ajs2,w_ajn2,w_akb2,w_akt2)
        call w_amgcri(nj,nk,nj2,nk2,j22,k22,m22,n22,jb2,je2,kb2,ke2,ff,w_ff2)
        goto 100
201     format(' *        total residual satisfied,   return        *')
202     format(' *        residual not satisfied, continuing        *')
203     format(' *        declev not reached,decrat exceeded        *')
204     format(' *        declev satisfied,  go out of L-B-L        *')
1000    return
        end


!c================================c
        subroutine w_solve3(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode,eps1,eps2,eps3,eps4,ressum0,ressum, &
                          ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!c================================c
        include  'table.prc'
        include  'table.w_gd3'
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
100     call w_solve0(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode1,eps1,eps2,eps3,eps4,ressum,ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!        write(*,12) ressum
12      format(' *','  l-3   ','postLBL residual ressum=',1pe10.3,8x,'*')

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

200     j13=1
        k13=2
        j23=j13+1
        k23=k13+1
        m23=m13-1
        n23=n13-1
        call w_amgcoi(nj,nk,nj3,nk3,j23,k23,m23,n23,jb3,je3,kb3,ke3,app,res,aiw,aie,ajs,ajn,akb,akt, &
                    w_app3,w_bpp3,w_aiw3,w_aie3,w_ajs3,w_ajn3,w_akb3,w_akt3)
        call w_solve2(j13,k13,m13,n13,jb3,je3,kb3,ke3,ni3,nj3,nk3,mode2,eps13,eps23,eps33,eps43,ressum03,ressum3, &
                    w_ff3,w_res3,w_app3,w_bpp3,w_aiw3,w_aie3,w_ajs3,w_ajn3,w_akb3,w_akt3)
        call w_amgcri(nj,nk,nj3,nk3,j23,k23,m23,n23,jb3,je3,kb3,ke3,ff,w_ff3)
        goto 100
201     format(' *        total residual satisfied,   return        *')
202     format(' *        residual not satisfied, continuing        *')
203     format(' *        declev not reached,decrat exceeded        *')
204     format(' *        declev satisfied,  go out of L-B-L        *')
1000    return
        end


!c================================c
        subroutine w_solve4(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode,eps1,eps2,eps3,eps4,ressum0,ressum, &
                          ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!c================================c
        include  'table.prc'
        include  'table.w_gd4'
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
100     call w_solve0(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode1,eps1,eps2,eps3,eps4,ressum,ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
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

200     j14=1
        k14=2
        j24=j14+1
        k24=k14+1
        m24=m14-1
        n24=n14-1
        call w_amgcoi(nj,nk,nj4,nk4,j24,k24,m24,n24,jb4,je4,kb4,ke4,app,res,aiw,aie,ajs,ajn,akb,akt, &
                    w_app4,w_bpp4,w_aiw4,w_aie4,w_ajs4,w_ajn4,w_akb4,w_akt4)
        call w_solve3(j14,k14,m14,n14,jb4,je4,kb4,ke4,ni4,nj4,nk4,mode2,eps14,eps24,eps34,eps44,ressum04,ressum4, &
                    w_ff4,w_res4,w_app4,w_bpp4,w_aiw4,w_aie4,w_ajs4,w_ajn4,w_akb4,w_akt4)
        call w_amgcri(nj,nk,nj4,nk4,j24,k24,m24,n24,jb4,je4,kb4,ke4,ff,w_ff4)
        goto 100
201     format(' *        total residual satisfied,   return        *')
202     format(' *        residual not satisfied, continuing        *')
203     format(' *        declev not reached,decrat exceeded        *')
204     format(' *        declev satisfied,  go out of L-B-L        *')
1000    return
        end


!c================================c
        subroutine w_solve5(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode,eps1,eps2,eps3,eps4,ressum0,ressum, &
                          ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!c================================c
        include  'table.prc'
        include  'table.w_gd5'
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
100     call w_solve0(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode1,eps1,eps2,eps3,eps4,ressum,ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
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

200     j15=1
        k15=2
        j25=j15+1
        k25=k15+1
        m25=m15-1
        n25=n15-1
        call w_amgcoi(nj,nk,nj5,nk5,j25,k25,m25,n25,jb5,je5,kb5,ke5,app,res,aiw,aie,ajs,ajn,akb,akt, &
                    w_app5,w_bpp5,w_aiw5,w_aie5,w_ajs5,w_ajn5,w_akb5,w_akt5)
        call w_solve4(j15,k15,m15,n15,jb5,je5,kb5,ke5,ni5,nj5,nk5,mode2,eps15,eps25,eps35,eps45,ressum05,ressum5, &
                    w_ff5,w_res5,w_app5,w_bpp5,w_aiw5,w_aie5,w_ajs5,w_ajn5,w_akb5,w_akt5)
        call w_amgcri(nj,nk,nj5,nk5,j25,k25,m25,n25,jb5,je5,kb5,ke5,ff,w_ff5)
        goto 100
201     format(' *        total residual satisfied,   return        *')
202     format(' *        residual not satisfied, continuing        *')
203     format(' *        declev not reached,decrat exceeded        *')
204     format(' *        declev satisfied,  go out of L-B-L        *')
1000    return
        end


!c================================c
        subroutine w_solve6(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode,eps1,eps2,eps3,eps4,ressum0,ressum, &
                          ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!c================================c
        include  'table.prc'
        include  'table.w_gd6'
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
11      format(' *','  l-6   ','initial residual ressum=',1pe10.3,8x,'*')
!c--------------------------------c
100     call w_solve0(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode1,eps1,eps2,eps3,eps4,ressum,ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
        write(*,12) ressum
12      format(' *','  l-6   ','postLBL residual ressum=',1pe10.3,8x,'*')

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

200     j16=1
        k16=2
        j26=j16+1
        k26=k16+1
        m26=m16-1
        n26=n16-1
        call w_amgcoi(nj,nk,nj6,nk6,j26,k26,m26,n26,jb6,je6,kb6,ke6,app,res,aiw,aie,ajs,ajn,akb,akt, &
                    w_app6,w_bpp6,w_aiw6,w_aie6,w_ajs6,w_ajn6,w_akb6,w_akt6)
        call w_solve5(j16,k16,m16,n16,jb6,je6,kb6,ke6,ni6,nj6,nk6,mode2,eps16,eps26,eps36,eps46,ressum06,ressum6, &
                    w_ff6,w_res6,w_app6,w_bpp6,w_aiw6,w_aie6,w_ajs6,w_ajn6,w_akb6,w_akt6)
        call w_amgcri(nj,nk,nj6,nk6,j26,k26,m26,n26,jb6,je6,kb6,ke6,ff,w_ff6)
        write(*,*)'*        The end of accelerating process           *'
        write(*,*)'*                                                  *'
        goto 100
201     format(' *        total residual satisfied,   return        *')
202     format(' *        residual not satisfied, continuing        *')
203     format(' *        declev not satisfied,   continuing        *')
204     format(' *        declev satisfied,  go out of L-B-L        *')
1000    return
        end

!c================================c
        subroutine w_solve7(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode,eps1,eps2,eps3,eps4,ressum0,ressum, &
                          ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!c================================c
        include  'table.prc'
        include  'table.w_gd7'
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
100     call w_solve0(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode1,eps1,eps2,eps3,eps4,ressum,ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
        write(*,12) ressum
12      format(' *','  l-7   ','postLBL residual ressum=',1pe10.3,8x,'*')

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

200     j17=1
        k17=2
        j27=j17+1
        k27=k17+1
        m27=m17-1
        n27=n17-1
        call w_amgcoi(nj,nk,nj7,nk7,j27,k27,m27,n27,jb7,je7,kb7,ke7,app,res,aiw,aie,ajs,ajn,akb,akt, &
                    w_app7,w_bpp7,w_aiw7,w_aie7,w_ajs7,w_ajn7,w_akb7,w_akt7)
        call w_solve6(j17,k17,m17,n17,jb7,je7,kb7,ke7,ni7,nj7,nk7,mode2,eps17,eps27,eps37,eps47,ressum07,ressum7, &
                    w_ff7,w_res7,w_app7,w_bpp7,w_aiw7,w_aie7,w_ajs7,w_ajn7,w_akb7,w_akt7)
        call w_amgcri(nj,nk,nj7,nk7,j27,k27,m27,n27,jb7,je7,kb7,ke7,ff,w_ff7)
        write(*,*)'*        The end of accelerating process           *'
        write(*,*)'*                                                  *'
        goto 100
201     format(' *        total residual satisfied,   return        *')
202     format(' *        residual not satisfied, continuing        *')
203     format(' *        declev not satisfied,   continuing        *')
204     format(' *        declev satisfied,  go out of L-B-L        *')
1000    return
        end

!c================================c
        subroutine w_solve8(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode,eps1,eps2,eps3,eps4,ressum0,ressum, &
                          ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
!c================================c
        include  'table.prc'
        include  'table.w_gd8'
        dimension app(nj,nk),bpp(nj,nk),aiw(nj,nk),aie(nj,nk),ajs(nj,nk),ajn(nj,nk),akb(nj,nk),akt(nj,nk)
        dimension ff(nj,nk),res(nj,nk)
!c--------------------------------c
        j2=j1+1
        k2=k1+1
        m2=m1-1
        n2=n1-1
        mode1=3
        mode2=2
        m18=m1
        call residup(j2,k2,m2,n2,jb,je,kb,ke,nj,nk,ressum0,ff,app,bpp,aiw,aie,ajs,ajn,akb,akt,res)
        write(*,11) ressum0
11      format(' *','  l-8   ','initial residual ressum=',1pe10.3,8x,'*')
!c--------------------------------c
100     call w_solve0(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,mode1,eps1,eps2,eps3,eps4,ressum,ff,res,app,bpp,aiw,aie,ajs,ajn,akb,akt)
        write(*,12) ressum
12      format(' *','  l-8   ','postLBL residual ressum=',1pe10.3,8x,'*')

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

200     j18=1
        k18=2
        j28=j18+1
        k28=k18+1
        m28=m18-1
        n28=n18-1
        call w_amgcoi(nj,nk,nj8,nk8,j28,k28,m28,n28,jb8,je8,kb8,ke8,app,res,aiw,aie,ajs,ajn,akb,akt, &
                    w_app8,w_bpp8,w_aiw8,w_aie8,w_ajs8,w_ajn8,w_akb8,w_akt8)
        call w_solve7(j18,k18,m18,n18,jb8,je8,kb8,ke8,ni8,nj8,nk8,mode2,eps18,eps28,eps38,eps48,ressum08,ressum8, &
                    w_ff8,w_res8,w_app8,w_bpp8,w_aiw8,w_aie8,w_ajs8,w_ajn8,w_akb8,w_akt8)
        call w_amgcri(nj,nk,nj8,nk8,j28,k28,m28,n28,jb8,je8,kb8,ke8,ff,w_ff8)
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
        subroutine w_amgcoi(nj1,nk1,nj2,nk2,j22,k22,m22,n22,jb2,je2,kb2,ke2,w_app1,w_res1,w_aiw1,w_aie1,&
                           w_ajs1,w_ajn1,w_akb1,w_akt1, w_app2,w_bpp2,w_aiw2,w_aie2,w_ajs2,w_ajn2,w_akb2,w_akt2)
!c===============================c
        include  'table.prc'
        dimension w_app1(nj1,nk1),w_res1(nj1,nk1),w_aiw1(nj1,nk1),w_aie1(nj1,nk1),w_ajs1(nj1,nk1),w_ajn1(nj1,nk1), &
                  w_akb1(nj1,nk1),w_akt1(nj1,nk1)
        dimension w_app2(nj2,nk2),w_bpp2(nj2,nk2),w_aiw2(nj2,nk2),w_aie2(nj2,nk2),w_ajs2(nj2,nk2),w_ajn2(nj2,nk2), &
                  w_akb2(nj2,nk2),w_akt2(nj2,nk2)
!c-------------------------------c
        do 100 j=j22,m22
        do 100 k=k22,n22
        if((j.ge.jb2.and.j.lt.je2).and.(k.ge.kb2.and.k.lt.ke2)) goto 100
        w_app2(j,k)=w_app1(j,2*k-1)+w_app1(j,2*k-2)-w_akb1(j,2*k-1)-w_akt1(j,2*k-2)
        w_aiw2(j,k)=0
        w_aie2(j,k)=0
        w_ajs2(j,k)=w_ajs1(j,2*k-2)+w_ajs1(j,2*k-1)
        w_ajn2(j,k)=w_ajn1(j,2*k-2)+w_ajn1(j,2*k-1)
        w_akb2(j,k)=w_akb1(j,2*k-2)
        w_akt2(j,k)=w_akt1(j,2*k-1)

        w_bpp2(j,k)=w_res1(2*j-2,k)+w_res1(2*j-1,k)
100     continue
        return
        end

!c===============================c
        subroutine w_amgcri(nj1,nk1,nj2,nk2,j22,k22,m22,n22,jb2,je2,kb2,ke2,ff1,ff2)
!c===============================c
        include  'table.prc'
        dimension ff1(nj1,nk1),ff2(nj2,nk2)
!c-------------------------------c
        do 100 j=j22,m22
        do 100 k=k22,n22
        if((j.ge.jb2.and.j.lt.je2).and.(k.ge.kb2.and.k.lt.ke2)) goto 100
        ff1(2*j-2,k)=ff1(2*j-2,k)+ff2(j,k)
        ff1(2*j-1,k)=ff1(2*j-1,k)+ff2(j,k)
100     continue
        return
        end
