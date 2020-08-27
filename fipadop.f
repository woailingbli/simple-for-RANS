c===============================c
        subroutine opd0fm
c===============================c
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
c----------------------------------------------------------------------c
20      format(' *',6x,'stralh    number    strh=',1pd10.3,6x,'*')
21      format(' *',6x,'renold    number    reno=',1pd10.3,6x,'*')
22      format(' *',6x,'prandt    number    prtl=',1pd10.3,6x,'*')
23      format(' *',6x,'eckert    number    eckt=',1pd10.3,6x,'*')

100     format(' *',4x,'the time marching step dt: dt=',1pd9.2,4x,'*')
101     format(' *',10x,' ip0=',i6,' jp0=',i6,' kp0=',i6,11x,'*')
102     format(' *',10x,' imx=',i6,' jmx=',i6,' kmx=',i6,11x,'*')
103     format(' *',10x,' iu0=',i6,' ju0=',i6,' ku0=',i6,11x,'*')
104     format(' *',10x,' iv0=',i6,' jv0=',i6,' kv0=',i6,11x,'*')
105     format(' *',10x,' iw0=',i6,' jw0=',i6,' kw0=',i6,11x,'*')
107     format(' *',4x,'   it0=',i6,'   jt0=',i6,'   kt0=',i6,4x,'*')
106     format
     *  (' *',8x,' jb1=',i6,' je1=',i6,' kb1=',i6,' ke1=',i6,8x,'*')
108     format(' *',2x,' ubo,vbo,wbo,tbo=',4d12.3,2x,' *')
109     format(' *',2x,' ubi,vbi,wbi,tbi=',4d12.3,2x,' *')
50      format(' *                                               *')
52      format(' *             The Criterion Numbers             *')
53      format(' *           The Coordinate Parameters           *')
54      format(' *-----------------------------------------------*')
c----------------------------------------------------------------------c
        open(unit=ioc,file='dat0.fmt',form='formatted')
        write(ioc,54)
        write(ioc,50)
        write(ioc,52)
        write(ioc,50)
        write(ioc,20) strh
        write(ioc,21) reno
        write(ioc,22) prtl
        write(ioc,23) eckt
        write(ioc,50)
        write(ioc,53)
        write(ioc,50)
        write(ioc,100) dt
        write(ioc,101) ip0,jp0,kp0
        write(ioc,102) imx,jmx,kmx
        write(ioc,103) iu0,ju0,ku0
        write(ioc,104) iv0,jv0,kv0
        write(ioc,105) iw0,jw0,kw0
        write(ioc,107) it0,jt0,kt0
        write(ioc,106) jb1,je1,kb1,ke1
        write(ioc,108) ubo,vbo,wbo,tbo
        write(ioc,109) ubi,vbi,wbi,tbi
1001    write(ioc,50)
        write(ioc,54)
        close(unit=ioc)
        return
        end
c===============================c
        subroutine opd0uf
c===============================c
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
c-------------------------------c
        open(unit=ioc,file='dat0.ufm',form='unformatted')
        write(ioc) strh,reno,prtl,eckt,dt
        write(ioc) jb1,je1,kb1,ke1
        write(ioc) ip0,jp0,kp0,imx,jmx,kmx
        write(ioc) iu0,ju0,ku0,iv0,jv0,kv0,iw0,jw0,kw0,it0,jt0,kt0
        write(ioc) ubo,vbo,wbo,tbo,ubi,vbi,wbi,tbi
1001    close(unit=ioc)
        return
        end

c===============================c
        subroutine opd1fm
c===============================c
        use mdugrd
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
c-------------------------------c
        open(unit=ioc,file='dat1.fmt',form='formatted')
        write(ioc,101)
101     format('*',25(1h-),'momentum grid on level one',25(1h-),'*')
        call wtcfmt(xu1,' xu1= ',ioc,iu0,imx,imx)
        call wtcfmt(yv1,' yv1= ',ioc,jv0,jmx,jmx)
        call wtcfmt(zw1,' zw1= ',ioc,kw0,kmx,kmx)
        close(unit=ioc)
        return
        end
c===============================c
        subroutine opd1uf
c===============================c
        use mdugrd
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
c-------------------------------c
        open(unit=ioc,file='dat1.ufm',form='unformatted')
        write(ioc)(xu1(i),i=iu0,imx)
        write(ioc)(yv1(j),j=jv0,jmx)
        write(ioc)(zw1(k),k=kw0,kmx)
        write(ioc)((ibp(j,k),j=jp0,jmx),k=kp0,kmx)
        write(ioc)((ibu(j,k),j=ju0,jmx),k=ku0,kmx)
        write(ioc)((ibv(j,k),j=jv0,jmx),k=kv0,kmx)
        write(ioc)((ibw(j,k),j=jw0,jmx),k=kw0,kmx)
        close(unit=ioc)
        return
        end
c===============================c
        subroutine opd2uf
c===============================c
        include 'table.prc'
        include 'table.cre'
        include 'table.gdm'
c-------------------------------c
        open(unit=ioc,file='dat2.ufm',form='unformatted')
        do i=0,levmgd
        write(ioc) j1m(i),m1m(i),k1m(i),n1m(i)
        write(ioc) jbm(i),jem(i),kbm(i),kem(i)
        write(ioc) epsp1(i),epsp2(i),epsp3(i),epsp4(i)
        end do
        write(ioc) epsm1,epsm2,epsm3,epsm4
        close(unit=ioc)
        return
        end


c===============================c
        subroutine opdp3d
c===============================c
        use mdugrd;use mduvar
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
c-------------------------------c
        open(unit=ioc,file='grid-t.p3d',form='unformatted')
        write(ioc) imx-ip0+1,jmx-jp0+1,kmx-kp0+1
        write(ioc)
     *  (((sngl(xm1(i)),i=ip0,imx),j=jp0,jmx),k=kp0,kmx),
     *  (((sngl(ym1(j)),i=ip0,imx),j=jp0,jmx),k=kp0,kmx),
     *  (((sngl(zm1(k)),i=ip0,imx),j=jp0,jmx),k=kp0,kmx),
     *  (((ibp(j,k),i=ip0,imx),j=jp0,jmx),k=kp0,kmx)
        close(unit=ioc)
        open(unit=ioc,file='func-t.p3d',form='unformatted')
        write(ioc) imx-it0+1,jmx-jt0+1,kmx-kt0+1,1
        write(ioc)
     *  (((sngl(ti(i,j,k)),i=it0,imx),j=jt0,jmx),k=kt0,kmx)
        close(unit=ioc)

        open(unit=ioc,file='grid-u.p3d',form='unformatted')
        write(ioc) imx-iu0+1,jmx-ju0+1,kmx-ku0+1
        write(ioc)
     *  (((sngl(xu1(i)),i=iu0,imx),j=ju0,jmx),k=ku0,kmx),
     *  (((sngl(ym1(j)),i=iu0,imx),j=ju0,jmx),k=ku0,kmx),
     *  (((sngl(zm1(k)),i=iu0,imx),j=ju0,jmx),k=ku0,kmx),
     *  (((ibu(j,k),i=iu0,imx),j=ju0,jmx),k=ku0,kmx)
        close(unit=ioc)
        open(unit=ioc,file='func-u.p3d',form='unformatted')
        write(ioc) imx-iu0+1,jmx-ju0+1,kmx-ku0+1,1
        write(ioc)
     *  (((sngl(ui(i,j,k)),i=iu0,imx),j=ju0,jmx),k=ku0,kmx)
        close(unit=ioc)

c        open(unit=ioc,file='grid-v.p3d',form='unformatted')
c        write(ioc) imx-iv0+1,jmx-jv0+1,kmx-kv0+1
c        write(ioc)
c     *  (((sngl(xm1(i)),i=iv0,imx),j=jv0,jmx),k=kv0,kmx),
c     *  (((sngl(yv1(j)),i=iv0,imx),j=jv0,jmx),k=kv0,kmx),
c     *  (((sngl(zm1(k)),i=iv0,imx),j=jv0,jmx),k=kv0,kmx),
c     *  (((ibv(j,k),i=iv0,imx),j=jv0,jmx),k=kv0,kmx)
c        close(unit=ioc)
c        open(unit=ioc,file='func-v.p3d',form='unformatted')
c        write(ioc) imx-iv0+1,jmx-jv0+1,kmx-kv0+1,1
c        write(ioc)
c     *  (((sngl(vi1(i,j,k)),i=iv0,imx),j=jv0,jmx),k=kv0,kmx)
c        close(unit=ioc)

c        open(unit=ioc,file='grid-w.p3d',form='unformatted')
c        write(ioc) imx-iw0+1,jmx-jw0+1,kmx-kw0+1
c        write(ioc)
c     *  (((sngl(xm1(i)),i=iw0,imx),j=jw0,jmx),k=kw0,kmx),
c     *  (((sngl(ym1(j)),i=iw0,imx),j=jw0,jmx),k=kw0,kmx),
c     *  (((sngl(zw1(k)),i=iw0,imx),j=jw0,jmx),k=kw0,kmx),
c     *  (((ibw(j,k),i=iw0,imx),j=jw0,jmx),k=kw0,kmx)
c        close(unit=ioc)
c        open(unit=ioc,file='func-w.p3d',form='unformatted')
c        write(ioc) imx-iw0+1,jmx-jw0+1,kmx-kw0+1,1
c        write(ioc)
c     *  (((sngl(wi1(i,j,k)),i=iw0,imx),j=jw0,jmx),k=kw0,kmx)
c        close(unit=ioc)
        return;end
        

c===============================c
        subroutine opdufm
c===============================c
        use mduvar
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
c-------------------------------c
100     format()
120     format(7x,23(1h*),6x,'INT of uu1',6x,23(1h*))
c-------------------------------c
        open(unit=ioc,file='dat2.fmt',form='formatted')
c-------the initial conditions ui1-----------c
        write(ioc,100)
        write(ioc,120)
        write(ioc,100)
        call wtffmt(ui,modeui,ioc,iw0,jw0,kw0,imx,jmx,kmx,imx,jmx,kmx)
        close(unit=ioc)
        return
        end
c===============================c
        subroutine opdvfm
c===============================c
        use mduvar
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
c-------------------------------c
100     format()
120     format(7x,23(1h*),6x,'INT of vv1',6x,23(1h*))
c-------------------------------c
        open(unit=ioc,file='dat3.fmt',form='formatted')
c-------the initial conditions vi1-------c
        write(ioc,100)
        write(ioc,120)
        write(ioc,100)
        call wtffmt(vi,modevi,ioc,iw0,jw0,kw0,imx,jmx,kmx,imx,jmx,kmx)
        close(unit=ioc)
        return
        end
c===============================c
        subroutine opdwfm
c===============================c
        use mduvar
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
c-------------------------------c
100     format()
120     format(7x,23(1h*),6x,'INT of ww1',6x,23(1h*))
c-------------------------------c
        open(unit=ioc,file='dat4.fmt',form='formatted')
c-------the initial conditions of wi1-------c
        write(ioc,100)
        write(ioc,120)
        write(ioc,100)
        call wtffmt(wi,modewi,ioc,iw0,jw0,kw0,imx,jmx,kmx,imx,jmx,kmx)
        close(unit=ioc)
        return
        end

c===============================c
        subroutine opdff0
c===============================c
        use mduvar
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
c-------------------------------c
        open(unit=ioc,file='data.ff0',form='unformatted')
        write(ioc) (((ui(i,j,k),i=iu0,imx),j=ju0,jmx),k=ku0,kmx)
        write(ioc) (((vi(i,j,k),i=iv0,imx),j=jv0,jmx),k=kv0,kmx)
        write(ioc) (((wi(i,j,k),i=iw0,imx),j=jw0,jmx),k=kw0,kmx)
        write(ioc) (((ti(i,j,k),i=it0,imx),j=jt0,jmx),k=kt0,kmx)
        close(unit=ioc)
        return
        end

        subroutine wtcfmt(x,cx,ioc,i0,i1,ni)
        double precision x(ni)
        character cx*6
c-------------------------------c
100     format()
101     format('*',76(1h-),'*')
111     format(2x,'i =',4x,5(i4,8x),i4)
112     format(a6,1p6d12.4)
c-------------------------------c
        write(ioc,100)
        write(ioc,101)
        iend=i0-1
10      if(iend.eq.i1) goto 11
        ibeg=iend+1
        iend=iend+6
        iend=min0(iend,i1)
        write(ioc,100)
        write(ioc,111)(i,i=ibeg,iend)
        write(ioc,112) cx,(x(i),i=ibeg,iend)
        goto 10
11      write(ioc,101)
        return
        end
        subroutine wtcufm(x,ioc,i0,i1,ni)
        double precision x(ni)
        write(ioc)(x(i),i=i0,i1)
        return
        end

        subroutine wtffmt(f,mode,ioc,i0,j0,k0,i1,j1,k1,ni,nj,nk)
        double precision f(ni,nj,nk)
c------------------------------------------c
100     format()
121     format(7x,'i =',i4)
122     format(7h    j =,i5,1x,6i10)
123     format(4h   k)
124     format(i4,1x,1p7d10.2)

131     format(7x,'j =',i4)
132     format(7h    k =,i5,1x,6i10)
133     format(4h   i)
134     format(i4,1x,1p7d10.2)

141     format(7x,'k =',i4)
142     format(7h    i =,i5,1x,6i10)
143     format(4h   j)
144     format(i4,1x,1p7d10.2)
c------------------------------------------c
        ifl=i0+i1
        jfl=j0+j1
        kfl=k0+k1
        if(mode.eq.1) then
        do 11 i=i0,i1
        write(ioc,121) i
        write(ioc,100)
        jbeg=j0-7
12      jbeg=jbeg+7
        jend=jbeg+6
        jend=min0(jend,j1)
        do 13 kk=k0,k1
        k=kfl-kk
        write(ioc,124) k,(f(i,j,k),j=jbeg,jend)
13      continue
        write(ioc,123)
        write(ioc,122)(j,j=jbeg,jend)
        write(ioc,100)
        if(jend.lt.j1) goto 12
11      continue
        else if(mode.eq.2) then
        do 14 j=j0,j1
        write(ioc,131) j
        write(ioc,100)
        kbeg=k0-7
15      kbeg=kbeg+7
        kend=kbeg+6
        kend=min0(kend,k1)
        do 16 ii=i0,i1
        i=ifl-ii
        write(ioc,134) i,(f(i,j,k),k=kbeg,kend)
16      continue
        write(ioc,133)
        write(ioc,132)(k,k=kbeg,kend)
        write(ioc,100)
        if(kend.lt.k1) goto 15
14      continue
        else if(mode.eq.3) then
        do 17 k=k0,k1
        write(ioc,141) k
        write(ioc,100)
        ibeg=i0-7
18      ibeg=ibeg+7
        iend=ibeg+6
        iend=min0(iend,i1)
        do 19 jj=j0,j1
        j=jfl-jj
        write(ioc,144) j,(f(i,j,k),i=ibeg,iend)
19      continue
        write(ioc,143)
        write(ioc,142)(i,i=ibeg,iend)
        write(ioc,100)
        if(iend.lt.i1) goto 18
17      continue
        end if
        return
        end
        subroutine wtfufm(f,ioc,i0,j0,k0,i1,j1,k1,ni,nj,nk)
        double precision f(ni,nj,nk)
        write(ioc) (((f(i,j,k),i=i0,i1),j=j0,j1),k=k0,k1)
        return
        end

        subroutine chk2d1(f,fnm,nch,j0,k0,j1,k1,nj,nk)
        include 'table.prc'
        dimension f(nj,nk)
        character fnm*8
c------------------------------------------c
100     format()
122     format(7h    j =,i5,1x,6i10)
123     format(4h   k)
124     format(i4,1x,1p7d10.2)
c------------------------------------------c
        jfl=j0+j1
        kfl=k0+k1
        open(unit=nch,file=fnm)
        write(nch,100)
        jbeg=j0-7
12      jbeg=jbeg+7
        jend=jbeg+6
        jend=min0(jend,j1)
        do 13 kk=k0,k1
        k=kfl-kk
        write(nch,124) k,(f(j,k),j=jbeg,jend)
13      continue
        write(nch,123)
        write(nch,122)(j,j=jbeg,jend)
        write(nch,100)
        if(jend.lt.j1) goto 12
        close(nch)
        return
        end
        subroutine chk2d2(f,j0,k0,j1,k1,nj,nk)
        include 'table.prc'
        dimension f(nj,nk)
c------------------------------------------c
100     format()
122     format(7h    j =,i5,1x,6i10)
123     format(4h   k)
124     format(i4,1x,1p7d10.2)
c------------------------------------------c
        jfl=j0+j1
        kfl=k0+k1
        write(*,100)
        jbeg=j0-7
12      jbeg=jbeg+7
        jend=jbeg+6
        jend=min0(jend,j1)
        do 13 kk=k0,k1
        k=kfl-kk
        write(*,124) k,(f(j,k),j=jbeg,jend)
13      continue
        write(*,123)
        write(*,122)(j,j=jbeg,jend)
        write(*,100)
        if(jend.lt.j1) goto 12
        return
        end

        subroutine chk3d1(ichk,f,fnm,nch,i,j0,k0,j1,k1,ni,nj,nk)
        include 'table.prc'
        dimension f(ni,nj,nk)
        character fnm*8
c------------------------------------------c
100     format()
121     format(7x,'i =',i4)
122     format(7h    j =,i5,1x,6i10)
123     format(4h   k)
124     format(i4,1x,1p7d10.2)
c------------------------------------------c
        jfl=j0+j1
        kfl=k0+k1
        if(i.eq.ichk) then
        open(unit=nch,file=fnm)
        write(nch,121) i
        write(nch,100)
        jbeg=j0-7
12      jbeg=jbeg+7
        jend=jbeg+6
        jend=min0(jend,j1)
        do 13 kk=k0,k1
        k=kfl-kk
        write(nch,124) k,(f(i,j,k),j=jbeg,jend)
13      continue
        write(nch,123)
        write(nch,122)(j,j=jbeg,jend)
        write(nch,100)
        if(jend.lt.j1) goto 12
        end if
        return
        end
        subroutine chk3d2(ichk,f,i,j0,k0,j1,k1,ni,nj,nk)
        include 'table.prc'
        dimension f(ni,nj,nk)
c------------------------------------------c
100     format()
121     format(7x,'i =',i4)
122     format(7h    j =,i5,1x,6i10)
123     format(4h   k)
124     format(i4,1x,1p7d10.2)
c------------------------------------------c
        jfl=j0+j1
        kfl=k0+k1
        if(i.eq.ichk) then
        write(*,121) i
        write(*,100)
        jbeg=j0-7
12      jbeg=jbeg+7
        jend=jbeg+6
        jend=min0(jend,j1)
        do 13 kk=k0,k1
        k=kfl-kk
        write(*,124) k,(f(i,j,k),j=jbeg,jend)
13      continue
        write(*,123)
        write(*,122)(j,j=jbeg,jend)
        write(*,100)
        if(jend.lt.j1) goto 12
        end if
        close(10)
        return
        end
        
        
        subroutine output1
c===============================c
        use mdugrd;use mduvar;use mdumgv;use mdupnt
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
        open(9,file='simple_re400_pre_correct_alpha_loop100.dat',
     *           form='formatted')
	write(9,*) 'TITLE    ="Plot3D DataSet"'
        write(9,*) 'VARIABLES = "y" "z" "delta_p" "p_corrected" "v" "w"'
        write(9,*) 'DATASETAUXDATA Common.SpeedOfSound="1.0"'
        write(9,*) 'DATASETAUXDATA Common.VectorVarsAreVelocity="FALSE"'
        write(9,*) 'ZONE T="Zone-original grid"'
        write(9,*) 'STRANDID=0, SOLUTIONTIME=0'
        write(9,*) 'I=514, J=514, K=1, ZONETYPE=Ordered'
        write(9,*) 'DATAPACKING=POINT'
        write(9,*) 'DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE)'
	do k=1,kmx
	do j=1,jmx
        write(9,*) ym1(j),zm1(k),pre1(j,k),pi1(j,k),vi1(j,k),wi1(j,k)
        enddo
        enddo
        close(9)
        end
        
        
