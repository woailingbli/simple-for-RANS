c===============================c
        subroutine wlamin
c===============================c
        use mdugrd;use mduvar;use mdumgv;use mdupnt
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
c-------------------------------c
        dimension Fjn(jmx,kmx),Fjs(jmx,kmx),Fkt(jmx,kmx),Fkb(jmx,kmx)
        dimension Djn(jmx,kmx),Djs(jmx,kmx),Dkt(jmx,kmx),Dkb(jmx,kmx)
        dimension vpwpn(jmx,kmx),vpwps(jmx,kmx)
        
        if(levmgd.eq.9) then
        ibp9_w=ibw;ibpm0=>ibp9_w
         ffm0=>ff9_w ;appm0=>app9_w;bppm0=>bpp9_w;resm0=>res9_w
        ajsm0=>ajs9_w;ajnm0=>ajn9_w;akbm0=>akb9_w;aktm0=>akt9_w
        end if
        if(levmgd.eq.8) then
        ibp8_w=ibw;ibpm0=>ibp8_w
         ffm0=>ff8_w ;appm0=>app8_w;bppm0=>bpp8_w;resm0=>res8_w
        ajsm0=>ajs8_w;ajnm0=>ajn8_w;akbm0=>akb8_w;aktm0=>akt8_w
        end if
c-------------------------------c
        j2=jw0+1;k2=kw0+1;m2=jmx-1;n2=kmx-1

        do j=j2,m2;do k=k2,n2
        if(ibpm0(j,k).eq.1) then
                               
        Fjn(j,k) = dzw(k)*(vi1(j+1,k)*fz2(k)+vi1(j+1,k-1)*fz1(k))
        Fjs(j,k) = dzw(k)*(vi1(j,k)*fz2(k)+vi1(j,k-1)*fz1(k))
        Fkt(j,k) = dym(j)*(wi1(j,k)+wi1(j,k+1))/2
        Fkb(j,k) = dym(j)*(wi1(j,k)+wi1(j,k-1))/2
        
        
        Djn(j,k) = dzw(k)/(reno*dyv(j+1))
        Djs(j,k) = dzw(k)/(reno*dyv(j))
        Dkt(j,k) = dym(j)/(reno*dzm(k))
        Dkb(j,k) = dym(j)/(reno*dzm(k-1))
        
        ajsm0(j,k) = Djs(j,k) + max(0.0d-00,Fjs(j,k))
        ajnm0(j,k) = Djn(j,k) + max(0.0d-00,-Fjn(j,k))
        akbm0(j,k) = Dkb(j,k) + max(0.0d-00,Fkb(j,k))
        aktm0(j,k) = Dkt(j,k) + max(0.0d-00,-Fkt(j,k))

!        ajsm0(j,k) = Djs(j,k) + 5.0d-01*Fjs(j,k)
!        ajnm0(j,k) = Djn(j,k) - 5.0d-01*Fjn(j,k)
!        akbm0(j,k) = Dkb(j,k) + 5.0d-01*Fkb(j,k)
!        aktm0(j,k) = Dkt(j,k) - 5.0d-01*Fkt(j,k)
        
        else
        ajsm0(j,k)=0.0d+00;ajnm0(j,k)=0.0d+00
        akbm0(j,k)=0.0d+00;aktm0(j,k)=0.0d+00
        end if
        end do;end do

        do j=j2,m2;do k=k2,n2
        if(ibpm0(j,k).eq.1) then
        vpwpn(j,k) = (vpwp(j,k-1)*fz1(k)+vpwp(j,k)*fz2(k))*fy1(j+1)
     *             +(vpwp(j+1,k-1)*fz1(k)+vpwp(j+1,k)*fz2(k))*fy2(j+1)
        vpwps(j,k) = (vpwp(j-1,k-1)*fz1(k)+vpwp(j-1,k)*fz2(k))*fy1(j)
     *             +(vpwp(j,k-1)*fz1(k)+vpwp(j,k)*fz2(k))*fy2(j)
        
        appm0(j,k)=ajsm0(j,k)+ajnm0(j,k)+akbm0(j,k)+aktm0(j,k)
        bppm0(j,k)=-(pi1(j,k)-pi1(j,k-1))*dym(j) /1.0d-05
     *            -(wpwp(j,k)-wpwp(j,k-1))*dym(j)
     *            -(vpwpn(j,k)-vpwps(j,k))*dzw(k)
     
        appm0(j,k)=appm0(j,k)/alpha_w
        bppm0(j,k)=bppm0(j,k)+(1.0d+00-alpha_w)*appm0(j,k)*wi1(j,k)
        else
        appm0(j,k)=1.0d+00;bppm0(j,k)=0.0d+00
        end if
        end do;end do
        !$omp end sections
        !$omp end parallel
c-------outer wall boundary-------c
        j=j2
        do k=k2,n2
        appm0(j,k)=appm0(j,k)+ajsm0(j,k)
        bppm0(j,k)=bppm0(j,k)-ajsm0(j,k)*wi1(j,k)

        enddo
        ffm0(j-1,k2:n2)=wbo;ajsm0(j,k2:n2)=0.0d+00
        j=m2
        do k=k2,n2
        appm0(j,k)=appm0(j,k)+ajnm0(j,k)
        bppm0(j,k)=bppm0(j,k)-ajnm0(j,k)*wi1(j,k)
        enddo
        ffm0(j+1,k2:n2)=wbo;ajnm0(j,k2:n2)=0.0d+00
!        
!        k=k2
!        do j=j2,m2
!        appm0(j,k)=appm0(j,k)+akbm0(j,k)
!        bppm0(j,k)=bppm0(j,k)-akbm0(j,k)*wi1(j,k)
!        enddo
!        ffm0(j2:m2,k-1)=vbo;akbm0(j2:m2,k)=0.0d+00
!        k=n2
!        do j=j2,m2
!        appm0(j,k)=appm0(j,k)+aktm0(j,k)
!        bppm0(j,k)=bppm0(j,k)-aktm0(j,k)*wi1(j,k)
!        
!!        enddo
!        ffm0(j2:m2,k+1)=vbo;aktm0(j2:m2,k)=0.0d+00
c-------inner wall boundary-------c
        j=je1
        do k=kb1,ke1-1
        appm0(j,k)=appm0(j,k)+ajsm0(j,k)
        bppm0(j,k)=bppm0(j,k)-ajsm0(j,k)*wi1(j,k)
        enddo
        ffm0(j-1,kb1:ke1-1)=wbi;ajsm0(j,kb1:ke1-1)=0.0d+00
        
        j=jb1-1
        do k=kb1,ke1-1
        appm0(j,k)=appm0(j,k)+ajnm0(j,k)
        bppm0(j,k)=bppm0(j,k)-ajnm0(j,k)*wi1(j,k)
        enddo
        ffm0(j,kb1:ke1-1)=wbi;ajnm0(j,kb1:ke1-1)=0.0d+00
!
!        k=ke1
!        do j=jb1,je1-1
!        appm0(j,k)=appm0(j,k)+akbm0(j,k)
!        bppm0(j,k)=bppm0(j,k)-akbm0(j,k)*wi1(j,k)
!        enddo
!        ffm0(jb1:je1-1,k-1)=vbi;akbm0(jb1:je1-1,k)=0.0d+00
!        
!        k=kb1-1
!        do j=jb1,je1-1
!        appm0(j,k)=appm0(j,k)+aktm0(j,k)
!        bppm0(j,k)=bppm0(j,k)-akbm0(j,k)*wi1(j,k)        
!        enddo
!        ffm0(jb1:je1-1,k)=vbi;aktm0(jb1:je1-1,k)=0.0d+00
c-------solving laminar profile-------c
        do j=j2,m2;do k=k2,n2
        apc2(j,k)=appm0(j,k)
        end do;end do

        call amgcof_w(levmgd)
        call solvew(levmgd,modeui,ressumb,ressume)

        do j=ju0,jmx;do k=ku0,kmx
        ww1(j,k)=ffm0(j,k)
        end do;end do
        
!        open(9,file='pointer_simple_w_re400_ml_stress.dat',
!     *           form='formatted')
!	write(9,*) 'TITLE    ="Plot3D DataSet"'
!        write(9,*) 'VARIABLES = "y" "z" "w"'
!        write(9,*) 'DATASETAUXDATA Common.SpeedOfSound="1.0"'
!        write(9,*) 'DATASETAUXDATA Common.VectorVarsAreVelocity="FALSE"'
!        write(9,*) 'ZONE T="Zone-original grid"'
!        write(9,*) 'STRANDID=0, SOLUTIONTIME=0'
!        write(9,*) 'I=514, J=514, K=1, ZONETYPE=Ordered'
!        write(9,*) 'DATAPACKING=POINT'
!        write(9,*) 'DT=(SINGLE SINGLE SINGLE SINGLE)'
!	do k=1,kmx
!	    do j=1,jmx
!                write(9,*) ym1(j),zw1(k),ww1(j,k)
!            enddo
!        enddo
!        close(9)
!        
        return
        end
        
c================================c
        subroutine solve0_w(lmgd,mode,ressumb,ressume)
c================================c
        use mduvar;use mdumgv;use mdupnt;use mduwrk
        use omp_lib

        include 'table.prc'
        include 'table.gd1'
        include 'table.gdm'
c--------------------------------c
        allocate(fw(m1mw(lmgd),n1mw(lmgd)))
        allocate(pk(n1mw(lmgd)),qk(n1mw(lmgd)))
        allocate(pj(m1mw(lmgd)),qj(m1mw(lmgd)))
c--------------------------------c
        nullify(ibpm0,ffm0,appm0,bppm0,resm0,ajsm0,ajnm0,akbm0,aktm0)
        if(lmgd.eq.9) then
        ibpm0=>ibp9_w
         ffm0=>ff9_w ;appm0=>app9_w;bppm0=>bpp9_w;resm0=>res9_w
        ajsm0=>ajs9_w;ajnm0=>ajn9_w;akbm0=>akb9_w;aktm0=>akt9_w
        end if
        if(lmgd.eq.8) then
        ibpm0=>ibp8_w
         ffm0=>ff8_w ;appm0=>app8_w;bppm0=>bpp8_w;resm0=>res8_w
        ajsm0=>ajs8_w;ajnm0=>ajn8_w;akbm0=>akb8_w;aktm0=>akt8_w
        end if
        if(lmgd.eq.7) then
        ibpm0=>ibp7_w
         ffm0=>ff7_w ;appm0=>app7_w;bppm0=>bpp7_w;resm0=>res7_w
        ajsm0=>ajs7_w;ajnm0=>ajn7_w;akbm0=>akb7_w;aktm0=>akt7_w
        end if
        if(lmgd.eq.6) then
        ibpm0=>ibp6_w
         ffm0=>ff6_w ;appm0=>app6_w;bppm0=>bpp6_w;resm0=>res6_w
        ajsm0=>ajs6_w;ajnm0=>ajn6_w;akbm0=>akb6_w;aktm0=>akt6_w
        end if
        if(lmgd.eq.5) then
        ibpm0=>ibp5_w
         ffm0=>ff5_w ;appm0=>app5_w;bppm0=>bpp5_w;resm0=>res5_w
        ajsm0=>ajs5_w;ajnm0=>ajn5_w;akbm0=>akb5_w;aktm0=>akt5_w
        end if
        if(lmgd.eq.4) then
        ibpm0=>ibp4_w
         ffm0=>ff4_w ;appm0=>app4_w;bppm0=>bpp4_w;resm0=>res4_w
        ajsm0=>ajs4_w;ajnm0=>ajn4_w;akbm0=>akb4_w;aktm0=>akt4_w
        end if
        if(lmgd.eq.3) then
        ibpm0=>ibp3_w
         ffm0=>ff3_w ;appm0=>app3_w;bppm0=>bpp3_w;resm0=>res3_w
        ajsm0=>ajs3_w;ajnm0=>ajn3_w;akbm0=>akb3_w;aktm0=>akt3_w
        end if
        if(lmgd.eq.2) then
        ibpm0=>ibp2_w
         ffm0=>ff2_w ;appm0=>app2_w;bppm0=>bpp2_w;resm0=>res2_w
        ajsm0=>ajs2_w;ajnm0=>ajn2_w;akbm0=>akb2_w;aktm0=>akt2_w
        end if
        if(lmgd.eq.1) then
        ibpm0=>ibp1_w
         ffm0=>ff1_w ;appm0=>app1_w;bppm0=>bpp1_w;resm0=>res1_w
        ajsm0=>ajs1_w;ajnm0=>ajn1_w;akbm0=>akb1_w;aktm0=>akt1_w
        end if
c--------------------------------c
        j2=j1mw(lmgd)+1;m2=m1mw(lmgd)-1
        k2=k1mw(lmgd)+1;n2=n1mw(lmgd)-1
        j3=j2+1;k3=k2+1;m3=m2-1;n3=n2-1
        iter=0
c--------------------------------c
        ressumb=0.0d+00
        do j=j2,m2;do k=k2,n2
        if(ibpm0(j,k).ne.0) then
        resm0(j,k)=bppm0(j,k)-appm0(j,k)*ffm0(j,k)+
     *  ajsm0(j,k)*ffm0(j-1,k)+ajnm0(j,k)*ffm0(j+1,k)+
     *  akbm0(j,k)*ffm0(j,k-1)+aktm0(j,k)*ffm0(j,k+1)
        ressumb=ressumb+dabs(resm0(j,k))
        end if
        end do;end do
        ressumw=ressumb
c        write(*,10) ressumb
10      format(' *',8x,'initial residual ressum=',1pe10.3,8x,'*')
c--------------------------------c
100     iter=iter+1;fw(:,:)=ffm0(:,:)
c--------------------------------c
        !$omp parallel
        !$omp sections private(j)
        !$omp do
        do j=j2,m2
        denom=1.0d+00/appm0(j,k2)
        pk(k2)=aktm0(j,k2)*denom
        temp=bppm0(j,k2)+
     *  ajsm0(j,k2)*ffm0(j-1,k2)+ajnm0(j,k2)*ffm0(j+1,k2)
        qk(k2)=temp*denom
        !$omp end do
        !$omp end parallel   
        
        !$omp parallel
        !$omp sections private(k)
        !$omp do
        do k=k3,n2
        denom=1.0d+00/(appm0(j,k)-pk(k-1)*akbm0(j,k))
        pk(k)=aktm0(j,k)*denom
        temp=bppm0(j,k)+ajsm0(j,k)*ffm0(j-1,k)+ajnm0(j,k)*ffm0(j+1,k)
        qk(k)=(temp+akbm0(j,k)*qk(k-1))*denom
        end do
        ffm0(j,n2)=qk(n2)
        do k=n3,k2,-1;ffm0(j,k)=ffm0(j,k+1)*pk(k)+qk(k);end do
        end do
        !$omp end do
        !$omp end parallel   
        
c--------------------------------c
        !$omp parallel
        !$omp sections private(k)
        !$omp do
        do k=k2,n2
        denom=1.0d+00/appm0(j2,k)
        pj(j2)=ajnm0(j2,k)*denom
        temp=bppm0(j2,k)+
     *  akbm0(j2,k)*ffm0(j2,k-1)+aktm0(j2,k)*ffm0(j2,k+1)
        qj(j2)=temp*denom
        !$omp end do
        !$omp end parallel 
        
        !$omp parallel
        !$omp sections private(j)
        !$omp do
        do j=j3,m2
        denom=1.0d+00/(appm0(j,k)-pj(j-1)*ajsm0(j,k))
        pj(j)=ajnm0(j,k)*denom
        temp=bppm0(j,k)+akbm0(j,k)*ffm0(j,k-1)+aktm0(j,k)*ffm0(j,k+1)
        qj(j)=(temp+ajsm0(j,k)*qj(j-1))*denom
        end do
        !$omp end do
        !$omp end parallel   
        
        ffm0(m2,k)=qj(m2)
        do j=m3,j2,-1;ffm0(j,k)=ffm0(j+1,k)*pj(j)+qj(j);end do
        end do
c--------------------------------c
        ressume=0.0d+00
        do j=j2,m2;do k=k2,n2
        if(ibpm0(j,k).ne.0) then
        resm0(j,k)=bppm0(j,k)-appm0(j,k)*ffm0(j,k)+
     *  ajsm0(j,k)*ffm0(j-1,k)+ajnm0(j,k)*ffm0(j+1,k)+
     *  akbm0(j,k)*ffm0(j,k-1)+aktm0(j,k)*ffm0(j,k+1)
        ressume=ressume+dabs(resm0(j,k))
        end if
        end do;end do
c        write(*,200) ressume

        if(ressume.ge.ressumw) then
        ffm0(:,:)=fw(:,:)
c        write(*,201)
        goto 1000
        end if

c-------mode=1: control residual precision
        if(mode.eq.1) then
        if(ressume.le.epsp2(lmgd)) then
c        write(*,202)
        goto 1000
        else
        goto 100
        end if
        end if

c-------mode=2: control residual decreasing level
        if(mode.eq.2) then
        if(ressume.le.epsp2(lmgd)) then
c        write(*,202)
        goto 1000
        end if
        declev=ressume/ressumb
        if(declev.gt.epsp3(lmgd)) then
c        write(*,206)
        goto 100
        else
c        write(*,207)
        goto 1000
        end if
        end if

c-------mode=3: control residual decreasing rate
        if(mode.eq.3) then
        if(ressume.le.epsp2(lmgd)) then
c        write(*,202)
        goto 1000
        end if
        decrat=ressume/ressumw
        if(decrat.le.epsp4(lmgd)) then
        ressumw=ressume
c        write(*,204)
        goto 100
        else
c        write(*,205)
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
1000    deallocate(fw,pj,qj,pk,qk)
        return
        end
                
c================================c
        recursive subroutine solvew(lmgd1,mode,ressumb,ressume)
c================================c
        use mduvar;use mduwrk;use mdumgv;use mdupnt
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
        include 'table.gdm'
c--------------------------------c
        nullify(ibpm0,ffm0,appm0,bppm0,resm0,ajsm0,ajnm0,akbm0,aktm0)
        nullify(ibpml,ffml,appml,bppml,resml,ajsml,ajnml,akbml,aktml)
        if(lmgd1.eq.9) then
        ibpm0=>ibp9_w; ffm0=>ff9_w;bppm0=>bpp9_w;resm0=>res9_w
        appm0=>app9_w;ajsm0=>ajs9_w
        ajnm0=>ajn9_w;akbm0=>akb9_w;aktm0=>akt9_w
        ibpml=>ibp8_w; ffml=>ff8_w;bppml=>bpp8_w;resml=>res8_w
        lmgd2=lmgd1-1
        end if
        if(lmgd1.eq.8) then
        ibpm0=>ibp8_w; ffm0=>ff8_w;bppm0=>bpp8_w;resm0=>res8_w
        appm0=>app8_w;ajsm0=>ajs8_w;ajnm0=>ajn8_w
        akbm0=>akb8_w;aktm0=>akt8_w
        ibpml=>ibp7_w; ffml=>ff7_w;bppml=>bpp7_w
        resml=>res7_w;lmgd2=lmgd1-1
        end if
        if(lmgd1.eq.7) then
        ibpm0=>ibp7_w; ffm0=>ff7_w;bppm0=>bpp7_w;resm0=>res7_w
        appm0=>app7_w;ajsm0=>ajs7_w;ajnm0=>ajn7_w
        akbm0=>akb7_w;aktm0=>akt7_w
        ibpml=>ibp6_w; ffml=>ff6_w;bppml=>bpp6_w
        resml=>res6_w;lmgd2=lmgd1-1
        end if
        if(lmgd1.eq.6) then
        ibpm0=>ibp6_w; ffm0=>ff6_w;bppm0=>bpp6_w;resm0=>res6_w;
        appm0=>app6_w;ajsm0=>ajs6_w;ajnm0=>ajn6_w
        akbm0=>akb6_w;aktm0=>akt6_w
        ibpml=>ibp5_w; ffml=>ff5_w;bppml=>bpp5_w
        resml=>res5_w;lmgd2=lmgd1-1
        end if
        if(lmgd1.eq.5) then
        ibpm0=>ibp5_w; ffm0=>ff5_w;bppm0=>bpp5_w;resm0=>res5_w;
        appm0=>app5_w;ajsm0=>ajs5_w;ajnm0=>ajn5_w
        akbm0=>akb5_w;aktm0=>akt5_w
        ibpml=>ibp4_w; ffml=>ff4_w;bppml=>bpp4_w
        resml=>res4_w;lmgd2=lmgd1-1
        end if
        if(lmgd1.eq.4) then
        ibpm0=>ibp4_w; ffm0=>ff4_w;bppm0=>bpp4_w;resm0=>res4_w;
        appm0=>app4_w;ajsm0=>ajs4_w;ajnm0=>ajn4_w
        akbm0=>akb4_w;aktm0=>akt4_w
        ibpml=>ibp3_w; ffml=>ff3_w;bppml=>bpp3_w
        resml=>res3_w;lmgd2=lmgd1-1
        end if
        if(lmgd1.eq.3) then
        ibpm0=>ibp3_w; ffm0=>ff3_w;bppm0=>bpp3_w;resm0=>res3_w;
        appm0=>app3_w;ajsm0=>ajs3_w;ajnm0=>ajn3_w
        akbm0=>akb3_w;aktm0=>akt3_w
        ibpml=>ibp2_w; ffml=> ff2_w;bppml=>bpp2_w
        resml=>res2_w;lmgd2=lmgd1-1
        end if
        if(lmgd1.eq.2) then
        ibpm0=>ibp2_w; ffm0=> ff2_w;bppm0=>bpp2_w;resm0=>res2_w;
        appm0=>app2_w;ajsm0=>ajs2_w;ajnm0=>ajn2_w
        akbm0=>akb2_w;aktm0=>akt2_w
        ibpml=>ibp1_w; ffml=> ff1_w;bppml=>bpp1_w
        resml=>res1_w;lmgd2=lmgd1-1
        end if        
        if(lmgd1.eq.1) then
        ibpm0=>ibp1_w; ffm0=> ff1_w;bppm0=>bpp1_w;resm0=>res1_w;
        appm0=>app1_w;ajsm0=>ajs1_w;ajnm0=>ajn1_w
        akbm0=>akb1_w;aktm0=>akt1_w
        ibpml=>ibp0_w; ffml=> ff0_w;bppml=>bpp0_w
        resml=>res0_w;lmgd2=lmgd1-1
        end if
c--------------------------------c
        j2=j1mw(lmgd1)+1;m2=m1mw(lmgd1)-1;
        k2=k1mw(lmgd1)+1;n2=n1mw(lmgd1)-1
c--------------------------------c
        ressumb=0.0d+00
        do j=j2,m2;do k=k2,n2
        if(ibpm0(j,k).ne.0) then
        resm0(j,k)=bppm0(j,k)-appm0(j,k)*ffm0(j,k)+
     *  ajsm0(j,k)*ffm0(j-1,k)+ajnm0(j,k)*ffm0(j+1,k)+
     *  akbm0(j,k)*ffm0(j,k-1)+aktm0(j,k)*ffm0(j,k+1)
        ressumb=ressumb+dabs(resm0(j,k))
        end if
        end do;end do
        ressum0(lmgd1)=ressumb
        if(lmgd1.eq.levmgd) write(*,11) lmgd1,ressumb
11      format
     *  (' *',2x,'L-',i2,2x,'initial residual ressum=',1pe10.3,8x,'*')
c--------------------------------c
        mode1=3
100     call solve0_w(lmgd1,mode1,ressumb,ressume)
        ressum1(lmgd1)=ressume
        if(lmgd1.eq.levmgd) write(*,12) lmgd1,ressume
12      format
     *  (' *',2x,'L-',i2,2x,'postLBL residual ressum=',1pe10.3,8x,'*')

c-------mode=1: control residual precision
        if(mode.eq.1) then
        if(ressum1(lmgd1).le.epsp2(lmgd1)) then
c        write(*,201)
c        goto 1000
        return
        else
c        write(*,202)
        goto 200
        end if
        end if

c-------mode=2: control residual decreasing level
        if(mode.eq.2) then
        if(ressum1(lmgd1).le.epsp2(lmgd1)) then
c        write(*,201)
c        goto 1000
        return
        end if
        declev=ressum1(lmgd1)/ressum0(lmgd1)
        if(declev.gt.epsp3(lmgd1)) then
c        write(*,203)
        goto 200
        else
c        write(*,204)
c        goto 1000
        if(lmgd1.lt.levmgd) return
        if(lmgd1.eq.levmgd) goto 1000
        end if
        end if

200     call amgres_w(lmgd1)
        mode2=2
        call solvew(lmgd2,mode2,ressumb,ressume)
!        if(lmgd2.eq.1) call solveb(lmgd2,mode2,ressumb,ressume)
        call amgcor_w(lmgd1)
        goto 100
201     format(' *        total residual satisfied,   return        *')
202     format(' *        residual not satisfied, continuing        *')
203     format(' *        declev not reached,decrat exceeded        *')
204     format(' *        declev satisfied,  go out of L-B-L        *')
1000    ressumb=ressum0(lmgd1)
        ressume=ressum1(lmgd1)
        return
        end


c===============================c
        subroutine amgcof_w(levmgd)
c===============================c
        use mdumgv;use mdupnt
        include 'table.prc'
        include 'table.gdm'
c-------------------------------c
        do lmgd=levmgd,0,-1
        nullify(ibpm0,ffm0,appm0,bppm0,resm0,ajsm0,ajnm0,akbm0,aktm0)
        nullify(ibpml,ffml,appml,bppml,resml,ajsml,ajnml,akbml,aktml)
        if(lmgd.eq.9) then
         ffm0=>ff9_w ;appm0=>app9_w;bppm0=>bpp9_w;ibpm0=>ibp9_w
        ajsm0=>ajs9_w;ajnm0=>ajn9_w;akbm0=>akb9_w;aktm0=>akt9_w
         ffml=>ff8_w;appml=>app8_w;bppml=>bpp8_w;ibpml=>ibp8_w
        ajsml=>ajs8_w;ajnml=>ajn8_w;aktml=>akt8_w;akbml=>akb8_w
        end if
        if(lmgd.eq.8) then
         ffm0=>ff8_w;appm0=>app8_w;bppm0=>bpp8_w;ibpm0=>ibp8_w
        ajsm0=>ajs8_w;ajnm0=>ajn8_w;akbm0=>akb8_w;aktm0=>akt8_w
         ffml=>ff7_w;appml=>app7_w;bppml=>bpp7_w;ibpml=>ibp7_w
        ajsml=>ajs7_w;ajnml=>ajn7_w;aktml=>akt7_w;akbml=>akb7_w
        end if
        if(lmgd.eq.7) then
         ffm0=>ff7_w;appm0=>app7_w;bppm0=>bpp7_w;ibpm0=>ibp7_w
        ajsm0=>ajs7_w;ajnm0=>ajn7_w;akbm0=>akb7_w;aktm0=>akt7_w
         ffml=>ff6_w;appml=>app6_w;bppml=>bpp6_w;ibpml=>ibp6_w
        ajsml=>ajs6_w;ajnml=>ajn6_w;akbml=>akb6_w;aktml=>akt6_w
        end if
        if(lmgd.eq.6) then
         ffm0=>ff6_w;appm0=>app6_w;bppm0=>bpp6_w;ibpm0=>ibp6_w
        ajsm0=>ajs6_w;ajnm0=>ajn6_w;akbm0=>akb6_w;aktm0=>akt6_w
         ffml=>ff5_w;appml=>app5_w;bppml=>bpp5_w;ibpml=>ibp5_w
        ajsml=>ajs5_w;ajnml=>ajn5_w;akbml=>akb5_w;aktml=>akt5_w
        end if
        if(lmgd.eq.5) then
         ffm0=>ff5_w;appm0=>app5_w;bppm0=>bpp5_w;ibpm0=>ibp5_w
        ajsm0=>ajs5_w;ajnm0=>ajn5_w;akbm0=>akb5_w;aktm0=>akt5_w
         ffml=>ff4_w;appml=>app4_w;bppml=>bpp4_w;ibpml=>ibp4_w
        ajsml=>ajs4_w;ajnml=>ajn4_w;akbml=>akb4_w;aktml=>akt4_w
        end if
        if(lmgd.eq.4) then
         ffm0=>ff4_w;appm0=>app4_w;bppm0=>bpp4_w;ibpm0=>ibp4_w
        ajsm0=>ajs4_w;ajnm0=>ajn4_w;akbm0=>akb4_w;aktm0=>akt4_w
         ffml=>ff3_w;appml=>app3_w;bppml=>bpp3_w;ibpml=>ibp3_w
        ajsml=>ajs3_w;ajnml=>ajn3_w;akbml=>akb3_w;aktml=>akt3_w
        end if
        if(lmgd.eq.3) then
         ffm0=>ff3_w;appm0=>app3_w;bppm0=>bpp3_w;ibpm0=>ibp3_w
        ajsm0=>ajs3_w;ajnm0=>ajn3_w;akbm0=>akb3_w;aktm0=>akt3_w
         ffml=>ff2_w;appml=>app2_w;bppml=>bpp2_w;ibpml=>ibp2_w
        ajsml=>ajs2_w;ajnml=>ajn2_w;akbml=>akb2_w;aktml=>akt2_w
        end if
        if(lmgd.eq.2) then
         ffm0=>ff2_w;appm0=>app2_w;bppm0=>bpp2_w;ibpm0=>ibp2_w
        ajsm0=>ajs2_w;ajnm0=>ajn2_w;akbm0=>akb2_w;aktm0=>akt2_w
         ffml=>ff1_w;appml=>app1_w;bppml=>bpp1_w;ibpml=>ibp1_w
        ajsml=>ajs1_w;ajnml=>ajn1_w;akbml=>akb1_w;aktml=>akt1_w
        end if
        if(lmgd.eq.1) then
         ffm0=>ff1_w;appm0=>app1_w;bppm0=>bpp1_w;ibpm0=>ibp1_w
        ajsm0=>ajs1_w;ajnm0=>ajn1_w;akbm0=>akb1_w;aktm0=>akt1_w
         ffml=>ff0_w;appml=>app0_w;bppml=>bpp0_w;ibpml=>ibp0_w
        ajsml=>ajs0_w;ajnml=>ajn0_w;akbml=>akb0_w;aktml=>akt0_w
        end if
        call coff12_w(lmgd,lmgd-1)
        end do
        nullify(ibpm0,ffm0,appm0,bppm0,resm0,ajsm0,ajnm0,akbm0,aktm0)
        nullify(ibpml,ffml,appml,bppml,resml,ajsml,ajnml,akbml,aktml)
        return;end   
     
c-------------------------------c
        subroutine coff12_w(mgd1,mgd2)
c-------------------------------c
        use mdumgv;use mdupnt
        include 'table.prc'
        include 'table.gdm'
c-------------------------------c
        j2=j1mw(mgd2)+1;m2=m1mw(mgd2)-1;k2=k1mw(mgd2)+1;n2=n1mw(mgd2)-1
                
        do j=j2,m2;do k=k2,n2
        jfs=2*j-j2;jfn=jfs+1
        if(ibpm0(jfs,k).eq.0.and.ibpm0(jfn,k).eq.0) then
        ibpml(j,k)=0
        elseif(ibpm0(jfs,k).eq.2.and.ibpm0(jfn,k).eq.2) then
        ibpml(j,k)=2
        else
        ibpml(j,k)=1
        end if
        end do;end do

        do j=j2,m2;do k=k2,n2
        if(ibpml(j,k).eq.0.or.ibpml(j,k).eq.2) then
        appml(j,k)=1.0d+00
        ajsml(j,k)=0.0d+00;ajnml(j,k)=0.0d+00
        akbml(j,k)=0.0d+00;aktml(j,k)=0.0d+00
        else
        jfs=2*j-j2;jfn=jfs+1
        fjfs=0.0d+00;if(ibpm0(jfs,k).ne.0) fjfs=1.0d+00
        fjfn=0.0d+00;if(ibpm0(jfn,k).ne.0) fjfn=1.0d+00
        appml(j,k)=fjfs*appm0(jfs,k)+fjfn*appm0(jfn,k)
     *            -fjfn*ajsm0(jfn,k)-fjfs*ajnm0(jfs,k)
        ajsml(j,k)=fjfs*ajsm0(jfs,k)
        ajnml(j,k)=fjfn*ajnm0(jfn,k)
        akbml(j,k)=fjfs*akbm0(jfs,k)+fjfn*akbm0(jfn,k)
        aktml(j,k)=fjfs*aktm0(jfs,k)+fjfn*aktm0(jfn,k)
        end if
        end do;end do

        return;end
c===============================c

        subroutine amgres_w(mgd1)
c===============================c
        use mdumgv;use mdupnt
        include 'table.prc'
        include 'table.gdm'
c-------------------------------c
        nullify(ibpm0,resm0,ibpml,bppml)
        if(mgd1.eq.9) then
        mgd2=mgd1-1;ibpml=>ibp8_w;bppml=>bpp8_w;
        ibpm0=>ibp9_w;resm0=>res9_w
        end if        
        if(mgd1.eq.8) then
        mgd2=mgd1-1;ibpml=>ibp7_w;bppml=>bpp7_w;
        ibpm0=>ibp8_w;resm0=>res8_w
        end if
        if(mgd1.eq.7) then
        mgd2=mgd1-1;ibpml=>ibp6_w;bppml=>bpp6_w;
        ibpm0=>ibp7_w;resm0=>res7_w
        end if
        if(mgd1.eq.6) then
        mgd2=mgd1-1;ibpml=>ibp5_w;bppml=>bpp5_w;
        ibpm0=>ibp6_w;resm0=>res6_w
        end if
        if(mgd1.eq.5) then
        mgd2=mgd1-1;ibpml=>ibp4_w;bppml=>bpp4_w;
        ibpm0=>ibp5_w;resm0=>res5_w
        end if
        if(mgd1.eq.4) then
        mgd2=mgd1-1;ibpml=>ibp3_w;bppml=>bpp3_w;
        ibpm0=>ibp4_w;resm0=>res4_w
        end if
        if(mgd1.eq.3) then
        mgd2=mgd1-1;ibpml=>ibp2_w;bppml=>bpp2_w;
        ibpm0=>ibp3_w;resm0=>res3_w
        end if
        if(mgd1.eq.2) then
        mgd2=mgd1-1;ibpml=>ibp1_w;bppml=>bpp1_w;
        ibpm0=>ibp2_w;resm0=>res2_w
        end if
        if(mgd1.eq.1) then
        mgd2=mgd1-1;ibpml=>ibp0_w;bppml=>bpp0_w;
        ibpm0=>ibp1_w;resm0=>res1_w
        end if
c-------------------------------c
        j2=j1mw(mgd2)+1;m2=m1mw(mgd2)-1;k2=k1mw(mgd2)+1;n2=n1mw(mgd2)-1
        
        do j=j2,m2;do k=k2,n2
        if(ibpml(j,k).eq.0) then
        bppml(j,k)=0.0d+00
        else
        jfs=2*j-j2;jfn=jfs+1
        fjfs=0.0d+00;if(ibpm0(jfs,k).ne.0) fjfs=1.0d+00
        fjfn=0.0d+00;if(ibpm0(jfn,k).ne.0) fjfn=1.0d+00
        bppml(j,k)=fjfs*resm0(jfs,k)+fjfn*resm0(jfn,k)
        end if
        end do;end do
        return;end
c===============================c
        subroutine amgcor_w(mgd1)
c===============================c
        use mdumgv;use mdupnt
        include 'table.prc'
        include 'table.gdm'
c-------------------------------c
        nullify(ibpm0,ffm0,ibpml,ffml)
        if(mgd1.eq.9) then
        mgd2=mgd1-1;ibpml=>ibp8_w;ffml=>ff8_w;ibpm0=>ibp9_w;ffm0=>ff9_w
        end if        
        if(mgd1.eq.8) then
        mgd2=mgd1-1;ibpml=>ibp7_w;ffml=>ff7_w;ibpm0=>ibp8_w;ffm0=>ff8_w
        end if
        if(mgd1.eq.7) then
        mgd2=mgd1-1;ibpml=>ibp6_w;ffml=>ff6_w;ibpm0=>ibp7_w;ffm0=>ff7_w
        end if
        if(mgd1.eq.6) then
        mgd2=mgd1-1;ibpml=>ibp5_w;ffml=>ff5_w;ibpm0=>ibp6_w;ffm0=>ff6_w
        end if
        if(mgd1.eq.5) then
        mgd2=mgd1-1;ibpml=>ibp4_w;ffml=>ff4_w;ibpm0=>ibp5_w;ffm0=>ff5_w
        end if
        if(mgd1.eq.4) then
        mgd2=mgd1-1;ibpml=>ibp3_w;ffml=>ff3_w;ibpm0=>ibp4_w;ffm0=>ff4_w
        end if
        if(mgd1.eq.3) then
        mgd2=mgd1-1;ibpml=>ibp2_w;ffml=>ff2_w;ibpm0=>ibp3_w;ffm0=>ff3_w
        end if
        if(mgd1.eq.2) then
        mgd2=mgd1-1;ibpml=>ibp1_w;ffml=>ff1_w;ibpm0=>ibp2_w;ffm0=>ff2_w
        end if
        if(mgd1.eq.1) then
        mgd2=mgd1-1;ibpml=>ibp0_w;ffml=>ff0_w;ibpm0=>ibp1_w;ffm0=>ff1_w
        end if
c-------------------------------c
        j2=j1mw(mgd2)+1;m2=m1mw(mgd2)-1;k2=k1mw(mgd2)+1;n2=n1mw(mgd2)-1
        do j=j2,m2;do k=k2,n2
        if(ibpml(j,k).eq.1) then
        jfs=2*j-j2;jfn=jfs+1
        if(ibpm0(jfs,k).ne.0) ffm0(jfs,k)=ffm0(jfs,k)+ffml(j,k)
        if(ibpm0(jfn,k).ne.0) ffm0(jfn,k)=ffm0(jfn,k)+ffml(j,k)
        end if
        end do;end do
        return;end

