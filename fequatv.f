c===============================c
        subroutine vlamin
c===============================c
        use mdugrd;use mduvar;use mdumgv;use mdupnt
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
c-------------------------------c
        dimension Fjn(jmx,kmx),Fjs(jmx,kmx),Fkt(jmx,kmx),Fkb(jmx,kmx)
        dimension Djn(jmx,kmx),Djs(jmx,kmx),Dkt(jmx,kmx),Dkb(jmx,kmx)
        dimension vpwpt(jmx,kmx),vpwpb(jmx,kmx)
        
        if(levmgd.eq.9) then
        ibp9_v=ibv;ibpm0=>ibp9_v
         ffm0=>ff9_v ;appm0=>app9_v;bppm0=>bpp9_v;resm0=>res9_v
        ajsm0=>ajs9_v;ajnm0=>ajn9_v;akbm0=>akb9_v;aktm0=>akt9_v
        end if
        if(levmgd.eq.8) then
        ibp8_v=ibv;ibpm0=>ibp8_v
         ffm0=>ff8_v ;appm0=>app8_v;bppm0=>bpp8_v;resm0=>res8_v
        ajsm0=>ajs8_v;ajnm0=>ajn8_v;akbm0=>akb8_v;aktm0=>akt8_v
        end if
c-------------------------------c
        j2=jv0+1;k2=kv0+1;m2=jmx-1;n2=kmx-1

!        ss = 0.0d-01        
!        do j=2,513
!        do k=2,513                       
!        delta_mass(j,k) = - (vi1(j+1,k)-vi1(j,k)) * dzm(k)
!     *    - (wi1(j,k+1)-wi1(j,k)) * dym(j)
!        ss = ss + dabs(delta_mass(j,k))
!        enddo
!        enddo
        
        do j=j2,m2;do k=k2,n2
        if(ibpm0(j,k).eq.1) then
                               
        Fjn(j,k) = dzm(k)*(vi1(j+1,k)+vi1(j,k))/2
        Fjs(j,k) = dzm(k)*(vi1(j-1,k)+vi1(j,k))/2
        Fkt(j,k) = dyv(j)*(wi1(j,k+1)*fy2(j)+wi1(j-1,k+1)*fy1(j))
        Fkb(j,k) = dyv(j)*(wi1(j,k)*fy2(j)+wi1(j-1,k)*fy1(j))                      

        Djn(j,k) = dzm(k)/(reno*dym(j))
        Djs(j,k) = dzm(k)/(reno*dym(j-1))
        Dkt(j,k) = dyv(j)/(reno*dzw(k+1))
        Dkb(j,k) = dyv(j)/(reno*dzw(k))
  
        ajsm0(j,k) = Djs(j,k) + max(0.0d-00,Fjs(j,k))
        ajnm0(j,k) = Djn(j,k) + max(0.0d-00,-Fjn(j,k))
        akbm0(j,k) = Dkb(j,k) + max(0.0d-00,Fkb(j,k))
        aktm0(j,k) = Dkt(j,k) + max(0.0d-00,-Fkt(j,k))
!        
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
        vpwpt(j,k)=(vpwp(j-1,k)*fy1(j)+vpwp(j,k)*fy2(j))*fz1(k+1)
     &            +(vpwp(j-1,k+1)*fy1(j)+vpwp(j,k+1)*fy2(j))*fz2(k+1)
                        
        vpwpb(j,k)=(vpwp(j-1,k-1)*fy1(j)+vpwp(j,k-1)*fy2(j))*fz1(k)
     &             +(vpwp(j-1,k)*fy1(j)+vpwp(j,k)*fy2(j))*fz2(k)
        
        appm0(j,k)=ajsm0(j,k)+ajnm0(j,k)+akbm0(j,k)+aktm0(j,k)
        bppm0(j,k)=-(pi1(j,k)-pi1(j-1,k))*dzm(k)/1.0d-05
     *                -(vpvp(j,k)-vpvp(j-1,k))*dzm(k)
     *                -(vpwpt(j,k)-vpwpb(j,k))*dyv(j)
        
        appm0(j,k)=appm0(j,k)/alpha_v
        bppm0(j,k)=bppm0(j,k)+(1-alpha_v)*appm0(j,k)*vi1(j,k)
        else
        appm0(j,k)=1.0d+00;bppm0(j,k)=0.0d+00
        end if
        end do;end do

c-------outer wall boundary-------c
!        j=j2
!        do k=k2,n2;appm0(j,k)=appm0(j,k)+ajsm0(j,k);enddo
!        ffm0(j-1,k2:n2)=vbo;ajsm0(j,k2:n2)=0.0d+00
!        j=m2
!        do k=k2,n2;appm0(j,k)=appm0(j,k)+ajnm0(j,k);enddo
!        ffm0(j+1,k2:n2)=vbo;ajnm0(j,k2:n2)=0.0d+00
        k=k2
        do j=j2,m2
        appm0(j,k)=appm0(j,k)+akbm0(j,k)
        bppm0(j,k)=bppm0(j,k)-akbm0(j,k)*vi1(j,k)
        enddo
        ffm0(j2:m2,k-1)=vbo;akbm0(j2:m2,k)=0.0d+00
        k=n2
        do j=j2,m2
        appm0(j,k)=appm0(j,k)+aktm0(j,k)
        bppm0(j,k)=bppm0(j,k)-aktm0(j,k)*vi1(j,k)
        
        enddo
        ffm0(j2:m2,k+1)=vbo;aktm0(j2:m2,k)=0.0d+00
c-------inner wall boundary-------c
!        j=je1
!        do k=kb1,ke1-1
!        appm0(j,k)=appm0(j,k)+ajsm0(j,k)
!        bppm0(j,k)=bppm0(j,k)-ajsm0(j,k)*vi1(j,k)
!        enddo
!        ffm0(j-1,kb1:ke1-1)=vbi;ajsm0(j,kb1:ke1-1)=0.0d+00
!        
!        j=jb1-1
!        do k=kb1,ke1-1
!        appm0(j,k)=appm0(j,k)+ajnm0(j,k)
!        bppm0(j,k)=bppm0(j,k)-ajnm0(j,k)*vi1(j,k)
!        enddo
!        ffm0(j,kb1:ke1-1)=vbi;ajnm0(j,kb1:ke1-1)=0.0d+00

        k=ke1
        do j=jb1,je1-1
        appm0(j,k)=appm0(j,k)+akbm0(j,k)
        bppm0(j,k)=bppm0(j,k)-akbm0(j,k)*vi1(j,k)
        enddo
        ffm0(jb1:je1-1,k-1)=vbi;akbm0(jb1:je1-1,k)=0.0d+00
        
        k=kb1-1
        do j=jb1,je1-1
        appm0(j,k)=appm0(j,k)+aktm0(j,k)
        bppm0(j,k)=bppm0(j,k)-akbm0(j,k)*vi1(j,k)        
        enddo
        ffm0(jb1:je1-1,k)=vbi;aktm0(jb1:je1-1,k)=0.0d+00
c-------solving laminar profile-------c
        do j=j2,m2;do k=k2,n2
        apc1(j,k)=appm0(j,k)
        end do;end do

        call amgcof_v(levmgd)
        call solvev(levmgd,modeui,ressumb,ressume)

        do j=ju0,jmx;do k=ku0,kmx
        vv1(j,k)=ffm0(j,k)
        end do;end do
        
!        open(9,file='pointer_simple_v_re400_ml_stress.dat',
!     *           form='formatted')
!	write(9,*) 'TITLE    ="Plot3D DataSet"'
!        write(9,*) 'VARIABLES = "y" "z" "v"'
!        write(9,*) 'DATASETAUXDATA Common.SpeedOfSound="1.0"'
!        write(9,*) 'DATASETAUXDATA Common.VectorVarsAreVelocity="FALSE"'
!        write(9,*) 'ZONE T="Zone-original grid"'
!        write(9,*) 'STRANDID=0, SOLUTIONTIME=0'
!        write(9,*) 'I=514, J=514, K=1, ZONETYPE=Ordered'
!        write(9,*) 'DATAPACKING=POINT'
!        write(9,*) 'DT=(SINGLE SINGLE SINGLE SINGLE)'
!	do k=1,kmx
!	    do j=1,jmx
!                write(9,*) yv1(j),zm1(k),vv1(j,k)
!            enddo
!        enddo
!        close(9)
        
        return
        end
        
c================================c
        subroutine solve0_v(lmgd,mode,ressumb,ressume)
c================================c
        use mduvar;use mdumgv;use mdupnt;use mduwrk
        use omp_lib

        include 'table.prc'
        include 'table.gd1'
        include 'table.gdm'
c--------------------------------c
        allocate(fw(m1mv(lmgd),n1mv(lmgd)))
        allocate(pk(n1mv(lmgd)),qk(n1mv(lmgd)))
        allocate(pj(m1mv(lmgd)),qj(m1mv(lmgd)))
c--------------------------------c
        nullify(ibpm0,ffm0,appm0,bppm0,resm0,ajsm0,ajnm0,akbm0,aktm0)
        if(lmgd.eq.9) then
        ibpm0=>ibp9_v
         ffm0=>ff9_v ;appm0=>app9_v;bppm0=>bpp9_v;resm0=>res9_v
        ajsm0=>ajs9_v;ajnm0=>ajn9_v;akbm0=>akb9_v;aktm0=>akt9_v
        end if
        if(lmgd.eq.8) then
        ibpm0=>ibp8_v
         ffm0=>ff8_v ;appm0=>app8_v;bppm0=>bpp8_v;resm0=>res8_v
        ajsm0=>ajs8_v;ajnm0=>ajn8_v;akbm0=>akb8_v;aktm0=>akt8_v
        end if
        if(lmgd.eq.7) then
        ibpm0=>ibp7_v
         ffm0=>ff7_v ;appm0=>app7_v;bppm0=>bpp7_v;resm0=>res7_v
        ajsm0=>ajs7_v;ajnm0=>ajn7_v;akbm0=>akb7_v;aktm0=>akt7_v
        end if
        if(lmgd.eq.6) then
        ibpm0=>ibp6_v
         ffm0=>ff6_v ;appm0=>app6_v;bppm0=>bpp6_v;resm0=>res6_v
        ajsm0=>ajs6_v;ajnm0=>ajn6_v;akbm0=>akb6_v;aktm0=>akt6_v
        end if
        if(lmgd.eq.5) then
        ibpm0=>ibp5_v
         ffm0=>ff5_v ;appm0=>app5_v;bppm0=>bpp5_v;resm0=>res5_v
        ajsm0=>ajs5_v;ajnm0=>ajn5_v;akbm0=>akb5_v;aktm0=>akt5_v
        end if
        if(lmgd.eq.4) then
        ibpm0=>ibp4_v
         ffm0=>ff4_v ;appm0=>app4_v;bppm0=>bpp4_v;resm0=>res4_v
        ajsm0=>ajs4_v;ajnm0=>ajn4_v;akbm0=>akb4_v;aktm0=>akt4_v
        end if
        if(lmgd.eq.3) then
        ibpm0=>ibp3_v
         ffm0=>ff3_v ;appm0=>app3_v;bppm0=>bpp3_v;resm0=>res3_v
        ajsm0=>ajs3_v;ajnm0=>ajn3_v;akbm0=>akb3_v;aktm0=>akt3_v
        end if
        if(lmgd.eq.2) then
        ibpm0=>ibp2_v
         ffm0=>ff2_v ;appm0=>app2_v;bppm0=>bpp2_v;resm0=>res2_v
        ajsm0=>ajs2_v;ajnm0=>ajn2_v;akbm0=>akb2_v;aktm0=>akt2_v
        end if
        if(lmgd.eq.1) then
        ibpm0=>ibp1_v
         ffm0=>ff1_v ;appm0=>app1_v;bppm0=>bpp1_v;resm0=>res1_v
        ajsm0=>ajs1_v;ajnm0=>ajn1_v;akbm0=>akb1_v;aktm0=>akt1_v
        end if
c--------------------------------c
        j2=j1mv(lmgd)+1;m2=m1mv(lmgd)-1
        k2=k1mv(lmgd)+1;n2=n1mv(lmgd)-1
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

        do j=j2,m2
        denom=1.0d+00/appm0(j,k2)
        pk(k2)=aktm0(j,k2)*denom
        temp=bppm0(j,k2)+
     *  ajsm0(j,k2)*ffm0(j-1,k2)+ajnm0(j,k2)*ffm0(j+1,k2)
        qk(k2)=temp*denom

        do k=k3,n2
        denom=1.0d+00/(appm0(j,k)-pk(k-1)*akbm0(j,k))
        pk(k)=aktm0(j,k)*denom
        temp=bppm0(j,k)+ajsm0(j,k)*ffm0(j-1,k)+ajnm0(j,k)*ffm0(j+1,k)
        qk(k)=(temp+akbm0(j,k)*qk(k-1))*denom
        end do
        ffm0(j,n2)=qk(n2)
        do k=n3,k2,-1;ffm0(j,k)=ffm0(j,k+1)*pk(k)+qk(k);end do
        end do

c--------------------------------c

        do k=k2,n2
        denom=1.0d+00/appm0(j2,k)
        pj(j2)=ajnm0(j2,k)*denom
        temp=bppm0(j2,k)+
     *  akbm0(j2,k)*ffm0(j2,k-1)+aktm0(j2,k)*ffm0(j2,k+1)
        qj(j2)=temp*denom

        do j=j3,m2
        denom=1.0d+00/(appm0(j,k)-pj(j-1)*ajsm0(j,k))
        pj(j)=ajnm0(j,k)*denom
        temp=bppm0(j,k)+akbm0(j,k)*ffm0(j,k-1)+aktm0(j,k)*ffm0(j,k+1)
        qj(j)=(temp+ajsm0(j,k)*qj(j-1))*denom
        end do

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
        recursive subroutine solvev(lmgd1,mode,ressumb,ressume)
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
        ibpm0=>ibp9_v; ffm0=>ff9_v;bppm0=>bpp9_v;resm0=>res9_v
        appm0=>app9_v;ajsm0=>ajs9_v
        ajnm0=>ajn9_v;akbm0=>akb9_v;aktm0=>akt9_v
        ibpml=>ibp8_v; ffml=>ff8_v;bppml=>bpp8_v;resml=>res8_v
        lmgd2=lmgd1-1
        end if
        if(lmgd1.eq.8) then
        ibpm0=>ibp8_v; ffm0=>ff8_v;bppm0=>bpp8_v;resm0=>res8_v
        appm0=>app8_v;ajsm0=>ajs8_v;ajnm0=>ajn8_v
        akbm0=>akb8_v;aktm0=>akt8_v
        ibpml=>ibp7_v; ffml=>ff7_v;bppml=>bpp7_v
        resml=>res7_v;lmgd2=lmgd1-1
        end if
        if(lmgd1.eq.7) then
        ibpm0=>ibp7_v; ffm0=>ff7_v;bppm0=>bpp7_v;resm0=>res7_v
        appm0=>app7_v;ajsm0=>ajs7_v;ajnm0=>ajn7_v
        akbm0=>akb7_v;aktm0=>akt7_v
        ibpml=>ibp6_v; ffml=>ff6_v;bppml=>bpp6_v
        resml=>res6_v;lmgd2=lmgd1-1
        end if
        if(lmgd1.eq.6) then
        ibpm0=>ibp6_v; ffm0=>ff6_v;bppm0=>bpp6_v;resm0=>res6_v;
        appm0=>app6_v;ajsm0=>ajs6_v;ajnm0=>ajn6_v
        akbm0=>akb6_v;aktm0=>akt6_v
        ibpml=>ibp5_v; ffml=>ff5_v;bppml=>bpp5_v
        resml=>res5_v;lmgd2=lmgd1-1
        end if
        if(lmgd1.eq.5) then
        ibpm0=>ibp5_v; ffm0=>ff5_v;bppm0=>bpp5_v;resm0=>res5_v;
        appm0=>app5_v;ajsm0=>ajs5_v;ajnm0=>ajn5_v
        akbm0=>akb5_v;aktm0=>akt5_v
        ibpml=>ibp4_v; ffml=>ff4_v;bppml=>bpp4_v
        resml=>res4_v;lmgd2=lmgd1-1
        end if
        if(lmgd1.eq.4) then
        ibpm0=>ibp4_v; ffm0=>ff4_v;bppm0=>bpp4_v;resm0=>res4_v;
        appm0=>app4_v;ajsm0=>ajs4_v;ajnm0=>ajn4_v
        akbm0=>akb4_v;aktm0=>akt4_v
        ibpml=>ibp3_v; ffml=>ff3_v;bppml=>bpp3_v
        resml=>res3_v;lmgd2=lmgd1-1
        end if
        if(lmgd1.eq.3) then
        ibpm0=>ibp3_v; ffm0=>ff3_v;bppm0=>bpp3_v;resm0=>res3_v;
        appm0=>app3_v;ajsm0=>ajs3_v;ajnm0=>ajn3_v
        akbm0=>akb3_v;aktm0=>akt3_v
        ibpml=>ibp2_v; ffml=> ff2_v;bppml=>bpp2_v
        resml=>res2_v;lmgd2=lmgd1-1
        end if
        if(lmgd1.eq.2) then
        ibpm0=>ibp2_v; ffm0=> ff2_v;bppm0=>bpp2_v;resm0=>res2_v;
        appm0=>app2_v;ajsm0=>ajs2_v;ajnm0=>ajn2_v
        akbm0=>akb2_v;aktm0=>akt2_v
        ibpml=>ibp1_v; ffml=> ff1_v;bppml=>bpp1_v
        resml=>res1_v;lmgd2=lmgd1-1
        end if        
        if(lmgd1.eq.1) then
        ibpm0=>ibp1_v; ffm0=> ff1_v;bppm0=>bpp1_v;resm0=>res1_v;
        appm0=>app1_v;ajsm0=>ajs1_v;ajnm0=>ajn1_v
        akbm0=>akb1_v;aktm0=>akt1_v
        ibpml=>ibp0_v; ffml=> ff0_v;bppml=>bpp0_v
        resml=>res0_v;lmgd2=lmgd1-1
        end if
c--------------------------------c
        j2=j1mv(lmgd1)+1;m2=m1mv(lmgd1)-1;
        k2=k1mv(lmgd1)+1;n2=n1mv(lmgd1)-1
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
100     call solve0_v(lmgd1,mode1,ressumb,ressume)
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

200     call amgres_v(lmgd1)
        mode2=2
        call solvev(lmgd2,mode2,ressumb,ressume)
!        if(lmgd2.eq.1) call solveb(lmgd2,mode2,ressumb,ressume)
        call amgcor_v(lmgd1)
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
        subroutine amgcof_v(levmgd)
c===============================c
        use mdumgv;use mdupnt
        include 'table.prc'
        include 'table.gdm'
c-------------------------------c
        do lmgd=levmgd,0,-1
        nullify(ibpm0,ffm0,appm0,bppm0,resm0,ajsm0,ajnm0,akbm0,aktm0)
        nullify(ibpml,ffml,appml,bppml,resml,ajsml,ajnml,akbml,aktml)
        if(lmgd.eq.9) then
         ffm0=>ff9_v ;appm0=>app9_v;bppm0=>bpp9_v;ibpm0=>ibp9_v
        ajsm0=>ajs9_v;ajnm0=>ajn9_v;akbm0=>akb9_v;aktm0=>akt9_v
         ffml=>ff8_v;appml=>app8_v;bppml=>bpp8_v;ibpml=>ibp8_v
        ajsml=>ajs8_v;ajnml=>ajn8_v;aktml=>akt8_v;akbml=>akb8_v
        end if
        if(lmgd.eq.8) then
         ffm0=>ff8_v;appm0=>app8_v;bppm0=>bpp8_v;ibpm0=>ibp8_v
        ajsm0=>ajs8_v;ajnm0=>ajn8_v;akbm0=>akb8_v;aktm0=>akt8_v
         ffml=>ff7_v;appml=>app7_v;bppml=>bpp7_v;ibpml=>ibp7_v
        ajsml=>ajs7_v;ajnml=>ajn7_v;aktml=>akt7_v;akbml=>akb7_v
        end if
        if(lmgd.eq.7) then
         ffm0=>ff7_v;appm0=>app7_v;bppm0=>bpp7_v;ibpm0=>ibp7_v
        ajsm0=>ajs7_v;ajnm0=>ajn7_v;akbm0=>akb7_v;aktm0=>akt7_v
         ffml=>ff6_v;appml=>app6_v;bppml=>bpp6_v;ibpml=>ibp6_v
        ajsml=>ajs6_v;ajnml=>ajn6_v;akbml=>akb6_v;aktml=>akt6_v
        end if
        if(lmgd.eq.6) then
         ffm0=>ff6_v;appm0=>app6_v;bppm0=>bpp6_v;ibpm0=>ibp6_v
        ajsm0=>ajs6_v;ajnm0=>ajn6_v;akbm0=>akb6_v;aktm0=>akt6_v
         ffml=>ff5_v;appml=>app5_v;bppml=>bpp5_v;ibpml=>ibp5_v
        ajsml=>ajs5_v;ajnml=>ajn5_v;akbml=>akb5_v;aktml=>akt5_v
        end if
        if(lmgd.eq.5) then
         ffm0=>ff5_v;appm0=>app5_v;bppm0=>bpp5_v;ibpm0=>ibp5_v
        ajsm0=>ajs5_v;ajnm0=>ajn5_v;akbm0=>akb5_v;aktm0=>akt5_v
         ffml=>ff4_v;appml=>app4_v;bppml=>bpp4_v;ibpml=>ibp4_v
        ajsml=>ajs4_v;ajnml=>ajn4_v;akbml=>akb4_v;aktml=>akt4_v
        end if
        if(lmgd.eq.4) then
         ffm0=>ff4_v;appm0=>app4_v;bppm0=>bpp4_v;ibpm0=>ibp4_v
        ajsm0=>ajs4_v;ajnm0=>ajn4_v;akbm0=>akb4_v;aktm0=>akt4_v
         ffml=>ff3_v;appml=>app3_v;bppml=>bpp3_v;ibpml=>ibp3_v
        ajsml=>ajs3_v;ajnml=>ajn3_v;akbml=>akb3_v;aktml=>akt3_v
        end if
        if(lmgd.eq.3) then
         ffm0=>ff3_v;appm0=>app3_v;bppm0=>bpp3_v;ibpm0=>ibp3_v
        ajsm0=>ajs3_v;ajnm0=>ajn3_v;akbm0=>akb3_v;aktm0=>akt3_v
         ffml=>ff2_v;appml=>app2_v;bppml=>bpp2_v;ibpml=>ibp2_v
        ajsml=>ajs2_v;ajnml=>ajn2_v;akbml=>akb2_v;aktml=>akt2_v
        end if
        if(lmgd.eq.2) then
         ffm0=>ff2_v;appm0=>app2_v;bppm0=>bpp2_v;ibpm0=>ibp2_v
        ajsm0=>ajs2_v;ajnm0=>ajn2_v;akbm0=>akb2_v;aktm0=>akt2_v
         ffml=>ff1_v;appml=>app1_v;bppml=>bpp1_v;ibpml=>ibp1_v
        ajsml=>ajs1_v;ajnml=>ajn1_v;akbml=>akb1_v;aktml=>akt1_v
        end if
        if(lmgd.eq.1) then
         ffm0=>ff1_v;appm0=>app1_v;bppm0=>bpp1_v;ibpm0=>ibp1_v
        ajsm0=>ajs1_v;ajnm0=>ajn1_v;akbm0=>akb1_v;aktm0=>akt1_v
         ffml=>ff0_v;appml=>app0_v;bppml=>bpp0_v;ibpml=>ibp0_v
        ajsml=>ajs0_v;ajnml=>ajn0_v;akbml=>akb0_v;aktml=>akt0_v
        end if
        call coff12_v(lmgd,lmgd-1)
        end do
        nullify(ibpm0,ffm0,appm0,bppm0,resm0,ajsm0,ajnm0,akbm0,aktm0)
        nullify(ibpml,ffml,appml,bppml,resml,ajsml,ajnml,akbml,aktml)
        return;end   
     
c-------------------------------c
        subroutine coff12_v(mgd1,mgd2)
c-------------------------------c
        use mdumgv;use mdupnt
        include 'table.prc'
        include 'table.gdm'
c-------------------------------c
        j2=j1mv(mgd2)+1;m2=m1mv(mgd2)-1;k2=k1mv(mgd2)+1;n2=n1mv(mgd2)-1
        do j=j2,m2;do k=k2,n2
        kfb=2*k-k2;kft=kfb+1
        if(ibpm0(j,kfb).eq.0.and.ibpm0(j,kft).eq.0) then
        ibpml(j,k)=0
        elseif(ibpm0(j,kfb).eq.2.and.ibpm0(j,kft).eq.2) then
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
        kfb=2*k-k2;kft=kfb+1
        fkfb=0.0d+00;if(ibpm0(j,kfb).ne.0) fkfb=1.0d+00
        fkft=0.0d+00;if(ibpm0(j,kft).ne.0) fkft=1.0d+00
        appml(j,k)=fkfb*appm0(j,kfb)+fkft*appm0(j,kft)
     *            -fkft*akbm0(j,kft)-fkfb*aktm0(j,kfb)
        ajsml(j,k)=fkfb*ajsm0(j,kfb)+fkft*ajsm0(j,kft)
        ajnml(j,k)=fkfb*ajnm0(j,kfb)+fkft*ajnm0(j,kft)
        akbml(j,k)=fkfb*akbm0(j,kfb)
        aktml(j,k)=fkft*aktm0(j,kft)
        end if
        end do;end do
        return;end
c===============================c
        subroutine amgres_v(mgd1)
c===============================c
        use mdumgv;use mdupnt
        include 'table.prc'
        include 'table.gdm'
c-------------------------------c
        nullify(ibpm0,resm0,ibpml,bppml)
        if(mgd1.eq.9) then
        mgd2=mgd1-1;ibpml=>ibp8_v;bppml=>bpp8_v;
        ibpm0=>ibp9_v;resm0=>res9_v
        end if        
        if(mgd1.eq.8) then
        mgd2=mgd1-1;ibpml=>ibp7_v;bppml=>bpp7_v;
        ibpm0=>ibp8_v;resm0=>res8_v
        end if
        if(mgd1.eq.7) then
        mgd2=mgd1-1;ibpml=>ibp6_v;bppml=>bpp6_v;
        ibpm0=>ibp7_v;resm0=>res7_v
        end if
        if(mgd1.eq.6) then
        mgd2=mgd1-1;ibpml=>ibp5_v;bppml=>bpp5_v;
        ibpm0=>ibp6_v;resm0=>res6_v
        end if
        if(mgd1.eq.5) then
        mgd2=mgd1-1;ibpml=>ibp4_v;bppml=>bpp4_v;
        ibpm0=>ibp5_v;resm0=>res5_v
        end if
        if(mgd1.eq.4) then
        mgd2=mgd1-1;ibpml=>ibp3_v;bppml=>bpp3_v;
        ibpm0=>ibp4_v;resm0=>res4_v
        end if
        if(mgd1.eq.3) then
        mgd2=mgd1-1;ibpml=>ibp2_v;bppml=>bpp2_v;
        ibpm0=>ibp3_v;resm0=>res3_v
        end if
        if(mgd1.eq.2) then
        mgd2=mgd1-1;ibpml=>ibp1_v;bppml=>bpp1_v;
        ibpm0=>ibp2_v;resm0=>res2_v
        end if
        if(mgd1.eq.1) then
        mgd2=mgd1-1;ibpml=>ibp0_v;bppml=>bpp0_v;
        ibpm0=>ibp1_v;resm0=>res1_v
        end if
c-------------------------------c
        j2=j1mv(mgd2)+1;m2=m1mv(mgd2)-1;k2=k1mv(mgd2)+1;n2=n1mv(mgd2)-1
        do j=j2,m2;do k=k2,n2
        if(ibpml(j,k).eq.0) then
        bppml(j,k)=0.0d+00
        else
        kfb=2*k-k2;kft=kfb+1
        fkfb=0.0d+00;if(ibpm0(j,kfb).ne.0) fkfb=1.0d+00
        fkft=0.0d+00;if(ibpm0(j,kft).ne.0) fkft=1.0d+00
        bppml(j,k)=fkfb*resm0(j,kfb)+fkft*resm0(j,kft)
        end if
        end do;end do
        return;end
c===============================c
        subroutine amgcor_v(mgd1)
c===============================c
        use mdumgv;use mdupnt
        include 'table.prc'
        include 'table.gdm'
c-------------------------------c
        nullify(ibpm0,ffm0,ibpml,ffml)
        if(mgd1.eq.9) then
        mgd2=mgd1-1;ibpml=>ibp8_v;ffml=>ff8_v;ibpm0=>ibp9_v;ffm0=>ff9_v
        end if        
        if(mgd1.eq.8) then
        mgd2=mgd1-1;ibpml=>ibp7_v;ffml=>ff7_v;ibpm0=>ibp8_v;ffm0=>ff8_v
        end if
        if(mgd1.eq.7) then
        mgd2=mgd1-1;ibpml=>ibp6_v;ffml=>ff6_v;ibpm0=>ibp7_v;ffm0=>ff7_v
        end if
        if(mgd1.eq.6) then
        mgd2=mgd1-1;ibpml=>ibp5_v;ffml=>ff5_v;ibpm0=>ibp6_v;ffm0=>ff6_v
        end if
        if(mgd1.eq.5) then
        mgd2=mgd1-1;ibpml=>ibp4_v;ffml=>ff4_v;ibpm0=>ibp5_v;ffm0=>ff5_v
        end if
        if(mgd1.eq.4) then
        mgd2=mgd1-1;ibpml=>ibp3_v;ffml=>ff3_v;ibpm0=>ibp4_v;ffm0=>ff4_v
        end if
        if(mgd1.eq.3) then
        mgd2=mgd1-1;ibpml=>ibp2_v;ffml=>ff2_v;ibpm0=>ibp3_v;ffm0=>ff3_v
        end if
        if(mgd1.eq.2) then
        mgd2=mgd1-1;ibpml=>ibp1_v;ffml=>ff1_v;ibpm0=>ibp2_v;ffm0=>ff2_v
        end if
        if(mgd1.eq.1) then
        mgd2=mgd1-1;ibpml=>ibp0_v;ffml=>ff0_v;ibpm0=>ibp1_v;ffm0=>ff1_v
        end if
c-------------------------------c
        j2=j1mv(mgd2)+1;m2=m1mv(mgd2)-1;k2=k1mv(mgd2)+1;n2=n1mv(mgd2)-1
        do j=j2,m2;do k=k2,n2
        if(ibpml(j,k).eq.1) then
        kfb=2*k-k2;kft=kfb+1
        if(ibpm0(j,kfb).ne.0) ffm0(j,kfb)=ffm0(j,kfb)+ffml(j,k)
        if(ibpm0(j,kft).ne.0) ffm0(j,kft)=ffm0(j,kft)+ffml(j,k)
        end if
        end do;end do
        return;end

