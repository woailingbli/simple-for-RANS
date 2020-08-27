c================================c
        subroutine solve0(lmgd,mode,ressumb,ressume)
c================================c
        use mduvar;use mdumgv;use mdupnt;use mduwrk
        use omp_lib
        include 'table.prc'
        include 'table.gd1'
        include 'table.gdm'
c--------------------------------c
        allocate(fw(m1m(lmgd),n1m(lmgd)))
        allocate(pk(n1m(lmgd)),qk(n1m(lmgd)))
        allocate(pj(m1m(lmgd)),qj(m1m(lmgd)))
c--------------------------------c
        nullify(ibpm0,ffm0,appm0,bppm0,resm0,ajsm0,ajnm0,akbm0,aktm0)
        if(lmgd.eq.9) then
        ibpm0=>ibp9
         ffm0=>ff9 ;appm0=>app9;bppm0=>bpp9;resm0=>res9
        ajsm0=>ajs9;ajnm0=>ajn9;akbm0=>akb9;aktm0=>akt9
        end if
        if(lmgd.eq.8) then
        ibpm0=>ibp8
         ffm0=>ff8 ;appm0=>app8;bppm0=>bpp8;resm0=>res8
        ajsm0=>ajs8;ajnm0=>ajn8;akbm0=>akb8;aktm0=>akt8
        end if
        if(lmgd.eq.7) then
        ibpm0=>ibp7
         ffm0=>ff7 ;appm0=>app7;bppm0=>bpp7;resm0=>res7
        ajsm0=>ajs7;ajnm0=>ajn7;akbm0=>akb7;aktm0=>akt7
        end if
        if(lmgd.eq.6) then
        ibpm0=>ibp6
         ffm0=>ff6 ;appm0=>app6;bppm0=>bpp6;resm0=>res6
        ajsm0=>ajs6;ajnm0=>ajn6;akbm0=>akb6;aktm0=>akt6
        end if
        if(lmgd.eq.5) then
        ibpm0=>ibp5
         ffm0=>ff5 ;appm0=>app5;bppm0=>bpp5;resm0=>res5
        ajsm0=>ajs5;ajnm0=>ajn5;akbm0=>akb5;aktm0=>akt5
        end if
        if(lmgd.eq.4) then
        ibpm0=>ibp4
         ffm0=>ff4 ;appm0=>app4;bppm0=>bpp4;resm0=>res4
        ajsm0=>ajs4;ajnm0=>ajn4;akbm0=>akb4;aktm0=>akt4
        end if
        if(lmgd.eq.3) then
        ibpm0=>ibp3
         ffm0=>ff3 ;appm0=>app3;bppm0=>bpp3;resm0=>res3
        ajsm0=>ajs3;ajnm0=>ajn3;akbm0=>akb3;aktm0=>akt3
        end if
        if(lmgd.eq.2) then
        ibpm0=>ibp2
         ffm0=>ff2 ;appm0=>app2;bppm0=>bpp2;resm0=>res2
        ajsm0=>ajs2;ajnm0=>ajn2;akbm0=>akb2;aktm0=>akt2
        end if
        if(lmgd.eq.1) then
        ibpm0=>ibp1
         ffm0=>ff1 ;appm0=>app1;bppm0=>bpp1;resm0=>res1
        ajsm0=>ajs1;ajnm0=>ajn1;akbm0=>akb1;aktm0=>akt1
        end if
c--------------------------------c
        j2=j1m(lmgd)+1;m2=m1m(lmgd)-1
        k2=k1m(lmgd)+1;n2=n1m(lmgd)-1
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
        recursive subroutine solvem(lmgd1,mode,ressumb,ressume)
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
        ibpm0=>ibp9; ffm0=>ff9 ;bppm0=>bpp9;resm0=>res9
        appm0=>app9;ajsm0=>ajs9;ajnm0=>ajn9;akbm0=>akb9;aktm0=>akt9
        ibpml=>ibp8; ffml=>ff8 ;bppml=>bpp8;resml=>res8;lmgd2=lmgd1-1
        end if
        if(lmgd1.eq.8) then
        ibpm0=>ibp8; ffm0=>ff8 ;bppm0=>bpp8;resm0=>res8
        appm0=>app8;ajsm0=>ajs8;ajnm0=>ajn8;akbm0=>akb8;aktm0=>akt8
        ibpml=>ibp7; ffml=>ff7 ;bppml=>bpp7;resml=>res7;lmgd2=lmgd1-1
        end if
        if(lmgd1.eq.7) then
        ibpm0=>ibp7; ffm0=>ff7 ;bppm0=>bpp7;resm0=>res7
        appm0=>app7;ajsm0=>ajs7;ajnm0=>ajn7;akbm0=>akb7;aktm0=>akt7
        ibpml=>ibp6; ffml=>ff6 ;bppml=>bpp6;resml=>res6;lmgd2=lmgd1-1
        end if
        if(lmgd1.eq.6) then
        ibpm0=>ibp6; ffm0=>ff6 ;bppm0=>bpp6;resm0=>res6;
        appm0=>app6;ajsm0=>ajs6;ajnm0=>ajn6;akbm0=>akb6;aktm0=>akt6
        ibpml=>ibp5; ffml=>ff5 ;bppml=>bpp5;resml=>res5;lmgd2=lmgd1-1
        end if
        if(lmgd1.eq.5) then
        ibpm0=>ibp5; ffm0=>ff5 ;bppm0=>bpp5;resm0=>res5;
        appm0=>app5;ajsm0=>ajs5;ajnm0=>ajn5;akbm0=>akb5;aktm0=>akt5
        ibpml=>ibp4; ffml=>ff4 ;bppml=>bpp4;resml=>res4;lmgd2=lmgd1-1
        end if
        if(lmgd1.eq.4) then
        ibpm0=>ibp4; ffm0=>ff4 ;bppm0=>bpp4;resm0=>res4;
        appm0=>app4;ajsm0=>ajs4;ajnm0=>ajn4;akbm0=>akb4;aktm0=>akt4
        ibpml=>ibp3; ffml=>ff3 ;bppml=>bpp3;resml=>res3;lmgd2=lmgd1-1
        end if
        if(lmgd1.eq.3) then
        ibpm0=>ibp3; ffm0=> ff3;bppm0=>bpp3;resm0=>res3;
        appm0=>app3;ajsm0=>ajs3;ajnm0=>ajn3;akbm0=>akb3;aktm0=>akt3
        ibpml=>ibp2; ffml=> ff2;bppml=>bpp2;resml=>res2;lmgd2=lmgd1-1
        end if
        if(lmgd1.eq.2) then
        ibpm0=>ibp2; ffm0=> ff2;bppm0=>bpp2;resm0=>res2;
        appm0=>app2;ajsm0=>ajs2;ajnm0=>ajn2;akbm0=>akb2;aktm0=>akt2
        ibpml=>ibp1; ffml=> ff1;bppml=>bpp1;resml=>res1;lmgd2=lmgd1-1
        end if        
        if(lmgd1.eq.1) then
        ibpm0=>ibp1; ffm0=> ff1;bppm0=>bpp1;resm0=>res1;
        appm0=>app1;ajsm0=>ajs1;ajnm0=>ajn1;akbm0=>akb1;aktm0=>akt1
        ibpml=>ibp0; ffml=> ff0;bppml=>bpp0;resml=>res0;lmgd2=lmgd1-1
        end if
c--------------------------------c
        j2=j1m(lmgd1)+1;m2=m1m(lmgd1)-1;k2=k1m(lmgd1)+1;n2=n1m(lmgd1)-1
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
100     call solve0(lmgd1,mode1,ressumb,ressume)
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

200     call amgres(lmgd1)
        mode2=2
        if(lmgd2.gt.1) call solvem(lmgd2,mode2,ressumb,ressume)
        if(lmgd2.eq.1) call solveb(lmgd2,mode2,ressumb,ressume)
        call amgcor(lmgd1)
        goto 100
201     format(' *        total residual satisfied,   return        *')
202     format(' *        residual not satisfied, continuing        *')
203     format(' *        declev not reached,decrat exceeded        *')
204     format(' *        declev satisfied,  go out of L-B-L        *')
1000    ressumb=ressum0(lmgd1)
        ressume=ressum1(lmgd1)
        return
        end

c================================c
        subroutine solveb(lmgd,mode,ressumb,ressume)
c================================c
        use mdumgv;use mdupnt
        include  'table.prc'
        include  'table.gd1'
        include  'table.gdm'
c--------------------------------c
        nullify(ibpm0,ffm0,appm0,bppm0,resm0,ajsm0,ajnm0,akbm0,aktm0)
        nullify(ibpml,ffml,appml,bppml,resml,ajsml,ajnml,akbml,aktml)
        if(lmgd.eq.1) then
        ibpml=>ibp0; ffml=>ff0 ;bppml=>bpp0;resml=>res0
        ibpm0=>ibp1; ffm0=>ff1 ;bppm0=>bpp1;resm0=>res1
        appm0=>app1;ajsm0=>ajs1;ajnm0=>ajn1;akbm0=>akb1;aktm0=>akt1
        end if
c--------------------------------c
        j2=j1m(lmgd)+1;m2=m1m(lmgd)-1;k2=k1m(lmgd)+1;n2=n1m(lmgd)-1;
        jw=j1m(0)+1;m0=m1m(0);kw=k1m(0)+1;n0=n1m(0)
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
        ressum0(lmgd)=ressumb
c        write(*,11) lmgd,ressumb
11      format
     *  (' *',2x,'L-',i2,2x,'initial residual ressum=',1pe10.3,8x,'*')
c--------------------------------c
100     mode1=3;call solve0(lmgd,mode1,ressumb,ressume)
        ressum1(lmgd)=ressume
c        write(*,12) ressume
12      format
     *  (' *','  l-2   ','postLBL residual ressum=',1pe10.3,8x,'*')
c-------mode=1: control residual precision
        if(mode.eq.1) then
        if(ressum1(lmgd).le.epsp2(lmgd)) then
c        write(*,201)
        goto 1000
        else
c        write(*,202)
        goto 200
        end if
        end if

c-------mode=2: control residual decreasing level
        if(mode.eq.2) then
        if(ressum1(lmgd).le.epsp2(lmgd)) then
c        write(*,201)
        goto 1000
        end if
        declev=ressum1(lmgd)/ressum0(lmgd)
        if(declev.gt.epsp3(lmgd)) then
c        write(*,203)
        goto 200
        else
c        write(*,204)
        goto 1000
        end if
        end if

200     call amgres(lmgd)
        ff0(jw,kw)=bpp0(jw,kw)/app0(jw,kw)
        call amgcor(lmgd)
        goto 100
201     format(' *        total residual satisfied,   return        *')
202     format(' *        residual not satisfied, continuing        *')
203     format(' *        declev not reached,decrat exceeded        *')
204     format(' *        declev satisfied,  go out of L-B-L        *')
205     format(' *',2x, 'nwg=',i7,1x,'<',1x,'nwm=',i7,1x,'stop',2x,' *')
c1000    deallocate(a,b,z,ma,mv,v)
1000    return
        end

c===============================c
        subroutine amgcof(levmgd)
c===============================c
        use mdumgv;use mdupnt
        include 'table.prc'
        include 'table.gdm'
c-------------------------------c
        do lmgd=levmgd,0,-1
        nullify(ibpm0,ffm0,appm0,bppm0,resm0,ajsm0,ajnm0,akbm0,aktm0)
        nullify(ibpml,ffml,appml,bppml,resml,ajsml,ajnml,akbml,aktml)
        if(lmgd.eq.9) then
         ffm0=>ff9 ;appm0=>app9 ;bppm0=>bpp9 ;ibpm0=>ibp9
        ajsm0=>ajs9;ajnm0=>ajn9 ;akbm0=>akb9 ;aktm0=>akt9
         ffml=>ff8 ;appml=>app8 ;bppml=>bpp8 ;ibpml=>ibp8
        ajsml=>ajs8;ajnml=>ajn8 ;aktml=>akt8 ;akbml=>akb8
        end if
        if(lmgd.eq.8) then
         ffm0=>ff8 ;appm0=>app8 ;bppm0=>bpp8 ;ibpm0=>ibp8
        ajsm0=>ajs8;ajnm0=>ajn8 ;akbm0=>akb8 ;aktm0=>akt8
         ffml=>ff7 ;appml=>app7 ;bppml=>bpp7 ;ibpml=>ibp7
        ajsml=>ajs7;ajnml=>ajn7 ;aktml=>akt7 ;akbml=>akb7
        end if
        if(lmgd.eq.7) then
         ffm0=>ff7 ;appm0=>app7;bppm0=>bpp7;ibpm0=>ibp7
        ajsm0=>ajs7;ajnm0=>ajn7;akbm0=>akb7;aktm0=>akt7
         ffml=>ff6 ;appml=>app6;bppml=>bpp6;ibpml=>ibp6
        ajsml=>ajs6;ajnml=>ajn6;akbml=>akb6;aktml=>akt6
        end if
        if(lmgd.eq.6) then
         ffm0=>ff6 ;appm0=>app6;bppm0=>bpp6;ibpm0=>ibp6
        ajsm0=>ajs6;ajnm0=>ajn6;akbm0=>akb6;aktm0=>akt6
         ffml=>ff5 ;appml=>app5;bppml=>bpp5;ibpml=>ibp5
        ajsml=>ajs5;ajnml=>ajn5;akbml=>akb5;aktml=>akt5
        end if
        if(lmgd.eq.5) then
         ffm0=>ff5 ;appm0=>app5;bppm0=>bpp5;ibpm0=>ibp5
        ajsm0=>ajs5;ajnm0=>ajn5;akbm0=>akb5;aktm0=>akt5
         ffml=>ff4 ;appml=>app4;bppml=>bpp4;ibpml=>ibp4
        ajsml=>ajs4;ajnml=>ajn4;akbml=>akb4;aktml=>akt4
        end if
        if(lmgd.eq.4) then
         ffm0=>ff4 ;appm0=>app4;bppm0=>bpp4;ibpm0=>ibp4
        ajsm0=>ajs4;ajnm0=>ajn4;akbm0=>akb4;aktm0=>akt4
         ffml=>ff3 ;appml=>app3;bppml=>bpp3;ibpml=>ibp3
        ajsml=>ajs3;ajnml=>ajn3;akbml=>akb3;aktml=>akt3
        end if
        if(lmgd.eq.3) then
         ffm0=>ff3 ;appm0=>app3;bppm0=>bpp3;ibpm0=>ibp3
        ajsm0=>ajs3;ajnm0=>ajn3;akbm0=>akb3;aktm0=>akt3
         ffml=>ff2 ;appml=>app2;bppml=>bpp2;ibpml=>ibp2
        ajsml=>ajs2;ajnml=>ajn2;akbml=>akb2;aktml=>akt2
        end if
        if(lmgd.eq.2) then
         ffm0=>ff2 ;appm0=>app2;bppm0=>bpp2;ibpm0=>ibp2
        ajsm0=>ajs2;ajnm0=>ajn2;akbm0=>akb2;aktm0=>akt2
         ffml=>ff1 ;appml=>app1;bppml=>bpp1;ibpml=>ibp1
        ajsml=>ajs1;ajnml=>ajn1;akbml=>akb1;aktml=>akt1
        end if
        if(lmgd.eq.1) then
         ffm0=>ff1 ;appm0=>app1;bppm0=>bpp1;ibpm0=>ibp1
        ajsm0=>ajs1;ajnm0=>ajn1;akbm0=>akb1;aktm0=>akt1
         ffml=>ff0 ;appml=>app0;bppml=>bpp0;ibpml=>ibp0
        ajsml=>ajs0;ajnml=>ajn0;akbml=>akb0;aktml=>akt0
        end if
        call coff12(lmgd,lmgd-1)
        end do
        nullify(ibpm0,ffm0,appm0,bppm0,resm0,ajsm0,ajnm0,akbm0,aktm0)
        nullify(ibpml,ffml,appml,bppml,resml,ajsml,ajnml,akbml,aktml)
        return;end
c-------------------------------c
        subroutine coff12(mgd1,mgd2)
c-------------------------------c
        use mdumgv;use mdupnt
        include 'table.prc'
        include 'table.gdm'
c-------------------------------c
        j2=j1m(mgd2)+1;m2=m1m(mgd2)-1;k2=k1m(mgd2)+1;n2=n1m(mgd2)-1
        do j=j2,m2;do k=k2,n2
        jfs=2*j-j2;jfn=jfs+1;kfb=2*k-k2;kft=kfb+1
        if(ibpm0(jfs,kfb).eq.0.and.ibpm0(jfn,kfb).eq.0.and.
     *     ibpm0(jfs,kft).eq.0.and.ibpm0(jfn,kft).eq.0) then
        ibpml(j,k)=0
        elseif(ibpm0(jfs,kfb).eq.2.and.ibpm0(jfn,kfb).eq.2.and.
     *     ibpm0(jfs,kft).eq.2.and.ibpm0(jfn,kft).eq.2) then
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
        jfs=2*j-j2;jfn=jfs+1;kfb=2*k-k2;kft=kfb+1
        fjfskfb=0.0d+00;if(ibpm0(jfs,kfb).ne.0) fjfskfb=1.0d+00
        fjfnkfb=0.0d+00;if(ibpm0(jfn,kfb).ne.0) fjfnkfb=1.0d+00
        fjfskft=0.0d+00;if(ibpm0(jfs,kft).ne.0) fjfskft=1.0d+00
        fjfnkft=0.0d+00;if(ibpm0(jfn,kft).ne.0) fjfnkft=1.0d+00
        appml(j,k)=fjfskfb*appm0(jfs,kfb)+fjfnkfb*appm0(jfn,kfb)
     *            +fjfskft*appm0(jfs,kft)+fjfnkft*appm0(jfn,kft)
     *            -fjfnkfb*ajsm0(jfn,kfb)-fjfnkft*ajsm0(jfn,kft)
     *            -fjfskfb*ajnm0(jfs,kfb)-fjfskft*ajnm0(jfs,kft)
     *            -fjfskft*akbm0(jfs,kft)-fjfnkft*akbm0(jfn,kft)
     *            -fjfskfb*aktm0(jfs,kfb)-fjfnkfb*aktm0(jfn,kfb)
        ajsml(j,k)=fjfskfb*ajsm0(jfs,kfb)+fjfskft*ajsm0(jfs,kft)
        ajnml(j,k)=fjfnkfb*ajnm0(jfn,kfb)+fjfnkft*ajnm0(jfn,kft)
        akbml(j,k)=fjfskfb*akbm0(jfs,kfb)+fjfnkfb*akbm0(jfn,kfb)
        aktml(j,k)=fjfskft*aktm0(jfs,kft)+fjfnkft*aktm0(jfn,kft)
        end if
        end do;end do
        return;end
        

        
c===============================c
        subroutine amgres(mgd1)
c===============================c
        use mdumgv;use mdupnt
        include 'table.prc'
        include 'table.gdm'
c-------------------------------c
        nullify(ibpm0,resm0,ibpml,bppml)
        if(mgd1.eq.9) then
        mgd2=mgd1-1;ibpml=>ibp8;bppml=>bpp8;ibpm0=>ibp9;resm0=>res9
        end if        
        if(mgd1.eq.8) then
        mgd2=mgd1-1;ibpml=>ibp7;bppml=>bpp7;ibpm0=>ibp8;resm0=>res8
        end if
        if(mgd1.eq.7) then
        mgd2=mgd1-1;ibpml=>ibp6;bppml=>bpp6;ibpm0=>ibp7;resm0=>res7
        end if
        if(mgd1.eq.6) then
        mgd2=mgd1-1;ibpml=>ibp5;bppml=>bpp5;ibpm0=>ibp6;resm0=>res6
        end if
        if(mgd1.eq.5) then
        mgd2=mgd1-1;ibpml=>ibp4;bppml=>bpp4;ibpm0=>ibp5;resm0=>res5
        end if
        if(mgd1.eq.4) then
        mgd2=mgd1-1;ibpml=>ibp3;bppml=>bpp3;ibpm0=>ibp4;resm0=>res4
        end if
        if(mgd1.eq.3) then
        mgd2=mgd1-1;ibpml=>ibp2;bppml=>bpp2;ibpm0=>ibp3;resm0=>res3
        end if
        if(mgd1.eq.2) then
        mgd2=mgd1-1;ibpml=>ibp1;bppml=>bpp1;ibpm0=>ibp2;resm0=>res2
        end if
        if(mgd1.eq.1) then
        mgd2=mgd1-1;ibpml=>ibp0;bppml=>bpp0;ibpm0=>ibp1;resm0=>res1
        end if
c-------------------------------c
        j2=j1m(mgd2)+1;m2=m1m(mgd2)-1;k2=k1m(mgd2)+1;n2=n1m(mgd2)-1
        do j=j2,m2;do k=k2,n2
        if(ibpml(j,k).eq.0) then
        bppml(j,k)=0.0d+00
        else
        jfs=2*j-j2;jfn=jfs+1;kfb=2*k-k2;kft=kfb+1
        fjfskfb=0.0d+00;if(ibpm0(jfs,kfb).ne.0) fjfskfb=1.0d+00
        fjfskft=0.0d+00;if(ibpm0(jfs,kft).ne.0) fjfskft=1.0d+00
        fjfnkfb=0.0d+00;if(ibpm0(jfn,kfb).ne.0) fjfnkfb=1.0d+00
        fjfnkft=0.0d+00;if(ibpm0(jfn,kft).ne.0) fjfnkft=1.0d+00
        bppml(j,k)=fjfskfb*resm0(jfs,kfb)+fjfskft*resm0(jfs,kft)
     *            +fjfnkfb*resm0(jfn,kfb)+fjfnkft*resm0(jfn,kft)
        end if
        end do;end do
        return;end
        

        
c===============================c
        subroutine amgcor(mgd1)
c===============================c
        use mdumgv;use mdupnt
        include 'table.prc'
        include 'table.gdm'
c-------------------------------c
        nullify(ibpm0,ffm0,ibpml,ffml)
        if(mgd1.eq.9) then
        mgd2=mgd1-1;ibpml=>ibp8;ffml=>ff8;ibpm0=>ibp9;ffm0=>ff9
        end if        
        if(mgd1.eq.8) then
        mgd2=mgd1-1;ibpml=>ibp7;ffml=>ff7;ibpm0=>ibp8;ffm0=>ff8
        end if
        if(mgd1.eq.7) then
        mgd2=mgd1-1;ibpml=>ibp6;ffml=>ff6;ibpm0=>ibp7;ffm0=>ff7
        end if
        if(mgd1.eq.6) then
        mgd2=mgd1-1;ibpml=>ibp5;ffml=>ff5;ibpm0=>ibp6;ffm0=>ff6
        end if
        if(mgd1.eq.5) then
        mgd2=mgd1-1;ibpml=>ibp4;ffml=>ff4;ibpm0=>ibp5;ffm0=>ff5
        end if
        if(mgd1.eq.4) then
        mgd2=mgd1-1;ibpml=>ibp3;ffml=>ff3;ibpm0=>ibp4;ffm0=>ff4
        end if
        if(mgd1.eq.3) then
        mgd2=mgd1-1;ibpml=>ibp2;ffml=>ff2;ibpm0=>ibp3;ffm0=>ff3
        end if
        if(mgd1.eq.2) then
        mgd2=mgd1-1;ibpml=>ibp1;ffml=>ff1;ibpm0=>ibp2;ffm0=>ff2
        end if
        if(mgd1.eq.1) then
        mgd2=mgd1-1;ibpml=>ibp0;ffml=>ff0;ibpm0=>ibp1;ffm0=>ff1
        end if
c-------------------------------c
        j2=j1m(mgd2)+1;m2=m1m(mgd2)-1;k2=k1m(mgd2)+1;n2=n1m(mgd2)-1
        do j=j2,m2;do k=k2,n2
        if(ibpml(j,k).eq.1) then
        jfs=2*j-j2;jfn=jfs+1;kfb=2*k-k2;kft=kfb+1
        if(ibpm0(jfs,kfb).ne.0) ffm0(jfs,kfb)=ffm0(jfs,kfb)+ffml(j,k)
        if(ibpm0(jfn,kfb).ne.0) ffm0(jfn,kfb)=ffm0(jfn,kfb)+ffml(j,k)
        if(ibpm0(jfs,kft).ne.0) ffm0(jfs,kft)=ffm0(jfs,kft)+ffml(j,k)
        if(ibpm0(jfn,kft).ne.0) ffm0(jfn,kft)=ffm0(jfn,kft)+ffml(j,k)
        end if
        end do;end do
        return;end

