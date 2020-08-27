
c===============================c
        subroutine alcgrd
c===============================c
        use mdugrd
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
c-------grid one-----------c
        allocate(xu1(imx),yv1(jmx),zw1(kmx))
        allocate(xm1(imx),dxm(imx),dxu(imx),fx1(imx),fx2(imx))
        allocate(ym1(jmx),dym(jmx),dyv(jmx),fy1(jmx),fy2(jmx))
        allocate(zm1(kmx),dzm(kmx),dzw(kmx),fz1(kmx),fz2(kmx))
        allocate(ibu(jmx,kmx),ibv(jmx,kmx),ibw(jmx,kmx),ibp(jmx,kmx))
        return;end

c===============================c
        subroutine alcvar
c===============================c
        use mduvar;use mdumgv
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
        include 'table.gdm'
c-------physical variables------c
        allocate
     *  (ulm(imx,jmx,kmx),vlm(imx,jmx,kmx),
     *   wlm(imx,jmx,kmx),tlm(imx,jmx,kmx))
        allocate
     *  (ui(imx,jmx,kmx),vi(imx,jmx,kmx),
     *   wi(imx,jmx,kmx),ti(imx,jmx,kmx))
        allocate
     *  (ui1(jmx,kmx),vi1(jmx,kmx),
     *   wi1(jmx,kmx),pi1(jmx,kmx))
        allocate
     *  (vpvp(jmx,kmx),vpwp(jmx,kmx),
     *   wpwp(jmx,kmx),upvp(jmx,kmx),
     *   upwp(jmx,kmx))
        allocate
     *  (vv1(jmx,kmx),ww1(jmx,kmx),
     *   pre1(jmx,kmx),apc1(jmx,kmx),apc2(jmx,kmx))
        allocate
     *  (um(jmx,kmx),vm(jmx,kmx),wm(jmx,kmx),pm(jmx,kmx))
        allocate(delta(jmx,kmx))

        allocate(disp(jmx,kmx))
c-------multigrid variables------c
        do lmgd=0,levmgd
        if(lmgd.eq.9) then
        allocate(ibp9(m1m(lmgd),n1m(lmgd)),
     *   ff9(m1m(lmgd),n1m(lmgd)),res9(m1m(lmgd),n1m(lmgd)),
     *  app9(m1m(lmgd),n1m(lmgd)),bpp9(m1m(lmgd),n1m(lmgd)),
     *  ajs9(m1m(lmgd),n1m(lmgd)),akb9(m1m(lmgd),n1m(lmgd)),
     *  ajn9(m1m(lmgd),n1m(lmgd)),akt9(m1m(lmgd),n1m(lmgd)))
        ibp9=0; ff9=0.0d+00; res9=0.0d+00; app9=0.0d+00 ;bpp9=0.0d+00
        akb9=0.0d+00; akt9=0.0d+00; ajs9=0.0d+00; ajn9=0.0d+00
        end if
        if(lmgd.eq.8) then
        allocate(ibp8(m1m(lmgd),n1m(lmgd)),
     *   ff8(m1m(lmgd),n1m(lmgd)),res8(m1m(lmgd),n1m(lmgd)),
     *  app8(m1m(lmgd),n1m(lmgd)),bpp8(m1m(lmgd),n1m(lmgd)),
     *  ajs8(m1m(lmgd),n1m(lmgd)),akb8(m1m(lmgd),n1m(lmgd)),
     *  ajn8(m1m(lmgd),n1m(lmgd)),akt8(m1m(lmgd),n1m(lmgd)))
        ibp8=0; ff8=0.0d+00; res8=0.0d+00; app8=0.0d+00 ;bpp8=0.0d+00
        akb8=0.0d+00; akt8=0.0d+00; ajs8=0.0d+00; ajn8=0.0d+00
        end if
        if(lmgd.eq.7) then
        allocate(ibp7(m1m(lmgd),n1m(lmgd)),
     *   ff7(m1m(lmgd),n1m(lmgd)),res7(m1m(lmgd),n1m(lmgd)),
     *  app7(m1m(lmgd),n1m(lmgd)),bpp7(m1m(lmgd),n1m(lmgd)),
     *  ajs7(m1m(lmgd),n1m(lmgd)),akb7(m1m(lmgd),n1m(lmgd)),
     *  ajn7(m1m(lmgd),n1m(lmgd)),akt7(m1m(lmgd),n1m(lmgd)))
        ibp7=0; ff7=0.0d+00; res7=0.0d+00; app7=0.0d+00; bpp7=0.0d+00
        akb7=0.0d+00; akt7=0.0d+00; ajs7=0.0d+00; ajn7=0.0d+00
        end if
        if(lmgd.eq.6) then
        allocate(ibp6(m1m(lmgd),n1m(lmgd)),
     *   ff6(m1m(lmgd),n1m(lmgd)),res6(m1m(lmgd),n1m(lmgd)),
     *  app6(m1m(lmgd),n1m(lmgd)),bpp6(m1m(lmgd),n1m(lmgd)),
     *  ajs6(m1m(lmgd),n1m(lmgd)),akb6(m1m(lmgd),n1m(lmgd)),
     *  ajn6(m1m(lmgd),n1m(lmgd)),akt6(m1m(lmgd),n1m(lmgd)))
        ibp6=0; ff6=0.0d+00; res6=0.0d+00; app6=0.0d+00; bpp6=0.0d+00
        akb6=0.0d+00; akt6=0.0d+00; ajs6=0.0d+00; ajn6=0.0d+00
        end if
        if(lmgd.eq.5) then
        allocate(ibp5(m1m(lmgd),n1m(lmgd)),
     *   ff5(m1m(lmgd),n1m(lmgd)),res5(m1m(lmgd),n1m(lmgd)),
     *  app5(m1m(lmgd),n1m(lmgd)),bpp5(m1m(lmgd),n1m(lmgd)),
     *  ajs5(m1m(lmgd),n1m(lmgd)),akb5(m1m(lmgd),n1m(lmgd)),
     *  ajn5(m1m(lmgd),n1m(lmgd)),akt5(m1m(lmgd),n1m(lmgd)))
        ibp5=0; ff5=0.0d+00; res5=0.0d+00; app5=0.0d+00; bpp5=0.0d+00
        akb5=0.0d+00; akt5=0.0d+00; ajs5=0.0d+00; ajn5=0.0d+00
        end if
        if(lmgd.eq.4) then
        allocate(ibp4(m1m(lmgd),n1m(lmgd)),
     *   ff4(m1m(lmgd),n1m(lmgd)),res4(m1m(lmgd),n1m(lmgd)),
     *  app4(m1m(lmgd),n1m(lmgd)),bpp4(m1m(lmgd),n1m(lmgd)),
     *  ajs4(m1m(lmgd),n1m(lmgd)),akb4(m1m(lmgd),n1m(lmgd)),
     *  ajn4(m1m(lmgd),n1m(lmgd)),akt4(m1m(lmgd),n1m(lmgd)))
        ibp4=0; ff4=0.0d+00; res4=0.0d+00; app4=0.0d+00; bpp4=0.0d+00
        akb4=0.0d+00; akt4=0.0d+00; ajs4=0.0d+00; ajn4=0.0d+00
        end if
        if(lmgd.eq.3) then
        allocate(ibp3(m1m(lmgd),n1m(lmgd)),
     *   ff3(m1m(lmgd),n1m(lmgd)),res3(m1m(lmgd),n1m(lmgd)),
     *  app3(m1m(lmgd),n1m(lmgd)),bpp3(m1m(lmgd),n1m(lmgd)),
     *  ajs3(m1m(lmgd),n1m(lmgd)),akb3(m1m(lmgd),n1m(lmgd)),
     *  ajn3(m1m(lmgd),n1m(lmgd)),akt3(m1m(lmgd),n1m(lmgd)))
        ibp3=0; ff3=0.0d+00; res3=0.0d+00; app3=0.0d+00; bpp3=0.0d+00
        akb3=0.0d+00; akt3=0.0d+00; ajs3=0.0d+00; ajn3=0.0d+00
        end if
        if(lmgd.eq.2) then
        allocate(ibp2(m1m(lmgd),n1m(lmgd)),
     *   ff2(m1m(lmgd),n1m(lmgd)),res2(m1m(lmgd),n1m(lmgd)),
     *  app2(m1m(lmgd),n1m(lmgd)),bpp2(m1m(lmgd),n1m(lmgd)),
     *  ajs2(m1m(lmgd),n1m(lmgd)),akb2(m1m(lmgd),n1m(lmgd)),
     *  ajn2(m1m(lmgd),n1m(lmgd)),akt2(m1m(lmgd),n1m(lmgd)))
        ibp2=0; ff2=0.0d+00; res2=0.0d+00; app2=0.0d+00; bpp2=0.0d+00
        akb2=0.0d+00; akt2=0.0d+00; ajs2=0.0d+00; ajn2=0.0d+00
        end if
        if(lmgd.eq.1) then
        allocate(ibp1(m1m(lmgd),n1m(lmgd)),
     *   ff1(m1m(lmgd),n1m(lmgd)),res1(m1m(lmgd),n1m(lmgd)),
     *  app1(m1m(lmgd),n1m(lmgd)),bpp1(m1m(lmgd),n1m(lmgd)),
     *  ajs1(m1m(lmgd),n1m(lmgd)),akb1(m1m(lmgd),n1m(lmgd)),
     *  ajn1(m1m(lmgd),n1m(lmgd)),akt1(m1m(lmgd),n1m(lmgd)))
        ibp1=0; ff1=0.0d+00; res1=0.0d+00; app1=0.0d+00; bpp1=0.0d+00
        akb1=0.0d+00; akt1=0.0d+00; ajs1=0.0d+00; ajn1=0.0d+00
        end if
        if(lmgd.eq.0) then
        allocate(ibp0(m1m(lmgd),n1m(lmgd)),
     *   ff0(m1m(lmgd),n1m(lmgd)),res0(m1m(lmgd),n1m(lmgd)),
     *  app0(m1m(lmgd),n1m(lmgd)),bpp0(m1m(lmgd),n1m(lmgd)),
     *  ajs0(m1m(lmgd),n1m(lmgd)),akb0(m1m(lmgd),n1m(lmgd)),
     *  ajn0(m1m(lmgd),n1m(lmgd)),akt0(m1m(lmgd),n1m(lmgd)))
        ibp0=0; ff0=0.0d+00; res0=0.0d+00; app0=0.0d+00; bpp0=0.0d+00
        akb0=0.0d+00; akt0=0.0d+00; ajs0=0.0d+00; ajn0=0.0d+00
        end if
        end do

        do lmgd=0,levmgd
        if(lmgd.eq.9) then
        allocate(ibp9_v(m1mv(lmgd),n1mv(lmgd)),
     *   ff9_v(m1mv(lmgd),n1mv(lmgd)),res9_v(m1mv(lmgd),n1mv(lmgd)),
     *  app9_v(m1mv(lmgd),n1mv(lmgd)),bpp9_v(m1mv(lmgd),n1mv(lmgd)),
     *  ajs9_v(m1mv(lmgd),n1mv(lmgd)),akb9_v(m1mv(lmgd),n1mv(lmgd)),
     *  ajn9_v(m1mv(lmgd),n1mv(lmgd)),akt9_v(m1mv(lmgd),n1mv(lmgd)))
        ibp9_v=0;ff9_v=0.0d+00;res9_v=0.0d+00;
        app9_v=0.0d+00;bpp9_v=0.0d+00
        akb9_v=0.0d+00; akt9_v=0.0d+00; ajs9_v=0.0d+00; ajn9_v=0.0d+00
        end if
        if(lmgd.eq.8) then
        allocate(ibp8_v(m1mv(lmgd),n1mv(lmgd)),
     *   ff8_v(m1mv(lmgd),n1mv(lmgd)),res8_v(m1mv(lmgd),n1mv(lmgd)),
     *  app8_v(m1mv(lmgd),n1mv(lmgd)),bpp8_v(m1mv(lmgd),n1mv(lmgd)),
     *  ajs8_v(m1mv(lmgd),n1mv(lmgd)),akb8_v(m1mv(lmgd),n1mv(lmgd)),
     *  ajn8_v(m1mv(lmgd),n1mv(lmgd)),akt8_v(m1mv(lmgd),n1mv(lmgd)))
        ibp8_v=0;ff8_v=0.0d+00;res8_v=0.0d+00;
        app8_v=0.0d+00;bpp8_v=0.0d+00
        akb8_v=0.0d+00; akt8_v=0.0d+00; ajs8_v=0.0d+00; ajn8_v=0.0d+00
        end if
        if(lmgd.eq.7) then
        allocate(ibp7_v(m1mv(lmgd),n1mv(lmgd)),
     *   ff7_v(m1mv(lmgd),n1mv(lmgd)),res7_v(m1mv(lmgd),n1mv(lmgd)),
     *  app7_v(m1mv(lmgd),n1mv(lmgd)),bpp7_v(m1mv(lmgd),n1mv(lmgd)),
     *  ajs7_v(m1mv(lmgd),n1mv(lmgd)),akb7_v(m1mv(lmgd),n1mv(lmgd)),
     *  ajn7_v(m1mv(lmgd),n1mv(lmgd)),akt7_v(m1mv(lmgd),n1mv(lmgd)))
        ibp7_v=0;ff7_v=0.0d+00;res7_v=0.0d+00;
        app7_v=0.0d+00;bpp7_v=0.0d+00
        akb7_v=0.0d+00; akt7_v=0.0d+00; ajs7_v=0.0d+00; ajn7_v=0.0d+00
        end if
        if(lmgd.eq.6) then
        allocate(ibp6_v(m1mv(lmgd),n1mv(lmgd)),
     *   ff6_v(m1mv(lmgd),n1mv(lmgd)),res6_v(m1mv(lmgd),n1mv(lmgd)),
     *  app6_v(m1mv(lmgd),n1mv(lmgd)),bpp6_v(m1mv(lmgd),n1mv(lmgd)),
     *  ajs6_v(m1mv(lmgd),n1mv(lmgd)),akb6_v(m1mv(lmgd),n1mv(lmgd)),
     *  ajn6_v(m1mv(lmgd),n1mv(lmgd)),akt6_v(m1mv(lmgd),n1mv(lmgd)))
        ibp6_v=0;ff6_v=0.0d+00;res6_v=0.0d+00;
        app6_v=0.0d+00;bpp6_v=0.0d+00
        akb6_v=0.0d+00; akt6_v=0.0d+00; ajs6_v=0.0d+00; ajn6_v=0.0d+00
        end if
        if(lmgd.eq.5) then
        allocate(ibp5_v(m1mv(lmgd),n1mv(lmgd)),
     *   ff5_v(m1mv(lmgd),n1mv(lmgd)),res5_v(m1mv(lmgd),n1mv(lmgd)),
     *  app5_v(m1mv(lmgd),n1mv(lmgd)),bpp5_v(m1mv(lmgd),n1mv(lmgd)),
     *  ajs5_v(m1mv(lmgd),n1mv(lmgd)),akb5_v(m1mv(lmgd),n1mv(lmgd)),
     *  ajn5_v(m1mv(lmgd),n1mv(lmgd)),akt5_v(m1mv(lmgd),n1mv(lmgd)))
        ibp5_v=0;ff5_v=0.0d+00;res5_v=0.0d+00;
        app5_v=0.0d+00;bpp_v=0.0d+00
        akb5_v=0.0d+00; akt5_v=0.0d+00; ajs5_v=0.0d+00; ajn5_v=0.0d+00
        end if
        if(lmgd.eq.4) then
        allocate(ibp4_v(m1mv(lmgd),n1mv(lmgd)),
     *   ff4_v(m1mv(lmgd),n1mv(lmgd)),res4_v(m1mv(lmgd),n1mv(lmgd)),
     *  app4_v(m1mv(lmgd),n1mv(lmgd)),bpp4_v(m1mv(lmgd),n1mv(lmgd)),
     *  ajs4_v(m1mv(lmgd),n1mv(lmgd)),akb4_v(m1mv(lmgd),n1mv(lmgd)),
     *  ajn4_v(m1mv(lmgd),n1mv(lmgd)),akt4_v(m1mv(lmgd),n1mv(lmgd)))
        ibp4_v=0;ff4_v=0.0d+00;res4_v=0.0d+00;
        app4_v=0.0d+00;bpp4_v=0.0d+00
        akb4_v=0.0d+00; akt4_v=0.0d+00; ajs4_v=0.0d+00; ajn4_v=0.0d+00
        end if
        if(lmgd.eq.3) then
        allocate(ibp3_v(m1mv(lmgd),n1mv(lmgd)),
     *   ff3_v(m1mv(lmgd),n1mv(lmgd)),res3_v(m1mv(lmgd),n1mv(lmgd)),
     *  app3_v(m1mv(lmgd),n1mv(lmgd)),bpp3_v(m1mv(lmgd),n1mv(lmgd)),
     *  ajs3_v(m1mv(lmgd),n1mv(lmgd)),akb3_v(m1mv(lmgd),n1mv(lmgd)),
     *  ajn3_v(m1mv(lmgd),n1mv(lmgd)),akt3_v(m1mv(lmgd),n1mv(lmgd)))
        ibp3_v=0;ff3_v=0.0d+00;res3_v=0.0d+00;
        app3_v=0.0d+00;bpp3_v=0.0d+00
        akb3_v=0.0d+00; akt3_v=0.0d+00; ajs3_v=0.0d+00; ajn3_v=0.0d+00
        end if
        if(lmgd.eq.2) then
        allocate(ibp2_v(m1mv(lmgd),n1mv(lmgd)),
     *   ff2_v(m1mv(lmgd),n1mv(lmgd)),res2_v(m1mv(lmgd),n1mv(lmgd)),
     *  app2_v(m1mv(lmgd),n1mv(lmgd)),bpp2_v(m1mv(lmgd),n1mv(lmgd)),
     *  ajs2_v(m1mv(lmgd),n1mv(lmgd)),akb2_v(m1mv(lmgd),n1mv(lmgd)),
     *  ajn2_v(m1mv(lmgd),n1mv(lmgd)),akt2_v(m1mv(lmgd),n1mv(lmgd)))
        ibp2_v=0;ff2_v=0.0d+00;res2_v=0.0d+00;
        app2_v=0.0d+00;bpp2_v=0.0d+00
        akb2_v=0.0d+00; akt2_v=0.0d+00; ajs2_v=0.0d+00; ajn2_v=0.0d+00
        end if
        if(lmgd.eq.1) then
        allocate(ibp1_v(m1mv(lmgd),n1mv(lmgd)),
     *   ff1_v(m1mv(lmgd),n1mv(lmgd)),res1_v(m1mv(lmgd),n1mv(lmgd)),
     *  app1_v(m1mv(lmgd),n1mv(lmgd)),bpp1_v(m1mv(lmgd),n1mv(lmgd)),
     *  ajs1_v(m1mv(lmgd),n1mv(lmgd)),akb1_v(m1mv(lmgd),n1mv(lmgd)),
     *  ajn1_v(m1mv(lmgd),n1mv(lmgd)),akt1_v(m1mv(lmgd),n1mv(lmgd)))
        ibp1_v=0;ff1_v=0.0d+00;res1_v=0.0d+00;
        app1_v=0.0d+00;bpp1_v=0.0d+00
        akb1_v=0.0d+00; akt1_v=0.0d+00; ajs1_v=0.0d+00; ajn1_v=0.0d+00
        end if
        if(lmgd.eq.0) then
        allocate(ibp0_v(m1mv(lmgd),n1mv(lmgd)),
     *   ff0_v(m1mv(lmgd),n1mv(lmgd)),res0_v(m1mv(lmgd),n1mv(lmgd)),
     *  app0_v(m1mv(lmgd),n1mv(lmgd)),bpp0_v(m1mv(lmgd),n1mv(lmgd)),
     *  ajs0_v(m1mv(lmgd),n1mv(lmgd)),akb0_v(m1mv(lmgd),n1mv(lmgd)),
     *  ajn0_v(m1mv(lmgd),n1mv(lmgd)),akt0_v(m1mv(lmgd),n1mv(lmgd)))
        ibp0_v=0;ff0_v=0.0d+00;res0_v=0.0d+00;
        app0_v=0.0d+00;bpp0_v=0.0d+00
        akb0_v=0.0d+00; akt0_v=0.0d+00; ajs0_v=0.0d+00; ajn0_v=0.0d+00
        end if
        end do
          
        do lmgd=0,levmgd
        if(lmgd.eq.9) then
        allocate(ibp9_w(m1mw(lmgd),n1mw(lmgd)))
        allocate(ff9_w(m1mw(lmgd),n1mw(lmgd)))
        allocate(res9_w(m1mw(lmgd),n1mw(lmgd)))
        allocate(app9_w(m1mw(lmgd),n1mw(lmgd)))
        allocate(bpp9_w(m1mw(lmgd),n1mw(lmgd)))
        allocate(ajs9_w(m1mw(lmgd),n1mw(lmgd)))
        allocate(akb9_w(m1mw(lmgd),n1mw(lmgd)))
        allocate(ajn9_w(m1mw(lmgd),n1mw(lmgd)))
        allocate(akt9_w(m1mw(lmgd),n1mw(lmgd)))
        ibp9_w=0;ff9_w=0.0d+00;res9_w=0.0d+00;
        app9_w=0.0d+00;bpp9_w=0.0d+00
        akb9_w=0.0d+00; akt9_w=0.0d+00; ajs9_w=0.0d+00; ajn9_w=0.0d+00
        end if
        if(lmgd.eq.8) then
        allocate(ibp8_w(m1mw(lmgd),n1mw(lmgd)),
     *   ff8_w(m1mw(lmgd),n1mw(lmgd)),res8_w(m1mw(lmgd),n1mw(lmgd)),
     *  app8_w(m1mw(lmgd),n1mw(lmgd)),bpp8_w(m1mw(lmgd),n1mw(lmgd)),
     *  ajs8_w(m1mw(lmgd),n1mw(lmgd)),akb8_w(m1mw(lmgd),n1mw(lmgd)),
     *  ajn8_w(m1mw(lmgd),n1mw(lmgd)),akt8_w(m1mw(lmgd),n1mw(lmgd)))
        ibp8_w=0;ff8_w=0.0d+00;res8_w=0.0d+00;
        app8_w=0.0d+00;bpp8_w=0.0d+00
        akb8_w=0.0d+00; akt8_w=0.0d+00; ajs8_w=0.0d+00; ajn8_w=0.0d+00
        end if
        if(lmgd.eq.7) then
        allocate(ibp7_w(m1mw(lmgd),n1mw(lmgd)),
     *   ff7_w(m1mw(lmgd),n1mw(lmgd)),res7_w(m1mw(lmgd),n1mw(lmgd)),
     *  app7_w(m1mw(lmgd),n1mw(lmgd)),bpp7_w(m1mw(lmgd),n1mw(lmgd)),
     *  ajs7_w(m1mw(lmgd),n1mw(lmgd)),akb7_w(m1mw(lmgd),n1mw(lmgd)),
     *  ajn7_w(m1mw(lmgd),n1mw(lmgd)),akt7_w(m1mw(lmgd),n1mw(lmgd)))
        ibp7_w=0;ff7_w=0.0d+00;res7_w=0.0d+00;
        app7_w=0.0d+00;bpp7_w=0.0d+00
        akb7_w=0.0d+00; akt7_w=0.0d+00; ajs7_w=0.0d+00; ajn7_w=0.0d+00
        end if
        if(lmgd.eq.6) then
        allocate(ibp6_w(m1mw(lmgd),n1mw(lmgd)),
     *   ff6_w(m1mw(lmgd),n1mw(lmgd)),res6_w(m1mw(lmgd),n1mw(lmgd)),
     *  app6_w(m1mw(lmgd),n1mw(lmgd)),bpp6_w(m1mw(lmgd),n1mw(lmgd)),
     *  ajs6_w(m1mw(lmgd),n1mw(lmgd)),akb6_w(m1mw(lmgd),n1mw(lmgd)),
     *  ajn6_w(m1mw(lmgd),n1mw(lmgd)),akt6_w(m1mw(lmgd),n1mw(lmgd)))
        ibp6_w=0;ff6_w=0.0d+00;res6_w=0.0d+00;
        app6_w=0.0d+00;bpp6_w=0.0d+00
        akb6_w=0.0d+00; akt6_w=0.0d+00; ajs6_w=0.0d+00; ajn6_w=0.0d+00
        end if
        if(lmgd.eq.5) then
        allocate(ibp5_w(m1mw(lmgd),n1mw(lmgd)),
     *   ff5_w(m1mw(lmgd),n1mw(lmgd)),res5_w(m1mw(lmgd),n1mw(lmgd)),
     *  app5_w(m1mw(lmgd),n1mw(lmgd)),bpp5_w(m1mw(lmgd),n1mw(lmgd)),
     *  ajs5_w(m1mw(lmgd),n1mw(lmgd)),akb5_w(m1mw(lmgd),n1mw(lmgd)),
     *  ajn5_w(m1mw(lmgd),n1mw(lmgd)),akt5_w(m1mw(lmgd),n1mw(lmgd)))
        ibp5_w=0;ff5_w=0.0d+00;res5_w=0.0d+00;
        app5_w=0.0d+00;bpp5_w=0.0d+00
        akb5_w=0.0d+00; akt5_w=0.0d+00; ajs5_w=0.0d+00; ajn5_w=0.0d+00
        end if
        if(lmgd.eq.4) then
        allocate(ibp4_w(m1mw(lmgd),n1mw(lmgd)),
     *   ff4_w(m1mw(lmgd),n1mw(lmgd)),res4_w(m1mw(lmgd),n1mw(lmgd)),
     *  app4_w(m1mw(lmgd),n1mw(lmgd)),bpp4_w(m1mw(lmgd),n1mw(lmgd)),
     *  ajs4_w(m1mw(lmgd),n1mw(lmgd)),akb4_w(m1mw(lmgd),n1mw(lmgd)),
     *  ajn4_w(m1mw(lmgd),n1mw(lmgd)),akt4_w(m1mw(lmgd),n1mw(lmgd)))
        ibp4_w=0;ff4_w=0.0d+00;res4_w=0.0d+00;
        app4_w=0.0d+00;bpp4_w=0.0d+00
        akb4_w=0.0d+00; akt4_w=0.0d+00; ajs4_w=0.0d+00; ajn4_w=0.0d+00
        end if
        if(lmgd.eq.3) then
        allocate(ibp3_w(m1mw(lmgd),n1mw(lmgd)),
     *   ff3_w(m1mw(lmgd),n1mw(lmgd)),res3_w(m1mw(lmgd),n1mw(lmgd)),
     *  app3_w(m1mw(lmgd),n1mw(lmgd)),bpp3_w(m1mw(lmgd),n1mw(lmgd)),
     *  ajs3_w(m1mw(lmgd),n1mw(lmgd)),akb3_w(m1mw(lmgd),n1mw(lmgd)),
     *  ajn3_w(m1mw(lmgd),n1mw(lmgd)),akt3_w(m1mw(lmgd),n1mw(lmgd)))
        ibp3_w=0;ff3_w=0.0d+00;res3_w=0.0d+00;
        app3_w=0.0d+00;bpp3_w=0.0d+00
        akb3_w=0.0d+00; akt3_w=0.0d+00; ajs3_w=0.0d+00; ajn3_w=0.0d+00
        end if
        if(lmgd.eq.2) then
        allocate(ibp2_w(m1mw(lmgd),n1mw(lmgd)),
     *   ff2_w(m1mw(lmgd),n1mw(lmgd)),res2_w(m1mw(lmgd),n1mw(lmgd)),
     *  app2_w(m1mw(lmgd),n1mw(lmgd)),bpp2_w(m1mw(lmgd),n1mw(lmgd)),
     *  ajs2_w(m1mw(lmgd),n1mw(lmgd)),akb2_w(m1mw(lmgd),n1mw(lmgd)),
     *  ajn2_w(m1mw(lmgd),n1mw(lmgd)),akt2_w(m1mw(lmgd),n1mw(lmgd)))
        ibp2_w=0;ff2_w=0.0d+00;res2_w=0.0d+00;
        app2_w=0.0d+00;bpp2_w=0.0d+00
        akb2_w=0.0d+00; akt2_w=0.0d+00; ajs2_w=0.0d+00; ajn2_w=0.0d+00
        end if
        if(lmgd.eq.1) then
        allocate(ibp1_w(m1mw(lmgd),n1mw(lmgd)),
     *   ff1_w(m1mw(lmgd),n1mw(lmgd)),res1_w(m1mw(lmgd),n1mw(lmgd)),
     *  app1_w(m1mw(lmgd),n1mw(lmgd)),bpp1_w(m1mw(lmgd),n1mw(lmgd)),
     *  ajs1_w(m1mw(lmgd),n1mw(lmgd)),akb1_w(m1mw(lmgd),n1mw(lmgd)),
     *  ajn1_w(m1mw(lmgd),n1mw(lmgd)),akt1_w(m1mw(lmgd),n1mw(lmgd)))
        ibp1_w=0;ff1_w=0.0d+00;res1_w=0.0d+00;
        app1_w=0.0d+00;bpp1_w=0.0d+00
        akb1_w=0.0d+00; akt1_w=0.0d+00; ajs1_w=0.0d+00; ajn1_w=0.0d+00
        end if
        if(lmgd.eq.0) then
        allocate(ibp0_w(m1mw(lmgd),n1mw(lmgd)),
     *   ff0_w(m1mw(lmgd),n1mw(lmgd)),res0_w(m1mw(lmgd),n1mw(lmgd)),
     *  app0_w(m1mw(lmgd),n1mw(lmgd)),bpp0_w(m1mw(lmgd),n1mw(lmgd)),
     *  ajs0_w(m1mw(lmgd),n1mw(lmgd)),akb0_w(m1mw(lmgd),n1mw(lmgd)),
     *  ajn0_w(m1mw(lmgd),n1mw(lmgd)),akt0_w(m1mw(lmgd),n1mw(lmgd)))
        ibp0_w=0;ff0_w=0.0d+00;res0_w=0.0d+00;
        app0_w=0.0d+00;bpp0_w=0.0d+00
        akb0_w=0.0d+00; akt0_w=0.0d+00; ajs0_w=0.0d+00; ajn0_w=0.0d+00
        end if
        end do
        return;end

c===============================c
        subroutine dalmgv
c===============================c
        use mdumgv;use mdupnt
        include 'table.prc'
        include 'table.gdm'
c------------------------------c
        deallocate(ff0,app0,bpp0,res0,akb0,ajs0,akt0,ajn0,ibp0)
        deallocate(ff1,app1,bpp1,res1,akb1,ajs1,akt1,ajn1,ibp1)
        deallocate(ff2,app2,bpp2,res2,akb2,ajs2,akt2,ajn2,ibp2)
        deallocate(ff3,app3,bpp3,res3,akb3,ajs3,akt3,ajn3,ibp3)
        deallocate(ff4,app4,bpp4,res4,akb4,ajs4,akt4,ajn4,ibp4)
        deallocate(ff5,app5,bpp5,res5,akb5,ajs5,akt5,ajn5,ibp5)
        deallocate(ff6,app6,bpp6,res6,akb6,ajs6,akt6,ajn6,ibp6)
        deallocate(ff7,app7,bpp7,res7,akb7,ajs7,akt7,ajn7,ibp7)
        deallocate(ff8,app8,bpp8,res8,akb8,ajs8,akt8,ajn8,ibp8)

        deallocate(ff0_v,app0_v,bpp0_v,res0_v,akb0_v,ajs0_v,
     *           akt0_v,ajn0_v,ibp0_v)
        deallocate(ff1_v,app1_v,bpp1_v,res1_v,akb1_v,ajs1_v,
     *           akt1_v,ajn1_v,ibp1_v)
        deallocate(ff2_v,app2_v,bpp2_v,res2_v,akb2_v,ajs2_v,akt2_v,
     *           ajn2_v,ibp2_v)
        deallocate(ff3_v,app3_v,bpp3_v,res3_v,akb3_v,ajs3_v,akt3_v,
     *           ajn3_v,ibp3_v)
        deallocate(ff4_v,app4_v,bpp4_v,res4_v,akb4_v,ajs4_v,akt4_v,
     *           ajn4_v,ibp4_v)
        deallocate(ff5_v,app5_v,bpp5_v,res5_v,akb5_v,ajs5_v,akt5_v,
     *           ajn5_v,ibp5_v)
        deallocate(ff6_v,app6_v,bpp6_v,res6_v,akb6_v,ajs6_v,akt6_v,
     *           ajn6_v,ibp6_v)
        deallocate(ff7_v,app7_v,bpp7_v,res7_v,akb7_v,ajs7_v,akt7_v,
     *           ajn7_v,ibp7_v)
        deallocate(ff8_v,app8_v,bpp8_v,res8_v,akb8_v,ajs8_v,akt8_v,
     *           ajn8_v,ibp8_v)

        deallocate(ff0_w,app0_w,bpp0_w,res0_w,akb0_w,ajs0_w,akt0_w,
     *           ajn0_w,ibp0_w)
        deallocate(ff1_w,app1_w,bpp1_w,res1_w,akb1_w,ajs1_w,akt1_w,
     *           ajn1_w,ibp1_w)
        deallocate(ff2_w,app2_w,bpp2_w,res2_w,akb2_w,ajs2_w,akt2_w,
     *           ajn2_w,ibp2_w)
        deallocate(ff3_w,app3_w,bpp3_w,res3_w,akb3_w,ajs3_w,akt3_w,
     *           ajn3_w,ibp3_w)
        deallocate(ff4_w,app4_w,bpp4_w,res4_w,akb4_w,ajs4_w,akt4_w,
     *           ajn4_w,ibp4_w)
        deallocate(ff5_w,app5_w,bpp5_w,res5_w,akb5_w,ajs5_w,akt5_w,
     *           ajn5_w,ibp5_w)
        deallocate(ff6_w,app6_w,bpp6_w,res6_w,akb6_w,ajs6_w,akt6_w,
     *           ajn6_w,ibp6_w)
        deallocate(ff7_w,app7_w,bpp7_w,res7_w,akb7_w,ajs7_w,akt7_w,
     *           ajn7_w,ibp7_w)
        deallocate(ff8_w,app8_w,bpp8_w,res8_w,akb8_w,ajs8_w,akt8_w,
     *           ajn8_w,ibp8_w)
        return;end
