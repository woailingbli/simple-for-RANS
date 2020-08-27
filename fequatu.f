c===============================c
        subroutine ulamin
c===============================c
        use mdugrd;use mduvar;use mdumgv;use mdupnt
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
c-------------------------------c
        dimension Fjn(jmx,kmx),Fjs(jmx,kmx),Fkt(jmx,kmx),Fkb(jmx,kmx)
        dimension Djn(jmx,kmx),Djs(jmx,kmx),Dkt(jmx,kmx),Dkb(jmx,kmx)
        dimension upvpn(jmx,kmx),upvps(jmx,kmx),vel(jmx,kmx)
        dimension upwpt(jmx,kmx),upwpb(jmx,kmx)
        
        if(levmgd.eq.9) then
        ibp9=ibp;ibpm0=>ibp9
         ffm0=>ff9;appm0=>app9;bppm0=>bpp9;resm0=>res9
        ajsm0=>ajs9;ajnm0=>ajn9;akbm0=>akb9;aktm0=>akt9
        end if
        if(levmgd.eq.8) then
        ibp8=ibp;ibpm0=>ibp8
         ffm0=>ff8;appm0=>app8;bppm0=>bpp8;resm0=>res8
        ajsm0=>ajs8;ajnm0=>ajn8;akbm0=>akb8;aktm0=>akt8
        end if
c-------------------------------c
        j2=ju0+1;k2=ku0+1;m2=jmx-1;n2=kmx-1

        
        do j=j2,m2;do k=k2,n2
        if(ibpm0(j,k).eq.1) then
                               
        Fjn(j,k) = dzm(k)*vi1(j+1,k)
        Fjs(j,k) = dzm(k)*vi1(j,k)
        Fkt(j,k) = dym(j)*wi1(j,k+1)
        Fkb(j,k) = dym(j)*wi1(j,k)                      

        Djn(j,k) = dzm(k)/(reno*dyv(j+1))
        Djs(j,k) = dzm(k)/(reno*dyv(j))
        Dkt(j,k) = dym(j)/(reno*dzw(k+1))
        Dkb(j,k) = dym(j)/(reno*dzw(k))
  
        ajsm0(j,k) = Djs(j,k) + max(0.0d-00,Fjs(j,k))
        ajnm0(j,k) = Djn(j,k) + max(0.0d-00,-Fjn(j,k))
        akbm0(j,k) = Dkb(j,k) + max(0.0d-00,Fkb(j,k))
        aktm0(j,k) = Dkt(j,k) + max(0.0d-00,-Fkt(j,k))

        else
        ajsm0(j,k)=0.0d+00;ajnm0(j,k)=0.0d+00
        akbm0(j,k)=0.0d+00;aktm0(j,k)=0.0d+00
        end if
        end do;end do
        
        do j=j2,m2;do k=k2,n2
        if(ibpm0(j,k).eq.1) then
        upvps(j,k)=upvp(j-1,k)*fy2(j)+upvp(j,k)*fy1(j)
        upvpn(j,k)=upvp(j,k)*fy2(j+1)+upvp(j+1,k)*fy1(j+1)
        
        upwpb(j,k)=upwp(j,k-1)*fz2(k)+upwp(j,k)*fz1(k)
        upwpt(j,k)=upwp(j,k)*fz2(k+1)+upwp(j,k+1)*fz1(k+1)

        appm0(j,k)=ajsm0(j,k)+ajnm0(j,k)+akbm0(j,k)+aktm0(j,k)
        bppm0(j,k)=-(upvpn(j,k)-upvps(j,k))*dzm(k)
     *             -(upwpt(j,k)-upwpb(j,k))*dym(j)
     *             +dym(k)*dzm(j)*2.0d-00
        else
        appm0(j,k)=1.0d+00;bppm0(j,k)=0.0d+00
        end if
        end do;end do
        
c-------outer wall boundary-------c
        j=j2
        do k=k2,n2;appm0(j,k)=appm0(j,k)+ajsm0(j,k);enddo
        ffm0(j-1,k2:n2)=ubo;ajsm0(j,k2:n2)=0.0d+00
        j=m2
        do k=k2,n2;appm0(j,k)=appm0(j,k)+ajnm0(j,k);enddo
        ffm0(j+1,k2:n2)=ubo;ajnm0(j,k2:n2)=0.0d+00
        k=k2
        do j=j2,m2;appm0(j,k)=appm0(j,k)+akbm0(j,k);enddo
        ffm0(j2:m2,k-1)=ubo;akbm0(j2:m2,k)=0.0d+00
        k=n2
        do j=j2,m2;appm0(j,k)=appm0(j,k)+aktm0(j,k);enddo
        ffm0(j2:m2,k+1)=ubo;aktm0(j2:m2,k)=0.0d+00
c-------inner wall boundary-------c
        j=je1
        do k=kb1,ke1-1;appm0(j,k)=appm0(j,k)+ajsm0(j,k);enddo
        ffm0(j-1,kb1:ke1-1)=ubi;ajsm0(j,kb1:ke1-1)=0.0d+00
        j=jb1-1
        do k=kb1,ke1-1;appm0(j,k)=appm0(j,k)+ajnm0(j,k);enddo
        ffm0(j,kb1:ke1-1)=ubi;ajnm0(j,kb1:ke1-1)=0.0d+00

        k=ke1
        do j=jb1,je1-1;appm0(j,k)=appm0(j,k)+akbm0(j,k);enddo
        ffm0(jb1:je1-1,k-1)=ubi;akbm0(jb1:je1-1,k)=0.0d+00
        k=kb1-1
        do j=jb1,je1-1;appm0(j,k)=appm0(j,k)+aktm0(j,k);enddo
        ffm0(jb1:je1-1,k)=ubi;aktm0(jb1:je1-1,k)=0.0d+00
c-------solving laminar profile-------c
        appm0(2,2) = 1.0d-00
        bppm0(2,2) = 0.0d-00
        ajsm0(2,2)=0.0d-00
        ajnm0(2,2)=0.0d-00
        akbm0(2,2)=0.0d-00
        aktm0(2,2)=0.0d-00

        call amgcof(levmgd)
        call solvem(levmgd,modeui,ressumb,ressume)


        do j=ju0,jmx;do k=ku0,kmx
        ui1(j,k)=ffm0(j,k)
        end do;end do
        
        pm = pi1
        do j=j2,m2;do k=k2,n2
        um(j,k)=ui1(j,k)
        vm(j,k)=5.0d-01*(vi1(j,k)+vi1(j+1,k))
        wm(j,k)=5.0d-01*(wi1(j,k)+wi1(j,k+1))
        vel(j,k)=sqrt(vm(j,k)*vm(j,k)+wm(j,k)*wm(j,k))
        end do;end do
                
        open(9,file='simple_re200_ml_stress_alpha_loop2000.dat',
     *           form='formatted')
	write(9,*) 'TITLE    ="Plot3D DataSet"'
        write(9,*) 'VARIABLES = "y" "z" "p" "u" "v" "w" "vel"'
        write(9,*) 'DATASETAUXDATA Common.SpeedOfSound="1.0"'
        write(9,*) 'DATASETAUXDATA Common.VectorVarsAreVelocity="FALSE"'
        write(9,*) 'ZONE T="Zone-original grid"'
        write(9,*) 'STRANDID=0, SOLUTIONTIME=0'
        write(9,*) 'I=514, J=514, K=1, ZONETYPE=Ordered'
        write(9,*) 'DATAPACKING=POINT'
        write(9,*) 'DT=(SINGLE SINGLE SINGLE SINGLE 
     *                 SINGLE SINGLE SINGLE)'
	do k=1,kmx
	do j=1,jmx
        write(9,*) ym1(j),zm1(k),pm(j,k),um(j,k),
     *             vm(j,k),wm(j,k),vel(j,k)
        enddo
        enddo
        close(9)
        
        return
        end
