c===============================c
        subroutine plamin
c===============================c
        use mdugrd;use mduvar;use mdumgv;use mdupnt
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
c-------------------------------c
        
        if(levmgd.eq.9) then
        ibp9=ibp;ibpm0=>ibp9
         ffm0=>ff9 ;appm0=>app9;bppm0=>bpp9;resm0=>res9
        ajsm0=>ajs9;ajnm0=>ajn9;akbm0=>akb9;aktm0=>akt9
        end if
        if(levmgd.eq.8) then
        ibp8=ibw;ibpm0=>ibp8
         ffm0=>ff8 ;appm0=>app8;bppm0=>bpp8;resm0=>res8
        ajsm0=>ajs8;ajnm0=>ajn8;akbm0=>akb8;aktm0=>akt8
        end if
c-------------------------------c
        j2=ju0+1;k2=ku0+1;m2=jmx-1;n2=kmx-1
        j1=ju0;k1=ku0;m1=jmx;n1=kmx
        
        do j=j2,m2;do k=k2,n2
        if(ibpm0(j,k).eq.1) then
                               
        ajsm0(j,k)=dzm(k)*dzm(k)/apc1(j,k)*alpha_v
        ajnm0(j,k)=dzm(k)*dzm(k)/apc1(j+1,k)*alpha_v
        akbm0(j,k)=dym(j)*dym(j)/apc2(j,k)*alpha_w
        aktm0(j,k)=dym(j)*dym(j)/apc2(j,k+1)*alpha_w
        
        else
        ajsm0(j,k)=0.0d+00;ajnm0(j,k)=0.0d+00
        akbm0(j,k)=0.0d+00;aktm0(j,k)=0.0d+00
        end if
        end do;end do



c-------outer wall boundary-------c
        j=j2
        ajsm0(j,k2:n2)=0.0d+00
        j=m2
        ajnm0(j,k2:n2)=0.0d+00
        k=k2
        akbm0(j2:m2,k)=0.0d+00
        k=n2
        aktm0(j2:m2,k)=0.0d+00
c-------inner wall boundary-------c
        j=je1
        do k=kb1,ke1-1
        appm0(j,k)=appm0(j,k)+ajsm0(j,k)
        enddo
        ajsm0(j,kb1:ke1-1)=0.0d+00
        j=jb1-1
        do k=kb1,ke1-1
        appm0(j,k)=appm0(j,k)+ajnm0(j,k)
        enddo
        ajnm0(j,kb1:ke1-1)=0.0d+00

        k=ke1
        do j=jb1,je1-1
        appm0(j,k)=appm0(j,k)+akbm0(j,k)
        enddo
        akbm0(jb1:je1-1,k)=0.0d+00
        k=kb1-1
        do j=jb1,je1-1
        appm0(j,k)=appm0(j,k)+aktm0(j,k)
        enddo
        aktm0(jb1:je1-1,k)=0.0d+00
                
        do j=j2,m2;do k=k2,n2
        if(ibpm0(j,k).eq.1) then
        appm0(j,k)=ajsm0(j,k)+ajnm0(j,k)+akbm0(j,k)+aktm0(j,k)
        bppm0(j,k)=-(vv1(j+1,k)-vv1(j,k))*dzm(k)
     &             -(ww1(j,k+1)-ww1(j,k))*dym(j)
        else
        appm0(j,k)=1.0d+00;bppm0(j,k)=0.0d+00
        end if
        end do;end do
        appm0(513,513) = 1.0d-00
        bppm0(513,513) = 0.0d-00
        ajsm0(513,513)=0.0d-00
        ajnm0(513,513)=0.0d-00
        akbm0(513,513)=0.0d-00
        aktm0(513,513)=0.0d-00

c-------solving laminar profile-------c

        call amgcof(levmgd)
        call solvem(levmgd,modeui,ressumb,ressume)

        do j=ju0,jmx;do k=ku0,kmx
        pre1(j,k)=ffm0(j,k)
        end do;end do

        do j=j2,m2
        pre1(j,k1)=pre1(j,k2)
     *            +dzw(k2)/dzw(k2+1)*(pre1(j,k2)-pre1(j,k2+1))
        pre1(j,n1)=pre1(j,n2)
     *            +dzw(n2+1)/dzw(n2)*(pre1(j,n2)-pre1(j,n2-1))
        enddo

        do k=k2,n2
        pre1(j1,k)=pre1(j2,k)
     *            +dyv(j2)/dyv(j2+1)*(pre1(j2,k)-pre1(j2+1,k))
        pre1(m1,k)=pre1(m2,k)
     *            +dyv(m2+1)/dyv(m2)*(pre1(m2,k)-pre1(m2-1,k))
        enddo
        
        prej=pre1(j2,k1)
     *      +dyv(j2)/dyv(j2+1)*(pre1(j2,k1)-pre1(j2+1,k1))
        prek=pre1(j1,k2)
     *      +dzw(k2)/dzw(k2+1)*(pre1(j1,k2)-pre1(j1,k2+1))
        pre1(j1,k1)=5.0d-01*(prej+prek)

        prej=pre1(j2,n1)
     *      +dyv(j2)/dyv(j2+1)*(pre1(j2,n1)-pre1(j2+1,n1))
        prek=pre1(j1,n2)
     *      +dzw(n2+1)/dzw(n2)*(pre1(j1,n2)-pre1(j1,n2-1))
        pre1(j1,n1)=5.0d-01*(prej+prek)

        prej=pre1(m2,k1)+dyv(m2+1)/dyv(m2)*(pre1(m2,k1)-pre1(m2-1,k1))
        prek=pre1(m1,k2)+dzw(k2)/dzw(k2+1)*(pre1(m1,k2)-pre1(m1,k2+1))
        pre1(m1,k1)=5.0d-01*(prej+prek)

        prej=pre1(m2,n1)+dyv(m2+1)/dyv(m2)*(pre1(m2,n1)-pre1(m2-1,n1))
        prek=pre1(m1,n2)+dzw(n2+1)/dzw(n2)*(pre1(m1,n2)-pre1(m1,n2-1))
        pre1(m1,n1)=5.0d-01*(prej+prek)

        do j=jb1+1,je1-2
        pre1(j,ke1-1)=pre1(j,ke1)
     *               +dzw(ke1)/dzw(ke1+1)*(pre1(j,ke1)-pre1(j,ke1+1))
        pre1(j,kb1)=pre1(j,kb1-1)
     *             +dzw(kb1)/dzw(kb1-1)*(pre1(j,kb1-1)-pre1(j,kb1-2))
        enddo

        do k=kb1+1,ke1-2
        pre1(je1-1,k)=pre1(je1,k)
     *               +dyv(je1)/dyv(je1+1)*(pre1(je1,k)-pre1(je1+1,k))
        pre1(jb1,k)=pre1(jb1-1,k)
     *             +dyv(jb1)/dyv(jb1-1)*(pre1(jb1-1,k)-pre1(jb1-2,k))
        enddo

        prej=pre1(je1,ke1-1)
     *      +dyv(je1)/dyv(je1+1)*(pre1(je1,ke1-1)-pre1(je1+1,ke1-1))
        prek=pre1(je1-1,ke1)
     *      +dzw(ke1)/dzw(ke1+1)*(pre1(je1-1,ke1)-pre1(je1-1,ke1+1))
        pre1(je1-1,ke1-1)=5.0d-01*(prej+prek)

        prej=pre1(je1,kb1)
     *      +dyv(je1)/dyv(je1+1)*(pre1(je1,kb1)-pre1(je1+1,kb1))
        prek=pre1(je1,kb1-1)
     *      +dzw(kb1)/dzw(kb1-1)*(pre1(je1,kb1-1)-pre1(je1,kb1-2))
        pre1(je1-1,kb1)=5.0d-01*(prej+prek)

        prej=pre1(jb1-1,ke1)
     *      +dyv(jb1)/dyv(jb1-1)*(pre1(jb1-1,ke1)-pre1(jb1-2,ke1))
        prek=pre1(jb1,ke1)
     *      +dzw(ke1)/dzw(ke1+1)*(pre1(jb1,ke1)-pre1(jb1,ke1+1))
        pre1(jb1,ke1-1)=5.0d-01*(prej+prek)

        prej=pre1(jb1-1,kb1)
     *      +dyv(jb1)/dyv(jb1-1)*(pre1(jb1-1,kb1)-pre1(jb1-2,kb1))
        prek=pre1(jb1,kb1-1)
     *      +dzw(kb1)/dzw(kb1-1)*(pre1(jb1,kb1-1)-pre1(jb1,kb1-2))
        pre1(jb1,kb1)=5.0d-01*(prej+prek)


!c-------correction of v by p-------c
        j1=2
        k1=1
        j2=j1+1
        k2=k1+1
                
        aa = 0
        jj = 0
        kk = 0
        do j=j2,m2
        do k=k2,n2
        if(ibpm0(j,k).eq.1) then
                if (abs(vi1(j,k))>aa) then
                        aa=abs(vi1(j,k))
                        jj=j
                        kk=k
                endif
        endif
        enddo
        enddo
        
        do j=j2,m2
        do k=k2,n2
        if(ibpm0(j,k).eq.1) then
        vi1(j,k)=vv1(j,k)-
     *           (pre1(j,k)-pre1(j-1,k))*dzm(k)/apc1(j,k)*alpha_v
        endif
        enddo
        enddo
        

        

        kkk= abs(abs(vi1(jj,kk))-aa)/aa * 1.0d+02
!c-------correction of w by p-------c
        j1=1
        k1=2
        j2=j1+1
        k2=k1+1
        do j=j2,m2
        do k=k2,n2
        if(ibpm0(j,k).eq.1) then
        wi1(j,k)=ww1(j,k)-
     *           (pre1(j,k)-pre1(j,k-1))*dym(j)/apc2(j,k)*alpha_w
        endif
        enddo
        enddo

        j1=1
        k1=1
        j2=j1+1
        k2=k1+1
        do j=j1,m1
        do k=k1,n1
        if(ibpm0(j,k).eq.1) then
        pi1(j,k)=pi1(j,k)+pre1(j,k)*1.0d-05*alpha_p
        endif
        enddo
        enddo        

        
!        open(9,file='mid_result_alpha.dat',
!     *           form='formatted')
!	write(9,*) 'TITLE    ="Plot3D DataSet"'
!        write(9,*) 'VARIABLES = "y""z""ap""bp0""an""as""ab""at"'
!        write(9,*) 'DATASETAUXDATA Common.SpeedOfSound="1.0"'
!        write(9,*) 'DATASETAUXDATA Common.VectorVarsAreVelocity="FALSE"'
!        write(9,*) 'ZONE T="Zone-original grid"'
!        write(9,*) 'STRANDID=0, SOLUTIONTIME=0'
!        write(9,*) 'I=514, J=514, K=1, ZONETYPE=Ordered'
!        write(9,*) 'DATAPACKING=POINT'
!        write(9,*) 'DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE)'
!	do k=1,kmx
!	do j=1,jmx
!        write(9,*) ym1(j),zm1(k),appm0(j,k),bppm0(j,k),
!     *           ajsm0(j,k),ajnm0(j,k),akbm0(j,k),aktm0(j,k)
!        enddo
!        enddo
!        close(9)
!        
!        open(9,file='simple_re400_pre_correct_alpha_loop10.dat',
!     *           form='formatted')
!	write(9,*) 'TITLE    ="Plot3D DataSet"'
!        write(9,*) 'VARIABLES = "y" "z" "delta_p" "p_corrected" "v" "w"'
!        write(9,*) 'DATASETAUXDATA Common.SpeedOfSound="1.0"'
!        write(9,*) 'DATASETAUXDATA Common.VectorVarsAreVelocity="FALSE"'
!        write(9,*) 'ZONE T="Zone-original grid"'
!        write(9,*) 'STRANDID=0, SOLUTIONTIME=0'
!        write(9,*) 'I=514, J=514, K=1, ZONETYPE=Ordered'
!        write(9,*) 'DATAPACKING=POINT'
!        write(9,*) 'DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE)'
!	do k=1,kmx
!	do j=1,jmx
!        write(9,*) ym1(j),zm1(k),pre1(j,k),pi1(j,k),vi1(j,k),wi1(j,k)
!        enddo
!        enddo
!        close(9)
                

        
        return
        end
                
