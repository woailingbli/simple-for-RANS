!c===============================c
        subroutine equatv
!c===============================c
        use varalc
        include 'table.prc'
        include 'table.cre'
        include 'table.v_gd11'
        include 'table.lei'

        dimension Fn(nj1,nk1),Fs(nj1,nk1),Ft(nj1,nk1),Fb(nj1,nk1),vpwpt(nj1,nk1),vpwpb(nj1,nk1)
        dimension Dn(nj1,nk1),Ds(nj1,nk1),Dt(nj1,nk1),Db(nj1,nk1),deltaF(nj1,nk1),deltaA(nj1,nk1)

!c-------------------------------c
        j11=2;k11=1
        j21=j11+1;k21=k11+1
        m21=m11-1;n21=n11-1
!c-------------------------------c
!c-------doing vyy,ajs,ajn-------c

        do j=j21,m21
            do k=k21,n21
                if((j.ge.jb1.and.j.lt.je1).and.(k.ge.kb1.and.k.lt.ke1)) cycle
                       
                        Fn(j,k) = dzm1(k)*(vi1(j+1,k)+vi1(j,k))/2
                        Fs(j,k) = dzm1(k)*(vi1(j-1,k)+vi1(j,k))/2
                        Ft(j,k) = dyv1(j)*(wi1(j,k+1)*fy21(j)+wi1(j-1,k+1)*fy11(j))
                        Fb(j,k) = dyv1(j)*(wi1(j,k)*fy21(j)+wi1(j-1,k)*fy11(j)) 
                        
                        Dn(j,k) = dzm1(k)/(reno*dym1(j))
                        Ds(j,k) = dzm1(k)/(reno*dym1(j-1))
                        Dt(j,k) = dyv1(j)/(reno*dzw1(k+1))
                        Db(j,k) = dyv1(j)/(reno*dzw1(k))
                        
                        ajs1(j,k) = Ds(j,k) + max(0.0d-01 , Fs(j,k))
                        ajn1(j,k) = Dn(j,k) + max(0.0d-01 , -Fn(j,k))
                        akb1(j,k) = Db(j,k) + max(0.0d-01 , Fb(j,k))
                        akt1(j,k) = Dt(j,k) + max(0.0d-01 , -Ft(j,k))
                        deltaF(j,k) = Fn(j,k) - Fs(j,k) + Ft(j,k) - Fb(j,k)
                      
                        apc1(j,k) = ajs1(j,k) + ajn1(j,k) + akb1(j,k) + akt1(j,k)
                        
                        deltaA(j,k) = (ajs1(j,k)+ajn1(j,k)+akb1(j,k)+akt1(j,k))-apc1(j,k)
                        vpwpt(j,k) = (vpwp(j-1,k)*fy11(j)+vpwp(j,k)*fy21(j))*fz11(k+1)+(vpwp(j-1,k+1)*fy11(j)+vpwp(j,k+1)*fy21(j))*fz21(k+1)
                        
                        vpwpb(j,k) = (vpwp(j-1,k-1)*fy11(j)+vpwp(j,k-1)*fy21(j))*fz11(k)+(vpwp(j-1,k)*fy11(j)+vpwp(j,k)*fy21(j))*fz21(k)

                        if (k==n21.or.((k==kb1-1).and.(j.ge.jb1.and.j.lt.je1))) vpwpt(j,k)=0.0d+00
                        if (k==k21.or.((k==ke1).and.(j.ge.jb1.and.j.lt.je1))) vpwpb(j,k)=0.0d+00
                        
                        aiw1(j,k)=0
                        aie1(j,k)=0
                        bpp1(j,k) = - (pre1(j,k)-pre1(j-1,k)) * dzm1(k)  * 1.0d-4
                        bpp1(j,k) = bpp1(j,k) - (vpvp(j,k) - vpvp(j-1,k)) * dzm1(k) - (vpwpt(j,k) - vpwpb(j,k)) * dyv1(j)
            enddo
        enddo
!!!c-------------------------------c
!!!c********************************c
        call bondv2(j11,k11,m11,n11,jb1,je1,kb1,ke1,ni1,nj1,nk1,vv1,apc1,bpp1,aiw1,ajs1,akb1,aie1,ajn1,akt1)
!!!c********************************c
        mode=1

        if(levgrd.eq.1) then
            call v_solve1(j11,k11,m11,n11,jb1,je1,kb1,ke1,ni1,nj1,nk1,mode,epsp11,epsp21,epsp31,epsp41,ressum0,ressum, &
                        vv1,res1,apc1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1)
        end if
        if(levgrd.eq.2) then
            call v_solve2(j11,k11,m11,n11,jb1,je1,kb1,ke1,ni1,nj1,nk1,mode,epsp11,epsp21,epsp31,epsp41,ressum0,ressum, &
                        vv1,res1,apc1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1)
        end if
        if(levgrd.eq.3) then
            call v_solve3(j11,k11,m11,n11,jb1,je1,kb1,ke1,ni1,nj1,nk1,mode,epsp11,epsp21,epsp31,epsp41,ressum0,ressum, &
                        vv1,res1,apc1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1)
        end if
        if(levgrd.eq.4) then
            call v_solve4(j11,k11,m11,n11,jb1,je1,kb1,ke1,ni1,nj1,nk1,mode,epsp11,epsp21,epsp31,epsp41,ressum0,ressum, &
                        vv1,res1,apc1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1)
        end if
        if(levgrd.eq.5) then
            call v_solve5(j11,k11,m11,n11,jb1,je1,kb1,ke1,ni1,nj1,nk1,mode,epsp11,epsp21,epsp31,epsp41,ressum0,ressum, &
                        vv1,res1,apc1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1)
        end if
        if(levgrd.eq.6) then
            call v_solve6(j11,k11,m11,n11,jb1,je1,kb1,ke1,ni1,nj1,nk1,mode,epsp11,epsp21,epsp31,epsp41,ressum0,ressum, &
                        vv1,res1,apc1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1)
        end if
        if(levgrd.eq.7) then
            call v_solve7(j11,k11,m11,n11,jb1,je1,kb1,ke1,ni1,nj1,nk1,mode,epsp11,epsp21,epsp31,epsp41,ressum0,ressum,&
                        vv1,res1,apc1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1)
        end if
        if(levgrd.eq.8) then
            call v_solve8(j11,k11,m11,n11,jb1,je1,kb1,ke1,ni1,nj1,nk1,mode,epsp11,epsp21,epsp31,epsp41,ressum0,ressum,&
                        vv1,res1,apc1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1)
        end if
        if(levgrd.eq.9) then
            call v_solve9(j11,k11,m11,n11,jb1,je1,kb1,ke1,ni1,nj1,nk1,mode,epsp11,epsp21,epsp31,epsp41,ressum0,ressum,&
                        vv1,res1,apc1,bpp1,aiw1,aie1,ajs1,ajn1,akb1,akt1)
        end if

	open(91,file='..\output\RANS-vw\vv_simple_test.dat',form='formatted')
	write(91,*) 'TITLE    ="Plot3D DataSet"'
        write(91,*) 'VARIABLES = "y" "z" "v"'
        write(91,*) 'DATASETAUXDATA Common.SpeedOfSound="1.0"'
        write(91,*) 'DATASETAUXDATA Common.VectorVarsAreVelocity="FALSE"'
        write(91,*) 'ZONE T="Zone-original grid"'
        write(91,*) 'STRANDID=0, SOLUTIONTIME=0'
        write(91,*) 'I=514, J=514, K=1, ZONETYPE=Ordered'
        write(91,*) 'DATAPACKING=POINT'
        write(91,*) 'DT=(SINGLE SINGLE SINGLE)'
	do k=1,514
	    do j=1,514
                write(91,*) ym1(j),zm1(k),vv1(j,k)
            enddo
        enddo
        close(91)
1000    return
end
