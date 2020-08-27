
c===============================c
        subroutine gddata_re200
c===============================c
        use mdugrd
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
c-------grid one-----------c
        call xugrid(iu0,imx,imx,x0,xl,xu1)
        open(199,file='re200_griddata.p2d',form='unformatted')
        read(199)(yv1(j),j=1,258)
        close(199)
        open(199,file='re200_griddata.p2d',form='unformatted')
        read(199)(zw1(j),j=1,258)
        close(199)
        
!        call yzgrid
!     *  (jv0,jmx,jmx,y0,yl,yv1,ymb,yme,jb1,je1,jmd1,cs1,jmd2,cs2)
!        call yzgrid
!     *  (kw0,kmx,kmx,z0,zl,zw1,zmb,zme,kb1,ke1,kmd1,cs1,kmd2,cs2)
        return
        end
        
c===============================c
        subroutine gddata_re400
c===============================c
        use mdugrd
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
c-------grid one-----------c
        call xugrid(iu0,imx,imx,x0,xl,xu1)
        open(199,file='griddata.txt')
        read(199,*)yv1
        close(199)
        open(190,file='griddata.txt')
        read(190,*)zw1
        close(190)
        
!        call yzgrid
!     *  (jv0,jmx,jmx,y0,yl,yv1,ymb,yme,jb1,je1,jmd1,cs1,jmd2,cs2)
!        call yzgrid
!     *  (kw0,kmx,kmx,z0,zl,zw1,zmb,zme,kb1,ke1,kmd1,cs1,kmd2,cs2)
        return
        end
        
        subroutine xugrid(i1,l1,ni,x0,xl,xu)
        implicit double precision (a-h,o-z)
        dimension xu(ni)
        xu(i1)=x0;dx=(xl-x0)/float(l1-i1)
        do i=i1+1,l1;xu(i)=xu(i-1)+dx;end do
        return
        end

        subroutine yzgrid
     *  (j1,m1,nj,y0,yl,yv,yb,ye,jmb,jme,jmd1,cs1,jmd2,cs2)
        implicit double precision (a-h,o-z)
        dimension yv(nj)
        yv(j1)=-1.0d+00;dy=1.0d+00/float(jmd1-j1)
        do j=j1+1,jmd1;yv(j)=yv(j-1)+dy;enddo

        dh=5.0d-01*(yb-y0)
        do j=j1,jmd1
        yv(j)=dh*(1.0d+00+dtanh(cs1*yv(j))/dtanh(cs1))
        end do

        yv(jmd1)=0.0d+00;dy=1.0d+00/float(jmb-jmd1)
        do j=jmd1+1,jmb;yv(j)=yv(j-1)-dy;enddo

        dh=5.0d-01*(yb-y0)
        do j=jmd1,jmb
        yv(j)=dh*(1.0d+00-dtanh(cs2*yv(j))/dtanh(cs2))
        end do

        yv(jmb)=-1.0d+00;dy=2.0d+00/float(jme-jmb)
        do j=jmb+1,jme;yv(j)=yv(j-1)+dy;enddo

        dh=5.0d-01*(ye-yb)
        do j=jmb,jme
        yv(j)=yb+dh*(1.0d+00+dtanh(cs2*yv(j))/dtanh(cs2))
        end do

        yv(jme)=-1.0d+00;dy=1.0d+00/float(jmd2-jme)
        do j=jme+1,jmd2;yv(j)=yv(j-1)+dy;enddo

        dh=5.0d-01*(yl-ye)
        do j=jme,jmd2
        yv(j)=ye+dh*(1.0d+00+dtanh(cs2*yv(j))/dtanh(cs2))
        end do

        yv(jmd2)=0.0d+00;dy=1.0d+00/float(m1-jmd2)
        do j=jmd2+1,m1;yv(j)=yv(j-1)-dy;enddo

        do j=jmd2,m1
        yv(j)=ye+dh*(1.0d+00-dtanh(cs1*yv(j))/dtanh(cs1))
        end do
        return
        end

c===============================c
        subroutine stdata_re400
c===============================c
        use mduvar;use mdumgv
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
        include 'table.gdm'
c-------stress and vel----------c

c-----------read pressure-----------c

        open(11,file='re400_p.p2d',form='unformatted')
        read(11) jmx,kmx,nvar
        read(11)((pi1(j,k),j=1,jmx),k=1,kmx)
        close(11)
c---------------------------------c
c----------read velocity-----------------c
c-------------Re=400,DNS------------------------c
        open(11,file='re400_u.p2d',form='unformatted')
        read(11) jmx,kmx,nvar
        read(11)((ui1(j,k),j=1,jmx),k=1,kmx)
        close(11)
!!
        open(11,file='re400_v.p2d',form='unformatted')
        read(11) jmx,kmx,nvar
        read(11)((vi1(j,k),j=1,jmx),k=1,kmx)
        close(11)

        open(11,file='re400_w.p2d',form='unformatted')
        read(11) jmx,kmx,nvar
        read(11)((wi1(j,k),j=1,jmx),k=1,kmx)
        close(11)
c--------------------------------------------------------------c
c-------------------Re=200,DNS interp-----------------------------c
!        open(191,file='v258to514ac.dat',form='formatted')
!        do j=1,514
!	    do k=1,514
!	        read(191,*)vi1(j,k)
!            enddo
!        enddo
!        close(191)
!        
!        open(192,file='w258to514ac.dat',form='formatted')
!        do j=1,514
!	    do k=1,514
!	        read(192,*)wi1(j,k)
!            enddo
!        enddo
!        close(192)
        
c----------------read stress---------------c
c----------------DNS stress--------------c
!        open(11,file='re400_vpvp.p2d',form='unformatted')
!        read(11) jmx,kmx,nvar
!        read(11)((vpvp(j,k),j=1,jmx),k=1,kmx)
!        close(11)
!                
!        open(11,file='re400_vpwp.p2d',form='unformatted')
!        read(11) jmx,kmx,nvar
!        read(11)((vpwp(j,k),j=1,jmx),k=1,kmx)
!        close(11)
!        
!        open(11,file='re400_wpwp.p2d',form='unformatted')
!        read(11) jmx,kmx,nvar
!        read(11)((wpwp(j,k),j=1,jmx),k=1,kmx)
!        close(11)
!        
        open(11,file='re400_upvp.p2d',form='unformatted')
        read(11) jmx,kmx,nvar
        read(11)((upvp(j,k),j=1,jmx),k=1,kmx)
        close(11)
                
        open(11,file='re400_upwp.p2d',form='unformatted')
        read(11) jmx,kmx,nvar
        read(11)((upwp(j,k),j=1,jmx),k=1,kmx)
        close(11)
       
c---------------ml stress-----------------c
        open(11,file='sym_400_vpvp.txt',form='formatted')
        read(11,*) vpvp  
        close(11)
        
        open(192,file='sym_400_vpwp.txt',form='formatted')
        read(192,*) vpwp
        close(192)
        
        open(192,file='sym_400_wpwp.txt',form='formatted')
        read(192,*) wpwp
        close(192)
        end
        
c===============================c
        subroutine stdata_re200
c===============================c
        use mduvar;use mdumgv
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
        include 'table.gdm'
c-------stress and vel----------c
        open(11,file='re200_p.p2d',form='unformatted')
        read(11) jmx,kmx,nvar
        read(11)((pi1(j,k),j=1,jmx),k=1,kmx)
        close(11)
        
        open(11,file='re200_v.p2d',form='unformatted')
        read(11) jmx,kmx,nvar
        read(11)((vi1(j,k),j=1,jmx),k=1,kmx)
        close(11)

        open(11,file='re200_w.p2d',form='unformatted')
        read(11) jmx,kmx,nvar
        read(11)((wi1(j,k),j=1,jmx),k=1,kmx)
        close(11)

        open(11,file='re200_vpvp.p2d',form='unformatted')
        read(11) jmx,kmx,nvar
        read(11)((vpvp(j,k),j=1,jmx),k=1,kmx)
        close(11)
                
        open(11,file='re200_vpwp.p2d',form='unformatted')
        read(11) jmx,kmx,nvar
        read(11)((vpwp(j,k),j=1,jmx),k=1,kmx)
        close(11)
        
        open(11,file='re200_wpwp.p2d',form='unformatted')
        read(11) jmx,kmx,nvar
        read(11)((wpwp(j,k),j=1,jmx),k=1,kmx)
        close(11)
        
        open(11,file='re200_upvp.p2d',form='unformatted')
        read(11) jmx,kmx,nvar
        read(11)((upvp(j,k),j=1,jmx),k=1,kmx)
        close(11)
                
        open(11,file='re200_upwp.p2d',form='unformatted')
        read(11) jmx,kmx,nvar
        read(11)((upwp(j,k),j=1,jmx),k=1,kmx)
        close(11)
        
        end
        
c===============================c
        subroutine augrid
c===============================c
        use mdugrd
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
c-------grid level one----------c
        call mngrid(ip0,imx,imx,xm1,xu1,dxm,dxu,fx1,fx2)
        call mngrid(ip0,jmx,jmx,ym1,yv1,dym,dyv,fy1,fy2)
        call mngrid(ip0,kmx,kmx,zm1,zw1,dzm,dzw,fz1,fz2)
        return
        end

        subroutine mngrid(i1,l1,ni,xm,xu,dxm,dxu,fx1,fx2)
        implicit double precision (a-h,o-z)
        dimension xm(ni),xu(ni),fx1(ni),fx2(ni),dxm(ni),dxu(ni)
c-------------------------------c
        i2=i1+1;l2=l1-1
        xm(i1)=xu(i2)-5.0d-01*(xu(i2+1)-xu(i2))
        do i=i2,l2;xm(i)=5.0d-01*(xu(i+1)+xu(i));enddo
        xm(l1)=xu(l1)+5.0d-01*(xu(l1)-xu(l1-1))
        do i=i2,l1;dxu(i)=xm(i)-xm(i-1);enddo
        do i=i2,l2;dxm(i)=xu(i+1)-xu(i);enddo
        do i=i2,l1
        fx1(i)=(xm(i)-xu(i))/dxu(i)
        fx2(i)=(xu(i)-xm(i-1))/dxu(i)
        end do
        return
        end

c===============================c
        subroutine ibdata
c===============================c
        use mdugrd
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
c-------------------------------c
        ibp(:,:)=1
        do j=jp0,jmx;do k=kp0,kmx
        if((j.ge.jb1.and.j.le.je1-1).and.
     *     (k.ge.kb1.and.k.le.ke1-1)) ibp(j,k)=0
        end do;end do
        ibp(jp0,:)=2;ibp(jmx,:)=2
        ibp(:,kp0)=2;ibp(:,kmx)=2
        ibp(jb1:je1-1,kb1)=2;ibp(jb1:je1-1,ke1-1)=2
        ibp(jb1,kb1:ke1-1)=2;ibp(je1-1,kb1:ke1-1)=2

        ibu(:,:)=1
        do j=ju0,jmx;do k=ku0,kmx
        if((j.ge.jb1.and.j.le.je1-1).and.
     *     (k.ge.kb1.and.k.le.ke1-1)) ibu(j,k)=0
        end do;end do
        ibu(ju0,:)=2;ibu(jmx,:)=2
        ibu(:,ku0)=2;ibu(:,kmx)=2
        ibu(jb1:je1-1,kb1)=2;ibu(jb1:je1-1,ke1-1)=2
        ibu(jb1,kb1:ke1-1)=2;ibu(je1-1,kb1:ke1-1)=2

        ibv(:,:)=1;
        do j=jv0,jmx;do k=kv0,kmx
        if((j.ge.jb1.and.j.le.je1).and.
     *     (k.ge.kb1.and.k.lt.ke1)) ibv(j,k)=0
        end do;end do
        ibv(jv0,:)=2;ibv(jmx,:)=2
        ibv(:,kv0)=2;ibv(:,kmx)=2
        ibv(jb1:je1,kb1)=2;ibv(jb1:je1,ke1-1)=2
        ibv(jb1,kb1:ke1-1)=2;ibv(je1,kb1:ke1-1)=2

        ibw(:,:)=1
        do j=jw0,jmx;do k=kw0,kmx
        if((j.ge.jb1.and.j.lt.je1).and.
     *     (k.ge.kb1.and.k.le.ke1)) ibw(j,k)=0
        end do;end do
        ibw(jw0,:)=2;ibw(jmx,:)=2
        ibw(:,kw0)=2;ibw(:,kmx)=2
        ibw(jb1:je1-1,kb1)=2;ibw(jb1:je1-1,ke1)=2
        ibw(jb1,kb1:ke1)=2;ibw(je1-1,kb1:ke1)=2
        return
        end
