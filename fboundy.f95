!c*******************************c
        subroutine bondv2(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,vv1,app,con,aim,ajm,akm,aip,ajp,akp)
!c*******************************c
        include 'table.prc'
        dimension app(nj,nk),con(nj,nk),aim(nj,nk),aip(nj,nk),ajm(nj,nk),ajp(nj,nk),akm(nj,nk),akp(nj,nk)
        dimension vv1(nj,nk)
!c*******************************c
        j2=j1+1
        k2=k1+1
        m2=m1-1
        n2=n1-1
        
!c-------------------------------c
!        do k=k2,n2
!            app(j,k)=app(j,k)+ajp(j,k)
!            con(j,k)=con(j,k)-ajp(j,k)*vv1(j,k)
!            ajp(j,k)=0.0d+00
!        enddo
!c-------------------------------c
        do j=j2,m2
            k=k2
            app(j,k)=app(j,k)+akm(j,k)
            con(j,k)=con(j,k)-akm(j,k)*vv1(j,k)
            akm(j,k)=0.0d+00
            k=n2
            app(j,k)=app(j,k)+akp(j,k)
            con(j,k)=con(j,k)-akp(j,k)*vv1(j,k)
            akp(j,k)=0.0d+00
        enddo
!c-------------------------------c 
!        do k=kb,ke-1
!            do j=jb,je
!                    app(j,k)=1.0
!            enddo
!        enddo
!        do k=kb,ke-1
!            j=jb
!            app(j,k)=1.0
!            j=je
!            app(j,k)=1.0
!        enddo
!c-------------------------------c
        do j=jb,je
            k=ke
!            j=j2
!            app(j,k)=app(j,k)+ajm(j,k)
!            con(j,k)=con(j,k)-ajm(j,k)*vv1(j,k)
!            ajm(j,k)=0.0d+00
!            j=m2
            app(j,k)=app(j,k)+akm(j,k)
            con(j,k)=con(j,k)-akm(j,k)*vv1(j,k)
            akm(j,k)=0.0d+00
            k=kb-1
            app(j,k)=app(j,k)+akp(j,k)
            con(j,k)=con(j,k)-akp(j,k)*vv1(j,k)
            akp(j,k)=0.0d+00
        enddo
!c-------------------------------c
        return
        end

        
!c*******************************c
        subroutine bondw2(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,ww1,app,con,aim,ajm,akm,aip,ajp,akp)
!c*******************************c
        include 'table.prc'
        dimension app(nj,nk),con(nj,nk),aim(nj,nk),aip(nj,nk),ajm(nj,nk),ajp(nj,nk),akm(nj,nk),akp(nj,nk)
        dimension ww1(nj,nk)
!c*******************************c
        j2=j1+1
        k2=k1+1
        m2=m1-1
        n2=n1-1
!c-------------------------------c
        do k=k2,n2 
            j=j2
            app(j,k)=app(j,k)+ajm(j,k)
            con(j,k)=con(j,k)-ajm(j,k)*ww1(j,k)
            ajm(j,k)=0.0d+00
            j=m2
            app(j,k)=app(j,k)+ajp(j,k)
            con(j,k)=con(j,k)-ajp(j,k)*ww1(j,k)
            ajp(j,k)=0.0d+00
        enddo
!c-------------------------------c
!        do j=j2,m2
!            k=k2
!            app(j,k)=app(j,k)+akm(j,k)
!            con(j,k)=con(j,k)-akm(j,k)*ww1(j,k)
!            akm(j,k)=0.0d+00
!            k=n2
!            app(j,k)=app(j,k)+akp(j,k)
!            con(j,k)=con(j,k)-akp(j,k)*ww1(j,k)
!            akp(j,k)=0.0d+00
!        enddo
!c-------------------------------c
        do k=kb,ke
            j=je
            app(j,k)=app(j,k)+ajm(j,k)
            con(j,k)=con(j,k)-ajm(j,k)*ww1(j,k)
            ajm(j,k)=0.0d+00
            j=jb-1
            app(j,k)=app(j,k)+ajp(j,k)
            con(j,k)=con(j,k)-ajp(j,k)*ww1(j,k)
            ajp(j,k)=0.0d+00
        enddo
!c-------------------------------c
!        do j=jb,je-1
!            k=ke+1
!            app(j,k)=app(j,k)+akm(j,k)
!            con(j,k)=con(j,k)-akm(j,k)*ww1(j,k)
!            akm(j,k)=0.0d+00
!            k=kb-1
!            app(j,k)=app(j,k)+akp(j,k)
!            con(j,k)=con(j,k)-akp(j,k)*ww1(j,k)
!            akp(j,k)=0.0d+00
!        enddo
!c-------------------------------c
        return
        end

        subroutine bondp2(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,ww1,app,con,aim,ajm,akm,aip,ajp,akp)
!c*******************************c
        include 'table.prc'
        dimension app(nj,nk),con(nj,nk),aim(nj,nk),aip(nj,nk),ajm(nj,nk),ajp(nj,nk),akm(nj,nk),akp(nj,nk)
        dimension ww1(nj,nk)
!c*******************************c
        j2=j1+1
        k2=k1+1
        m2=m1-1
        n2=n1-1
!c-------------------------------c
        do k=k2,n2 
            j=j2
            app(j,k)=app(j,k)-ajm(j,k)
            con(j,k)=con(j,k)-ajm(j,k)*ww1(j,k)
            ajm(j,k)=0.0d+00
            j=m2
            app(j,k)=app(j,k)-ajp(j,k)
            con(j,k)=con(j,k)-ajp(j,k)*ww1(j,k)
            ajp(j,k)=0.0d+00
        enddo
!c-------------------------------c
!        do j=j2,m2
!            k=k2
!            app(j,k)=app(j,k)+akm(j,k)
!            con(j,k)=con(j,k)-akm(j,k)*ww1(j,k)
!            akm(j,k)=0.0d+00
!            k=n2
!            app(j,k)=app(j,k)+akp(j,k)
!            con(j,k)=con(j,k)-akp(j,k)*ww1(j,k)
!            akp(j,k)=0.0d+00
!        enddo
!c-------------------------------c
        do k=kb,ke
            j=je
            app(j,k)=app(j,k)-ajm(j,k)
            con(j,k)=con(j,k)-ajm(j,k)*ww1(j,k)
            ajm(j,k)=0.0d+00
            j=jb-1
            app(j,k)=app(j,k)-ajp(j,k)
            con(j,k)=con(j,k)-ajp(j,k)*ww1(j,k)
            ajp(j,k)=0.0d+00
        enddo
!c-------------------------------c
!        do j=jb,je-1
!            k=ke+1
!            app(j,k)=app(j,k)+akm(j,k)
!            con(j,k)=con(j,k)-akm(j,k)*ww1(j,k)
!            akm(j,k)=0.0d+00
!            k=kb-1
!            app(j,k)=app(j,k)+akp(j,k)
!            con(j,k)=con(j,k)-akp(j,k)*ww1(j,k)
!            akp(j,k)=0.0d+00
!        enddo
!c-------------------------------c
        return
        end

!c*******************************c
        subroutine bondu2(j1,k1,m1,n1,jb,je,kb,ke,ni,nj,nk,uu1,app,con,aim,ajm,akm,aip,ajp,akp)
!c*******************************c
        include 'table.prc'
        dimension app(nj,nk),con(nj,nk),aim(nj,nk),aip(nj,nk),ajm(nj,nk),ajp(nj,nk),akm(nj,nk),akp(nj,nk)
        dimension uu1(nj,nk)
!c*******************************c
        j2=j1+1
        k2=k1+1
        m2=m1-1
        n2=n1-1
!c-------------------------------c
        do k=k2,n2
            j=j2
            con(j,k)=con(j,k)+ajm(j,k)*uu1(j-1,k)
            ajm(j,k)=0.0d+00
            j=m2
            con(j,k)=con(j,k)+ajp(j,k)*uu1(j+1,k)
            ajp(j,k)=0.0d+00
        enddo
!c-------------------------------c
        do j=j2,m2
            k=k2
            con(j,k)=con(j,k)+akm(j,k)*uu1(j,k-1)
            akm(j,k)=0.0d+00
            k=n2
            con(j,k)=con(j,k)+akp(j,k)*uu1(j,k+1)
            akp(j,k)=0.0d+00
        enddo
!c-------------------------------c 
        do k=kb,ke-1
            j=je
            con(j,k)=con(j,k)+ajm(j,k)*uu1(j-1,k)
            ajm(j,k)=0.0d+00
            j=jb-1
            con(j,k)=con(j,k)+ajp(j,k)*uu1(j+1,k)
            ajp(j,k)=0.0d+00
        enddo
!c-------------------------------c
        do j=jb,je-1
            k=ke
            con(j,k)=con(j,k)+akm(j,k)*uu1(j,k-1)
            akm(j,k)=0.0d+00
            k=kb-1
            con(j,k)=con(j,k)+akp(j,k)*uu1(j,k+1)
            akp(j,k)=0.0d+00
        enddo
!c-------------------------------c
        return
        end
