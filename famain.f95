!c===============================c
        include 'table.mdu'
!c===============================c
        program simple
!c-------------------------------c
        call ipd0fm
        call ipd1uf
        call ipdff0
        call augrid
	
        write(*,*)'*--------------------------------------------------*'
        write(*,*)'*                                                  *'
        write(*,*)'*    The program is to calculate 2D-rans           *'
        write(*,*)'*                                                  *'
        write(*,*)'*--------------------------------------------------*'
        write(*,*)'*                                                  *'
        call begins
        call rans_simple
!!        call update
!!        call initia
        write(*,*)'*                                                  *'
        write(*,*)'*--------------------------------------------------*'
        write(*,*)'*                                                  *'
        write(*,*)'*         output the divergence-free field         *'
!!c        call opd0fm
!!c        call opd1fm
!!c        call opdff0
!!c        call opdhf0
        call shuchu
        write(*,*)'*                                                  *'
        write(*,*)'*         The----End----Of---The---Program         *'
        write(*,*)'*                                                  *'
        write(*,*)'*--------------------------------------------------*'
        end

!c======================================================================c
        subroutine rans_simple 
!c======================================================================c
        use varalc
        include 'table.prc'
        write(*,*)'*--------------------------------------------------*'
        write(*,*)'*                                                  *'
        write(*,*)'*        solving  process for momentum Eqs.        *'

        write(*,*)'*                                                  *'
        call reno_stress

100 write(*,*)'*                                                  *'
            write(*,*)'*        solving  momentum  equation  of  v        *'
!            call ipd1fm
            call initiv
            call equatv

            write(*,*)'*                                                  *'
            write(*,*)'*        solving  momentum  equation  of  w        *'
            call initiw
            call equatw
            write(*,*)'*                                                  *'
            write(*,*)'*--------end of handling momentum equations--------*'
            write(*,*)'*--------------------------------------------------*'
            write(*,*)'*                                                  *'
            write(*,*)'*        solving pressure-correct  equation        *'
            call initip
            call equatp
            write(*,*)'*                                                  *'
            write(*,*)'*--------end of handling pressure-correct equations--------*'
            write(*,*)'*--------------------------------------------------*'
            call residual_all(j2,k2,m2,n2,jb,je,kb,ke,nj,nk,ressum,pre2,pre0,res4)
            write(*,*)'total  resudual:',resume
            if (resume < epsm21) then
            goto 666
            else
            call update
            goto 100
            endif
           write(*,*)'*                                                  *'
            write(*,*)'*--------------------------------------------------*'

666        return
        end

!!c======================================================================c
!        subroutine initia
!!c======================================================================c
!        call initiv
!        call initiw
!        call initip
!        return
!        end
