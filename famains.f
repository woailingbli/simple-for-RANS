c===============================c
        include 'table.md1'
        include 'table.md2'
        include 'table.md3'
        include 'table.md4'
        include 'table.md5'
        include 'table.md6'
c===============================c
        program initia
        USE omp_lib
        use mdugrd;use mduvar;use mdumgv;use mdupnt

        call alcgrd
        call alcvar

        call gddata_re400
        call stdata_re400        
!        call gddata_re200
!        call stdata_re200

        call augrid
        call ibdata
!        call opd0uf
!        call opd0fm
!        call opd1uf
!        call opd1fm
!        call opd2uf
        
        do i=1,500
100     write(*,*) '           solving v-momentum equation           '
        call vlamin
        write(*,*) '           solving w-momentum equation           '
        call wlamin
        write(*,*) '       solving pressure-correction equation       '
        call plamin
        enddo
        
1000    write(*,*)'          solving u-momentum equation           '
        call ulamin

        end
