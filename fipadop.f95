!c================================c
        subroutine ipd0fm
!c================================c
        include 'table.prc'
        include 'table.v_cre'
        include 'table.v_gd1'
        include 'table.v_gd2'
        include 'table.v_gd3'
        include 'table.v_gd4'
        include 'table.v_gd5'
        include 'table.v_gd6'
        include 'table.v_gd7'
        include 'table.v_gd8'
        include 'table.v_gd9'
!c----------------------------------------------------------------------c
        strh=1.0;reno=400;dt=0.00005
        m11=514;jb1=130;je1=386;n11=514;kb1=130;ke1=386  
        m12=514;jb2=130;je2=386;n12=4;kb2=100;ke2=1!!!!!    
        m13=514;jb3=130;je3=386;n13=6;kb3=3;ke3=5     
        m14=514;jb4=130;je4=386;n14=10;kb4=4;ke4=8     
        m15=514;jb5=130;je5=386;n15=18;kb5=6;ke5=14     
        m16=514;jb6=130;je6=386;n16=34;kb6=10;ke6=26
        m17=514;jb7=130;je7=386;n17=66;kb7=18;ke7=50
        m18=514;jb8=130;je8=386;n18=130;kb8=34;ke8=98  
        m19=514;jb9=130;je9=386;n19=258;kb9=66;ke9=194  
        return
        end
        
!c================================c
        subroutine ipd1fm
!c================================c
        include 'table.prc'
        include 'table.w_cre'
        include 'table.w_gd1'
        include 'table.w_gd2'
        include 'table.w_gd3'
        include 'table.w_gd4'
        include 'table.w_gd5'
        include 'table.w_gd6'
        include 'table.w_gd7'
        include 'table.w_gd8'
!c----------------------------------------------------------------------c
        strh=1.0;reno=400;dt=0.00005
        m11=514;jb1=130;je1=386;n11=514;kb1=130;ke1=386  
        m12=514;jb2=130;je2=386;n12=6;kb2=3;ke2=5     
        m13=514;jb3=130;je3=386;n13=10;kb3=4;ke3=8     
        m14=514;jb4=130;je4=386;n14=18;kb4=6;ke4=14     
        m15=514;jb5=130;je5=386;n15=34;kb5=10;ke5=26     
        m16=514;jb6=130;je6=386;n16=66;kb6=18;ke6=50
        m17=514;jb7=130;je7=386;n17=130;kb7=34;ke7=98
        m18=514;jb8=130;je8=386;n18=258;kb8=66;ke8=194  
      
        return
        end
        
!c===============================c
        subroutine ipd1uf
!c===============================c
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
!c-------------------------------c
!c--------------------------c
        open(199,file='..\input\RANS-vw\griddata.txt')
        read(199,*)yv1
        close(199)
        open(190,file='..\input\RANS-vw\griddata.txt')
        read(190,*)zw1
        close(190)
        return
        end
        
!c===============================c
        subroutine ipdff0
!c===============================c
        use varalc 
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
!c-------------------------------c
        include 'table.alc'
!c-------------------------------c
        open(191,file='..\input\RANS-vw\v258to514ac.dat',form='formatted')
        do j=1,514
	    do k=1,514
	        read(191,*)vi1(j,k)
            enddo
        enddo
        close(191)
        
        open(192,file='..\input\RANS-vw\w258to514ac.dat',form='formatted')
        do j=1,514
	    do k=1,514
	        read(192,*)wi1(j,k)
            enddo
        enddo
        close(192)
        
        open(193,file='..\input\RANS-vw\p_te.dat',form='formatted')
        do k=1,514
	    do j=1,514
	        read(193,*)pi1(j,k)
            enddo
        enddo
        close(193)        


        return
        end
        
!c*********************************c
!c*********************************c
!c*********************************c
        subroutine shuchu
        use varalc
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
        include 'table.lei'
	

	open(91,file='..\output\RANS-vw\shuchu2_test.dat',form='formatted')
	write(91,*) 'TITLE    ="Plot3D DataSet"'
        write(91,*) 'VARIABLES = "y" "z" "v" "w" "pre"'
        write(91,*) 'DATASETAUXDATA Common.SpeedOfSound="1.0"'
        write(91,*) 'DATASETAUXDATA Common.VectorVarsAreVelocity="FALSE"'
        write(91,*) 'ZONE T="Zone-original grid"'
        write(91,*) 'STRANDID=0, SOLUTIONTIME=0'
        write(91,*) 'I=514, J=514, K=1, ZONETYPE=Ordered'
        write(91,*) 'DATAPACKING=POINT'
        write(91,*) 'DT=(SINGLE SINGLE SINGLE SINGLE SINGLE)'
	do k=1,514
	    do j=1,514
                write(91,*) ym1(j),zm1(k),vv2(j,k),ww2(j,k),pre2(j,k)
            enddo
        enddo
end

