c----------------------------------------------------------------------c
        block data nsedat
c----------------------------------------------------------------------c
        include 'table.prc'
        include 'table.cre'
        include 'table.gd1'
        include 'table.gdm'
c-------the grid level-------c
        data levmgd/9/
c-------the grid sizes-------c
        data dt/2.0d-05/
        data jmd1,jmd2,kmd1,kmd2,cs1,cs2
     *  /48,458,48,458,2.0d+00,2.5d+00/
!        data jmd1,jmd2,kmd1,kmd2,cs1,cs2
!     *  /24,236,24,236,2.0d+00,2.5d+00/
        data iu0,ju0,ku0,iv0,jv0,kv0,iw0,jw0,kw0,ip0,jp0,kp0,it0,jt0,kt0
     *  /2,1,1,1,2,1,1,1,2,1,1,1,1,1,1/
        data imx/259/
        data jmx,jb1,je1/514,130,386/
        data kmx,kb1,ke1/514,130,386/
!        data jmx,jb1,je1/258,66,194/
!        data kmx,kb1,ke1/258,66,194/

        data j1m(0),m1m(0),k1m(0),n1m(0)/1,3,1,3/
        data j1m(1),m1m(1),k1m(1),n1m(1)/1,4,1,4/
        data j1m(2),m1m(2),k1m(2),n1m(2)/1,6,1,6/
        data j1m(3),m1m(3),k1m(3),n1m(3)/1,10,1,10/
        data j1m(4),m1m(4),k1m(4),n1m(4)/1,18,1,18/
        data j1m(5),m1m(5),k1m(5),n1m(5)/1,34,1,34/
        data j1m(6),m1m(6),k1m(6),n1m(6)/1,66,1,66/
        data j1m(7),m1m(7),k1m(7),n1m(7)/1,130,1,130/
        data j1m(8),m1m(8),k1m(8),n1m(8)/1,258,1,258/
        data j1m(9),m1m(9),k1m(9),n1m(9)/1,514,1,514/

        data jbm(0),jem(0),kbm(0),kem(0)/0,0,0,0/
        data jbm(1),jem(1),kbm(1),kem(1)/0,0,0,0/
        data jbm(2),jem(2),kbm(2),kem(2)/3,5,3,5/
        data jbm(3),jem(3),kbm(3),kem(3)/4,8,4,8/
        data jbm(4),jem(4),kbm(4),kem(4)/6,14,6,14/
        data jbm(5),jem(5),kbm(5),kem(5)/10,26,10,26/
        data jbm(6),jem(6),kbm(6),kem(6)/18,50,18,50/
        data jbm(7),jem(7),kbm(7),kem(7)/34,98,34,98/
        data jbm(8),jem(8),kbm(8),kem(8)/66,194,66,194/
        data jbm(9),jem(9),kbm(9),kem(9)/130,386,130,386/

        data j1mv(0),m1mv(0),k1mv(0),n1mv(0)/2,514,1,3/
        data j1mv(1),m1mv(1),k1mv(1),n1mv(1)/2,514,1,4/
        data j1mv(2),m1mv(2),k1mv(2),n1mv(2)/2,514,1,6/
        data j1mv(3),m1mv(3),k1mv(3),n1mv(3)/2,514,1,10/
        data j1mv(4),m1mv(4),k1mv(4),n1mv(4)/2,514,1,18/
        data j1mv(5),m1mv(5),k1mv(5),n1mv(5)/2,514,1,34/
        data j1mv(6),m1mv(6),k1mv(6),n1mv(6)/2,514,1,66/
        data j1mv(7),m1mv(7),k1mv(7),n1mv(7)/2,514,1,130/
        data j1mv(8),m1mv(8),k1mv(8),n1mv(8)/2,514,1,258/
        data j1mv(9),m1mv(9),k1mv(9),n1mv(9)/2,514,1,514/

        data jbmv(0),jemv(0),kbmv(0),kemv(0)/130,386,0,0/
        data jbmv(1),jemv(1),kbmv(1),kemv(1)/130,386,0,0/
        data jbmv(2),jemv(2),kbmv(2),kemv(2)/130,386,3,5/
        data jbmv(3),jemv(3),kbmv(3),kemv(3)/130,386,4,8/
        data jbmv(4),jemv(4),kbmv(4),kemv(4)/130,386,6,14/
        data jbmv(5),jemv(5),kbmv(5),kemv(5)/130,386,10,26/
        data jbmv(6),jemv(6),kbmv(6),kemv(6)/130,386,18,50/
        data jbmv(7),jemv(7),kbmv(7),kemv(7)/130,386,34,98/
        data jbmv(8),jemv(8),kbmv(8),kemv(8)/130,386,66,194/
        data jbmv(9),jemv(9),kbmv(9),kemv(9)/130,386,130,386/

        data j1mw(0),m1mw(0),k1mw(0),n1mw(0)/1,3,2,514/
        data j1mw(1),m1mw(1),k1mw(1),n1mw(1)/1,4,2,514/
        data j1mw(2),m1mw(2),k1mw(2),n1mw(2)/1,6,2,514/
        data j1mw(3),m1mw(3),k1mw(3),n1mw(3)/1,10,2,514/
        data j1mw(4),m1mw(4),k1mw(4),n1mw(4)/1,18,2,514/
        data j1mw(5),m1mw(5),k1mw(5),n1mw(5)/1,34,2,514/
        data j1mw(6),m1mw(6),k1mw(6),n1mw(6)/1,66,2,514/
        data j1mw(7),m1mw(7),k1mw(7),n1mw(7)/1,130,2,514/
        data j1mw(8),m1mw(8),k1mw(8),n1mw(8)/1,258,2,514/
        data j1mw(9),m1mw(9),k1mw(9),n1mw(9)/1,514,2,514/

        data jbmw(0),jemw(0),kbmw(0),kemw(0)/0,0,130,386/
        data jbmw(1),jemw(1),kbmw(1),kemw(1)/0,0,130,386/
        data jbmw(2),jemw(2),kbmw(2),kemw(2)/3,5,130,386/
        data jbmw(3),jemw(3),kbmw(3),kemw(3)/4,8,130,386/
        data jbmw(4),jemw(4),kbmw(4),kemw(4)/6,14,130,386/
        data jbmw(5),jemw(5),kbmw(5),kemw(5)/10,26,130,386/
        data jbmw(6),jemw(6),kbmw(6),kemw(6)/18,50,130,386/
        data jbmw(7),jemw(7),kbmw(7),kemw(7)/34,98,130,386/
        data jbmw(8),jemw(8),kbmw(8),kemw(8)/66,194,130,386/
        data jbmw(9),jemw(9),kbmw(9),kemw(9)/130,386,130,386/       
 
        data epsm1,epsm2,epsm3,epsm4
     *       /1.0d-15,1.0d-10,1.0d-08,5.0d-01/
        data epsp1(0),epsp2(0),epsp3(0),epsp4(0)
     *       /1.0d-15,1.0d-12,1.0d-01,5.0d-01/
        data epsp1(1),epsp2(1),epsp3(1),epsp4(1)
     *       /1.0d-15,1.0d-12,1.0d-01,5.0d-01/
        data epsp1(2),epsp2(2),epsp3(2),epsp4(2)
     *       /1.0d-15,1.0d-12,1.0d-01,5.0d-01/
        data epsp1(3),epsp2(3),epsp3(3),epsp4(3)
     *       /1.0d-15,1.0d-12,1.0d-01,5.0d-01/
        data epsp1(4),epsp2(4),epsp3(4),epsp4(4)
     *       /1.0d-15,1.0d-12,1.0d-01,5.0d-01/
        data epsp1(5),epsp2(5),epsp3(5),epsp4(5)
     *       /1.0d-15,1.0d-12,1.0d-01,5.0d-01/
        data epsp1(6),epsp2(6),epsp3(6),epsp4(6)
     *       /1.0d-15,1.0d-12,1.0d-01,5.0d-01/
        data epsp1(7),epsp2(7),epsp3(7),epsp4(7)
     *       /1.0d-15,1.0d-12,1.0d-01,5.0d-01/
        data epsp1(8),epsp2(8),epsp3(8),epsp4(8)
     *       /1.0d-15,1.0d-12,1.0d-06,5.0d-01/
        data epsp1(9),epsp2(9),epsp3(9),epsp4(9)
     *       /1.0d-13,1.0d-10,1.0d-10,5.0d-01/

        data alpha_v,alpha_w,alpha_p/9.0d-01,9.0d-01,1.0d-01/

        data x0,y0,z0/0.0d+00,0.0d+00,0.0d+00/
        data xl,yl,zl/1.2665318414472237d+01,4.0d+00,4.0d+00/
        data ymb,yme,zmb,zme/1.0d+00,3.0d+00,1.0d+00,3.0d+00/
c----------------------------c
        data modeui,modevi,modewi,modeti/2,2,2,2/
c-------criterion number-------c
        data strh,reno,prtl,eckt/1.0d+00,4.0d+02,0.7d+00,1.5d-10/
c-------bdc u,v,w,t-------c
        data ubo,vbo,wbo,tbo,ubi,vbi,wbi,tbi
     *  /0.0d+0,0.0d+0,0.0d+0,1.0d+0,0.0d+0,0.0d+0,0.0d+0,0.0d+0/
c-------turbulent control parameter-------c
        data csty,cnos,cwvn/5.0d-02,0.0d+00,1.0d+00/
c---------------------------------------------------------------------c
        end
