      subroutine umat(stress,statev,ddsdde,sse,spd,scd,
     1 rpl,ddsddt,drplde,drpldt,stran,dstran,
     2 time,dtime,temp,dtemp,predef,dpred,materl,ndi,nshr,ntens,
     3 nstatv,props,nprops,coords,drot,pnewdt,celent,
     4 dfgrd0,dfgrd1,noel,npt,kslay,kspt,kstep,kinc)

      include 'aba_param.inc'

      character*80 materl
      double precision a_kapa,a_kapa0,a_kxl,x,a_l,dkapa0,f,f1,beta
	  integer kewton,newton
	  
      dimension a(17),statev(nstatv),props(nprops),
     1 stran(ntens),dstran(ntens),dstress(ntens),
     1 stress(ntens),astress(ntens),s_dev(ntens),d_stres(ntens),
     1 d_eplas(ntens),astrain(ntens),
     1 eplas(ntens),e_elas(ntens),eplas0(ntens),eplasStari(ntens),
     1 ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     1 time(2),predef(1),dpred(1),
     1 coords(ndi),drot(ndi,ndi),e_elas_n1(ntens),
     1 dfgrd0(ndi,ndi),dfgrd1(ndi,ndi),
     1 a_mu(ntens),kroneker(ntens),
     1 Replas(ntens),kroneker2(ntens,ntens),deplas(ntens)

      parameter (one=1.0d0,two=2.0d0,three=3.0d0,six=6.0d0,zero=0.0d0)
      data toler,temp0,coef,yield0/1.d-6,273.,0.0d0,25.0d0/
	  ! inicijalizacija
	  a_kapa = 0
	  a_kapa0=0
	  kewton = 1
	  newton = 50000
	  a_kxl =0
	  x = 0
	  a_l = 0
	  dkapa0 = 0
	  f = 0
        f1 = 0
	  beta = 0

	  
	  
	  
      tol1 = 1.d-7
      tol2 = 1.d-7
      tolk = 1.d-4
      !C -----------------------------------------------------------
      !c
      !c     paneerselvam phd   (page 165 tens constants)
      !c
      a3 = -1.81e5
      a4 = -0.906e5
      a6 = 0.
      a7 = 0.
      a8 = 0.
      a9 = 0.
      x = 1.03e-3
      a_l = 2
      alfa = 14.2  !17.1
      h = 50
      beta = 8.23e4     !1.34e4
      a_m  = one
      gama = -5.e-3
      alpha_t  = 4.e-5
      !stari
      !a1 = 70.3e2
      !a2 = 0.036e5 !7.5e4
      !a5 = -0.6046e5

      !novi
      a1 =  0.036e5
      a2 = -0.6046e5
      a5 = 70.3e2
      !c -----------------------------------------------------------
       !c -----------------------------------------------------------
      a(1)=a1
      a(2)=a2
      a(3)=a3
      a(4)=a4
      a(5)=a5
      a(6)=a6
      a(7)=a7
      a(8)=a8
      a(9)=a9
      a(10)=x
      a(11)=a_l
      a(12)=alfa
      a(13)=h
      a(14)=beta
      a(15)=m
      a(16)=gama
      a(17)=alpha_t

      tol2 = 1.0d-6

      do k2=1, ndi
         kroneker(k2)     = one
         kroneker(k2+ndi) = zero
      end do 
      !  KONSTANTE
      !     NDI: Broj direkthih komponenti napona u datom trenutku!     NDI=3
      !     NSHR: Broj smicajnih komponenti napona u datom trenutku!      NSHR=3
      !     NTENS: velicina niza napona ili deformacija (NDI + NSHR)!     NTENS=6
      !     NOEL: Broj elementa
      !     NPT: Broj integracione ta?ke      
      !  UILAZNE VELICINE U UMAT
      !     STRAN(NTENS): Niz koji sadrzi ukupne deformacije na pocetku inkrenenta
      !     DSTRAN(NTENS): Niz inkremenata deformacija
      !  IZLAZNE VELICINE    
      !     STRESS(NTENS): TENZOR NAPONA OVO JE I ULAZ I IZLAZ
      !     NA POCETKU DOBIJEMO TENZOR NAPONA NA POCETKU INKREMENTA
      !     PA U UMATU TREBA DA IZRACUNAMO TENZOR NAPONA NA KRAJU INKREMENTA    
      !     DDSDDE(NTENS,NTENS)   PARCIJALNI IZVOD INKREMENTA NAPONA PO INKREMENTU DEFORMAICJE
      !     OVO JE CISTO IZLAZNA VELICINA KOJA SE RACUNA U UMATU
      !     STRESS JE MNOGO VAZNIJI OD DDSDDE
      !     DDSDDE SE KORISTI U ABAQUSU ZA PROVERU KONVERGENCIJE       

      !  OSTALE ULAZNE ILI IZLAZNE VELICINE NE KORISTIMO
    
      !  UNUTRASNJE PROMENLJIVE
      !     STATEV NIZ UNUTRASNJIH PROMENLJIVIH
      !     NSTATV BROJ UNUTRASNJIH PROMENLJIVIH
      !     OVO VISE NE KORISTIMO SADA SVE CUVAMO U COMMON STRUKTURAMA          

        !c********************************************************************
        !c
        !c               predictor corector algorithm
        !c
        !c********************************************************************     
        !cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        !c***************      1) predictor phase      ***********************  
        !cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      ! ucitava iz prethodnog koraka
        do k1=1,ntens
        stress(k1)=0
        eplas0(k1) = statev(k1)
        eplas(k1) =eplas0(k1)
      !write(6,*) 'eplas0',k1,'=',eplas0(k1)
        e_elas_n1(k1) = stran(k1)+dstran(k1)-eplas0(k1)
        enddo
        a_kapa0  = statev(17)  ! 6.20  ! ucitava iz prethodnog koraka
        if (a_kapa0.lt.tolk) a_kapa0=tolk
      ! ucitava iz prethodnog koraka

      !invarijante
        e_i1n1 = e_elas_n1(1)+e_elas_n1(2)+e_elas_n1(3)
        
        e_i2n1 = e_elas_n1(1)*e_elas_n1(2)+e_elas_n1(2)*e_elas_n1(3)+
     1   e_elas_n1(3)*e_elas_n1(1)-e_elas_n1(4)*e_elas_n1(4)-
     1   e_elas_n1(5)*e_elas_n1(5)-e_elas_n1(6)*e_elas_n1(6)
    
        e_i3n1 = e_elas_n1(1)*e_elas_n1(2)*e_elas_n1(3)-
     1   e_elas_n1(1)*e_elas_n1(6)*e_elas_n1(6)-
     1   e_elas_n1(2)*e_elas_n1(5)*e_elas_n1(5)-
     1   e_elas_n1(3)*e_elas_n1(4)*e_elas_n1(4)+
     1   2*e_elas_n1(4)*e_elas_n1(5)*e_elas_n1(6)
         !invarijante                                          
      
         astress(1)= (2*a5*e_i1n1+3*a3*e_i1n1*e_i1n1+a4*e_i2n1)+
     1     (a1+a4*e_i1n1)*e_elas_n1(1)+a2*(e_elas_n1(1)*e_elas_n1(1)
     1    +e_elas_n1(4)*e_elas_n1(4)+e_elas_n1(5)*e_elas_n1(5))
    
         astress(2)= (2*a5*e_i1n1+3*a3*e_i1n1*e_i1n1+a4*e_i2n1)+
     1     (a1+a4*e_i1n1)*e_elas_n1(2)+a2*(e_elas_n1(4)*e_elas_n1(4)
     1    +e_elas_n1(2)*e_elas_n1(2)+e_elas_n1(6)*e_elas_n1(6))
    
         astress(3)= (2*a5*e_i1n1+3*a3*e_i1n1*e_i1n1+a4*e_i2n1)+
     1     (a1+a4*e_i1n1)*e_elas_n1(3)+a2*(e_elas_n1(5)*e_elas_n1(5)
     1    +e_elas_n1(6)*e_elas_n1(6)+e_elas_n1(3)*e_elas_n1(3))
    
    
         astress(4)= (a1+a4*e_i1n1)*e_elas_n1(4)
     1     +a2*(e_elas_n1(1)*e_elas_n1(4)
     1    +e_elas_n1(4)*e_elas_n1(2)+e_elas_n1(5)*e_elas_n1(6))
    
         astress(5)= (a1+a4*e_i1n1)*e_elas_n1(5)
     1     +a2*(e_elas_n1(1)*e_elas_n1(5)
     1    +e_elas_n1(4)*e_elas_n1(6)+e_elas_n1(5)*e_elas_n1(3))
    
         astress(6)= (a1+a4*e_i1n1)*e_elas_n1(6)
     1     +a2*(e_elas_n1(4)*e_elas_n1(5)
     1    +e_elas_n1(2)*e_elas_n1(6)+e_elas_n1(6)*e_elas_n1(3))
    

        call noviddsdde(a,ddsdde,ntens,e_elas_n1)
        ! call    hyperconstitutive (a,ddsdde,ntens,e_elas_n1)
      do k1=1,ntens              
          stress(k1) =astress(k1)
      enddo  
         dkapa = 0
    
       !c------------------ compute  loading surface f-----------------------
        call loadingf(f1,stress,a,ntens,ndi,s_dev,a_j2,a_mu) ! 6.21
        
        f = abs(f1) - h*a_kapa0 
        a_kapa = a_kapa0
        !goto 52
        if (f.le.zero) then 
        write(6,*) 'elasticnost'
        goto 52   
        endif
        !c------------------  end of elastic predictor ----------------------
        !cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        !c***************       1) corector phase      ***********************  
        !cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
              
        
        a_kxl = x+a_kapa0**a_l    

        do kewton = 1,newton
            
            call loadingf(f1,stress,a,ntens,ndi,s_dev,a_j2,a_mu)
            f  = abs(f1) - h*a_kapa            
            
            
            dkapa0  =  (((f/beta)**1)/a_kxl)*dtime
            a_kapa = a_kapa0 + dkapa0
            a_kxl = x+a_kapa**a_l
            dkapa  =  (((f/beta)**1)/a_kxl)*dtime
            
            do k1 = 1,ntens
            
            eplasStari(k1)   = eplas(k1)
            d_eplas(k1)   =  (dkapa)*a_mu(k1)
            eplas(k1)   = eplas0(k1) + d_eplas(k1)
            Replas(k1) = eplas(k1) - eplasStari(k1)
            enddo
                deplas_int = (Replas(1)**2+Replas(2)**2+Replas(3)**2+
     1    two*Replas(4)**2+two*Replas(5)**2+two*Replas(6)**2)**0.5
     
            f  = abs(f1) - h*a_kapa
            a_kxl = x+a_kapa**a_l
       !write(6,*) 'a_kxl3',kewton,'=',a_kxl
            skonvergencija = abs(a_kapa - a_kapa0)
            a_kapa0 = a_kapa
       if ((skonvergencija.lt.tolk).and.(deplas_int.lt.tol2)) then
            write(6,*) 'konvergirao u',kewton,'iteraciji'
            write(6,*) 'S =',skonvergencija,'||R||=',deplas_int
            goto 33
       !else
          !write(6,*) 'sk=',skonvergencija,dkapa0,kewton 
       endif
            
        enddo
           
            
            
        !dkapa = a_kapa - a_kapa0
 33        continue
         do k1 = 1,ntens
        !eplas(k1)   =  dkapa*a_mu(k1) !6.14
        !eplas(k1)   = eplas(k1) + (a_kapa-a_kapa0)*a_mu(k1)*dtime
        !eplas(k1)   = eplas(k1) + dkapa*a_mu(k1)*dtime
        !eplas(k1)   = eplas(k1) + dkapa*a_mu(k1)*dtime
        e_elas_n1(k1) = stran(k1)+dstran(k1)-eplas(k1)
        statev(17)=a_kapa  ! 6.20  ! ucitava iz prethodnog koraka
        enddo
        

        !invarijante    
        e_i1n1 = e_elas_n1(1)+e_elas_n1(2)+e_elas_n1(3)
        
        e_i2n1 = e_elas_n1(1)*e_elas_n1(2)+e_elas_n1(2)*e_elas_n1(3)+
     1   e_elas_n1(3)*e_elas_n1(1)-e_elas_n1(4)*e_elas_n1(4)-
     1   e_elas_n1(5)*e_elas_n1(5)-e_elas_n1(6)*e_elas_n1(6)
    
        e_i3n1 = e_elas_n1(1)*e_elas_n1(2)*e_elas_n1(3)-
     1   e_elas_n1(1)*e_elas_n1(6)*e_elas_n1(6)-
     1   e_elas_n1(2)*e_elas_n1(5)*e_elas_n1(5)-
     1   e_elas_n1(3)*e_elas_n1(4)*e_elas_n1(4)+
     1   2*e_elas_n1(4)*e_elas_n1(5)*e_elas_n1(6)
         !invarijante                                          
      
         astress(1)= (2*a5*e_i1n1+3*a3*e_i1n1*e_i1n1+a4*e_i2n1)+
     1     (a1+a4*e_i1n1)*e_elas_n1(1)+a2*(e_elas_n1(1)*e_elas_n1(1)
     1    +e_elas_n1(4)*e_elas_n1(4)+e_elas_n1(5)*e_elas_n1(5))
    
         astress(2)= (2*a5*e_i1n1+3*a3*e_i1n1*e_i1n1+a4*e_i2n1)+
     1     (a1+a4*e_i1n1)*e_elas_n1(2)+a2*(e_elas_n1(4)*e_elas_n1(4)
     1    +e_elas_n1(2)*e_elas_n1(2)+e_elas_n1(6)*e_elas_n1(6))
    
         astress(3)= (2*a5*e_i1n1+3*a3*e_i1n1*e_i1n1+a4*e_i2n1)+
     1     (a1+a4*e_i1n1)*e_elas_n1(3)+a2*(e_elas_n1(5)*e_elas_n1(5)
     1    +e_elas_n1(6)*e_elas_n1(6)+e_elas_n1(3)*e_elas_n1(3))
    
    
         astress(4)= (a1+a4*e_i1n1)*e_elas_n1(4)
     1     +a2*(e_elas_n1(1)*e_elas_n1(4)
     1    +e_elas_n1(4)*e_elas_n1(2)+e_elas_n1(5)*e_elas_n1(6))
    
         astress(5)= (a1+a4*e_i1n1)*e_elas_n1(5)
     1     +a2*(e_elas_n1(1)*e_elas_n1(5)
     1    +e_elas_n1(4)*e_elas_n1(6)+e_elas_n1(5)*e_elas_n1(3))
    
         astress(6)= (a1+a4*e_i1n1)*e_elas_n1(6)
     1     +a2*(e_elas_n1(4)*e_elas_n1(5)
     1    +e_elas_n1(2)*e_elas_n1(6)+e_elas_n1(6)*e_elas_n1(3))
    
       !do k1=1,ntens            
       !stress(k1) =astress(k1)
       !enddo
        goto 52    
       ! call noviddsdde(a,ddsdde,ntens,e_elas_n1)

 
      !c       corrector phase -  kraj 
 52      continue
        
        do k1=1,ntens
         statev(k1)  =  eplas(k1)
         !eplas(k1)   = eplas(k1) + dkapa*a_mu(k1)*dtime
         enddo 
         statev(17)=a_kapa ! 6.20  ! ucitava iz prethodnog koraka
        
                if(NPT.eq.1) then
                  !write(6,*)stran(1)+stran(2)+stran(3)
                  !write(6,*)a_j2
                  !write(6,*)w
                  !write(6,*)time, w
                  !write(6,*)a_kapa, a_kapa0, dkapau
                  !write(6,*)stran(2)
                  !write(6,*)stress(2)
                  !write(6,*)stran(2),dstran(2),time
                  !write(6,*)'U', stran(2)+dstran(2),'P',eplas(2)
                  !write(6,*)stress(2)
                  !write(6,*) 'e(2)=',stran(2)
                  !write(6,*)'eP(2)=',eplas(2)
                  !write(6,*)eplas(2)
                  !write(6,*)'eE(2)=',e_elas_n(2)
                endif
        
         
      return    
      end

!C
!c-----------------------------------------  umat   kraj
!c-----------------------------------------  umat   kraj
      
      
      subroutine hyperconstitutive(a,ddsdde,ntens,e_elas)

      include 'aba_param.inc'
        !     brojne konstante   
        !     real zero,one,two,three,four,six
        !     brojne konstante 
        !          
        !     ulazne promenljive   
        !     real a(17), a1,a2,a3,a4,a5,a6,a7,a8,a9      
        !     integer ntens
        !     real  e_elas(ntens)
        !     ulazne promenljive  
        !
        !     unutrasnje promenljive
        !     integer k1,k2
        !     real ee_i,ee_ii,ee_iii
        !     unutrasnje promenljive 
        !
        !     izlazne promenljive
        !     real ddsdde(ntens,ntens)
        !     izlazne promenljive

      dimension ddsdde(ntens,ntens), a(17), 
     1 e_elas(ntens)
     
      parameter (one=1.0d0,two=2.0d0,three=3.0d0,six=6.0d0)
      a1=a(1)
      a2=a(2)
      a3=a(3)
      a4=a(4)
      a5=a(5)
      a6=a(6)
      a7=a(7)
      a8=a(8)
      a9=a(9)
      do k1=1,ntens
      do k2=2,ntens
      ddsdde(k1,k2)=0
      enddo
      enddo

       ddsdde(1,1)= 2*a1 + a2 + 2*e_elas(1)*a4 + 2*e_elas(1)*a5 + 
     1 3*a3*(2*e_elas(1) + 2*e_elas(2) + 2*e_elas(3)) + a4*(e_elas(1) 
     1 + e_elas(2) + e_elas(3))
 
       ddsdde(1,2)= 2*a1 + e_elas(1)*a4 + e_elas(2)*a4 + 
     1  3*a3*(2*e_elas(1) + 2*e_elas(2) + 2*e_elas(3))
      
       ddsdde(1,3)= 2*a1 + e_elas(1)*a4 + e_elas(3)*a4 + 
     1  3*a3*(2*e_elas(1) + 2*e_elas(2) + 2*e_elas(3))
     
       ddsdde(1,4)= 2*e_elas(4)*a4 + 2*e_elas(4)*a5
       ddsdde(1,5)= 2*e_elas(5)*a4 + 2*e_elas(5)*a5
       ddsdde(1,6)= 2*e_elas(6)*a4

       ddsdde(2,1)= 2*a1 + e_elas(1)*a4 + e_elas(2)*a4 + 
     1  3*a3*(2*e_elas(1) + 2*e_elas(2) + 2*e_elas(3))
     
       ddsdde(2,2)= 2*a1 + a2 + 2*e_elas(2)*a4 + 2*e_elas(2)*a5 + 
     1  3*a3*(2*e_elas(1) + 2*e_elas(2) + 2*e_elas(3)) + a4*(e_elas(1) 
     1  + e_elas(2) + e_elas(3))
     
       ddsdde(2,3)= 2*a1 + e_elas(2)*a4 + e_elas(3)*a4 + 
     1  3*a3*(2*e_elas(1) + 2*e_elas(2) + 2*e_elas(3))
     
       ddsdde(2,4)= 2*e_elas(4)*a4 + 2*e_elas(4)*a5
       ddsdde(2,5)= 2*e_elas(5)*a4
       ddsdde(2,6)= 2*e_elas(6)*a4 + 2*e_elas(6)*a5
             
       ddsdde(3,1)= 2*a1 + e_elas(1)*a4 + e_elas(3)*a4 + 
     1  3*a3*(2*e_elas(1) + 2*e_elas(2) + 2*e_elas(3))
     
       ddsdde(3,2)= 2*a1 + e_elas(2)*a4 + e_elas(3)*a4 + 
     1  3*a3*(2*e_elas(1) + 2*e_elas(2) + 2*e_elas(3))
     
       ddsdde(3,3)= 2*a1 + a2 + 2*e_elas(3)*a4 + 2*e_elas(3)*a5 + 
     1  3*a3*(2*e_elas(1) + 2*e_elas(2) + 2*e_elas(3)) + a4*(e_elas(1) 
     1 + e_elas(2) + e_elas(3))
     
       ddsdde(3,4)= 2*e_elas(4)*a4
       ddsdde(3,5)= 2*e_elas(5)*a4 + 2*e_elas(5)*a5
       ddsdde(3,6)= 2*e_elas(6)*a4 + 2*e_elas(6)*a5
       
       ddsdde(4,1)= 2*e_elas(4)*a4 + 2*e_elas(4)*a5
       ddsdde(4,2)= 2*e_elas(4)*a4 + 2*e_elas(4)*a5
       ddsdde(4,3)= 2*e_elas(4)*a4
       ddsdde(4,4)= 2*a2 + a5*(2*e_elas(1) + 2*e_elas(2)) + 
     1  2*a4*(e_elas(1) + e_elas(2) + e_elas(3))
       ddsdde(4,5)= 2*e_elas(6)*a5
       ddsdde(4,6)= 2*e_elas(5)*a5
       
       ddsdde(5,1)= 2*e_elas(5)*a4 + 2*e_elas(5)*a5
       ddsdde(5,2)= 2*e_elas(5)*a4
       ddsdde(5,3)= 2*e_elas(5)*a4 + 2*e_elas(5)*a5
       ddsdde(5,4)= 2*e_elas(6)*a5
       ddsdde(5,5)= 2*a2 + a5*(2*e_elas(1) + 2*e_elas(3)) + 
     1  2*a4*(e_elas(1) + e_elas(2) + e_elas(3))
       ddsdde(5,6)= 2*e_elas(4)*a5
       
       ddsdde(6,1)= 2*e_elas(6)*a4
       ddsdde(6,2)= 2*e_elas(6)*a4 + 2*e_elas(6)*a5
       ddsdde(6,3)= 2*e_elas(6)*a4 + 2*e_elas(6)*a5
       ddsdde(6,4)= 2*e_elas(5)*a5
       ddsdde(6,5)= 2*e_elas(4)*a5
       ddsdde(6,6)= 2*a2 + a5*(2*e_elas(2) + 2*e_elas(3)) + 
     1  2*a4*(e_elas(1) + e_elas(2) + e_elas(3))    
       return
      end
           
      ! c----------------------subroutine hyperconstitutive     
      ! c----------------------------------------------------                            
        
      subroutine loadingf(f1,stress,a,ntens,ndi,s_dev,a_j2,a_mu)
     
      include 'aba_param.inc' 
      !c-
      double precision f1
      !     brojne konstante   
      !     real zero,one,two,three
      !     brojne konstante 
      !          
      !     ulazne promenljive   
      !     real a(17), alfa,gama       
      !     integer ntens,ndi
      !     real  stress(ntens)
      !     ulazne promenljive  
      !
      !     unutrasnje promenljive
      !     integer k1,k2
      !     real kroneker(ntens)
      !     unutrasnje promenljive 
      !
      !     izlazne promenljive
      !     real  f1
      !     real a_j2
      !     real a_mu(ntens)
      !     real s_dev(ntens)
      !     izlazne promenljive     

        dimension stress(ntens),s_dev(ntens),a(17),a_mu(ntens),
     1   kroneker(ntens),st(ntens)
      
            parameter (one=1.0d0,two=2.0d0,three=3.0d0,six=6.0d0) 
            
            f1 = 0
            
            alfa =  a(12)
            gama = a(16)
                s_i1=0
                s_i2=0
                s_i3=0
             do k2=1, ndi
               kroneker(k2)     = one
               kroneker(k2+ndi) = zero
               st(k2) = stress(k2) ! unutrasnja inverzija znaka  za racunanje funkcije tecenja
               st(k2+ndi) = stress(k2+ndi)
             end do 
            
       !c     s_napon invarijante  i funkcija f 


       !  s_i1 = stress(1)+stress(2)+stress(3)
      
            s_i1 = st(1)+st(2)+st(3)
      
       !  s_i2 = (stress(1)**2+stress(2)**2+stress(3)**2)/two+
       ! 1       stress(4)**2+stress(5)**2+stress(6)**2
       ! 1       -(s_i1**2)/two
           
           
             s_i2 = st(1)*st(2)+st(2)*st(3) 
     1         +st(1)*st(3)-st(4)**2-st(5)**2-st(6)**2
             !      s_i3 = stress(1)*stress(2)*stress(3)+
             !     1 two*stress(4)*stress(5)*stress(6)-stress(1)*stress(6)**2- 
             !     1 stress(2)*stress(5)**2-stress(3)*stress(4)**2
          s_i3 = st(1)*st(2)*st(3)+
     1 two*st(4)*st(5)*st(6)-st(1)*st(6)**2- 
     1 st(2)*st(5)**2-st(3)*st(4)**2

       do  k1=1,ndi
           !s_dev(k1)     = stress(k1) - s_i1/three
           !s_dev(k1+ndi) = stress(k1+ndi)
           s_dev(k1)     = st(k1) - s_i1/three
           s_dev(k1+ndi) = st(k1+ndi)
       enddo
  
      a_j2 = ( s_dev(1)**2+s_dev(2)**2+s_dev(3)**2)/two+
     1          s_dev(4)**2+s_dev(5)**2+ s_dev(6)**2 
     
      f1 = s_i1*s_i2 + alfa*s_i3 
      !if(f1.ne.f1)then
      !    write (6,*) "f1 je nan=",f1,s_i1,s_i2,alfa,s_i3
      !endif
         !                write(6,77), s_i1,s_i2,s_i3,f1,a_j2
         ! 77     format ('s_i1=',E12.4,' s_i2=',E12.4,' s_i3=',E12.4,' f1=',E12.4,' a_j2=',E12.4)
       
      do k1=1,ntens
           if (a_j2.gt.0) then
           a_mu(k1) = gama*kroneker(k1)+( 0.5*s_dev(k1)/(a_j2**0.5) ) !5.21
           endif
      enddo   
      return
      end
                   !c        
            !c---------------------------------subroutine loadingf      
               !c---------------------------------------------------- 

      subroutine noviddsdde(a,ddsdde,ntens,e_elas)

      include 'aba_param.inc'
        !     brojne konstante   
        !     real zero,one,two,three,four,six
        !     brojne konstante 
        !          
        !     ulazne promenljive   
        !     real a(17), a1,a2,a3,a4,a5,a6,a7,a8,a9      
        !     integer ntens
        !     real  e_elas(ntens)
        !     ulazne promenljive  

        !     izlazne promenljive
        !     real ddsdde(ntens,ntens)
        !     izlazne promenljive

      dimension ddsdde(ntens,ntens), a(17), 
     1 e_elas(ntens)
     
      parameter (one=1.0d0,two=2.0d0,three=3.0d0,six=6.0d0)
      a1=a(1)
      a2=a(2)
      a3=a(3)
      a4=a(4)
      a5=a(5)
      a6=a(6)
      a7=a(7)
      a8=a(8)
      a9=a(9)
      do k1=1,ntens
      do k2=2,ntens
      ddsdde(k1,k2)=0
      enddo
      enddo

       ddsdde(1,1)= 2*a1 + 4*a5 + 4*e_elas(1)*a2 + 
     1   4*e_elas(1)*a4 + 6*a3*(2*e_elas(1) + 
     1 2*e_elas(2) + 2*e_elas(3)) + 2*a4*(e_elas(1) + 
     1 e_elas(2) + e_elas(3))
 
 
 
 
       ddsdde(1,2)= 4*a5 + 2*e_elas(1)*a4 + 2*e_elas(2)*a4 + 
     1   6*a3*(2*e_elas(1) + 2*e_elas(2) + 2*e_elas(3))
      
       ddsdde(1,3)= 4*a5 + 2*e_elas(1)*a4 + 2*e_elas(3)*a4 + 
     1   6*a3*(2*e_elas(1) + 2*e_elas(2) + 2*e_elas(3))
     
       ddsdde(1,4)= e_elas(4)*a2 + e_elas(4)*a4 + e_elas(4)*a2 + 
     1   e_elas(4)*a4
       ddsdde(1,5)= e_elas(5)*a2 + e_elas(5)*a4 + e_elas(5)*a2 + 
     1   e_elas(5)*a4
       ddsdde(1,6)= E23*a4 + e_elas(6)*a4

       ddsdde(2,1)= 4*a5 + 2*e_elas(1)*a4 + 2*e_elas(2)*a4 + 
     1   6*a3*(2*e_elas(1) + 2*e_elas(2) + 2*e_elas(3))
     
       ddsdde(2,2)= 2*a1 + 4*a5 + 4*e_elas(2)*a2 + 4*e_elas(2)*a4 + 
     1   6*a3*(2*e_elas(1) + 2*e_elas(2) + 2*e_elas(3)) + 
     1   2*a4*(e_elas(1) + e_elas(2) + e_elas(3))
     
       ddsdde(2,3)= 4*a5 + 2*e_elas(2)*a4 + 2*e_elas(3)*a4 + 
     1   6*a3*(2*e_elas(1) + 2*e_elas(2) + 2*e_elas(3))
     
       ddsdde(2,4)= e_elas(4)*a2 + e_elas(4)*a4 + e_elas(4)*a2 + 
     1   e_elas(4)*a4
       ddsdde(2,5)= e_elas(5)*a4 + e_elas(5)*a4
       ddsdde(2,6)= E23*a2 + E23*a4 + e_elas(6)*a2 + e_elas(6)*a4
             
       ddsdde(3,1)= 4*a5 + 2*e_elas(1)*a4 + 2*e_elas(3)*a4 + 
     1   6*a3*(2*e_elas(1) + 2*e_elas(2) + 2*e_elas(3))
     
       ddsdde(3,2)= 4*a5 + 2*e_elas(2)*a4 + 2*e_elas(3)*a4 + 
     1   6*a3*(2*e_elas(1) + 2*e_elas(2) + 2*e_elas(3))
     
       ddsdde(3,3)= 2*a1 + 4*a5 + 4*e_elas(3)*a2 + 4*e_elas(3)*a4 + 
     1   6*a3*(2*e_elas(1) + 2*e_elas(2) + 2*e_elas(3)) + 
     1   2*a4*(e_elas(1) + e_elas(2) + e_elas(3))
     
       ddsdde(3,4)= e_elas(4)*a4 + e_elas(4)*a4
       ddsdde(3,5)= e_elas(5)*a2 + e_elas(5)*a4 + e_elas(5)*a2 + 
     1   e_elas(5)*a4
       ddsdde(3,6)= E23*a2 + E23*a4 + e_elas(6)*a2 + e_elas(6)*a4
       
       ddsdde(4,1)= 2*e_elas(4)*a2 + 2*e_elas(4)*a4
       ddsdde(4,2)= 2*e_elas(4)*a2 + 2*e_elas(4)*a4
       ddsdde(4,3)= 2*e_elas(4)*a4
       ddsdde(4,4)= a1 + a2*(e_elas(1) + e_elas(2)) + a4*(e_elas(1) + 
     1   e_elas(2) + e_elas(3))
       ddsdde(4,5)= e_elas(6)*a2
       ddsdde(4,6)= e_elas(5)*a2
       
       ddsdde(5,1)= 2*e_elas(5)*a2 + 2*e_elas(5)*a4
       ddsdde(5,2)= 2*e_elas(5)*a4
       ddsdde(5,3)= 2*e_elas(5)*a2 + 2*e_elas(5)*a4
       ddsdde(5,4)= e_elas(6)*a2
       ddsdde(5,5)= a1 + a2*(e_elas(1) + e_elas(3)) + a4*(e_elas(1) + 
     1   e_elas(2) + e_elas(3))
       ddsdde(5,6)= e_elas(4)*a2
       
       ddsdde(6,1)= 2*e_elas(6)*a4
       ddsdde(6,2)= 2*e_elas(6)*a2 + 2*e_elas(6)*a4
       ddsdde(6,3)= 2*e_elas(6)*a2 + 2*e_elas(6)*a4
       ddsdde(6,4)= e_elas(5)*a2
       ddsdde(6,5)= e_elas(4)*a2
       ddsdde(6,6)= a1 + a2*(e_elas(2) + e_elas(3)) + a4*(e_elas(1)
     1   + e_elas(2) + e_elas(3))
       return
      end
           
           !c----------------------subroutine noviddsdde 
