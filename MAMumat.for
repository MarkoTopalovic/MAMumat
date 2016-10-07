      subroutine umat(stress,statev,ddsdde,sse,spd,scd,
     1 rpl,ddsddt,drplde,drpldt,stran,dstran,
     2 time,dtime,temp,dtemp,predef,dpred,materl,ndi,nshr,ntens,
     3 nstatv,props,nprops,coords,drot,pnewdt,celent,
     4 dfgrd0,dfgrd1,noel,npt,kslay,kspt,kstep,kinc)
c
      include 'aba_param.inc'
      common /glavni/ ipoziv
      common /ak/ akapa(100000,10)
      common /comeplast/  ceplast(100000,8,6)
      common /comdtplast/ dtplast(100000,8,6)

      character*80 materl

      dimension a(17),statev(nstatv),props(nprops),
     1 stran(ntens),dstran(ntens),
     1 stress(ntens),s_dev(ntens),d_stres(ntens),
     1 d_eplas(ntens),dt_plas(ntens),dt_plas0(ntens),
     1 eplas(ntens),eplas0(ntens),
     1 ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     1 time(2),predef(1),dpred(1),  
     1 coords(ndi),drot(ndi,ndi),
     1 dfgrd0(ndi,ndi),dfgrd1(ndi,ndi),  
     1 a_mu(ntens),kroneker(ntens), 
     1 replas(ntens),kroneker2(ntens,ntens)  !deplas(ntens) 

      parameter (one=1.0d0,two=2.0d0,three=3.0d0,six=6.0d0,zero=0.0d0)
      data newton,toler,temp0,coef,yield0/5,1.d-6,273.,0.0d0,25.0d0/    

C -----------------------------------------------------------
c
c     paneerselvam phd   (page 165 tens constants)
c
      a1 = 70.3e2                  
      a2 = 0.036e5 !7.5e4          
      a3 = -1.81e5
      a4 = -0.906e5
      a5 = -0.6046e5
      a6 = 0.
      a7 = 0.
      a8 = 0.
      a9 = 0.
      x = 1.03e-3
      a_l = 2          
      alfa = 14.2  !17.1   
      h = 50.
      beta = 1.34e4  !8.23e4
      a_m  = one
      gama = -5.
      
      alpha_t  = 4.e-5
c -----------------------------------------------------------
c -----------------------------------------------------------
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

!     NDI: Broj direkthih komponenti napona u datom trenutku
!     NDI=3
!     NSHR: Broj smicajnih komponenti napona u datom trenutku
!      NSHR=3
!     NTENS: veli?ina niza napona ili deformacija (NDI + NSHR)
!     NTENS=6

!     NOEL: Broj elementa
!     NPT: Broj integracione ta?ke 


      
!  UILAZNE VELICINE U UMAT

!     STRAN(NTENS): Niz koji sadr�i ukupne deformacije na po?etku inkrenenta
!     DSTRAN(NTENS): Niz inkremenata deformacija
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  IZLAZNE VELICINE
      
!     STRESS(NTENS): TENZOR NAPONA OVO JE I ULAZ I IZLAZ
!     NA POCETKU DOBIJEMO TENZOR NAPONA NA POCETKU INKREMENTA
!     PA U UMATU TREBA DA IYRACUNAMO TENZOR NAPONA NA KRAJU INKREMENTA    

!     DDSDDE(NTENS,NTENS)   PARCIJALNI IZVOD INKREMENTA NAPONA PO INKREMENTU DEFORMAICJE
!     OVO JE CISTO IZLAZNA VELICINA KOJA SE RACUNA U UMATU
!     STRESS JE MNOGO VAZNIJI OD DDSDDE
!     DDSDDE SE KORISTI U ABAQUSU ZA PROVERU KONVERGENCIJE       


!  OSTALE ULAZNE ILI IZLAZNE VELICINE NE KORISTIMO


      
!  UNUTRASNJE PROMENLJIVE
!     STATEV NIZ UNUTRASNJIH PROMENLJIVIH
!     NSTATV BROJ UNUTRASNJIH PROMENLJIVIH
!     OVO VISE NE KORISTIMO SADA SVE CUVAMO U COMMON STRUKTURAMA      
      

c********************************************************************
c
c               predictor corector algorithm
c
c********************************************************************     

cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c***************      1) predictor phase      ***********************  
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      call hyperconstitutive(a,ddsdde,ntens,stran)
           do k1 = 1,ntens
              do k2 = 1,ntens
                 stress(k2)=stress(k2)+ddsdde(k2,k1)*dstran(k1) ! 6.19 
              enddo                                        
          enddo


c------------------ compute  loading surface f-----------------------
      call loadingf(f1,stress,a,ntens,ndi,s_dev,a_j2,a_mu) ! 6.21

      a_kapa0  = akapa(noel,npt)  ! 6.20   
      if (a_kapa0.lt.toler) then                 
         a_kapa0 = zero           
      endif
      f   = f1 - h*a_kapa0     
    
      if ((f.le.zero).or.(a_j2.lt.yield0)) then   
       write(6,*) 'ELASTIC'              
        goto 52           
      endif    
c------------------  end of elastic predictor ----------------------
!********************************************************************
!********************************************************************
c------------------     2) corrector phase    ----------------------
c    eplas - eplas0    ->   epsilon_0_n+1 - epsilon_n
c    deplas            ->   korekcija  (iteracija k)
c    a_kapa - a_kapa0 ->   kapa_0_n+1 - kapa_n
c    dkapa            ->   korekcija  (iteracija k)
c-------------------------------------------------------------------

c      zapocinjanje (inicijalizacija neravnoteznih delova)
          do k1=1,ntens
             eplas0(k1)= ceplast(noel,npt,k1)   !staro
             eplas(k1) = eplas0(k1) !
             dt_plas0(k1)=dtplast(noel,npt,k1)!?
          end do 
             a_kapa =  a_kapa0 ! linija 
             d_kapa =  zero       !staro
!      zapocinanje (inicijalizacija neravnoteznih delova)

!-----------------------------------glavna petlja      
        do 22 kewton=1,newton
!-----------------------------------glavna petlja 
!          write(6,*) 'kewton-newton-',kewton
        call loadingf(f1,stress,a,ntens,ndi,s_dev,a_j2,a_mu)
     
          f  = f1 - h*a_kapa     
          f = (abs(f) + f)/two 
       
          if (a_kapa.gt.toler) then
          a_kxl = x+a_kapa**a_l
          else
          a_kxl = x
          endif
          
          dkapa  = (f/beta)*dtime/a_kxl  !novo
          
          skapa = a_kapa-a_kapa0-dkapa
          absskapa = abs(skapa)
            
          do k1=1,ntens        
             replas(k1) = eplas(k1)-eplas0(k1)-dkapa*a_mu(k1)
             
!                write(6,*) 'eplas(k1)', eplas(k1)
!                write(6,*) 'eplas0(k1)', eplas0(k1)
!                write(6,*) 'dkapa', dkapa
!                write(6,*) 'amu', a_mu(k1)
!             write(6,*) 'replas', replas(k1)
          enddo    
     
        replas_int = (replas(1)**2+replas(2)**2+replas(3)**2+
     1    two*replas(4)**2+two*replas(5)**2+two*replas(6)**2)**0.5

      if ((replas_int.gt.tol1).or.(f.gt.tol2)) then
	  
		if ((replas_int.gt.tol1).and.(NPT.eq.2)) then
           write(6,*) 'R', replas_int, 'I=', kewton,'NPT',NPT
		endif
!			if (f.gt.tol2) then
!           write(6,*) 'f', f, 'I=', kewton, 'NPT',NPT
!		endif
		   
		   
      call  obnovi_s_kapa(stress,d_stres,d_eplas,ddsdde,ntens,
     1      ndi,a,a_kapa,replas,skapa,d_kapa,
     2      dtime,npt,noel,toler)


       do k1=1,ntens
             eplas(k1)   = eplas(k1) + d_eplas(k1)  
!             write(6,*) 'd_eplas posle obnovi s drugi NAN', d_eplas(k1) 
             dt_plas(k1) = d_eplas(k1)/dtime  ! brz.plast.def. 
             
             do k2 = 1,ntens
               stress(k1)= stress(k1)-ddsdde(k1,k2)*d_eplas(k2)
             enddo  
       enddo  
                    
             a_kapa = a_kapa0 + d_kapa  
             end if                
22    continue        

c  corrector phase -  kraj       
      
 52   continue
 

      do k1=1,ntens
         ceplast(noel,npt,k1)  =  eplas(k1)
         dtplast(noel,npt,k1)  =  dt_plas(k1)
      enddo                   
         akapa(noel,npt)=a_kapa               
       
      return    
      end
C
c-----------------------------------------  umat   kraj
c-----------------------------------------  umat   kraj
      
      
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
           
c----------------------subroutine hyperconstitutive     
c----------------------------------------------------                            
        
      subroutine loadingf(f1,stress,a,ntens,ndi,s_dev,a_j2,a_mu)
     
      include 'aba_param.inc' 
      
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
     1   kroneker(ntens)
      
      parameter (one=1.0d0,two=2.0d0,three=3.0d0,six=6.0d0) 
      alfa =  a(12)
      gama = -5.0d0

      do k2=1, ndi
         kroneker(k2)     = one
         kroneker(k2+ndi) = zero
      end do 
            
c     s_napon invarijante  i funkcija f 
      s_i1 = stress(1)+stress(2)+stress(3)
      s_i2 = (stress(1)**2+stress(2)**2+stress(3)**2)/two+
     1       stress(4)**2+stress(5)**2+stress(6)**2
     1       -(s_i1**2)/two
      s_i3 = stress(1)*stress(2)*stress(3)+
     1 two*stress(4)*stress(5)*stress(6)-stress(1)*stress(6)**2- 
     1 stress(2)*stress(5)**2-stress(3)*stress(4)**2

       do  k1=1,ndi
           s_dev(k1)     = stress(k1) - s_i1/three
           s_dev(k1+ndi) = stress(k1+ndi)
       enddo
  
      a_j2 = ( s_dev(1)**2+s_dev(2)**2+s_dev(3)**2)/two+
     1          s_dev(4)**2+s_dev(5)**2+ s_dev(6)**2 
     
      f1 = s_i1*s_i2 + alfa*s_i3 
             
      do k1=1,ntens
           a_mu(k1) = gama*kroneker(k1)+( 0.5*s_dev(k1)/(a_j2**0.5) )
      enddo   
      return
      end
c        
c---------------------------------subroutine loadingf      
c----------------------------------------------------  
      
      subroutine  obnovi_s_kapa(s_napon,d_stres,d_eplas,
     1      ddsdde,ntens,ndi,a,a_kapa,replas,skapa,d_kapa,
     2      dtime,npt,noel,toler)
      
      include 'aba_param.inc'
      
!     brojne konstante   
!     real zero,one,two,three,four,six,omega
!     brojne konstante 
!         
!     ulazne promenljive   
!     real a(17), a2,x,a_l,alfa,h,beta,m,gama
!           
!     integer ntens,ndi,npt,noel
!     real toler
!     real ddsdde(ntens,ntens)
!     real  s_napon(ntens)
!     real a_kapa,dtime,skapa
!     real replas(ntens)
!     real d_stres(ntens)
!     ulazne promenljive  
!
!     unutrasnje promenljive
!     integer k1,k2,k3
!     real kroneker(ntens)
!     real kroneker2(ntens,ntens),p1(ntens,ntens)
!     real proizvod
!     real  f1, f2,f
!     real a_j2
!     real s_dev(ntens),a_mu(ntens),b(ntens)
!     real a_kxl,y1
!     real z(ntens),proizvod2(ntens)
!     real xx2(ntens,ntens),x_matrica(ntens,ntens)
!     real x2(ntens,ntens),x1(ntens),y2(ntens)      
!     unutrasnje promenljive 
!
!     izlazne promenljive
!     real d_eplas(ntens)
!     real d_kapa
!     izlazne promenljive
      
c
c----------------  ulaz    replas,  skapa, s_napon   
c----------------  izlaz   d_eplas,  d_kapa, s_napon 
c
      dimension b(ntens),s_napon(ntens),s_dev(ntens),d_stres(ntens), 
     1 kroneker2(ntens,ntens),ddsdde(ntens,ntens),z(ntens),
     1 p1(ntens,ntens),a(17),d_eplas(ntens),
     2 proizvod2(ntens),replas(ntens),
     3 xx2(ntens,ntens),x_matrica(ntens,ntens),a_mu(ntens),
     4 x2(ntens,ntens),x1(ntens),y2(ntens),kroneker(ntens) 
     
      parameter (one=1.0d0,two=2.0d0,three=3.0d0,six=6.0d0,zero=0.0d0)
          
      a2   = a(2)
      x    = a(10)
      a_l  = a(11)
      alfa = a(12)
      h    = a(13)
      beta = a(14)
      m    = a(15)
      gama = a(16)
      
      omega = zero
c -----------------------------------------------------------

        call loadingf(f1,s_napon,a,ntens,ndi,s_dev,a_j2,a_mu)

        if (a_j2.gt.six) then
     
        f2 = h*a_kapa
        f  = f1 - f2     
        f = (abs(f) + f)/two 
            
       do k1=1, ntens
           do k2=1, ntens
             p1(k1,k2)=zero
             kroneker2(k1,k2)=zero
           end do
       end do
           do k2=1, ntens
             kroneker2(k2,k2)=one
           end do 
                             
           do k2=1, ndi
             kroneker(k2)=one
             kroneker(k2+ndi)=zero
           end do  

       do k1=1, ndi ! 5.60
           do k2=1, ndi
             p1(k1,k2)=-one/three
           end do
           p1(k1,k1)=two/three
           p1(ndi+k1,ndi+k1)=one
       end do    
             ! gde su ovde ABC kako je to uprosceno? !6.32
       z(1)=(1./2.)*((s_napon(4)**2+s_napon(5)**2+s_napon(6)**2)
     1       -two*s_napon(1)*s_napon(2)-three*s_napon(2)*s_napon(3)
     2       -two*s_napon(1)*s_napon(3)-s_napon(2)**2-s_napon(3)**2)
     3       +alfa*(s_napon(2)*s_napon(3)-s_napon(6)**2)
     
       z(2)=(1./2.)*((s_napon(4)**2+s_napon(5)**2+s_napon(6)**2)
     1       -two*s_napon(1)*s_napon(2)-three*s_napon(1)*s_napon(3)
     2       -two*s_napon(2)*s_napon(3)-s_napon(1)**2-s_napon(3)**2)
     3       +alfa*(s_napon(1)*s_napon(3)-s_napon(5)**2)
       
       z(3)=(1./2.)*((s_napon(4)**2+s_napon(5)**2+s_napon(6)**2)
     1       -two*s_napon(2)*s_napon(3)-three*s_napon(1)*s_napon(3)
     2       -two*s_napon(2)*s_napon(3)-s_napon(1)**2-s_napon(2)**2)
     3       +alfa*(s_napon(1)*s_napon(2)-s_napon(4)**2)
       
       z(4)=(1./2.)*(two*s_napon(1)*s_napon(4)+two*s_napon(2)*s_napon(4)
     1       +two*s_napon(3)*s_napon(4))
     2       +two*alfa*(s_napon(5)*s_napon(6)-s_napon(3)*s_napon(4))
       
       z(5)=(1./2.)*(two*s_napon(1)*s_napon(5)+two*s_napon(2)*s_napon(5)
     1       +two*s_napon(3)*s_napon(5))
     2       +two*alfa*(s_napon(4)*s_napon(6)-s_napon(2)*s_napon(5))

       z(6)=(1./2.)*(two*s_napon(1)*s_napon(6)+two*s_napon(2)*s_napon(6)
     1       +two*s_napon(3)*s_napon(6))
     2       +two*alfa*(s_napon(5)*s_napon(4)-s_napon(1)*s_napon(6)) 
         
       proizvod =    s_dev(1)**2+s_dev(2)**2+s_dev(3)**2+
     1            two*(s_dev(4)**2+s_dev(5)**2+s_dev(6)**2)

       do k1=1, ntens
           proizvod2(k1)=zero
           do k2=1, ntens
             proizvod2(k1)=proizvod2(k1)+s_dev(k2)*p1(k2,k1)
           end do
       end do 
        
      
          if (a_kapa.gt.toler) then
          a_kxl = x+a_kapa**a_l
          else
          a_kxl = x
          endif
      
       do k1=1, ntens
           do k2=1, ntens   !6.36
             xx2(k1,k2)=p1(k1,k2)/((two*proizvod)**0.5)
     3       -(s_dev(k1)* proizvod2(k2))/(two*( (proizvod/2)**1.5 ) )

           end do
       end do
       
       do k1=1, ntens
           do k2=1, ntens
                x2(k1,k2)=- dtime !*((f/beta)**(m-1))  = 1
     1          *(z(k1)/beta)*( a_mu(k2)/a_kxl )
     2          +( xx2(k1,k2)*f/beta )
     3          /(  (x+a_kapa**a_l)*(one-omega)  )
           end do
       end do
       
       do k1=1, ntens
            x1(k1)= one          !*((f/beta)**(m-1))    = 1
     1          *(h*a_mu(k1)*dtime)/(beta*(x+a_kapa**a_l))
     2          +(f/beta)*(a_mu(k1)*a_l*(a_kapa**(a_l-one))
     3          *dtime)/((x+(a_kapa)**a_l)**2)
       end do     
       
       do k1=1, ntens
            y2(k1)=-one*(f/beta)*(z(k1)*dtime)/(beta*a_kxl)
       end do     
      
            y1  = one + one  !*((f/beta)**(m-1))    = 1
     1          *(h*dtime)/(beta*a_kxl)
     2          +dtime*(f/beta)* a_l*( a_kapa**(a_l-one) ) 
     3          /(a_kxl**2)
!         write(6,*) 'y1 =' ,y1
!         write(6,*) 'skapa =' ,skapa
       
c   formula (1.2) na algoritmu - pocetak

       do k1=1, ntens
           do k2=1, ntens
              x_matrica(k1,k2)=kroneker2(k1,k2)
              do k3=1, ntens
                x_matrica(k1,k2) = x_matrica(k1,k2) - 
     1            x2(k1,k3)*ddsdde(k3,k2) + 
     1            (one/y1)*x1(k1)*y2(k3)*ddsdde(k3,k2)                 
              end do
           end do
       end do
      
       do k1=1, ntens
             b(k1)=-replas(k1)+(skapa/y1)*x1(k1)
!             write(6,*) 'ogromno za NaN replas(k1) =' ,replas(k1)
!            write(6,*) 'ogromno za NaN b(k1) =' ,b(k1)

       end do   
        
c  inverzna matrica iz (6.42)       
         call gauss(x_matrica,ntens)                
c  19-07-2013    -   vazi za gauss

       do k1=1,ntens
          d_eplas(k1) = zero   ! treba 0.
          do k2=1,ntens
            d_eplas(k1) = d_eplas(k1) + x_matrica(k1,k2)*b(k2)
          end do
         
       end  do

       do k1=1,ndi
          d_eplas(k1)=d_eplas(k1)-
     1                (d_eplas(1)+d_eplas(2)+d_eplas(3))/three
     
!      write(6,*) 'prvi NAN d_eplas(k1) =' ,d_eplas(k1)
       end do  
c  19-07-2013    -   vazi za gauss
   
       
   
    
       d_kapa = -(skapa/y1)
       do k1=1, ntens
       d_kapa=d_kapa + (a2/y1)*y2(k1)*d_eplas(k1)
       end do   
       

       
       
       
           
       d_kapa = abs(d_kapa)
       
       do k1 = 1,ntens
            d_stres(k1)= zero
            do k2 = 1,ntens
               d_stres(k1)=d_stres(k1)-ddsdde(k1,k2)*d_eplas(k2)
            enddo  
c            s_napon(k1)=s_napon(k1)+d_stres(k1)   
       enddo    
       
       endif  !  if (a_j2.gt.six) then
c   formula (1.2) na algoritmu - kraj

      return 
      end
      
c-------------------------------------------subroutine  obnovi_s_kapa

! --------------------------------------------------------------------
        subroutine gauss (a,n)       ! invert matrix by gauss method
! --------------------------------------------------------------------
        implicit none
        integer :: n
        real :: a(n,n)
! - - - local variables - - -
        real :: b(n,n), c, d, temp(n)
        integer :: i, j, k, m, imax(1), ipvt(n)
! - - - - - - - - - - - - - -
        b = a
        ipvt = (/ (i, i = 1, n) /)
        do k = 1,n
           imax = maxloc(abs(b(k:n,k)))
           m = k-1+imax(1)
           if (m /= k) then
              ipvt( (/m,k/) ) = ipvt( (/k,m/) )
              b((/m,k/),:) = b((/k,m/),:)
           end if
           d = 1/b(k,k)
           temp = b(:,k)
           do j = 1, n
              c = b(k,j)*d
              b(:,j) = b(:,j)-temp*c
              b(k,j) = c
           end do
           b(:,k) = temp*(-d)
           b(k,k) = d
        end do
        a(:,ipvt) = b
        end 
        