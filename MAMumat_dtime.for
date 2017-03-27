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

      character*80 materl

      dimension a(17),statev(nstatv),props(nprops),
     1 stran(ntens),dstran(ntens),
     1 stress(ntens),astress(ntens),s_dev(ntens),d_stres(ntens),
     1 d_eplas(ntens),
     1 eplas(ntens),e_elas(ntens),eplas0(ntens),
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
      beta = 8.23e4     !1.34e4
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
					astress(k2)=stress(k2)
                 stress(k2)=stress(k2)+ddsdde(k2,k1)*dstran(k1) ! 6.19 
              enddo                                        
          enddo


c------------------ compute  loading surface f-----------------------
      call loadingf(f1,stress,a,ntens,ndi,s_dev,a_j2,a_mu) ! 6.21

      a_kapa0  = akapa(noel,npt)  ! 6.20   ! ucitava iz prethodnog koraka
      if (a_kapa0.lt.toler) then                 
         a_kapa0 = zero           
      endif
      f   = f1 - h*a_kapa0     
    
      if ((f.le.zero).or.(a_j2.lt.yield0)) then   
       write(6,*) 'ELASTIC'              
        goto 52           
      endif    
c------------------  end of elastic predictor ----------------------
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c***************      1) corector phase      ***********************  
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
			write(6,*) 'plastic'
      a_kapa = a_kapa0
			do k3 = 1,10		
! ucitava iz prethodnog koraka
		do k1=1,ntens
         eplas0(k1) = ceplast(noel,npt,k1)  
		enddo
		
		do k1 = 1,ntens
          do k2 = 1,ntens
          stress(k2)=astress(k2)+ddsdde(k2,k1)*dstran(k1)*(1.1-k3*0.1)
          enddo                                        
		enddo
		
		call loadingf(f1,stress,a,ntens,ndi,s_dev,a_j2,a_mu)
		
		          f  = f1 - h*a_kapa     
					f = (abs(f) + f)/two 
      write(6,*) 'f1=', f1
      if (a_kapa.gt.toler) then
      a_kxl = x+a_kapa0**a_l
      else
      a_kxl = x
      endif
!      write (6,*) 'dtime=',dtime
      !dkapa  = ((f/beta)**1)/a_kxl  !novo
	  dkapa  = dtime*((f/beta)**1)/a_kxl  !mart 2017
	    
       write (6,*) 'dkapa=',dkapa
      a_kapa = a_kapa + dkapa
      
        do k1=1,ntens
			d_eplas(k1) = dkapa*a_mu(k1)*dtime
			eplas(k1) = eplas0(k1)+d_eplas(k1)
			e_elas(k1) = dstran(k1)-eplas(k1)
!      write(6,*)  'a_mu(',k1,')=', a_mu(k1)
!      write(6,*)  'eplas(',k1,')=', eplas(k1)
!      write(6,*)  'dstran(',k1,')=', dstran(k1)
!  	  write(6,*)  'e_elas(',k1,')=', e_elas(k1)
!			write(6,*)  '-------------------------------------------------------'
		enddo
      call hyperconstitutive(a,ddsdde,ntens,e_elas)
			do k1 = 1,ntens
				do k2 = 1,ntens
             stress(k2)=stress(k2)+ddsdde(k2,k1)*e_elas(k1) 
				enddo                                        
			enddo
		  

		  skapa = a_kapa-a_kapa0-dkapa
            
          do k1=1,ntens        
             replas(k1) = eplas(k1)-eplas0(k1)-dkapa*a_mu(k1)
!              write(6,*)  'eplas(',k1,')=', eplas(k1)
!              write(6,*)  'eplas0(',k1,')=', eplas0(k1)
!              write(6,*)  'a_mu(',k1,')=', a_mu(k1)
!				write(6,*)  'dstran(',k1,')=', dstran(k1)
!              write(6,*) 'dkapa', dkapa
          enddo 
			
		replas_int = (replas(1)**2+replas(2)**2+replas(3)**2+
     1    two*replas(4)**2+two*replas(5)**2+two*replas(6)**2)**0.5	
			
		if ((replas_int.lt.tol1).or.(skapa.lt.tol2)) then
!		zadovoljena konvergencija
			write(6,*) 'konvergira', 'k3=', k3
      write(6,*)  'Ri=', replas_int, 'SK', skapa
!		call xit
		goto 52
		endif	
			
			enddo

c  corrector phase -  kraj       
      
 52   continue
 

      do k1=1,ntens
         ceplast(noel,npt,k1)  =  eplas(k1)
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
c-

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
      write(6,*) 's_i1=',s_i1,'s_i2=',s_i2,'s_i3=',s_i3       
      do k1=1,ntens
           a_mu(k1) = gama*kroneker(k1)+( 0.5*s_dev(k1)/(a_j2**0.5) )
      enddo   
      return
      end
c        
c---------------------------------subroutine loadingf      
c----------------------------------------------------  
 