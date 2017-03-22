      subroutine umat(stress,statev,ddsdde,sse,spd,scd,
     1 rpl,ddsddt,drplde,drpldt,stran,dstran,
     2 time,dtime,temp,dtemp,predef,dpred,materl,ndi,nshr,ntens,
     3 nstatv,props,nprops,coords,drot,pnewdt,celent,
     4 dfgrd0,dfgrd1,noel,npt,kslay,kspt,kstep,kinc)
c
      include 'aba_param.inc'
      character*80 materl

      dimension a(17),statev(nstatv),props(nprops),
     1 stran(ntens),dstran(ntens),
     1 stress(ntens),astress(ntens),s_dev(ntens),d_stres(ntens),
     1 d_eplas(ntens),astrain(ntens),
     1 eplas(ntens),e_elas(ntens),eplas0(ntens),
     1 ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     1 time(2),predef(1),dpred(1),  
     1 coords(ndi),drot(ndi,ndi),
     1 dfgrd0(ndi,ndi),dfgrd1(ndi,ndi),  
     1 a_mu(ntens),kroneker(ntens), 
     1 replas(ntens),kroneker2(ntens,ntens)  !deplas(ntens) 

      
        if(NPT.eq.1) then
		write(6,*)time
		endif
		
      return    
      end
C
c-----------------------------------------  umat   kraj
c-----------------------------------------  umat   kraj
      
      
  