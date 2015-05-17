c *** calculates the SHG radiated efficiency
c *** input X_{ijk} in fort.25 from KK
c *** info in fort.15
c *** eps(w) bulk in fort.14 from KK
c **** output fort.36 with Rpp, Rps, Rsp, Rss
      PROGRAM rad
       implicit real (a-b)
     &        , complex (c)
     &        , real (d-h,p-z)
       parameter(max=2000)
       dimension w(max)
     & ,radpp(max),radps(max),radsp(max),radss(max)
       complex cXzzz(max),cXzxx(max),cXxxz(max),cXxxx(max)
     & ,cE(2,3),cepsb(max),cepss(max)
     & ,RESULT,EPS1C,OMEGA
       DATA PI/3.141592654/
c **** inits
c **** constants
c          timei = mclock()
	  ci = cmplx(0.,1.)
	  cero = cmplx(0.,0.)
	  cuno = cmplx(1.,0.)
c *** units and prefactor
        hbar = 1.0545e-27 !erg s
        sc = 2.997925e10 !cm/s
	Ab = 0.543e-8 ! Angstroms
	a0 = 5.43e-8 !A for Si
        Ry = 13.606 !eV Rydberg
	Area = (1./(2.*sqrt(2.))**2)*2.*sqrt(3.)
	 gamma = (1.e-40)*(1.e7)*(1.6021e-12)/(hbar**2*sc**3)
c **
	 gamma = gamma*32.*pi**3  !32pi^3 \approx 1000
c **
	 gamma = gamma*a0**2*Ab**3*(2.*Ry)**3
c *** reads angels
c       write(*,*)'2-3layer=>1,2,teta,phi,f1,f2'
c       read(*,*)layer,teta,phi,f1,f2
c       write(*,*)'2-3layer=>1,2,teta,phi,TB=>1 PW=>2,bulk(1) vac(2) 3L(3)'
c	read(*,*)layer,teta,phi,itb,model
      layer = 1
      teta = 65.
      phi = 30.
      itb = 2
      model= 1
c ** or vs teta
c      write(*,*)'Maxt,phi'
c      read(*,*)Maxt,phi
c *** reads from X'_{ijk} output
	if(itb.eq.1) then
c       open(unit=41,file='/u/bms/roma/shg/si111/1x1/fort.61'
c    & ,form='formatted',status='old')
	read(3,*)aaa,Dw,aaa
	Emax = 12.95
	end if
c *** from X_{ijk} to \chi_{ijk}
       cpf = -ci*(Ab/a0)**2/(2.*Area)
       cpf = cpf * (2.*Ry)**2
c      cpf = 6.933 * a0**3 * cpf  ![\chi] in 10^{-22} m^2 V^{-1}
c        write(*,*)cpf
c** prefactor
c **  1/A -> eta=#atoms/cell/A(cell) = 1/A_{1x1} for any reconstruction
	   if(itb.eq.1) then
	 gamma = (gamma/(Area)**2)*1.e21  !in 10^{-18} cm^2/Watt
	   end if
c ** For PW we need a factor of 1/(10.26)^{4} since R = 1/A^2
c ** and a factor of 2Ry for each P => (2Ry)**3
	   if(itb.eq.2) then
	 fv1 = (2.*Ry)**9/(a0/Ab)**4
	 fv1 = fv1 * 1.e7 *  1.6021e-12
	 fv2 = Ab**5/(hbar**2 * sc**3)
	 gamma = fv1 * fv2
	 gamma = (gamma*32.*pi**3/Area**2)!cm^2/Watt
c *** in terms of \chi in esu.cm
	 hbeV = 6.582e-16 !hbar in eV.s
	 gammaesu = (32.*pi**3/(hbeV**2*sc**3))
c *** from X_{ijk} to \chi_{ijk}
	  cpf = ci * (2.*Ry)**5 * (ab/a0)**5/Area
c ***
c	  write(*,*)gammaesu
	   end if
c **** reads eps(w) NOTICE that it must have the same broadening
c **** in KK transform. Also eps(w)  .05 < w < 13 and w must be same
c **** intervals as those of X_{ijk}
	 if(itb.eq.2)Nint=1200
	 if(itb.eq.1)Nint=Emax/Dw
	 do i = 1,Nint
c	  read(14,*)w(i),epsr,epsi
 	  read(14,*)xx,epsr,epsi
	  cepsb(i) = epsr + ci * epsi
	 end do
c **** reads X_{ijk} Re and Im
c *** zzz,  zxx,  xxz, xxx,
	 do i = 1,Nint
 	  read(25,*)w(i),x1,y1,x2,y2,x3,y3,x4,y4
	 cXzzz(i) = x1 + ci * y1
	 cXzxx(i) = x2 + ci * y2
	 cXxxz(i) = x3 + ci * y3
	 cXxxx(i) = x4 + ci * y4
c ** 3-layer
	 if(layer.eq.2) then
	  read(13,*)xx,epsr,epsi
	  cepss(i) = epsr + ci * epsi
	 end if
c **
	 end do
c** cycle for teta
c      do it =1,Maxt+1
c      ang = float(it-1) !degrees
c       teta = ang
      teta = teta * pi/180. ! degrees to radians
      phi  = phi  * pi/180. ! degrees to radians
c ** wave vectors
      qz = cos(teta)
      Q  = sin(teta)
c ** prefactor
      sect2=1./qz**2
      pf = gamma*sect2
c      pfesu = gammaesu*sect2*1.e7
c      pfesu = pfesu*1.e21 ! R in 10^{-21} cm^2/W
      pfesu = (32.*(pi**3))
      pfesu = pfesu/((hbeV**2)*(sc**3)*(cos(teta)**2)*1.e-7*1.e-21)
c      write (*,*) pfesu
c ** cycle for radiation
c ** since we need w and 2w we only go in X_{ijk} up to w/2
	 maxi=Nint/2
c         write(*,*)maxi,Nint
	do iw = 1,maxi
c**   Fresnel factors and more ***
c** eps_s,b(w) and eps_s,b(2w)
       cepsbw = cepsb(iw)
       cepsb2w = cepsb(2*iw)
c       write(*,*) cepsb2w
c** bulk definition
	if(model.eq.1) then
       cepssw  = cepsbw
       cepss2w = cepsb2w
	end if
c** vacuum definition
	if(model.eq.2) then
       cepssw  = (1.,0.)
       cepss2w = (1.,0.)
	end if
c** three-layer
	if(model.eq.3) then
       cepssw = cepss(iw)
       cepss2w = cepss(2*iw)
	end if
c** k_z^s,b(w)
       ckzbw  = csqrt(cepsbw - Q**2)
       ckzbw  =   real(ckzbw) + ci * abs(aimag(ckzbw))   !Im > 0
       ckzb2w = csqrt(cepsb2w - Q**2)
       ckzb2w =   real(ckzb2w) + ci *abs(aimag(ckzb2w)) !Im > 0
c**
       ckzsw  = csqrt(cepssw - Q**2)
       ckzsw  =   real(ckzsw) + ci * abs(aimag(ckzsw))   !Im > 0
       ckzs2w = csqrt(cepss2w - Q**2)
       ckzs2w =   real(ckzs2w) + ci *abs(aimag(ckzs2w)) !Im > 0
c** t_s(w) vacuum->surface s-transmition
      ctsvsw = 2.*qz/(qz + ckzsw)
c** t_s(2w) vacuum->surface s-transmition
      ctsvs2w = 2.*qz/(qz + ckzs2w)
c** t_p(w) vacuum->surface p-transmition
      ctpvsw = 2.*qz/(qz*cepssw + ckzsw)
c** t_p(2w) vacuum->surface p-transmition
      ctpvs2w = 2.*qz/(qz*cepss2w + ckzs2w)
c** t_s(w) surface->bulk s-transmition
      ctssbw = 2.*ckzsw/(ckzsw + ckzbw)
c** t_s(2w) surface->bulk s-transmition
      ctssb2w = 2.*ckzs2w/(ckzs2w + ckzb2w)
c** t_p(w) surface->bulk p-transmition
      ctpsbw = 2.*ckzsw/(cepsbw*ckzsw + cepssw*ckzbw)
c** t_p(2w) surface->bulk p-transmition
      ctpsb2w = 2.*ckzs2w/(cepsb2w*ckzs2w + cepss2w*ckzb2w)
c**
c**  for chi in esu.cm=cm^3/e
      if(itb.eq.2)then
      esu = 2.08e-15 * (a0/1.e-8)**3
      cesu = cpf * esu / w(iw)**3
      cXzzz(iw)= cesu * cXzzz(iw)
      cXzxx(iw)= cesu * cXzxx(iw)
      cXxxz(iw)= cesu * cXxxz(iw)
      cXxxx(iw)= cesu * cXxxx(iw)
c      write(*,*)w(iw),cesu,cXzzz(iw)
      end if
c** r_{ij}
      cradpp = Q*cepsb2w*(  Q**2 * cepsbw**2 * cXzzz(iw)
     &                     + ckzbw**2 *cepssw**2 * cXzxx(iw) )
     &        + cepssw*cepss2w*ckzbw*ckzb2w
     &          * ( - 2.*Q*cepsbw* cXxxz(iw)
     &             + ckzbw * cepssw * cXxxx(iw) * cos(3.*phi) )
c ** s -> P
      cradsP = Q * cepsb2w * cXzxx(iw)
     &         - cepss2w * ckzb2w * cXxxx(iw) * cos(3.*phi)
c ** p -> S
      cradpS = -ckzbw**2 *cepssw**2 * cXxxx(iw) * sin(3.*phi)
c ** s -> S
      cradsS =  cXxxx(iw) * sin(3.*phi)
c **
c ** bare \chi_{ijk}
	czzzb = cXzzz(iw)
	czxxb = cXzxx(iw)
	cxxzb = cXxxz(iw)
	cxxxb = cXxxx(iw)
c **
c       czzzb = cpf*cXzzz(iw)/w(iw)**3
c       czxxb = cpf*cXzxx(iw)/w(iw)**3
c       cxxzb = cpf*cXxxz(iw)/w(iw)**3
c       cxxxb = cpf*cXxxx(iw)/w(iw)**3
c ** with fresnel
c ** Bulk
	 czzzfb = Q**2*czzzb
	 czxxfb = ckzbw**2*czxxb
	 cxxzfb = 2.*ckzbw*ckzb2w*cxxzb
c ** 3-layer
      czzzf = ctpsb2w*ctpsbw**2*cepsb2w* Q**2*cepsbw**2*czzzb
      czxxf = ctpsb2w*ctpsbw**2*cepsb2w* ckzbw**2 *cepssw**2 * czxxb
      cxxzf = 2.*ctpsb2w*ctpsbw**2*cepssw*cepss2w*ckzbw
     &           *ckzb2w*cepsbw*cxxzb
c ** R_pP with each X_ijk
      cradppzzz = Q*cepsb2w* Q**2*cepsbw**2*cXzzz(iw)
      cradppzxx = Q*cepsb2w* ckzbw**2 *cepssw**2 * cXzxx(iw)
      cradppxxz = -2.*cepssw*cepss2w*ckzbw*ckzb2w*Q*cepsbw*cXxxz(iw)
c**
c** radiated efficiency in 10^{-21} cm^2/watt
c** For us {\cal R}_{out,in}, usually {\cal R}_{in,out}, However there
c** is a much better way of knowing which polarization is 'in-w photon'
c** or 'out-2w photon', by usin small letters for 'in' and capital
c** letters for 'out'. Thus our {\cal R}_{ij} -> {\cal R}_{Ij}, i.e.
c** {\cal R}_{Pp}, {\cal R}_{Ps}, {\cal R}_{Sp} and {\cal R}_{Ss}, equal to
c** {\cal R}_{pP}, {\cal R}_{sP}, {\cal R}_{pS} and {\cal R}_{sS}
c **
	ctpP = ctpvs2w*ctpsb2w*(ctpvsw*ctpsbw)**2/(2.*ci)
	ctsP = ctpvs2w*ctpsb2w*(ctsvsw*ctssbw)**2/(2.*ci)
	ctpS = ctsvs2w*ctssb2w*(ctpvsw*ctpsbw)**2/(2.*ci)
	ctsS = ctsvs2w*ctssb2w*(ctsvsw*ctssbw)**2/(2.*ci)
c ** for chi in esu.cm
!!      write(*,*)w(iw),cradps
	radpP(iw) = pfesu*cabs(ctpP*cradpP)**2*w(iw)**2
	radsP(iw) = pfesu*cabs(ctsP*cradsP)**2*w(iw)**2
	radpS(iw) = pfesu*cabs(ctpS*cradpS)**2*w(iw)**2
	radsS(iw) = pfesu*cabs(ctsS*cradsS)**2*w(iw)**2
c **
c       radpP(iw) = pf*cabs(ctpP*cradpP)**2/w(iw)**4
c       radsP(iw) = pf*cabs(ctsP*cradsP)**2/w(iw)**4
c       radpS(iw) = pf*cabs(ctpS*cradpS)**2/w(iw)**4
c       radsS(iw) = pf*cabs(ctsS*cradsS)**2/w(iw)**4
c ** R_pP with each X_ijk
	radpPzzz = pf*cabs(ctpP*cradpPzzz)**2/w(iw)**4
	radpPzxx = pf*cabs(ctpP*cradpPzxx)**2/w(iw)**4
	radpPxxz = pf*cabs(ctpP*cradpPxxz)**2/w(iw)**4
c**  shifted energy
	 ws = 2.*w(iw)
       if(itb.eq.2)ws= 2.*w(iw)
       write(301,*)w(iw),real(ctsvsw),aimag(ctsvsw)
       write(302,*)w(iw),real(ctpvsw),aimag(ctpvsw)
       write(303,*)w(iw),real(ctssbw),aimag(ctssbw)
       write(304,*)w(iw),real(ctpsbw),aimag(ctpsbw)
!!       write(69,*)ws,cabs(cradpP),cabs(cradpS),cabs(cradsP),cabs(cradsS)
       write(69,*)ws,real(cradpP),aimag(cradpP)
       write(36,*)ws,radpP(iw),radsP(iw),radpS(iw),radsS(iw)
       write(38,*)ws,radpPxxz,radpPzxx,radpPzzz
	write(91,901)ws,cabs(czzzb),cabs(czxxb),cabs(cxxzb),cabs(cxxxb)
c       write(92,901)ws,czzzf,czxxf,cxxzf
	write(93,901)ws,cabs(czzzfb),cabs(czxxfb),cabs(cxxzfb)
c	write(92,901)w(iw),cXzzz(iw),cXzxx(iw),cXxxz(iw),cXxxx(iw)
c       write(92,901)ws,cabs(czzzf),cabs(czxxf),cabs(cxxzf)
c       write(93,901)ws,tpPzzz,tpPzxx,tpPxxz,tsP
c***                       Pp        Ps=sP     Sp=pS     Ss
c     write(36,901)ang,ws,radpp(iw),radps(iw),radsp(iw),radss(iw)
c     write(91,901)ang,ws,rzzz,rzxx,rxxz
c**
c      end if
c**
	end do !iw
c*****
c*****
c     write(36,902)
c     write(91,902)
c       end do !it
c*****
 901    format(f8.3,10e13.5)
 902    format()
	 end
c*****
