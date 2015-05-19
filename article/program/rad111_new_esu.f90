! calculates the SHG radiated efficiency
! input X_{ijk} in fort.25 from KK
! eps(w) bulk in fort.14 from KK
! eps(w) layer in fort.13 from KK
! output fort.36 with Rpp, Rps, Rsp, Rss
PROGRAM rad
implicit real (a-b), complex (c), real (d-h,p-z)
parameter(max=2000)
dimension w(max),radpp(max),radps(max),radsp(max),radss(max)
complex cXzzz(max),cXzxx(max),cXxxz(max),cXxxx(max),cE(2,3),cepsb(max),cepss(max),RESULT,EPS1C,OMEGA
DATA PI/3.141592654/
! **** inits
! **** constants
! timei = mclock()
ci = cmplx(0.,1.)
cero = cmplx(0.,0.)
cuno = cmplx(1.,0.)
! *** units and prefactor
hbar = 1.0545e-27 !erg s
hbeV = 6.582e-16 !erg s
sc = 2.997925e10 !cm/s
ab = 0.53e-8 ! cm
a0 = 5.34e-8 ! for Si
Ry = 13.606  ! Rydberg (eV)
! **  1/A -> eta=#atoms/cell/A(cell) = 1/A_{1x1} for any reconstruction
Area = (1./(2.*sqrt(2.))**2)*2.*sqrt(3.)
! ** to put \chi in esu.cm
esu = 2.08e-15*(a0/1.d-8)**3
gamma = (32.*pi**3/(hbeV**2*sc**3))*1.d7 !in cm^2/W
! **
! *** reads angels
!       write(*,*)'2-3layer=>1,2,teta,phi,f1,f2'
!       read(*,*)layer,teta,phi,f1,f2
!       write(*,*)'2-3layer=>1,2,teta,phi,TB=>1 PW=>2,bulk(1) vac(2) 3L(3)'
! read(*,*)layer,teta,phi,itb,model
layer = 2
teta = 65.
phi = 30.
itb = 2
model= 3
! ** or vs teta
!      write(*,*)'Maxt,phi'
!      read(*,*)Maxt,phi
! *** reads from X'_{ijk} output
if(itb.eq.1) then !SETB
  !  open(unit=41,file='/u/bms/roma/shg/si111/1x1/fort.61',form='formatted',status='old')
  read(3,*)aaa,Dw,aaa
  Emax = 12.95
  ! *** from X_{ijk} to \chi_{ijk}
  cz = cuno* (2.*Ry)**2 * (ab/a0)**2/Area
end if
! **
if(itb.eq.2) then !PW
  ! *** in terms of \chi in esu.cm
  ! *** from X_{ijk} to \chi_{ijk}
  cz = ci * (2.*Ry)**5 * (ab/a0)**5/Area
  ! ***
end if
! **** reads eps(w) NOTICE that it must have the same broadening
! **** in KK transform. Also eps(w)  .05 < w < 13 and w must be same
! **** intervals as those of X_{ijk}
if(itb.eq.2)Nint=2000
if(itb.eq.1)Nint=Emax/Dw-1
write(*,*)Nint
do i = 1,Nint
  ! read(14,*)w(i),epsr,epsi
  read(14,*)xx,epsr,epsi
  cepsb(i) = epsr + ci * epsi
end do
! **** reads X_{ijk} Re and Im
! *** zzz,  zxx,  xxz, xxx,
do i = 1,Nint
  read(25,*)w(i),x1,y1,x2,y2,x3,y3,x4,y4
  cXzzz(i) = x1 + ci * y1
  cXzxx(i) = x2 + ci * y2
  cXxxz(i) = x3 + ci * y3
  cXxxx(i) = x4 + ci * y4
  ! ** 3-layer
  if(layer.eq.2) then
    read(13,*)xx,epsr,epsi
    cepss(i) = epsr + ci * epsi
  end if
  ! **
end do
!** cycle for teta
! do it =1,Maxt+1
!   ang = float(it-1) !degrees
!   teta = ang
teta = teta * pi/180. ! degrees to radians
phi  = phi  * pi/180. ! degrees to radians
! ** wave vectors
qz = cos(teta)
Q  = sin(teta)
! ** prefactor
sect2=1./qz**2
pf = gamma*sect2
pf = pf*1.e21 ! R in 10^{-21} cm^2/W for \chi^s in esu.cm
pf = 1
write(*,*)'pf=',pf
! ** cycle for radiation
! ** since we need w and 2w we only go in X_{ijk} up to w/2
maxi=Nint/2
! write(*,*)maxi,Nint
do iw = 1,maxi          
  !**   Fresnel factors and more ***
  !** eps_s,b(w) and eps_s,b(2w)
  cepsbw = cepsb(iw)
  cepsb2w = cepsb(2*iw)
  !** bulk definition
  if(model.eq.1) then
    cepssw  = cepsbw
    cepss2w = cepsb2w
  end if
  !** vacuum definition
  if(model.eq.2) then
    cepssw  = (1.,0.)
    cepss2w = (1.,0.)
  end if
  !** three-layer
  if(model.eq.3) then
    cepssw = cepss(iw)
    cepss2w = cepss(2*iw)
  end if
  !** k_z^s,b(w)
  ckzbw  = csqrt(cepsbw - Q**2)
  ckzbw  = real(ckzbw) + ci * abs(aimag(ckzbw))   !Im > 0
  ckzb2w = csqrt(cepsb2w - Q**2)
  ckzb2w = real(ckzb2w) + ci *abs(aimag(ckzb2w)) !Im > 0
  !**
  ckzsw  = csqrt(cepssw - Q**2)
  ckzsw  = real(ckzsw) + ci * abs(aimag(ckzsw))   !Im > 0
  ckzs2w = csqrt(cepss2w - Q**2)
  ckzs2w = real(ckzs2w) + ci *abs(aimag(ckzs2w)) !Im > 0
  !** t_s(w) vacuum->surface s-transmition
  ctsvsw = 2.*qz/(qz + ckzsw)
  !** t_s(2w) vacuum->surface s-transmition
  ctsvs2w = 2.*qz/(qz + ckzs2w)
  !** t_p(w) vacuum->surface p-transmition
  ctpvsw = 2.*qz/(qz*cepssw + ckzsw)
  !** t_p(2w) vacuum->surface p-transmition
  ctpvs2w = 2.*qz/(qz*cepss2w + ckzs2w)
  !** t_s(w) surface->bulk s-transmition
  ctssbw = 2.*ckzsw/(ckzsw + ckzbw)
  !** t_s(2w) surface->bulk s-transmition
  ctssb2w = 2.*ckzs2w/(ckzs2w + ckzb2w)
  !** t_p(w) surface->bulk p-transmition
  ctpsbw = 2.*ckzsw/(cepsbw*ckzsw + cepssw*ckzbw)
  !** t_p(2w) surface->bulk p-transmition
  ctpsb2w = 2.*ckzs2w/(cepsb2w*ckzs2w + cepss2w*ckzb2w)
  !**
  !**  for chi in esu.cm=cm^3/e
  ! cesu = esu*cz / w(iw)**3
  ! cXzzz(iw)= cesu * cXzzz(iw)
  ! cXzxx(iw)= cesu * cXzxx(iw)
  ! cXxxz(iw)= cesu * cXxxz(iw)
  ! cXxxx(iw)= cesu * cXxxx(iw)
  !** r_{ij}
  ! ** p -> P
  cradpP = Q*cepsb2w*(  Q**2 * cepsbw**2 * cXzzz(iw) &
                      + ckzbw**2 *cepssw**2 * cXzxx(iw) ) &
                      + cepssw*cepss2w*ckzbw*ckzb2w &
                    * ( - 2.*Q*cepsbw* cXxxz(iw) &
                      + ckzbw * cepssw * cXxxx(iw) * cos(3.*phi) ) 
  !** p -> S
  cradpS = -ckzbw**2 *cepssw**2 * cXxxx(iw) * sin(3.*phi)
  !** s -> P
  cradsP = Q * cepsb2w * cXzxx(iw) - cepss2w * ckzb2w * cXxxx(iw) * cos(3.*phi)
  !** s -> S
  cradsS =  cXxxx(iw) * sin(3.*phi)
  !**
  !** bare \chi_{ijk} in 10^{-13} esu.cm
  fes = 1.e13
  fes = 1
  czzzb = fes*cXzzz(iw)
  czxxb = fes*cXzxx(iw)
  cxxzb = fes*cXxxz(iw)
  cxxxb = fes*cXxxx(iw)
  ! **
  ! czzzb = cpf*cXzzz(iw)/w(iw)**3
  ! czxxb = cpf*cXzxx(iw)/w(iw)**3
  ! cxxzb = cpf*cXxxz(iw)/w(iw)**3
  ! cxxxb = cpf*cXxxx(iw)/w(iw)**3
  ! ** with fresnel
  ! ** Bulk
  czzzfb = Q**2*czzzb
  czxxfb = ckzbw**2*czxxb
  cxxzfb = 2.*ckzbw*ckzb2w*cxxzb
  ! ** 3-layer
  czzzf = ctpsb2w*ctpsbw**2*cepsb2w* Q**2*cepsbw**2*czzzb
  czxxf = ctpsb2w*ctpsbw**2*cepsb2w* ckzbw**2 *cepssw**2 * czxxb
  cxxzf = 2.*ctpsb2w*ctpsbw**2*cepssw*cepss2w*ckzbw*ckzb2w*cepsbw*cxxzb
  ! ** R_pP with each X_ijk
  cradppzzz = Q*cepsb2w* Q**2*cepsbw**2*cXzzz(iw)
  cradppzxx = Q*cepsb2w* ckzbw**2 *cepssw**2 * cXzxx(iw)
  cradppxxz = -2.*cepssw*cepss2w*ckzbw*ckzb2w*Q*cepsbw*cXxxz(iw)
  !**
  !** radiated efficiency in 10^{-21} cm^2/watt
  !** For us {\cal R}_{out,in}, usually {\cal R}_{in,out}, However there
  !** is a much better way of knowing which polarization is 'in-w photon'
  !** or 'out-2w photon', by usin small letters for 'in' and capital
  !** letters for 'out'. Thus our {\cal R}_{ij} -> {\cal R}_{Ij}, i.e.
  !** {\cal R}_{Pp}, {\cal R}_{Ps}, {\cal R}_{Sp} and {\cal R}_{Ss}, equal to
  !** {\cal R}_{pP}, {\cal R}_{sP}, {\cal R}_{pS} and {\cal R}_{sS}
  ! **
  ctpP = ctpvs2w*ctpsb2w*(ctpvsw*ctpsbw)**2/(2.*ci)
  ctsP = ctpvs2w*ctpsb2w*(ctsvsw*ctssbw)**2/(2.*ci)
  ctpS = ctsvs2w*ctssb2w*(ctpvsw*ctpsbw)**2/(2.*ci)
  ctsS = ctsvs2w*ctssb2w*(ctsvsw*ctssbw)**2/(2.*ci)
  ! ** for chi in esu.cm
  radpP(iw) = pf*w(iw)**2*cabs(ctpP*cradpP)**2
  radsP(iw) = pf*w(iw)**2*cabs(ctsP*cradsP)**2
  radpS(iw) = pf*w(iw)**2*cabs(ctpS*cradpS)**2
  radsS(iw) = pf*w(iw)**2*cabs(ctsS*cradsS)**2
  ! **
  ! radpP(iw) = pf*cabs(ctpP*cradpP)**2/w(iw)**4
  ! radsP(iw) = pf*cabs(ctsP*cradsP)**2/w(iw)**4
  ! radpS(iw) = pf*cabs(ctpS*cradpS)**2/w(iw)**4
  ! radsS(iw) = pf*cabs(ctsS*cradsS)**2/w(iw)**4
  ! ** R_pP with each X_ijk
  radpPzzz = pf*cabs(ctpP*cradpPzzz)**2*w(iw)**2
  radpPzxx = pf*cabs(ctpP*cradpPzxx)**2*w(iw)**2
  radpPxxz = pf*cabs(ctpP*cradpPxxz)**2*w(iw)**2
  !**  shifted energy
  if(itb.eq.2)ws= 2.*w(iw)
  write(36,*)ws,radpP(iw),radpS(iw),radsP(iw),radsS(iw)
  write(301,*)w(iw),cabs(cepsbw),cabs(cepssw),cabs(cepsb2w),cabs(cepss2w)
  write(302,*)w(iw),cabs(ckzbw),cabs(ckzb2w),cabs(ckzsw),cabs(ckzs2w)
  write(303,*)w(iw),cabs(ctsvsw),cabs(ctpvsw),cabs(ctssbw),cabs(ctpsbw)
  write(304,*)w(iw),cabs(ctsvs2w),cabs(ctpvs2w),cabs(ctssb2w),cabs(ctpsb2w)
  write(305,*)w(iw),cabs(cradpP),cabs(cradpS),cabs(cradsP),cabs(cradsS)
  ! write(38,*)ws,radpPzzz,radpPzxx,radpPxxz
  ! write(91,901)ws,cabs(czzzb),cabs(czxxb),cabs(cxxzb),cabs(cxxxb)
  ! write(92,901)w(iw),cXzzz(iw),cXzxx(iw),cXxxz(iw)
  ! write(93,901)ws,cabs(czzzfb),cabs(czxxfb),cabs(cxxzfb)
  ! write(93,901)w(iw),czzzfb,czxxfb,cxxzfb
  ! write(94,901)w(iw),cradppzzz,cradppzxx,cradppxxz
  ! write(95,901)w(iw),ctpP*cradppzzz,ctpP*cradppzxx,ctpP*cradppxxz
  ! write(96,901)w(iw),cradppzzz,cradppzxx,cradppxxz
  ! write(94,901)w(iw),ctpP,cradppzzz,cradppzxx,cradppxxz
  ! write(94,901)w(iw),cepsb(iw),cepss(iw)
  ! write(92,901)w(iw),cXzzz(iw),cXzxx(iw),cXxxz(iw),cXxxx(iw)
  ! write(92,901)ws,cabs(czzzf),cabs(czxxf),cabs(cxxzf)
  ! write(93,901)ws,tpPzzz,tpPzxx,tpPxxz,tsP
  !***                       Pp        Ps=sP     Sp=pS     Ss
  !     write(36,901)ang,ws,radpp(iw),radps(iw),radsp(iw),radss(iw)
  !     write(91,901)ang,ws,rzzz,rzxx,rxxz
  !**
  !      end if
  !**
end do !iw
!*****
!*****
!     write(36,902)
!     write(91,902)
!       end do !it
!*****
901    format(f8.3,10e13.5)
902    format()
   end
!*****
