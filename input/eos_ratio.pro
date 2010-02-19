filename='hyades2d_1.3ns.out
transform='n
.r getpict
logrho= alog10(w(*,*,0)*1000.0)
logTe = alog10(w(*,*,2)*1.16e7)
iXe   = where(w(*,*,7) eq 0)
iBe   = where(w(*,*,7) eq 1)
iPl   = where(w(*,*,7) eq 2)
iAu   = where(w(*,*,7) eq 3)

transform='r'   ; transform when reading files
nxreg=[256,256]
plotmode='contbar
dotransform='n' ; do not transform .r plotfunc is run for the ratio

set_device,'eos_ratio_xe.eps',/eps,/land
filename='eos_crash.out eos_sesame_Xe.out'
xreglimits=[0.5,3.0,3.0,6.5]
.r getpict
wreg = wreg1/wreg0(*,*,0:1)
func='EintXe pXe'
plottitle='EintXeSESAME/EintXeCRASH;pXeSESAME/pXeCRASH
.r plotfunc
oplot,logrho(iXe),logTe(iXe),psym=1
close_device,/pdf

set_device,'eos_ratio_be.eps',/eps,/land
filename='eos_crash.out eos_sesame_Be.out'
xreglimits=[2.0,2.5,5.0,7.5]
.r getpict
wreg = wreg1/wreg0(*,*,2:3)
autorange='n
plotmode='contbar
fmin=[0.2,0.2]
fmax=[5.0,5.0]
func='{EintBe}>0.2 {pBe}>0.2'
plottitle='EintBeSESAME/EintBeCRASH;pBeSESAME/pBeCRASH
.r plotfunc
oplot,logrho(iBe),logTe(iBe),psym=1
close_device,/pdf

set_device,'eos_ratio_pl.eps',/eps,/land
filename='eos_crash.out eos_sesame_Pl.out'
xreglimits=[1.0,2.5,3.2,5.5]
.r getpict
wreg = wreg1/wreg0(*,*,4:5)
autorange='n
plotmode='contbar
fmin=[0.2,0.2]
fmax=[6.0,6.0]
func='{EintPl} {pPl}'
plottitle='EintPlSESAME/EintPlCRASH;pPlSESAME/pPlCRASH
.r plotfunc
oplot,logrho(iPl),logTe(iPl),psym=1
close_device,/pdf

set_device,'eos_ratio_au.eps',/eps,/land
filename='eos_crash.out eos_sesame_Au.out'
xreglimits=[2.0,2.5,5.0,7.5]
.r getpict
wreg = wreg1/wreg0(*,*,6:7)
autorange='n
plotmode='contbar
fmin=[0.2,0.2]
fmax=[6.0,6.0]
func='{EintAu} {pAu}'
plottitle='EintAuSESAME/EintAuCRASH;pAuSESAME/pAuCRASH
.r plotfunc
oplot,logrho(iAu),logTe(iAu),psym=1
close_device,/pdf
