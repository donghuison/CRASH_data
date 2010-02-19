filename='hyades2d_1.3ns.out
transform='n
.r getpict
logrho= alog10(w(*,*,0)*1000.0)
logTe = alog10(w(*,*,2)*1.16e7)
iXe   = where(w(*,*,7) eq 0)
iBe   = where(w(*,*,7) eq 1)
iPl   = where(w(*,*,7) eq 2)
iAu   = where(w(*,*,7) eq 3)

transform='r'
nxreg=[256,256]
plotmode='contbarlog
fixaspect=0

set_device,'eos_crash_sesame_xe.eps',/eps,/land
filename='eos_crash.out eos_sesame_Xe.out'
xreglimits=[0.5,3.0,3.0,6.5]
func='{pXe}>1e4 {EintXe}>1e4'
.r animate
oplot,logrho(iXe),logTe(iXe),psym=1
close_device,/pdf

set_device,'eos_crash_sesame_be.eps',/eps,/land
filename='eos_crash.out eos_sesame_Be.out'
xreglimits=[2.0,2.5,5.0,7.5]
func='{pBe}>1e4 {EintBe}>1e4'
.r animate
oplot,logrho(iBe),logTe(iBe),psym=1
close_device,/pdf

set_device,'eos_crash_sesame_pl.eps',/eps,/land
filename='eos_crash.out eos_sesame_Pl.out'
xreglimits=[1.0,2.5,3.2,5.5]
func='{pPl}>1e4 {EintPl}>1e4'
.r animate
oplot,logrho(iPl),logTe(iPl),psym=1
close_device,/pdf

set_device,'eos_crash_sesame_au.eps',/eps,/land
filename='eos_crash.out eos_sesame_Au.out'
xreglimits=[2.0,2.5,5.0,7.5]
func='{pAu}>1e4 {EintAu}>1e4'
.r animate
oplot,logrho(iAu),logTe(iAu),psym=1
close_device,/pdf
