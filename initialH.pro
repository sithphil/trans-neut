pro initialH
;; analytical solution for initial f(x)=x, Boundary condition f(x=0)=0,f(x=p)=0

safe_colors, /first
path="data"
;g=file_import("data/circle.grd.hl2a.nc")
;g=file_import("circle.grd.nc")
g=file_import("data/circular.grd.hl2a22535.831ms.68x256.0.65to1.1.nc")
;g=file_import("slab.grd.nc")
;g5=file_import("cbm18_dens5.grid.nc")
;g6=file_import("cbm18_dens6.grid.nc")
g5=g
g6=g

density_unit = 1.e19                          ; 1/m^3
ee = 1.6022e-19                               ; elementary charge 
Ni_x = 1
Te_x = 10     ; eV


NX=g.NX
;NY=g.NY
;NX=68

lbar=g.rmag
yy=fltarr(NX)
ni=fltarr(NX)
ti=fltarr(NX)
te=fltarr(NX)
pei=fltarr(NX)

Diff=fltarr(NX)
chii=fltarr(NX)
chie=fltarr(NX)

unit_psi=g.rmag*g.rmag*g.bmag
jy=30

rxy=g.Rxy[*,jy]/lbar
bpxy=g.Bpxy[*,jy]/g.bmag

x=(g.psixy[*,jy]-g.psi_axis)/(g.psi_bndry-g.psi_axis)
;xreal=g.psixy[*,jy]/unit_psi
;x=(g.rxy[*,jy]-min(g.rxy))/(max(g.rxy)-min(g.rxy))
xreal=g.rxy[*,jy]/lbar

x5=(g5.psixy[*,jy]-g5.psi_axis)/(g5.psi_bndry-g5.psi_axis)
x5r=g5.psixy[*,jy]/(g5.rmag*g5.rmag*g5.bmag)

p5=g5.pressure[*,jy]

x6=(g6.psixy[*,jy]-g6.psi_axis)/(g6.psi_bndry-g6.psi_axis)
x6r=g6.psixy[*,jy]/(g6.rmag*g6.rmag*g6.bmag)

p6=g6.pressure[*,jy]

dt=100


w_nn=0.1

p_core = 950.     ; pmin PW Xi

delta_ped=0.062 ; width of pedestal
p_a = 0.01      ; parameter for gradient at pedestal top
p_del=18.       ; parameter for enlarge of core value
x0_p = 0.90     ; center of pedestal relative to normalzied psi

p_bottom = p_core -1.
p_edge = 2.      ; pmin - pbottom in PW Xi

;*** smaller than dens6 H3

;ni_core = 0.065
;ni_edge = 0.048
;ti_core = 4.5
;ti_edge = 0.38
;te_core = 4.5
;te_edge = 0.38

;*** as dens5   H4
;ni_core = 0.075
;ni_edge = 0.048
;ti_core = 4.5
;ti_edge = 0.38
;te_core = 4.5
;te_edge = 0.38
  

;*** smaller than dens6 H5
;ni_core = 0.08
;ni_edge = 0.048

;ti_core = 4.5
;ti_edge = 0.38     

;te_core = 4.5
;te_edge = 0.38 

;*** relative to dens6  H6
ni_core = 0.1
ni_edge = 0.047

ti_core = 5.
ti_edge = 0.37     

te_core = 5.
te_edge = 0.37

;*** upper than dens6
;ni_core = 0.14
;ni_edge = 0.046

;ti_core = 6.
;ti_edge = 0.35     

;te_core = 6.
;te_edge = 0.35  

;*** much higher than dens6
;ni_core = 0.2             ; unit Ni_x*p_del, normalized
;ni_edge = 0.045

;ti_core = 7.              ; unit Te_x*p_del
;ti_edge = 0.32     

;te_core = 7.              ; unit Te_x*p_del
;te_edge = 0.32


Diff_xin = 1.       ; m^2/s
;Grad_xin = -7000.     ; d(var)/dpsi, unit(var)/(m^2 tessla )
;flux_xin = - Diff_xin*rxy[0]*bpxy[0]*Grad_xin 
;flux_xin = - Diff_xin*Grad_xin 

for jx = 0,NX-1  do begin
     

        xprm = (x0_p - x[jx])/delta_ped
        
        Hmode_tanh =(exp(-xprm)+(1.+p_a*xprm)*exp(xprm)*p_del)/(exp(-xprm)+exp(xprm))
        yy[jx]= p_edge +p_core* (Hmode_tanh-1.) 
        ;print,jx,x[jx],xprm,yy[jx]

        ni[jx]= ni_edge +ni_core* (Hmode_tanh-1.) 
        ti[jx]= ti_edge +ti_core* (Hmode_tanh-1.) 
        te[jx]= te_edge +te_core* (Hmode_tanh-1.) 

  endfor

pei=ni*(ti+te)
pei = pei*Ni_x*density_unit*ee*Te_x             ; SI unit Pascals

print,'Var', '    Min',    '      Max'
print,'ni : ', min(ni),'---', max(ni)

print,'ti : ', min(ti),'---', max(ti)
print,'te : ', min(te),'---', max(te)
print,'pei : ', min(pei),'---', max(pei)

dy5dx=deriv(x5r,p5)
dy6dx=deriv(x6r,p6)

dyydx=deriv(xreal,pei)
;print,dyydx

dnidx=deriv(xreal,ni)
dtidx=deriv(xreal,ti)
dtedx=deriv(xreal,te)
;print,dnidx

print,'Gradient', '    Core',    '      Edge'
print,'ni : ', dnidx[0:1],'---', dnidx[NX-2:*]

print,'ti : ', dtidx[0:1],'---', dtidx[NX-2:*]
print,'te : ', dtedx[0:1],'---', dtedx[NX-2:*]



flux_xin_ni = -Diff_xin*dnidx[0]*rxy[0]*bpxy[0]
flux_xin_ti = -Diff_xin*dtidx[0]*rxy[0]*bpxy[0]
flux_xin_te = -Diff_xin*dtedx[0]*rxy[0]*bpxy[0]

for jx = 1,NX-1  do begin

Diff[jx] = -flux_xin_ni/rxy[jx]/bpxy[jx]/dnidx[jx]
chii[jx] = -flux_xin_ti/rxy[jx]/bpxy[jx]/dtidx[jx]
chie[jx] = -flux_xin_te/rxy[jx]/bpxy[jx]/dtedx[jx]


;if Diff[jx] GT 10. then BEGIN
 ;Diff[jx]=10.

;ENDIF


;Diff[jx] = -flux_xin/dyydx[jx]

ENDFOR

Diff[0]=Diff_xin
chii[0]=Diff_xin
chie[0]=Diff_xin

print,'Diffusion', '    Min',    '      Max'
print,'D : ', min(Diff),'---', max(Diff)
print,'chi_i : ', min(chii),'---', max(chii)
print,'chi_e : ', min(chie),'---', max(chie)

ddxdiff=deriv(xreal,Diff)
ddxchii=deriv(xreal,chii)
ddxchie=deriv(xreal,chie)

print,'Gradient', '    Core',    '      Edge'
print,'Diff_ni : ', ddxdiff[0],'---', ddxdiff[NX-1]

print,'chii : ', ddxchii[0],'---', ddxchii[NX-1]
print,'chie : ', ddxchie[0],'---', ddxchie[NX-1]



window,0
plot,x,dyydx,/xst, title='dyydx',chars=2,charthick=2,thick=2
oplot,x5,dy5dx,color=2
oplot,x6,dy6dx,color=4

window,10
plot,x,pei,/xst, title='P[pa]',chars=2,charthick=2,thick=2


window,1
plot,x,yy,/xst,chars=2,charthick=2,thick=2,title='P fitting'
oplot,x6,g6.pressure[*,jy],color=2,thick=2

window,2
plot,x,Diff,/xst,title='Diffusion Coef',chars=2,charthick=2,thick=2

window,12
plot,x,chii,/xst,title='chii Coef',chars=2,charthick=2,thick=2

window,13
plot,x,chie,/xst,title='chie Coef',chars=2,charthick=2,thick=2


window,3
plot,x,ni,/xst,chars=2,charthick=2,thick=2,title='Ni'

window,4
plot,x,ti,/xst,chars=2,charthick=2,thick=2,title='Ti'

window,5
plot,x,te,/xst,chars=2,charthick=2,thick=2,title='Te'




end


