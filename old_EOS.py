import math as m
import numpy as nu
from scipy import optimize
from scipy import interpolate
from scipy.interpolate import griddata
from scipy import constants as cst
import seafreeze as sf
import warnings
warnings.filterwarnings("ignore", category=nu.VisibleDeprecationWarning)

Tc=647.096 #K, temperature at the critical point of water
pc=22.064e6 #Pa, pressure at the critical point of water
pt=611.655 #Pa, pressure at the triple point of water
rhoc=322./18.015268e-3 #mol.m-3
Rig=8.3143714 #ideal gas constant
M_h2o=18.015268e-3 #kg.mol-1


domega=50

#################################################################
#################################################################
#################################################################
################################################################# PRETRAITMENT FOR THE USE OF MAZEVET EOS:

eos2=nu.loadtxt("./dataIN/eosMazevet.dat") ### this file should be a SQUARE GRID IN T-rho

#eos2 goes this way:
#rho[g/cc], T [K], P [GPa], P/(n_i kT), F/(N_i kT), U/(N_i kT), C_V/(N_i k), chi_T, chi_r, S/(N_i k)
s_offset=-54.280936454849495 #J.mol-1.K-1
Nlgn=len(eos2)


#Grid 1 = square grid in T, rho
Rhog1=eos2[:,0]
pg1=eos2[:,2]
Tg1=eos2[:,1]
sg1=eos2[:,9]
ug1=eos2[:,5]
chiR1=eos2[:,8]


#let's do the interpolation on the p, T rectangular grid
Ng2=20 #1 side size of the new grid
logpg2=nu.linspace(m.log10(1e-1),m.log10(1e4), Ng2 )#GPa, constant step in logspace
pg2=10**logpg2 #out of the logspace
Tg2=nu.linspace(m.log10(nu.min(Tg1+5)), m.log10(nu.max(Tg1-5)), Ng2 ) #grid in temperature
Tg2=10**Tg2

rhog2=nu.linspace(m.log10(nu.min(eos2[:,0])), m.log10(nu.max(eos2[:,0])), Ng2 ) # grid in density
rhog2=rhog2*1000./M_h2o #mol.m-3

tt1,pp1=nu.meshgrid(Tg2,pg2)


S_forRBS_Tp    = griddata( (Tg1,pg1),   sg1, (tt1,pp1) , method="linear")
U_forRBS_Tp    = griddata( (Tg1,pg1),   ug1, (tt1,pp1) , method="linear")
RHO_forRBS_Tp  = griddata( (Tg1,pg1), Rhog1, (tt1,pp1) , method="linear")
CHIR_forRBS_Tp = griddata( (Tg1,pg1), chiR1, (tt1,pp1) , method="linear")
V_forRBS_Tp    = griddata( (Tg1,pg1), 1./(Rhog1*1000./18.015268e-3), (tt1,pp1) , method="linear", fill_value=0.)



P_forRBS_Trho    = nu.zeros( (int(Nlgn**0.5) , int(Nlgn**0.5)) )
S_forRBS_Trho    = nu.zeros( (int(Nlgn**0.5) , int(Nlgn**0.5)) )
U_forRBS_Trho    = nu.zeros( (int(Nlgn**0.5) , int(Nlgn**0.5)) )
CHIR_forRBS_Trho = nu.zeros( (int(Nlgn**0.5) , int(Nlgn**0.5)) )


cpt=0
for i in range(0,int(Nlgn**0.5)):#should be T
    for j in range(0,int(Nlgn**0.5)): #should be rho
        P_forRBS_Trho[i,j]=eos2[cpt,2]*1e9
        S_forRBS_Trho[j,i]=eos2[cpt,9]*cst.k*3.*cst.N_A + s_offset
        U_forRBS_Trho[i,j]=eos2[cpt,5]*cst.k*3.*cst.N_A*eos2[cpt,1]
        cpt+=1
        
s_mazevet_Trho    =interpolate.RectBivariateSpline(rhog2, Tg2, S_forRBS_Tp )

#just to put in the relevant units:
for i in range(0,Ng2):
    for j in range(0,Ng2):
        S_forRBS_Tp[i,j]=S_forRBS_Tp[i,j]*cst.k*3.*cst.N_A + s_offset
        
        


#RBS require 2 1D vectors for x and y, and a "meshgridded" 2D vector for Z
s_mazevet_Tp    =interpolate.RectBivariateSpline(pg2*1e9, Tg2, S_forRBS_Tp )
u_mazevet_Tp    =interpolate.RectBivariateSpline(pg2*1e9, Tg2, U_forRBS_Tp*cst.k*3.*cst.N_A*tt1) # T in K, p in Pa, gives U in J.mol-1
rho_mazevet_Tp  =interpolate.RectBivariateSpline(pg2*1e9, Tg2, RHO_forRBS_Tp*1000./18.015268e-3 ) # T in K, p in Pa, gives RHO in mol.m-3
chiR_mazevet_Tp =interpolate.RectBivariateSpline(pg2*1e9, Tg2, CHIR_forRBS_Tp )# T in K, p in Pa, ChiR dimensionless



###############################################################
###############################################################

hpi=nu.loadtxt("./dataIN/dataHPices.dat")


NlgnHPI=len(hpi)


THPI=nu.zeros(NlgnHPI)
PHPI=nu.zeros(NlgnHPI)
rhoHPI=nu.zeros(NlgnHPI)
alphaHPI=nu.zeros(NlgnHPI)
CpHPI=nu.zeros(NlgnHPI)
sHPI=nu.zeros(NlgnHPI)

aa_HPI =     nu.zeros((int(NlgnHPI**0.5),int(NlgnHPI**0.5)))
ss_HPI =     nu.zeros((int(NlgnHPI**0.5),int(NlgnHPI**0.5)))
cc_HPI =     nu.zeros((int(NlgnHPI**0.5),int(NlgnHPI**0.5)))
rrho_HPI =   nu.zeros((int(NlgnHPI**0.5),int(NlgnHPI**0.5)))

tt =   nu.zeros((int(NlgnHPI**0.5),int(NlgnHPI**0.5)))
pp =   nu.zeros((int(NlgnHPI**0.5),int(NlgnHPI**0.5)))


for i in range(0,NlgnHPI):
    
    THPI[i] =     hpi[i,0]
    PHPI[i] =     hpi[i,1]
    rhoHPI[i] =   hpi[i,2]
    alphaHPI[i] = hpi[i,3]
    CpHPI[i] =    hpi[i,4]
    sHPI[i] =     hpi[i,5]
    
    


thpi = nu.linspace( nu.min(THPI), nu.max(THPI), int(NlgnHPI**0.5) )
phpi = nu.linspace( m.log10(nu.min(PHPI)), m.log10(nu.max(PHPI)), int(NlgnHPI**0.5) )
phpi=10**phpi

tt,pp=nu.meshgrid(thpi,phpi)#to double-check the indices


cpt=0
for i in range(0,int(NlgnHPI**0.5)):
    for j in range(0,int(NlgnHPI**0.5)):
        aa_HPI[i,j]=alphaHPI[cpt]
        ss_HPI[i,j]=sHPI[cpt]
        cc_HPI[i,j]=CpHPI[cpt]
        rrho_HPI[i,j]=rhoHPI[cpt]
        
        #print(i,j,cpt, )
    
        cpt+=1
        
s_HPI_Tp     =interpolate.RectBivariateSpline(phpi, thpi, ss_HPI )
alpha_HPI_Tp =interpolate.RectBivariateSpline(phpi, thpi, aa_HPI )
rho_HPI_Tp   =interpolate.RectBivariateSpline(phpi, thpi, rrho_HPI)
Cp_HPI_Tp    =interpolate.RectBivariateSpline(phpi, thpi, cc_HPI)


#####################################


datarholiq=nu.loadtxt("./dataIN/liqgrid_500.dat")
datarhovap=nu.loadtxt("./dataIN/vapgrid_500.dat")
datarhosc=nu.loadtxt("./dataIN/scgrid_500.dat")

tliq=nu.zeros(len(datarholiq))
plgrid=nu.zeros(len(datarholiq))
rholiq=nu.zeros(len(datarholiq))

tvap=nu.zeros(len(datarhovap))
pvgrid=nu.zeros(len(datarhovap))
rhovap=nu.zeros(len(datarhovap))

tsc=nu.zeros(len(datarhosc))
psc=nu.zeros(len(datarhosc))
rhosc=nu.zeros(len(datarhosc))







for i in range(0,len(datarhovap)):
    tvap[i]= m.log10(datarhovap[i,0])
    pvgrid[i]= datarhovap[i,1]
    rhovap[i]= m.log10(datarhovap[i,2])
    


Grmax=nu.max(pvgrid)
Grmin=nu.min(pvgrid)
Tmax=10**nu.max(tvap)
Tmin=10**nu.min(tvap)
tvaptest=nu.linspace( m.log10(Tmin),m.log10(Tmax),int(len(datarhovap)**0.5)  )
lgridtest=nu.linspace( Grmin,Grmax,int(len(datarhovap)**0.5) )
rhovaptest=nu.zeros( (int(len(datarhovap)**0.5),int(len(datarhovap)**0.5)) )

cpt=0
for i in range(0,int(len(datarhovap)**0.5) ):
    for j in range(0,int(len(datarhovap)**0.5) ):
        rhovaptest[i,j]=m.log10(datarhovap[cpt,2])
        cpt+=1
        

testvap_RBS    =interpolate.RectBivariateSpline( tvaptest,lgridtest, rhovaptest )



##--------------------------------------------------------------------------------------------
##---------------Now the liquid:
##--------------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------------

for i in range(0,len(datarholiq)):
    tliq[i]= m.log10(datarholiq[i,0])
    plgrid[i]=datarholiq[i,1]
    rholiq[i]= m.log10(datarholiq[i,2])


Grmax=nu.max(plgrid)
Grmin=nu.min(plgrid)
Tmax=10**nu.max(tliq)
Tmin=10**nu.min(tliq)

tliqtest=nu.linspace( m.log10(Tmin),m.log10(Tmax),int(len(datarholiq)**0.5)  )
plgridtest=nu.linspace( Grmin,Grmax,int(len(datarholiq)**0.5) )
rholiqtest=nu.zeros( (int(len(datarholiq)**0.5),int(len(datarholiq)**0.5)) )

cpt=0
for i in range(0,int(len(datarholiq)**0.5) ):
    for j in range(0,int(len(datarholiq)**0.5) ):
        rholiqtest[i,j]=m.log10(datarholiq[cpt,2])
        cpt+=1
        
testliq_RBS    =interpolate.RectBivariateSpline( tliqtest,plgridtest, rholiqtest )

##--------------------------------------------------------------------------------------------
##---------------And finally the supercritical grid:
##--------------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------------

for i in range(0,len(datarhosc)):
    tsc[i]= m.log10(datarhosc[i,0])
    psc[i]=m.log10(datarhosc[i,1])
    rhosc[i]= m.log10(datarhosc[i,2])


Pmax=10**nu.max(psc)
Pmin=10**nu.min(psc)
Tmax=10**nu.max(tsc)
Tmin=10**nu.min(tsc)

tsctest=nu.linspace( m.log10(Tmin),m.log10(Tmax),int(len(datarhosc)**0.5)  )
psctest=nu.linspace( m.log10(Pmin),m.log10(Pmax),int(len(datarhosc)**0.5) )
rhosctest=nu.zeros( (int(len(datarhosc)**0.5),int(len(datarhosc)**0.5)) )

cpt=0
for i in range(0,int(len(datarhosc)**0.5) ):
    for j in range(0,int(len(datarhosc)**0.5) ):
        rhosctest[i,j]=m.log10(datarhosc[cpt,2])
        cpt+=1
        

testsc_RBS    =interpolate.RectBivariateSpline( tsctest,psctest, rhosctest )


#####################################
dataSVP=nu.loadtxt("./dataIN/SVP.dat")

tsvp=nu.zeros(len(dataSVP))
psvp=nu.zeros(len(dataSVP))
rhosvp_v=nu.zeros(len(dataSVP))
rhosvp_l=nu.zeros(len(dataSVP))

for i in range(0,len(dataSVP)):
    tsvp[i]= m.log10(dataSVP[i,0])
    psvp[i]=m.log10(dataSVP[i,1])
    rhosvp_v[i]= m.log10(dataSVP[i,2])
    rhosvp_l[i]= m.log10(dataSVP[i,3])


f_SVP_v=interpolate.interp1d(tsvp,rhosvp_v)
f_SVP_l=interpolate.interp1d(tsvp,rhosvp_l)

#############################################################
#############################################################
#############################################################







def w(T,rho,p): #transition function from Mazevet to IAPWS in the fluid phase

    delT=domega
    T0=1273.
    omega=1.
    
    if T<=T0-delT and rho>1000./M_h2o:
        omega=1.
    elif T>T0-delT and T<=T0+delT and rho>1000./M_h2o :
        a=-1./(2.*delT)
        b=0.5+T0/(2.*delT)
        omega=a*T+b
    elif T>T0+delT :
        omega=0.
        
    if (rho*M_h2o/1000.<=1.):
        omega=1.
        
    if (p>HPice(T)):
        omega=0.
  
  
    return omega
    
    
    
    


#triple point ice VI, VII, Liquid
pVIVIIL=2.2e9
TVIVIIL=353.5 #K
    
    
    
    
    
    
    #################################################################
    #################################################################
    #################################################################
    #################################################################
    

def WEOS(T,ps,phase):
    #print("In WEOS",T,ps,phase)
    
    ps_tot=ps
    
    #Ok. First let's determine where we are
    #Options:
    #numbers 1 > 6 are phases of ice from seafreeze


    whereIam=-1

   
    
    if (T<500. and T>190. and ps<pVIVIIL and ps>Psat(T) ) and (phase!="v_svp" and phase!="l_svp" and phase!="v" and phase!="vap"):
        #print("INSEAFREEZE", ps,Psat(T),T)
        
        PT = nu.empty((1,), nu.object)
        PT[0] = (ps/1e6,T)
        whereIam=sf.whichphase(PT)[0] #seafreeze business
        
                
        if (whereIam==1):
            material="Ih"
        elif (whereIam==2):
            material="II"
        elif (whereIam==3):
            material="III"
        elif (whereIam==5 ):
            material="V"
        elif (whereIam==6 ):
            material="VI"
        elif (whereIam==0):
            material="water1"#_IAPWS95"
        else:
            print("Problem picking an ice for", PT, whereIam)
            exit()
                
                
        if (material=="V") and ps>1e9 :
            material="VI"
                
        out=sf.seafreeze(PT,material)
            
        rhoi_tot = out.rho[0]/M_h2o
        cp_tot   = out.Cp[0]*M_h2o
        s_tot    = out.S[0]*M_h2o
        u_tot    = out.U[0]*M_h2o
                
                
                #derivatives:
                
        PT2 = nu.empty((1,), nu.object)
        PT2[0] = (ps/1e6+1e-5,T)
        outdsdp=sf.seafreeze(PT2,material)
                
                
        dsdp_tot = out.alpha[0]*1./rhoi_tot*(-1.) #(outdsdp.S[0]-out.S[0])/1e-5*M_h2o #.. J.kg-1.K-1.   kg.mol-1
             
        dpdrho_tot = 0. #
            

            
            
            
            
        
            
            
    else: #make the distiction between IAPWS, IAPWS>Mazevet transition region, Mazevet and LTHP mazevet, and the SVP line
        #print(T,ps,phase,T<1273., ps>pVIVIIL, ps>HPice(T))
        if ( T<1273. and ps>pVIVIIL and ps>=HPice(T)) :
            whereIam=999 #LTHP situation
                
                # then compute HPI stuff
            rhoi_tot =    rho_HPI_Tp(ps,T)[0,0]/M_h2o
            
            cp_tot   =    s_HPI_Tp(ps,T,dy=1)[0,0] *T #Cp_HPI_Tp(ps,T)[0,0]
            s_tot    =    s_HPI_Tp(ps,T)[0,0]
            dsdp_tot =    s_HPI_Tp(ps,T,dx=1)[0,0]#-1.*alpha_HPI_Tp(ps,T)[0,0]/rhoi_tot
            dpdrho_tot =  0.#1./(rho_HPI_Tp(ps,T,dx=1.)[0,0])*M_h2o

            u_tot    = -1027.45143885 #value of U at the triple point
                
            
        elif ((T<1273.-domega and T>500 and ps<=HPice(T)) or (T<500. and ps<Psat(T)) or (T<500. and ps<HPice(T) and ps>pVIVIIL) ): #doesn't matter, we deal with IAPWS in any case
            #(HP ices at LT have normally been taken care of in the parent if)
                whereIam=100 #IAPWS
               
            
        elif (T>=1273.-domega): 
                #need to compute the density from Mazevet:
            rhoi_eos2 = rho_mazevet_Tp(ps,T)[0,0] #mol.m-3
            
            
            if (rhoi_eos2 < (1000./M_h2o) ):
                whereIam=100 #IAPWS
               
                
            elif ( rhoi_eos2 >= (1000./M_h2o) and T<=1273.+domega and ps<=HPice(T)):
                #transition between the eos:
                whereIam=150 #need to compute both eos
                
            elif (rhoi_eos2 >= (1000./M_h2o) and T>1273.+domega or ps>HPice(T) ):
                whereIam=200 #Mazevet only
                
            else:
                print("Parameters not accounted for:")
                print(T,ps,rhoi_eos2*M_h2o)
                exit()
    
    
    
   
        if phase=="v_svp" or phase=="l_svp" :
            whereIam = 120
    
        if ps==Psat(T) and phase=="l":
            whereIam = 120
            phase="l_svp"
        elif ps==Psat(T) and phase=="v":
            whereIam = 120
            phase="v_svp"
            
        if phase=="v" and whereIam!=200 and whereIam!=150 and whereIam!=999:
            whereIam=100
        
     
        #let's compute stuff based on what we assigned previously



        #print("whereIam",whereIam,"\n")

       
        
        if (whereIam == 100 or whereIam==120):

          
            if whereIam==100 :
            
                if  (T<Tc and ps<Psat(T)) or (T<Tc and phase=="vap") :
                    Pmax=Psat(T)+Psat(T)*0.005
                    Pmin=1e-2
                    A=12./(m.log10(Pmax)-m.log10(Pmin))
                    pput=A*m.log10(ps)-A*m.log10(Pmax)
                    
                    rhoi_tot=testvap_RBS(m.log10(T),pput)[0,0]
                    rhoi_tot=10**rhoi_tot
                    whereIam=102
                    
                elif T<Tc and ps>Psat(T) and phase!="vap":
                    Pmin=Psat(T)-Psat(T)*0.005
                    Pmax=1e10
                        
                    A=12./(m.log10(Pmax)-m.log10(Pmin))
                    pput=A*m.log10(ps)-A*m.log10(Pmax)
                    
                    rhoi_tot=testliq_RBS(m.log10(T),pput)[0,0]
                    rhoi_tot=10**rhoi_tot
                    whereIam=103
                    
                elif T>=Tc:
                    rhoi_tot=testsc_RBS(m.log10(T),m.log10(ps))[0,0]
                    rhoi_tot=10**rhoi_tot
                    whereIam=104
                 
                    
                
            elif whereIam==120 :
                if phase=="v_svp":
                    rhoi_tot=10**f_SVP_v(m.log10(T))
                    whereIam=121
                elif phase=="l_svp" :
                    rhoi_tot=10**f_SVP_l(m.log10(T))
                    whereIam=122
                    
                    
            #print(T,ps,phase)
            delt=rhoi_tot/rhoc
            tau=Tc/T
                    
                
            phi0,dphi0t,dphi0tt=CALC_phi01(delt,tau)
            phir,dphird,dphirt,dphirtt,dphirdt,dphirdd=CALC_phir1(delt,tau)
                    
            cp_tot = (-1.*tau**2.*(dphi0tt+dphirtt) + (1.+delt*dphird-delt*tau*dphirdt)**2. / (1. + 2.*delt*dphird + delt**2.*dphirdd))*Rig
                
                
            dhdp_T= ( 1.-(1.+ delt*dphird - delt*tau*dphirdt)/(1.+2.*delt*dphird+delt**2*dphirdd) )*1./rhoi_tot
            dsdp_tot=(1./(T*rhoi_tot) - 1./T*dhdp_T)*(-1.)
                
                
            dpdrho_tot=0. #(1. + 2.*delt*dphird + delt**2*dphirdd - (1.+delt*dphird-delt*tau*dphirdt)**2./(tau**2*(dphi0tt+dphirtt)))*Rig*T
            s_tot = (tau*(dphi0t+dphirt)-phi0-phir)*Rig
            u_tot = (tau*(dphi0t+dphirt))*Rig*T
            
            
            
        elif (whereIam == 150):
        #IAPWS and mazevet transition
            if  T<Tc and ps<Psat(T)  :
                Pmax=Psat(T)+Psat(T)*0.005
                Pmin=1e-2
                A=12./(m.log10(Pmax)-m.log10(Pmin))
                pput=A*m.log10(ps)-A*m.log10(Pmax)
            
                rhoi=testvap_RBS(m.log10(T),pput)[0,0] #interpolate.bisplev(m.log10(T),m.log10(ps),tck_vap)
                rhoi=10**rhoi
                whereIam=152
            
            elif T<Tc and ps>Psat(T) :
                Pmin=Psat(T)-Psat(T)*0.005
                Pmax=1e10
                
                A=12./(m.log10(Pmax)-m.log10(Pmin))
                pput=A*m.log10(ps)-A*m.log10(Pmax)
            
                rhoi=testliq_RBS(m.log10(T),pput)[0,0] #interpolate.bisplev(m.log10(T),m.log10(ps),tck_liq)
                rhoi=10**rhoi
                whereIam=153
            
            elif T>=Tc:
                rhoi=testsc_RBS(m.log10(T),m.log10(ps))[0,0] #interpolate.bisplev(m.log10(T),m.log10(ps),tck_liq)
                rhoi=10**rhoi
                whereIam=154
        
        
                
            delt=rhoi/rhoc
            tau=Tc/T
                    
                
            phi0,dphi0t,dphi0tt=CALC_phi01(delt,tau)
            phir,dphird,dphirt,dphirtt,dphirdt,dphirdd=CALC_phir1(delt,tau)
                    
            cp = (-1.*tau**2.*(dphi0tt+dphirtt) + (1.+delt*dphird-delt*tau*dphirdt)**2. / (1. + 2.*delt*dphird + delt**2.*dphirdd))*Rig
                
                
            dhdp_T= ( 1.-(1.+ delt*dphird - delt*tau*dphirdt)/(1.+2.*delt*dphird+delt**2*dphirdd) )*1./rhoi
            dsdp=(1./(T*rhoi) - 1./T*dhdp_T)*(-1.)
                
                
            dpdrho=0.#(1. + 2.*delt*dphird + delt**2*dphirdd - (1.+delt*dphird-delt*tau*dphirdt)**2./(tau**2*(dphi0tt+dphirtt)))*Rig*T
            s = (tau*(dphi0t+dphirt)-phi0-phir)*Rig
            u = (tau*(dphi0t+dphirt))*Rig*T
            
            ##########################################--- EOS #2 i.e. Mazevet
                
                 
            rhoi_eos2 = rho_mazevet_Tp(ps,T)
            s_eos2    = s_mazevet_Tp  (ps,T)
            u_eos2    = u_mazevet_Tp  (ps,T)
            cp_eos2   = T*s_mazevet_Tp(ps,T,dx=0,dy=1)
                 
            ps_tot=ps #user input
            ps_eos2=ps
                
             
            ##########################################--COMBINE THE TWO EOS:
              
            rhoi_tot = w(T,rhoi,ps)*rhoi+(1.-w(T,rhoi,ps))*rhoi_eos2[0,0]
            dpdrho_eos2=0.#chiR_mazevet_Tp(ps,T)*ps/rhoi_tot
             
            #############################################################
            #############################################################
            #############################################################
            
            

            s_tot=w(T,rhoi,ps_tot)*s+(1.-w(T,rhoi,ps_tot))*s_eos2[0,0]
            u_tot=w(T,rhoi,ps_tot)*u+(1.-w(T,rhoi,ps_tot))*u_eos2[0,0]
            cp_tot=w(T,rhoi,ps_tot)*cp+(1.-w(T,rhoi,ps_tot))*cp_eos2[0,0]
             
            dsdp_eos2=s_mazevet_Tp(ps_tot,T,dx=1,dy=0)
             
            dsdp_tot=w(T,rhoi,ps_tot)*dsdp+(1.-w(T,rhoi,ps_tot))*dsdp_eos2[0,0]
            dpdrho_tot=0.#w(T,rhoi,ps_tot,phase)*dpdrho+(1.-w(T,rhoi,ps_tot,phase))*dpdrho_eos2[0,0]
             
             
             
             
        
        elif (whereIam == 200):
        #Mazevet

                
            rhoi_tot = rho_mazevet_Tp(ps,T)[0,0]
            cp_tot  = T*s_mazevet_Tp(ps,T,dx=0,dy=1)[0,0]
            s_tot    = s_mazevet_Tp  (ps,T)[0,0]
            u_tot    = u_mazevet_Tp  (ps,T)[0,0]
            dsdp_tot=s_mazevet_Tp(ps_tot,T,dx=1,dy=0)[0,0]
           
            dpdrho_tot=0.#chiR_mazevet_Tp(ps,T)*ps/rhoi_tot
            
            
            
    
    
    return  ps_tot, rhoi_tot, cp_tot, s_tot, u_tot, dsdp_tot, dpdrho_tot, whereIam








#################################################################
#################################################################
#################################################################
#################################################################

def search_p(rhoi,T,ps):

    delt=rhoi/rhoc
    tau=Tc/T
    phir,dphird,dphirt,dphirtt,dphirdt,dphirdd=CALC_phir1(delt,tau)

    return ((1. + delt*dphird)*rhoi*Rig*T - ps)/ps # everything in Pa


#################################################################
#################################################################
#################################################################
#################################################################

def pressure(rhoi,T):


    delt=rhoi/rhoc
    tau=Tc/T

    phir,dphird,dphirt,dphirtt,dphirdt,dphirdd=CALC_phir1(delt,tau)

    return (1. + delt*dphird)*rhoi*Rig*T


#################################################################
#################################################################
#################################################################
#################################################################


def CALC_phi01(delt,tau):

    n0=nu.zeros(8)
    gam=nu.zeros(8)

    n0[0]=-8.32044648201; gam[0]=0.
    n0[1]=6.6832105268;   gam[1]=0.
    n0[2]=3.00632;        gam[2]=0.
    n0[3]=0.012436;       gam[3]=1.28728967
    n0[4]=0.97315;        gam[4]=3.53734222
    n0[5]=1.27950;        gam[5]=7.74073708
    n0[6]=0.96956;        gam[6]=9.24437796
    n0[7]=0.24873;        gam[7]=27.5075105


    sum=0.
    dsum=0.
    dstt=0.
    for i in range(3,8):
        sum=sum + n0[i]*m.log( 1.-m.exp(-1.*gam[i]*tau) )
        dsum=dsum + n0[i]*gam[i]*( (1. - m.exp(-1.*gam[i]*tau))**(-1.) -1. )
        dstt=dstt + n0[i]*gam[i]**2.*m.exp(-1.*gam[i]*tau)*(1. - m.exp(-1.*gam[i]*tau))**(-2.)


    phi0=m.log(delt) + n0 [0] + n0 [1]*tau + n0 [2]*m.log(tau) + sum
    dphi0t=n0 [1] +n0 [2]/tau + dsum
    dphi0tt=-1.*n0[2]/tau**2. - dstt

    return phi0, dphi0t, dphi0tt



#################################################################
#################################################################
#################################################################
#################################################################




def CALC_phir1(delt,tau):


#---------------------INITIALIZATIONS

    c=nu.zeros(56)
    d=nu.zeros(56)
    ti=nu.zeros(56)
    n=nu.zeros(56)
    al=nu.zeros(56)
    be=nu.zeros(56)
    gam=nu.zeros(56)
    eps=nu.zeros(56)
    a=nu.zeros(56)
    b=nu.zeros(56)
    B=nu.zeros(56)
    C=nu.zeros(56)
    D=nu.zeros(56)
    A=nu.zeros(56)
    bet=nu.zeros(56)
#------------------------------------


    c[0]=0.; d[0]=1.; ti[0]=-0.5;  n[0]=0.12533547935523e-1;
    c[1]=0.; d[1]=1.; ti[1]=0.875; n[1]=0.78957634722828e1;
    c[2]=0.; d[2]=1.; ti[2]=1.;    n[2]=-0.87803203303561e1;
    c[3]=0.; d[3]=2.; ti[3]=0.5;   n[3]=0.31802509345418e0;
    c[4]=0.; d[4]=2.; ti[4]=0.75;  n[4]=-0.26145533859358e0;
    c[5]=0.; d[5]=3.; ti[5]=0.375; n[5]=-0.78199751687981e-2;
    c[6]=0.; d[6]=4.; ti[6]=1.;    n[6]=0.88089493102134e-2;


    sum1=0.
    ds1d=0.
    ds1t=0.
    ds1tt=0.
    ds1dd=0.
    ds1dt=0.

    for i in range(0,7):
        sum1=sum1 + n[i]*delt**d[i]*tau**ti[i]
        ds1d=ds1d + n[i]*d[i]*delt**(d[i]-1.)*tau**ti[i]
        ds1t=ds1t + n[i]*ti[i]*delt**d[i]*tau**(ti[i]-1.)
        ds1tt=ds1tt + n[i]*ti[i]*(ti[i]-1.)*delt**d[i]*tau**(ti[i]-2.)
        ds1dd=ds1dd + n[i]*d[i]*(d[i]-1.)*delt**(d[i]-2.)*tau**ti[i]
        ds1dt=ds1dt + n[i]*d[i]*ti[i]*tau**(ti[i]-1.)*delt**(d[i]-1.)




    c[7]=1.;  d[7]=1. ;  ti[7]=4. ;    n[7]=-0.66856572307965
    c[8]=1.;  d[8]=1. ;  ti[8]=6. ;    n[8]=0.20433810950965
    c[9]=1.;  d[9]=1. ;  ti[9]=12. ;   n[9]=-0.66212605039687e-4 #>>>> ti
    c[10]=1.; d[10]=2. ; ti[10]=1. ;   n[10]=-0.19232721156002
    c[11]=1.; d[11]=2. ; ti[11]=5. ;   n[11]=-0.25709043003438
    c[12]=1.; d[12]=3. ; ti[12]=4. ;   n[12]=0.16074868486251
    c[13]=1.; d[13]=4. ; ti[13]=2. ;   n[13]=-0.40092828925807e-1
    c[14]=1.; d[14]=4. ; ti[14]=13. ;  n[14]=0.39343422603254e-6
    c[15]=1.; d[15]=5. ; ti[15]=9. ;   n[15]=-0.75941377088144e-5
    c[16]=1.; d[16]=7. ; ti[16]=3. ;   n[16]=0.56250979351888e-3
    c[17]=1.; d[17]=9. ; ti[17]=4. ;   n[17]=-0.15608652257135e-4
    c[18]=1.; d[18]=10. ; ti[18]=11. ; n[18]=0.11537996422951e-8
    c[19]=1.; d[19]=11. ; ti[19]=4. ;  n[19]=0.36582165144204e-6
    c[20]=1.; d[20]=13. ; ti[20]=13. ; n[20]=-0.13251180074668e-11
    c[21]=1.; d[21]=15. ; ti[21]=1. ;  n[21]=-0.62639586912454e-9
    c[22]=2.; d[22]=1. ;  ti[22]=7. ;  n[22]=-0.10793600908932
    c[23]=2.; d[23]=2. ;  ti[23]=1. ;  n[23]=0.17611491008752e-1
    c[24]=2.; d[24]=2. ;  ti[24]=9. ;  n[24]=0.22132295167546
    c[25]=2.; d[25]=2. ;  ti[25]=10. ; n[25]=-0.40247669763528
    c[26]=2.; d[26]=3. ;  ti[26]=10. ; n[26]=0.58083399985759
    c[27]=2.; d[27]=4. ;  ti[27]=3. ;  n[27]=0.49969146990806e-2
    c[28]=2.; d[28]=4. ;  ti[28]=7. ;  n[28]=-0.31358700712549e-1
    c[29]=2.; d[29]=4. ;  ti[29]=10. ; n[29]=-0.74315929710341
    c[30]=2.; d[30]=5. ;  ti[30]=10. ; n[30]=0.47807329915480
    c[31]=2.; d[31]=6. ;  ti[31]=6. ;  n[31]=0.20527940895948e-1
    c[32]=2.; d[32]=6. ;  ti[32]=10. ; n[32]=-0.13636435110343
    c[33]=2.; d[33]=7. ;  ti[33]=10. ; n[33]=0.14180634400617e-1
    c[34]=2.; d[34]=9. ;  ti[34]=1. ;  n[34]=0.83326504880713e-2 #>>>>?
    c[35]=2.; d[35]=9. ;  ti[35]=2. ;  n[35]=-0.29052336009585e-1
    c[36]=2.; d[36]=9. ;  ti[36]=3. ;  n[36]=0.38615085574206e-1
    c[37]=2.; d[37]=9. ;  ti[37]=4. ;  n[37]=-0.20393486513704e-1
    c[38]=2.; d[38]=9. ;  ti[38]=8. ;  n[38]=-0.16554050063734e-2
    c[39]=2.; d[39]=10. ; ti[39]=6. ;  n[39]=0.19955571979541e-2
    c[40]=2.; d[40]=10. ; ti[40]=9. ;  n[40]=0.15870308324157e-3
    c[41]=2.; d[41]=12. ; ti[41]=8. ;  n[41]=-0.16388568342530e-4
    c[42]=3.; d[42]=3. ;  ti[42]=16. ; n[42]=0.43613615723811e-1
    c[43]=3.; d[43]=4. ;  ti[43]=22. ; n[43]=0.34994005463765e-1
    c[44]=3.; d[44]=4. ;  ti[44]=23. ; n[44]=-0.76788197844621e-1
    c[45]=3.; d[45]=5. ;  ti[45]=23. ; n[45]=0.22446277332006e-1
    c[46]=4.; d[46]=14. ; ti[46]=10. ; n[46]=-0.62689710414685e-4
    c[47]=6.; d[47]=3. ;  ti[47]=50. ; n[47]=-0.55711118565645e-9
    c[48]=6.; d[48]=6. ;  ti[48]=44. ; n[48]=-0.19905718354408
    c[49]=6.; d[49]=6. ;  ti[49]=46. ; n[49]=0.31777497330738
    c[50]=6.; d[50]=6. ;  ti[50]=50. ; n[50]=-0.11841182425981

    sum2=0.
    ds2d=0.
    ds2t=0.
    ds2tt=0.
    ds2dd=0.
    ds2dt=0.

    for i in range(7,51):
        sum2=sum2 + n[i] * delt**(d[i])*tau**ti[i]*m.exp(-1. *delt**c[i])
        ds2d=ds2d + n[i] * m.exp(-1.*delt**c[i]) * ( delt**(d[i]-1.) * tau**ti[i] * (d[i]-c[i]*delt**c[i]) )
        ds2t=ds2t + n[i] * ti[i]*delt**d[i]*tau**(ti[i]-1.) * m.exp(-1.*delt**c[i])
        ds2tt=ds2tt + n[i] * ti[i]*(ti[i]-1.)*delt**d[i]*tau**(ti[i]-2.) * m.exp(-1.*delt**c[i])
        ds2dd=ds2dd + n[i] * m.exp(-1.*delt**c[i])*( delt**(d[i]-2.)*tau**ti[i]*( (d[i]-c[i]*delt**c[i])*(d[i]-1.-c[i]*delt**c[i]) - c[i]**2.*delt**c[i] ) )
        ds2dt=ds2dt + n[i] * ti[i]* tau**(ti[i]-1.)*delt**(d[i]-1.)*(d[i]-c[i]*delt**c[i])*m.exp(-1.*delt**c[i])







    c[51]=0. ; d[51]=3. ;  ti[51]=0. ; n[51]=-0.31306260323435e2; al[51]=20. ; be[51]=150 ; gam[51]=1.21 ; eps[51]=1. ;
    c[52]=0. ; d[52]=3. ;  ti[52]=1. ; n[52]=0.31546140237781e2;  al[52]=20. ; be[52]=150 ; gam[52]=1.21 ; eps[52]=1. ;
    c[53]=0. ; d[53]=3. ;  ti[53]=4. ; n[53]=-0.25213154341695e4; al[53]=20. ; be[53]=250 ; gam[53]=1.25 ; eps[53]=1. ;

    sum3=0.
    ds3d=0.
    ds3t=0.
    ds3tt=0.
    ds3dd=0.
    ds3dt=0.
    for i in range(51,54):
        sum3=sum3 +     n[i]*delt**d[i] * tau**ti[i] * m.exp(-1.*al[i]*(delt-eps[i])**2. - be[i]*(tau-gam[i])**2.)
        ds3d=ds3d +     n[i]*delt**d[i] * tau**ti[i] * m.exp(-1.*al[i]*(delt-eps[i])**2. - be[i]*(tau-gam[i])**2.) * ( d[i]/delt - 2.*al[i]*(delt-eps[i]) )
        ds3t=ds3t +     n[i]*delt**d[i] * tau**ti[i] * m.exp(-1.*al[i]*(delt-eps[i])**2. - be[i]*(tau-gam[i])**2.) * ( ti[i]/tau-2. *be[i]*(tau-gam[i]) )
        ds3tt=ds3tt +   n[i]*delt**d[i] * tau**ti[i] * m.exp(-1.*al[i]*(delt-eps[i])**2. - be[i]*(tau-gam[i])**2.) * ( (ti[i]/tau-2.*be[i]*(tau-gam[i]))**2. - ti[i]/tau**2. - 2.*be[i] )

        ds3dd=ds3dd +  n[i]*tau**ti[i] * m.exp( -1.*al[i]*(delt-eps[i])**2. - be[i]*(tau-gam[i])**2. ) * ( -2.*al[i]*delt**d[i] + 4.*al[i]**2*delt**d[i]*(delt-eps[i])**2. - 4.*d[i]*al[i]*delt**(d[i]-1.)*(delt-eps[i]) + d[i]*(d[i]-1.)*delt**(d[i]-2.) )

        ds3dt=ds3dt +   n[i]*delt**d[i] * tau**ti[i] * m.exp(-1.*al[i]*(delt-eps[i])**2. - be[i]*(tau-gam[i])**2.) * (ti[i]/tau-2.*be[i]*(tau-gam[i])) * (d[i]/delt - 2.*al[i]*(delt-eps[i]))






    a[54]=3.5 ; b[54]=0.85 ;  B[54]=0.2 ; n[54]=-0.14874640856724 ; C[54]=28. ; D[54]=700 ; A[54]=0.32 ; bet[54]=0.3 ;
    a[55]=3.5 ; b[55]=0.95 ;  B[55]=0.2 ; n[55]=0.31806110878444 ;  C[55]=32. ; D[55]=800 ; A[55]=0.32 ; bet[55]=0.3 ;

    sum4=0.
    ds4d=0.
    ds4t=0.
    ds4tt=0.
    ds4dd=0.
    ds4dt=0.
    for i in range(54,56):
        psi=m.exp( -1.*C[i]*(delt-1.)**2.  - D[i]*(tau-1.)**2.  )# ok
        theta=(1.-tau)+A[i]*((delt-1.)**2. )**(1./(2.*bet[i]))# ok
        Delta=theta**2+B[i]*((delt-1.)**2.)**(a[i])# ok

        dDeltad=(delt-1.) * (A[i]*theta*2./bet[i] * ((delt-1.)**2. )**(1./(2.*bet[i])-1.) + 2.*B[i]*a[i]*((delt-1.)**2. )**(a[i]-1. )  )
        dDeltbid=b[i]*Delta**(b[i]-1.)*dDeltad
        dDeltbit=-2.*theta*b[i]*Delta**( b[i]-1.)
        dpsid=-2.*C[i]*(delt-1.)*psi
        dpsit=-2.*D[i]*(tau-1.)*psi

        dpsitt=( 2.*D[i]*(tau-1.)**2. - 1. )*2.*D[i]*psi
        dpsidd=(2.*C[i]*(delt-1)**2.-1.)*2.*C[i]*psi
        dpsidt=4.*C[i]*D[i]*(delt-1.)*(tau-1.)*psi

        dDeldd=1./(delt-1.)*dDeltad + (delt-1.)**2*( 4.*B[i]*a[i]*(a[i]-1.)*((delt-1.)**2.)**(a[i]-2.) + 2.*A[i]**2.*(1./bet[i])**2.*( ((delt-1.)**2.)**(1./(2.*bet[i])-1.) )**2 + A[i]*theta*4./bet[i]*(1./(2.*bet[i])-1.)*((delt-1.)**2.)**(1./(2.*bet[i])-2.) )
        dDelbitt=2.*b[i]*Delta**(b[i]-1.) + 4.*theta**2*b[i]*(b[i]-1.)*Delta**(b[i]-2.)
        dDelbidd=b[i]*( Delta**(b[i]-1.)*dDeldd + (b[i]-1.)*Delta**(b[i]-2.)*dDeltad**2. )
        dDelbidt=-1.*A[i]*b[i]*2./bet[i]*Delta**(b[i]-1.)*(delt-1.)*((delt-1.)**2. )**(1./(2.*bet[i])-1.) - 2.*theta*b[i]*(b[i]-1.)*Delta**(b[i]-2.)*dDeltad

        sum4=sum4 + n[i]*Delta**(b[i])*delt*psi
        ds4t=ds4t + n [i]*delt*( dDeltbit*psi + Delta**b[i]*dpsit )
        ds4d=ds4d + n[i]*( Delta**b[i]*(psi+delt*dpsid) + dDeltbid*delt*psi )
        ds4tt=ds4tt + n[i]*delt*( dDelbitt*psi + 2.*dDeltbit*dpsit + Delta**b[i]*dpsitt)

        ds4dd=ds4dd + n[i]*( Delta**b[i]*(2.*dpsid+delt*dpsidd) + 2.*dDeltbid*(psi+delt*dpsid) + dDelbidd*delt*psi )
        ds4dt=ds4dt + n[i]*( Delta**b[i]*(dpsit+delt*dpsidt) + delt*dDeltbid*dpsit + dDeltbit*(psi+delt*dpsid) + dDelbidt*delt*psi)




    phir=sum1+sum2+sum3+sum4
    dphird=ds1d+ds2d+ds3d+ds4d
    dphirt=ds1t+ds2t+ds3t+ds4t
    dphirtt=ds1tt+ds2tt+ds3tt+ds4tt
    dphirdd=ds1dd+ds2dd+ds3dd+ds4dd
    dphirdt=ds1dt+ds2dt+ds3dt+ds4dt

    return phir,dphird,dphirt,dphirtt,dphirdt,dphirdd







    #################################################################
    #################################################################
    #################################################################
    #################################################################
def HPice(T): # give Pmelt in PASCALS
    
    pmelt=0.
    
    
    if (T>=273.31 and T<353.5): #ice VI
        pmelt=( 1.-1.07476*(1.-(T/273.31)**4.6) )*632.4 #MPa, Wagner 1994
    elif (T>=353.5 and T<800.) : # ice VII
        theta=T/355.
        pmelt=0.85*((theta)**3.47 -1.)+2.17 #GPa, Lin+ 2004
        pmelt=pmelt*1000. #MPa

    elif (T>=800 and T<=1273.+domega) :
        pmelt= 3.3488+1.8696e-5*T**2 #fitted from Hernandez+ 2018
        pmelt=pmelt*1000.
        
#    elif (T>1419. and T<=5000.):
#        pmelt= 26.528+1.2261e-3*T+6.2997e-6*T**2 #fitted from Hernandez+ 2018
#        pmelt=pmelt*1000.
#
#    elif (T>5000.):
#        pmelt=190e3 #Millot+ 2018


    if pmelt>0.:
        pmelt=pmelt*1e6 #Pa
        
    return pmelt



#################################################################
#################################################################
#################################################################
#################################################################

def Psat(T) :#Sat VP in Pa !!
    
    th=1.-T/Tc
    #print th
    #    if (T>273.15) : Psat=23.6027 - 4097.33 / (-34.7719+T)
    #    if (T<=273.15) : Psat=28.8760 - 6135.40 / (0.0+T)
    #    Psat=m.exp(Psat) #presumably, it's in Pascals (comparison w the fig. 1 of Kuramoto and Matsui)
    
    
    if (T>273.15 and T<=Tc) :
        Psat=Tc/T*(-7.85951783*th+1.84408259*th**1.5-11.7866497*th**3+22.6807411*th**3.5-15.9618719*th**4+1.80122502*th**7.5) #ddbst, after Wagner 2002
        Psat=m.exp(Psat)*pc #Pascals
        #valid from 273.15K to 633.15K ???? >> in Wagner 2002, valid up to the critical point, with an incertitude <0.1%
        
    if (T<=273.15) :
        th=T/273.16
        #print(T,th)
        Psat=-13.928169*(1.-th**(-1.5))+34.7078238*(1.-th**(-1.25))
        Psat=m.exp(Psat)*611.657 #Pascals
        #valid from 273.16K to 190K
    
    
    if (T>Tc):
        Psat=1e30


    return Psat

