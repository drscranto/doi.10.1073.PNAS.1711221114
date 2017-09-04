# DDE for a host with temperature dependent demographic rates

from PyDDE import pydde
from numpy import *

def define_params(spp):
    if spp=='med':
        # Parameter values for Harlequin bug
        bTopt = 0.8921
        Toptb = 298.2617
        sb =  3.0850
        TRm = 297.
        mJTR = 0.26638
        AmJ = 12651.
        ALE = -100000
        AHE = 53338. 
        TLE = 288.1 
        THE = 305.
        dJTR = 0.0547
        AdJ = 11690. 
        TRdJ = 297.
        dATR = 0.00287 
        AdA = 16824.
        TRdA = 297.
        Aq = float(AdA)
        Toptq = float(Toptb)
        qs = float(sb)
        qTR=0.2
        TRq = float(TRdA)
        #Seasonal variation in temperature - Mediterranean
        meanT = 290.0955
        amplT = 4.879776
        shiftT = 4.1904668
        Tmax = meanT + amplT
        qTopt = qTR*exp(Aq*(1./TRq - 1./Tmax))
    if spp=='trop':
        # Parameter values for pod sucking bug
        bTopt = 8.9313 
        Toptb = 300.44
        sb = 3.6419 
        TRm = 298.
        mJTR = 0.037495313
        AmJ = 5831.3
        ALE = 0. 
        AHE = 0. 
        TLE = 0. 
        THE = 0. 
        dJTR = 0.012867748 
        AdJ = 23770.
        TRdJ = 298.
        dATR = 0.026525199 
        AdA = 9710. 
        TRdA = 298. 
        #Temperature dependence of competition
        Aq = float(AdA) 
        Toptq = float(Toptb) 
        qs = float(sb) 
        qTR=0.2
        TRq = float(TRdA)
        #Seasonal variation in temperature - Tropical
        meanT = 300.2086
        amplT = 1.375938
        shiftT = 0.5814971
        Tmax = meanT + amplT
        qTopt = qTR*exp(Aq*(1./TRq - 1./Tmax))
    if spp=='temp':
        # Parameters for the Green plant bug
        bTopt = 2.71727
        Toptb = 298.8156
        sb = 8.096 
        TRm = 293.
        mJTR = 0.0377
        AmJ = 4132.8 
        ALE = -100000 
        AHE = 39404.
        TLE = 273. 
        THE = 310.396 
        dJTR = 0.02265
        AdJ = 6268.
        TRdJ = 298.
        dATR = 0.0293
        AdA = 4366.
        TRdA = 298.
        Aq = float(AdA) 
        Toptq = float(Toptb) 
        qs = float(sb)
        qTR=0.2
        TRq = float(TRdA)
        #Seasonal variation in temperature - Temperate
        meanT = 285.1977
        amplT = 15.243049
        shiftT = 4.4732788
        Tmax = meanT + amplT
        qTopt = qTR*exp(Aq*(1./TRq - 1./Tmax))

    return [bTopt,Toptb,sb,TRm,mJTR,AmJ,ALE,AHE,TLE,THE,
            dJTR,AdJ,TRdJ,dATR,AdA,TRdA,
            Aq,Toptq,qs,qTR,TRq, meanT,amplT,shiftT,qTopt]


#################################
## CONSTANT PARAMETERS
yr = 365.00
delta_years = 100.
dde_tstep = 1.
dde_stsc=array([1e-20,0,0,0])
#################################
## PARAMETERS TO DEFINE
spp = 'med'
max_years = delta_years + 2. 
keep_years = 3. 
delta_mean = 0.
delta_ampl = 0.

delta_ind = int(delta_mean)

q_form = 'unimodal' #'monotonic' #'unimodal
dd = 'fec' #'fec' #'Amort' #'Jmort'
sub_name ='ampl'

## babysit the algorithm ##
tol=1e-8
dde_dt=0.1
dde_hbsize=1e6

#################################

mypars = define_params(spp)
bTopt = mypars[0]
Toptb = mypars[1]
sb = mypars[2]
TRm = mypars[3]
mJTR = mypars[4]
AmJ = mypars[5]
ALE = mypars[6]
AHE = mypars[7]
TLE = mypars[8]
THE = mypars[9]
dJTR = mypars[10]
AdJ = mypars[11]
TRdJ = mypars[12]
dATR = mypars[13]
AdA = mypars[14]
TRdA = mypars[15]
Aq = mypars[16]
Toptq = mypars[17]
qs = mypars[18]
qTR = mypars[19]
TRq = mypars[20]
meanT = mypars[21]
amplT = mypars[22]
shiftT = mypars[23]
qTopt = mypars[24]

#################################
## FUNCTIONS

def T(t):
    # temperature in K as a function of time in days
    # define slope: climate will warm delta_mean degrees in delta_years
    m_mean = (delta_mean/delta_years)*(1/yr)
    m_ampl = (delta_ampl/delta_years)*(1/yr)
    if t < 0:
        return meanT
    elif t < (delta_years*yr):
        return (meanT + m_mean*t) + (amplT+m_ampl*t)*sin(2*pi*t/yr + shiftT)
    else:
        return (meanT + delta_mean) + (amplT+delta_ampl)*sin(2*pi*t/yr + shiftT)

def boltzmann_arrhenius(kTr, Ak, Tr, T):
    return kTr * exp(Ak * (1./Tr - 1./T))

def sharpe_schoolfield(kTr, A, Tr, AL, TL, AH, TH, T):
    return kTr * T/Tr * exp(A * (1./Tr - 1./T)) / (1.+exp(AL*(1./TL-1./T))+exp(AH*(1./TH-1./T)))

def gaussian(kTr, Topt, s, T):
    return kTr * exp(-(T-Topt)**2/2./s**2)

# density dependence of competition
def unimodal(t):
    return gaussian(qTopt, Toptq, qs, T(t))
def monotonic(t):
    return boltzmann_arrhenius(qTR, Aq, TRq, T(t))
def constant(t):
    return qTR

## life history functions
# fecundity
def b(t):
    return gaussian(bTopt, Toptb, sb, T(t))

# death rates
def dJ(t):
    return boltzmann_arrhenius(dJTR, AdJ, TRdJ, T(t))
def dA(t):
    return boltzmann_arrhenius(dATR, AdA, TRdA, T(t))

# maturation rates
if (spp=='trop'):
    def mat_J(t):
        return boltzmann_arrhenius(mJTR, AmJ, TRm, T(t))
else:
    def mat_J(t):
        return sharpe_schoolfield(mJTR, AmJ, TRm, ALE, TLE, AHE, THE, T(t))

# density dependence
def define_allqs(dd,q_form):

    if q_form=='unimodal':
        q=unimodal
    elif q_form=='monotonic':
        q=monotonic
    elif q_form=='constant':
        q=constant

    def q1(x, t):
        if dd == 'fec':
            return exp(-q(t)*x)
        else:
            return 1
    def q2(x, t):
        if dd == 'Amort':
            return q(t)*x
        else:
            return 0
    #def q3(x, t):
        #if dd == 'Jg':
            #return exp(-q(t)*x)
        #else:
            #return 1
    def q4(x, t):
        if dd == 'Jmort':
            return q(t)*x
        else:
            return 0

    return [q1,q2,q4]

###################
# system of ddes

[q1,q2,q4] = define_allqs(dd,q_form)

def ddegrad(s, c, t):
    
    # constant past history
    J_lag = 0.
    A_lag = 0.
    MJ = 0.
    
    if t > s[3]:
        J_lag = pydde.pastvalue(0, t-s[3], 0)
        A_lag = pydde.pastvalue(1, t-s[3], 0)
        MJ = A_lag*b(t-s[3])*q1(A_lag, t-s[3])*(mat_J(t)/mat_J(t-s[3]))*s[2]

    if any(s < 0):
        print('spp: %s, +mean: %s, +ampl, %s' % (spp,delta_mean,delta_ampl))
        print('t: %e, J_lag: %e, s: %s' % (t, J_lag, s))
    
    dJdt = s[1]*b(t)*q1(s[1],t) - MJ - (1+q4(s[0],t))*dJ(t)*s[0]
    
    dAdt = MJ - (1+q2(s[1],t))*dA(t)*s[1]
    
    dSdt = s[2]*( (mat_J(t)/mat_J(t-s[3]))*(q4(J_lag, t-s[3])+1)*dJ(t-s[3]) - dJ(t)*(q4(s[0],t)+1) )
    
    dtaudt = 1. -  mat_J(t)/mat_J(t-s[3])

    return array([dJdt,dAdt,dSdt,dtaudt])

## EXEC DDE SOLVER

dde_times = arange(0.0, max_years*yr, dde_tstep)
dde_eg = pydde.dde()
ddeist = array([ 0., .1, exp(-dJ(-1e-3)/mJTR), 1./mJTR ])

dde_eg.dde(y=ddeist, times=dde_times,
       func=ddegrad, parms=(), hbsize=dde_hbsize,
       tol=tol, nlag=2, dt=dde_dt)

filename=sub_name+'/ddedata_'+spp+'_'+str(delta_ind)+'.txt'      
outFile=open(filename,'w')
outStr = str([0]+[x for x in dde_eg.data[((max_years-keep_years)*365./dde_tstep):-1,0]]+[0])
outFile.write(outStr)
outFile.write('\n')
outStr = str([0]+[x for x in dde_eg.data[((max_years-keep_years)*365./dde_tstep):-1,1]]+[0])
outFile.write(outStr)
outFile.write('\n')
outStr = str([0]+[x for x in dde_eg.data[((max_years-keep_years)*365./dde_tstep):-1,2]]+[0])
outFile.write(outStr)
outFile.write('\n')
outStr = str([0]+[x for x in dde_eg.data[((max_years-keep_years)*365./dde_tstep):-1,3]]+[0])
outFile.write(outStr)
outFile.write('\n')
outStr = str([0]+[x for x in dde_eg.data[((max_years-keep_years)*365./dde_tstep):-1,4]]+[0])
outFile.write(outStr)
outFile.write('\n')
outStr = str([0]+[T(x) for x in dde_eg.data[((max_years-keep_years)*365./dde_tstep):-1,0]]+[0])
outFile.write(outStr)
outFile.write('\n')
outFile.close()

