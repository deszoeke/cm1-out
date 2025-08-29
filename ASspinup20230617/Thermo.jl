# Thermo module
# ! note change of Exner, theta interface !

module Thermo
# thermodynamics of moist air functions

using ForwardDiff

export theta_equiv
export moistad
# export theta_equiv_sat

# constants
const Cp = 1005.7  # from my Davies-Jones function, was 1005.
const Cpv = 1870.0 # J/kg/K
const Cw  = 4190.0
const L0 = 2.501e6 # J/kg

const C = 273.15 # K

const Rd = 287.04
const Rv = 461.5
const RdoRv=Rd/Rv

"latent heat of water vapor"
LvK(TempK) = 2.501e6 + (Cpv-Cw) * (TempK-273.0)

"""
es(T,p) = is saturation vapor pressure based on Wexler's formula,
with enhancement factor for moist air rather than water vapor.
The enhancement factor requires a pressure.
T [degrees C], p [Pa] (note the reversed input order), es [Pa]
Calling with optional keywords changes the units and
ignores the positional arguments.
es(T,p; TK=tk[Kelvin], P=pr[hPa])
From A. L. Buck 1981: JAM, 20, 1527-1532.
SPdeS 7 July 2004
"""
function es(T, p=1e5)
    # P in hPa
    P=p*1e-2
    esat = 1e2 * 6.1121*(1.0007 + 3.46e-8*P)*exp((17.502*T)/(240.97 + T)) # es in Pa
end

"es in hPa, same as keyword P[hPa]"
function es(T; P)
    # P in hPa
    esat =       6.1121*(1.0007 + 3.46e-8*P)*exp((17.502*T)/(240.97 + T)) # es in hPa, same as P
end

# # supply TK [Kelvin] by keyword, ignores positional T!!
# function es(T,p=1e5; TK=T+C, P=p*1e-2)
#     # P in hPa
#     T = TK - C
#     esat = 1e2 * 6.1121*(1.0007 + 3.46e-8*P)*exp((17.502*T)/(240.97 + T)) # convert es to Pa
# end

"""
qs(p,T) is saturation specific humidity based on Wexler's formula for es
with enhancement factor (see es.m).
p [Pa], T [degrees C], qs [kg/kg]
From A. L. Buck 1981: JAM, 20, 1527-1532.
SPdeS 7 July 2004
"""
function qs(p,T)
    esat = es(T,p) # T[C], p[Pa] method
    qsat = RdoRv*esat / (p + (RdoRv-1)*esat)
end

"dqsdT(p,T[C]) derivative of qs with respect to T at p,T by autodiff of Bolton's qs"
dqsdT(p,T) = ForwardDiff.derivative(t -> qs(p,t), T)

# wet bulb temperature methods
# for approximating the evap process

"General single Newton iteration to update x toward f(x) = fhat for a univariate function f"
updatex(f, x, fhat) = x + (fhat-f(x)) / ForwardDiff.derivative(f, x)

"""
Twet_autodiff(T[K], q[kg/kg], p[Pa]; niter=2) wet bulb temperature using Newton's method
for target specific humidity q[kg/kg]. Uses automatic differntiation.
"""
function Twet_autodiff(T, q, p; niter=2)
    f(t) = (t - T) + LvK((T+t)/2)/Cp * (qs(p,t-C) - q)
    t=T
    for i in 1:niter
        t = updatex(f, t, 0)
    end
    t
end
# 2 iterations converges to ~0.001 K

# call as...
# q = rh*qs(pa, Ta)
# Twet_autodiff(Ta, rh*qs(pa, Ta-C), pa)


# thermo functions

"Exner(p/p0) = T/θ = (p/p0)^(Rd/Cp)"
Exner(p) = p^0.287
Exner(p, qv) = p ^ ((Rd/Cp) * (1 - 0.28*qv))

"potential temperature [K] = theta(T[K], p[hPa])"
theta(T, p) = T / Exner(p/1000.0)
theta(T, p, qv) = T / Exner(p/1000.0, qv)

"vapor pressure ev(p,qv); ev has units of p"
ev(p, qv) = p*qv / (RdoRv + qv)

"T[K] LCL from Bolton, Temp[K], ev[hPa]"
Tlcl(Temp, ev) = 2840/(3.5*log(Temp) - log(ev) - 4.805) + 55;
Tlcl(Temp, p, qv) = Tlcl(Temp, ev(p,qv))

"equivalent potential temperature"
function theta_equiv(T::Real, p::Real, qv::Real)
    Tl = Tlcl(T, p, qv)
    thetae = theta(T,p,qv) * exp((3376.0/Tl - 2.54) * qv * (1 + 0.81*qv))
end

# qs(1e2*1000, 290.0-C)
"equivalent saturated potential temperature theta_equiv_sat(T[K], p[hPa])"
theta_equiv_sat(T, p) = theta_equiv( T, p, qs(1e2*p, T-C) )

g = 9.8
"""
moistad(T_Kelvin, p_Pa)
saturated moist adiabatic lapse rate (dT/dz)_m < 0
"""
moistad(T, p) = -g / ( LvK(T)*dqsdT(p,T-C) + Cp )           # fast to calculate
#moistad(T, p) = -(g/Cp) / ( LvK(T)/Cp * dqsdT(p,T-C) + 1 ) # clear relationship to g/Cp

gamma(T,p) = LvK(T)/g * dqsdT(p,T-C)
"""
Nsq_moist_add(T,p) adds latent heating effect on buoyancy.
    Nsq_moist_add(T,p) = -g^2/(Cp*T) * gamma/(1+gamma)
    gamma = L/g * dqsdT(p,T-C)
Compute moist Brunt-Vaisala frequency squared 
    Nsq_moist = Nsq + Nsq_moist_add(T,p)
Avoids differencing, avoids computing vertical gradients.
"""
Nsq_moist_add(T,p) = - g*g/(Cp*T) * gamma(T,p)/(1+gamma(T,p))

end # module Thermo
