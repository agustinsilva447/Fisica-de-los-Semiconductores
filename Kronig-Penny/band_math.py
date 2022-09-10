import numpy as np
import matplotlib.pylab as plt
from kp_solution import KPM_Soln, Eband_KP

#System Variables
a = 2.7         #Angstroms
b = 2.7       #Angstroms
U0 = 5          #eV
eps_range = 6

epslist, f_eps = KPM_Soln(a,b,U0,eps_range)
k, eps = Eband_KP(epslist, f_eps)

#k values to plot for the RHS Solution
klist=np.linspace(0,eps_range,10)  
k_min=klist*0-1
k_max=klist*0+1

#Plot Kronig-Penny Solution
plt.figure(1,dpi=120)
plt.title("Kronig-Penny Solution")
plt.xlabel("$\epsilon$")
plt.ylabel("f($\epsilon$)")
plt.ylim(-5,5)
plt.plot(epslist,f_eps,label="LHS")
plt.plot(klist,k_min, linestyle="dashed",color="grey",label="RHS")
plt.plot(klist,k_max, linestyle="dashed",color="grey")
plt.legend()

plt.figure(2,dpi=120)
plt.title("Kronig-Penny Energy Bands")
plt.xlabel(r'k / $\frac{\pi}{(a+b)}$')
plt.ylabel("$\epsilon$")
for i in range(len(k)):
    plt.plot(k[i],eps[i])

plt.show()    

Ec, Ev = eps[2], eps[1]
Eg = (np.min(Ec)-np.max(Ev))*U0