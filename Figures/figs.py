#!/usr/bin/env python
# coding: utf-8

# Main Text

# Fig. 1: Left is BTA results and Right is PEA results

# In[2]:


from pylab import *
from numpy import random
import matplotlib.pyplot as plt
import numba
from scipy.optimize import curve_fit
from matplotlib.ticker import MaxNLocator
import matplotlib as mpl
from matplotlib import rcParams



f1 = {'fontname':'cmss10'}
f2 = {'fontname':'cmr10'}
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.it'] = 'serif:italic' 



fig, (ax1, ax2) = plt.subplots(1, 2,figsize=(12, 4.8))




n = (np.array([1,2,3,4,5]))







EnExp = np.array([467,251,177,157,125])




QC = [366.7,232.0,176.3,138.5,117.0]
QC_is = [273.2,164.3,119.2,93.6,79.4]
QC_quartz = [416.0]

err = [26,21,19,19,29]
err_qc=[np.sqrt(1.313),np.sqrt(1.458),np.sqrt(2.728),np.sqrt(1.261),np.sqrt(3.828)]
err_qc_is=[np.sqrt(6.574),np.sqrt(0.846),np.sqrt(0.877),np.sqrt(0.694),np.sqrt(1.079)]
err_qc_quartz = [np.sqrt(1.054)]

ax1.plot(n,EnExp,marker='o',color='m',label='Exp. (Exfoliated)',linestyle='--',linewidth=3.6,markersize=10)
ax1.plot(n,QC,marker='s',color='royalblue',label='Film',linestyle='--',linewidth=3.6,markersize=10)
ax1.plot(n,QC_is, marker='D',color='g',label='Crystal',linestyle='--',linewidth=3.6,markersize=10)
#ax.scatter([1],QC_quartz, marker='^',color='royalblue', label = 'Quartz')
err_qc_quartz = []

plt.subplot(1,2,1)

plt.errorbar(n, QC, yerr = err_qc, fmt ='s',color='royalblue',capsize=6,markersize=10)
plt.errorbar(n, EnExp, yerr = err, fmt ='o',color='m',capsize=6,markersize=10)
plt.errorbar(n, QC_is, yerr = err_qc_is, fmt ='D',color='g',capsize=6,markersize=10)






#ax.legend(loc='best', frameon=False,fontsize='20')
ax1.tick_params(axis='x', length=6.8)
ax1.tick_params(axis='y', length=6.8)

plt.xticks(ticks=[1,2,3,4,5],fontsize='24',**f1)

plt.yticks(fontsize='24',**f1)
plt.ylim([50,530])
plt.xlim([0.8,5.2])
ax1.set_ylabel(r"$\mathit{E}_{\mathrm{B}}$ (meV)",fontsize='26',**f2)
ax1.set_xlabel(r"$\it{n}$",fontsize='26',**f2)
#ax.set_title("Renormalized Exciton Binding (BA)",fontsize = '18')
plt.tight_layout()



#ax.plot(n[:-1],En[:-1],linestyle='--', marker='o', label='Struve-Neumann PIMC')
#plt.plot(T, y(T), linestyle='dashed', color='r')
#ax.set_title(r"Binding Energy vs. Number of Sublayers (BA)")
#ax.set_ylabel(r"E (meV)")
#ax.set_xlabel(r"n")


#QC = [399.2,278.2,221.7,187.4,163.8]
QC_is_pea = [198.2,139.0,108.1,91.3,78.4]
QC_pea = [236.1,165.7,128.6,105.0,90.3]

EnExp_pea = np.array([220,170,125,100,82])



err_qc_is_pea = [np.sqrt(6.998),np.sqrt(10.907),np.sqrt(3.365),np.sqrt(7.832),np.sqrt(15.207)]
err_qc_pea=[np.sqrt(9.669),np.sqrt(36.760),np.sqrt(26.619),np.sqrt(13.420),np.sqrt(15.296)]


ax2.plot(n,EnExp_pea,marker='o',color='m',label='Exp. (Exfoliated)',linestyle='--',linewidth=3.6,markersize=10)
ax2.plot(n,QC_pea,marker='s',color='royalblue',label='Film',linestyle='--',linewidth=3.6,markersize=10)
ax2.plot(n,QC_is_pea, marker='D',color='g',label='Crystal',linestyle='--',linewidth=3.6,markersize=10)

plt.subplot(1,2,2)

plt.errorbar(n, QC_pea, yerr = err_qc_pea, fmt ='s',color='royalblue',capsize=6,markersize=10)
plt.errorbar(n, QC_is_pea, yerr = err_qc_is_pea, fmt ='D',color='g',capsize=6,markersize=10)






#ax.legend(loc='best', frameon=False,fontsize='20')
ax2.tick_params(axis='x', length=6.8)
ax2.tick_params(axis='y', length=6.8)

plt.xticks(ticks=[1,2,3,4,5],fontsize='24',**f1)

plt.yticks(fontsize='24',**f1)
plt.ylim([50,530])
plt.xlim([0.8,5.2])
ax2.set_ylabel(r"$\it{E}_{\mathrm{B}}$ (meV)",fontsize='26',**f2)
ax2.set_xlabel(r"$\it{n}$",fontsize='26',**f2)

plt.tight_layout()
#print(plt.rcParams['mathtext.it'])
#plt.gcf().set_size_inches(18, 6)


# Fig. 2) Here, we get the $P(r)$ for both perovskites. The first is for the BTA, and the second is for the PEA

# In[16]:


from pylab import *
from numpy import random
import matplotlib.pyplot as plt
import numba
from scipy.optimize import curve_fit
from matplotlib.ticker import MaxNLocator
import matplotlib as mpl
import numpy as np
from matplotlib import rcParams
f1 = {'fontname':'cmss10'}
f2 = {'fontname':'cmr10'}
##Put mathtex font 
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.it'] = 'serif:italic'

r_m = np.loadtxt('./bta/n1_r_nl.txt')
r_m = r_m.flatten()
dumb = np.histogram(r_m[:], density=True, bins = 100)

M = dumb[1][:-1]
yVal = dumb[0]



r_m1 = np.loadtxt('./bta/n1_r.txt')
r_m1 = r_m1.flatten() 
dumb1 = np.histogram(r_m1[:], density=True, bins = 100)

M1 = dumb1[1][:-1]
yVal1 = dumb1[0]


r_m2 = np.loadtxt('./bta/n3_r_nl.txt')
r_m2 = r_m2.flatten() 
dumb2 = np.histogram(r_m2[:], density=True, bins = 100)

M2 = dumb2[1][:-1]
yVal2 = dumb2[0]

r_m3 = np.loadtxt('./bta/n3_r.txt')
r_m3 = r_m3.flatten() 
dumb3 = np.histogram(r_m3[:], density=True, bins = 100)

M3 = dumb3[1][:-1]
yVal3 = dumb3[0]

r_m4 = np.loadtxt('./bta/n5_r_nl.txt')
r_m4 = r_m4.flatten()
dumb4 = np.histogram(r_m4[:], density=True, bins = 100)

M4 = dumb4[1][:-1]
yVal4 = dumb4[0]  



r_m6 = np.loadtxt('./bta/n5_r.txt')
r_m6 = r_m6.flatten() 
dumb6 = np.histogram(r_m6[:], density=True, bins = 100)

M6 = dumb6[1][:-1]
yVal6 = dumb6[0]





fig,ax = plt.subplots()
plt.figsize=(10, 6)
ax.tick_params(axis='x', length=6.8)
ax.tick_params(axis='y', length=6.8)

plt.xticks(fontsize='20',**f1)

plt.yticks(fontsize='20',**f1)
plt.ylim([-0.001,0.032])
plt.xlim([0,60])
ax.set_ylabel(r"$\it{P(r)}$",fontsize='22',**f2)
ax.set_xlabel(r"$\it{r} \: (\it{\AA})$",fontsize='22',**f2)
ax.plot(M/1.8897259886, yVal, label = r"$\it{n = }1$",color='red', linewidth=2.8, linestyle='-')
ax.plot(M1/1.8897259886, yVal1, color='red', linewidth=2.8, linestyle='--')


ax.plot(M2/1.8897259886, yVal2, color='blue', label = r"$\it{n = }3$", linewidth=2.8, linestyle='-')
ax.plot(M3/1.8897259886, yVal3, linewidth=2.8, color = 'blue', linestyle='--')



ax.plot(M4/1.8897259886, yVal4, label = r"$\it{n = }5$",color='green', linewidth=2.8, linestyle='-') 
ax.plot(M6/1.8897259886, yVal6, linewidth=2.8, color = 'green', linestyle='--')

ax.legend(loc='best', frameon=False,fontsize='22')
#ax.set_title("P$(r)$ for Film (BA)",fontsize = '18')
plt.gcf().set_size_inches(8, 6)


# In[15]:


from pylab import *
from numpy import random
import matplotlib.pyplot as plt
import numba
from scipy.optimize import curve_fit
from matplotlib.ticker import MaxNLocator
import matplotlib as mpl
import numpy as np
from matplotlib import rcParams
f1 = {'fontname':'cmss10'}
f2 = {'fontname':'cmr10'}
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.it'] = 'serif:italic'


r_m = np.loadtxt('./pea/n1_r_nl.txt')
r_m = r_m.flatten()
dumb = np.histogram(r_m[:], density=True, bins = 100)

M = dumb[1][:-1]
yVal = dumb[0]


r_m1 = np.loadtxt('./pea/n1_r.txt')
r_m1 = r_m1.flatten() 
dumb1 = np.histogram(r_m1[:], density=True, bins = 100)

M1 = dumb1[1][:-1]
yVal1 = dumb1[0]


r_m2 = np.loadtxt('./pea/n3_r_nl.txt')
r_m2 = r_m2.flatten() 
dumb2 = np.histogram(r_m2[:], density=True, bins = 100)

M2 = dumb2[1][:-1]
yVal2 = dumb2[0]

r_m3 = np.loadtxt('./pea/n3_r.txt')
r_m3 = r_m3.flatten() 
dumb3 = np.histogram(r_m3[:], density=True, bins = 100)

M3 = dumb3[1][:-1]
yVal3 = dumb3[0]


r_m4 = np.loadtxt('./pea/n5_r_nl.txt')
r_m4 = r_m4.flatten()
dumb4 = np.histogram(r_m4[:], density=True, bins = 100)

M4 = dumb4[1][:-1]
yVal4 = dumb4[0]  



r_m6 = np.loadtxt('./pea/n5_r.txt')
r_m6 = r_m6.flatten() 
dumb6 = np.histogram(r_m6[:], density=True, bins = 100)

M6 = dumb6[1][:-1]
yVal6 = dumb6[0]






fig,ax = plt.subplots()
plt.figsize=(10, 6)
ax.tick_params(axis='x', length=6.8)
ax.tick_params(axis='y', length=6.8)

plt.xticks(fontsize='20',**f1)

plt.yticks(fontsize='20',**f1)
plt.ylim([-0.001,0.032])
plt.xlim([0,60])
ax.set_ylabel(r"$\it{P(r)}$",fontsize='22',**f2)
ax.set_xlabel(r"$\it{r} \: (\it{\AA})$",fontsize='22',**f2)
ax.plot(M/1.8897259886, yVal, label = r"$\it{n = }1$",color='red', linewidth=2.8, linestyle='-')
ax.plot(M1/1.8897259886, yVal1, color='red', linewidth=2.8, linestyle='--')


ax.plot(M2/1.8897259886, yVal2, color='blue', label = r"$\it{n = }3$", linewidth=2.8, linestyle='-')
ax.plot(M3/1.8897259886, yVal3, linewidth=2.8, color = 'blue', linestyle='--')



ax.plot(M4/1.8897259886, yVal4, label = r"$\it{n = }5$",color='green', linewidth=2.8, linestyle='-') 
ax.plot(M6/1.8897259886, yVal6, linewidth=2.8, color = 'green', linestyle='--')

ax.legend(loc='best', frameon=False,fontsize='22')
#ax.set_title("P$(r)$ for Film (PEA)",fontsize = '18')

plt.gcf().set_size_inches(8, 6)


# Fig. 3) Now, we generate the figures for $\Delta E$ and $\Delta E_{B}$

# In[34]:


from pylab import *
from numpy import random
import matplotlib.pyplot as plt
import numba
from scipy.optimize import curve_fit
from matplotlib.ticker import MaxNLocator
import matplotlib as mpl

f1 = {'fontname':'cmss10'}
f2 = {'fontname':'cmr10'}
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.it'] = 'serif:italic' 


n = (np.array([1,2,3,4,5]))

ax = plt.figure().gca()
ax.xaxis.set_major_locator(MaxNLocator(integer=True))





QC = [30.96,33.04,32.56,31.59,30.33]
QC_pea = [21.15,24.74,25.82,26.06,25.92]

err_qc= [np.sqrt(0.018),np.sqrt(0.008),np.sqrt(0.273),np.sqrt(0.172),np.sqrt(0.177)]
err_qc_pea = [np.sqrt(0.084),np.sqrt(0.026),np.sqrt(0.005),np.sqrt(0.005),np.sqrt(0.005)]

ax.plot(n,QC,marker='o',color='royalblue',label='BTA', linestyle = '--',linewidth=3.6,markersize=10)
ax.plot(n,QC_pea,marker='o',color='firebrick',label='PEA',linestyle = '--',linewidth=3.6,markersize=10)
plt.errorbar(n, QC, yerr = err_qc, fmt ='o',color='royalblue',capsize=5,markersize=10)
plt.errorbar(n, QC_pea, yerr = err_qc_pea, fmt ='o',color='firebrick',capsize=5,markersize=10)




plt.xlim([0,5.5])
#plt.ylim([20,34.5])
plt.figsize=(8, 15)
ax.legend(loc='best', frameon=False,fontsize='22')
ax.tick_params(axis='x', length=6.8)
ax.tick_params(axis='y', length=6.8)

plt.xticks(ticks=[1,2,3,4,5],fontsize='20',**f1)

plt.yticks(fontsize='20',**f1)
plt.xlim([0.8,5.2])
plt.ylim([20,34])
ax.set_ylabel(r"$\mathit{\Delta E}$ (meV)",fontsize='22',**f2)
ax.set_xlabel(r"$\mathit{n}$",fontsize='22',**f2)
#ax.set_title("E-Ph Stabilization",fontsize = '18')
plt.gcf().set_size_inches(8, 6)


# In[35]:


from pylab import *
from numpy import random
import matplotlib.pyplot as plt
import numba
from scipy.optimize import curve_fit
from matplotlib.ticker import MaxNLocator
import matplotlib as mpl
from matplotlib import rcParams
f1 = {'fontname':'cmss10'}
f2 = {'fontname':'cmr10'}
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.it'] = 'serif:italic' 


ax = plt.figure().gca()
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
#ax.plot(n[:-1],En[:-1],linestyle='--', marker='o', label='Struve-Neumann PIMC')
#plt.plot(T, y(T), linestyle='dashed', color='r')






QC = [404.0,279.9,221.7,187.4,163.8]
QC_is = [313.1,208.0,163.4,141.9,128.1]




err_qc=[1.2,0.7,1.2,1.5,2.6]
err_qc_is = [3.4,0.80,0.6,2.1,1.0]

QCp = [366.7,232.0,176.3,138.5,117.0]
err_qcp=[np.sqrt(1.313),np.sqrt(1.458),np.sqrt(2.728),np.sqrt(1.261),np.sqrt(3.828)]

err_qc_isp=[np.sqrt(6.574),np.sqrt(0.846),np.sqrt(0.877),np.sqrt(0.694),np.sqrt(1.079)]
QC_isp = [273.2,164.3,119.2,93.6,79.4]

err_qcp_pea=[np.sqrt(9.669),np.sqrt(36.760),np.sqrt(26.619),np.sqrt(13.420),np.sqrt(15.296)]
QCp_pea = [236.1,165.7,128.6,105.0,90.3]

err_qc_pea=[1.8,3.1,0.7,1.1,0.4]
QC_pea = [257.5,192.0,157.6,137.0,122.1]

ex_ba = []
ex_pea = []
err_ba = []
err_pea = []

n=[1,2,3,4,5]
for i in range(len(QC_isp)):
    ex_ba.append(-(QCp[i]-QC[i]))
    ex_pea.append(-(QCp_pea[i]-QC_pea[i]))
    err_ba.append(np.sqrt(err_qc[i]**2 + err_qcp[i]**2))
    err_pea.append(np.sqrt(err_qc_pea[i]**2 + err_qcp_pea[i]**2))
    
    

ax.plot(n,ex_ba,marker='o',color='royalblue',label='BTA', linestyle = '--',linewidth=3.6,markersize=10)
ax.plot(n,ex_pea,marker='o',color='firebrick',label='PEA', linestyle = '--',linewidth=3.6,markersize=10)
#ax.plot(n,QC1,marker='o',color='firebrick',label='PEA (RK)',linestyle = '--',linewidth=1.6)
plt.errorbar(n, ex_ba, yerr = err_ba, fmt ='o',color='royalblue',capsize = 5,markersize=10)
plt.errorbar(n, ex_pea, yerr = err_pea, fmt ='o',color='firebrick',capsize = 5,markersize=10)
#plt.errorbar(n, QC1, yerr = err_qc1, fmt ='o',color='firebrick',capsize = 5)



plt.figsize=(8, 15)
ax.legend(loc='best', frameon=False,fontsize='22')
ax.tick_params(axis='x', length=6.8)
ax.tick_params(axis='y', length=6.8)

plt.xticks(ticks=[1,2,3,4,5],fontsize='20',**f1)

plt.yticks(fontsize='20',**f1)
#plt.ylim([0,500])
plt.xlim([0.8,5.2])
ax.set_ylabel(r"$\it{\Delta E}_{\mathrm{B}}$ (meV)",fontsize='22',**f2)
ax.set_xlabel(r"$\it{n}$",fontsize='22',**f2)
#ax.set_title("Ex-Ph Stabilization (BA)",fontsize = '18')



plt.gcf().set_size_inches(8, 6)


# Supplementary

# Fig.S1) Here, we plot the variational bound as a function of $n$ for both perovskites

# In[39]:


from pylab import *
from numpy import random
import matplotlib.pyplot as plt
import numba
from scipy.optimize import curve_fit
from matplotlib.ticker import MaxNLocator
import matplotlib as mpl
f1 = {'fontname':'cmss10'}
f2 = {'fontname':'cmr10'}
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.it'] = 'serif:italic' 


ep_st = 13.0
ep_in = 6.1
ep = ep_in
Lz = 1*6.39*1.8897259886
def var(ep_s,epo,w,m,L):
    w = w/1000/27.211399 
    ep_o = epo
    rp = Lz*epo*( (1/3)*(epo/ep) + (0.5*ep/epo)*(1-  (epo/ep)**2  )  )
    r = rp/Lz
    r = r*L
    alpha = L*(ep_s-ep_in)/(16*(r**2)*w)
    
    a = (r/ep_o)/(np.sqrt(1/(2*m*w)))


    s = ( (a**4 - a**2)*np.log(a) +np.pi*(a**3) - a**2 - a**4  )/ ( (1+a**2)**2  )
    gam = 4*alpha*s
    En = gam*w*1000*27.211399
    # print("Gamma is: ",gam)
    
    return (En)


n = np.linspace(1,16,1000)
bta = []
pea = []

for i in range(len(n)):
    bta = append(bta, var(ep_st,2.1,12.4,0.2,n[i]*Lz))
    pea = append(pea, var(ep_st,3.3,14.0,0.2,n[i]*Lz))
    
ax = plt.figure().gca()
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.plot(n,bta,color='royalblue',label='BA (Var.)',linestyle = '-',linewidth=1.8)
ax.plot(n,pea,color='firebrick',label='PEA (Var.)',linestyle = '-',linewidth=1.8)



plt.figsize=(8, 5)
#ax.legend(loc='best', frameon=False,fontsize='18')
ax.tick_params(axis='x', length=6.8)
ax.tick_params(axis='y', length=6.8)

plt.xticks(ticks=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15],fontsize='20',**f1)
plt.xlim([0.8,10.2])
plt.yticks(fontsize='20',**f1)
#plt.ylim([12,33])
ax.set_ylabel(r"$\Delta$E (meV)",fontsize='22',**f2)
ax.set_xlabel(r"n",fontsize='22',**f2)


# In[ ]:


matplotlib.pyplot.show()


# In[ ]:




