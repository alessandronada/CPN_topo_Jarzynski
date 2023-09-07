import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import sys

import matplotlib.pyplot as plt
from collections import Counter
from matplotlib import gridspec


plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Computer Modern Roman']
plt.rcParams['text.usetex'] = True
plt.rcParams['lines.linewidth'] = 0.8
plt.rcParams['figure.figsize']=[14.54,7.34]
#colours = plt.rcParams['axes.color_cycle']
plt.rcParams['errorbar.capsize'] = 2
plt.rcParams['lines.markersize'] = 2
plt.matplotlib.rc('font', size=11)

@ticker.FuncFormatter
def major_formatter(x, pos):
    return "%d" % x

def autocorr(in_data, t):
    Nd=len(in_data)
    abst = np.abs(t)
    avg  = np.average(in_data)
    autc=0.
    for i in range(Nd-abst):
        autc = autc + (in_data[i]-avg)*(in_data[i+abst]-avg)
    out = autc / ((Nd-abst)*1.)
    return out

def autc_and_int(in_data, Num):

    ts = np.arange(0,Num)
    norm=autocorr(chi,0)
    aut = []
    for it in ts:
        aut.append( autocorr(in_data,it)/norm )
    aut_int = [ ]
    for it in ts:
        aut_int.append( np.sum(aut[:it]))
    return aut, aut_int

def read_topo(faddr):
    rawdata = np.genfromtxt(faddr,usecols=(0,1),dtype=[('t',np.int32),('T','f8')])
    print(rawdata)
    Tc = rawdata['T']
    return Tc

def read_chi(faddr):
    rawdata = np.genfromtxt(faddr,usecols=(0),dtype=[('chi','f8')])
    chi=rawdata['chi']
    return chi

def tau_int(in_data, M):
    aut_int = 0.5*np.sum([ in_data[i] for i in range(-M,M+1)])
    return aut_int



mode = sys.argv[1]
faddr = sys.argv[2]
Num=100

if( mode == "Top" ):
    Tc = read_topo(faddr)
    chi = Tc**2
    print(chi)
    aut, aut_int = autc_and_int(chi,100)
elif( mode == "chi" ):
    chi = read_chi(faddr)
    aut, aut_int = autc_and_int(chi,100)
else :
    print(sys.argv[1])
    exit("Error")


if( mode == 'chi'):
    fig, axs  = plt.subplots(2,1)

    axs[0].set_xlabel(r'Trajectory')
    axs[0].set_ylabel(r'$\hat{\hat{C}}(t)$')

    axs[0].plot(aut, marker='o')
    
    axs[1].set_xlabel(r'Trajectory')
    axs[1].set_ylabel(r'$\tau_{\mathrm{int}}$')
    M = np.argmin( np.abs([i-4.*tau_int(aut,i) for i in range(0,50)]))
    print(M)
    x_tauint = range(0,50)
    y_tauint = [tau_int(aut, i) for i in x_tauint]
    yerr_tauint = [np.sqrt(2.*(2.*i+1)*tau_int(aut,i)**2/len(aut)) for i in x_tauint]
    axs[1].errorbar(x_tauint, y_tauint, yerr=yerr_tauint, marker='o')
    axs[1].errorbar(M, tau_int(aut, M), yerr=np.sqrt(2.*(2.*M+1)*tau_int(aut,M)**2/len(aut)) , marker='o')

    plt.savefig(faddr+'.pdf')
    plt.show()

elif( mode == 'Top'):
    fig = plt.figure()
    gs = gridspec.GridSpec(3, 2, width_ratios=[2.3, 0.4])
    gs.update(wspace=0.05, hspace=0.1)
    
    axs01 = plt.subplot(gs[0,:])
    axs01.set_xlabel(r'Trajectory')
    axs01.set_ylabel(r'$\hat{\hat{C}}(t)$')
    axs01.plot(aut, marker='o')
    
    axs02 = plt.subplot(gs[1,:])
    axs02.set_xlabel(r'Trajectory')
    axs02.set_ylabel(r'$\tau_{\mathrm{int}}$')
    M = np.argmin( np.abs([i-4.*tau_int(aut,i) for i in range(0,20)]))
    print(M)
    axs02.errorbar(range(0,20), [tau_int(aut, i) for i in range(20)], yerr=[np.sqrt(2.*(2.*i+1)*tau_int(aut,i)**2/len(aut)) for i in range(20)], marker='o')
    axs02.errorbar(M, tau_int(aut, M), yerr=np.sqrt(2.*(2.*M+1)*tau_int(aut,M)**2/len(aut)) , marker='*')
    
    maxTC = np.max( np.abs(Tc))+2
    axs0 = plt.subplot(gs[2,0])
    axs0.yaxis.set_major_formatter(major_formatter)
    axs0.set_ylim( -maxTC, maxTC)
    axs0.set_ylabel(r'$Q$')
    #tc,stc=es.bs_avg_binned(TC['TC'])
    #lab=r'$N_c='+str(int(N))+'$, $\\beta='+str(beta)+'$, $L='+str(L)+'a$,\\\\ $\langle Q_L \\rangle = {:0.2uS}$'.format(ufloat(tc,stc)) 
    axs0.step( range(len(Tc)), Tc, lw=0.7, alpha=0.6)
    axs0.legend(frameon=False, fontsize=9, loc='upper right')
    axs0.set_xticks([])
    axs0.set_xlabel(r'Trajectory')
    
    axs1=plt.subplot(gs[2,1])
    axs1.set_ylim( -maxTC, maxTC)
    axs1.set_xticks([])
    axs1.set_yticklabels([])
    
    Q_bins = Counter(np.rint(Tc))
    range_min = min(min(Q_bins), -max(Q_bins)) - 1
    range_max = -range_min + 1
    Q_range = np.arange(range_min, range_max)
    
    Q_counts = [Q_bins[Q] for Q in Q_range]
    
    axs1.step(Q_counts, Q_range - 0.5)
    axs1.set_xlabel(r'Count')
    
    xdata=Q_range
    ydata=Q_counts
    ydata_err = np.sqrt(Q_counts)

    from scipy.optimize import curve_fit
    def f(x,A,m,s):
        return A*np.exp(-(x-m)**2 / (2.*s**2))
    popt, pcov = curve_fit(f, xdata, ydata)
    smooth_Q_range = np.linspace(range_min - 0.5, range_max - 0.5, 1000)
    axs1.plot(f(smooth_Q_range, *popt),smooth_Q_range)
    
    
    chi2 = np.sum( (ydata-f(xdata,*popt))**2/f(xdata,*popt))/( len(Q_counts)-len(popt))
    perr = np.sqrt(np.diag(pcov))
    
    plt.show()





