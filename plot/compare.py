import matplotlib.pyplot as plt
import matplotlib.font_manager
myfont = matplotlib.font_manager.FontProperties(fname='/System/Library/Fonts/PingFang.ttc')
from helper import get_data_from

ci = 1.747
with open('data/casida.dat') as f:
    omegaHartree, omegaFock, omegaStatic, omegaDynamic = list(map(float, f.read().strip().split()))

omega, spectrumHartree, spectrumFock, spectrumStatic, spectrumDynamic = get_data_from('data/spectrum.dat')
omegaspec = [item - .001 for item in omega]
plt.vlines(ci, 0, 2000, colors='black', linestyles='dashed', label='CI')
plt.plot(omega, spectrumHartree, label='Hartree', color='blue')
plt.vlines(omegaHartree, 0, 2000, colors='blue', linestyles='dashed')
plt.plot(omega, spectrumFock, label='Hartree-Fock', color='green')
plt.vlines(omegaFock, 0, 2000, colors='green', linestyles='dashed')
plt.plot(omega, spectrumStatic, label='BSE_static', color='orange')
plt.vlines(omegaStatic, 0, 2000, colors='orange', linestyles='dashed')
plt.plot(omegaspec, spectrumDynamic, label='BSE_dynamic', color='red')
plt.vlines(omegaDynamic, 0, 2000, colors='red', linestyles='dashed')
plt.xlim((1.67, 1.77))
plt.ylim((0, 2000))
plt.xlabel(r'$\omega$')
plt.ylabel(r'$I(\omega)$')
plt.legend()
plt.show()
