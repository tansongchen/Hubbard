import matplotlib.pyplot as plt
import matplotlib.font_manager
myfont = matplotlib.font_manager.FontProperties(fname='/System/Library/Fonts/PingFang.ttc')
from Helper import get_data_from

omega, spectrum_re, spectrm_im, spectrum_abs = get_data_from('output/bse/spectrum.dat')

plt.plot(omega, spectrum_re, label = r'$\Re \tilde d$')
plt.plot(omega, spectrm_im, label = r'$\Im \tilde d$')
plt.plot(omega, spectrum_abs, label = r'$|\tilde d|$')
plt.xlim((1.58, 1.62))
plt.ylim((-2000, 2000))
plt.xlabel(r'$\omega$')
plt.ylabel(r'$I(\omega)$')
plt.legend()
plt.show()