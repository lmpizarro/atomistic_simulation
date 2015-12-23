import numpy as np
import pyfftw as fft
import sys

N = 128 

# mas info https://hgomersall.github.io/pyFFTW/sphinx/tutorial.html#the-workhorse-pyfftw-fftw-class

r = np.zeros(N)

a = fft.n_byte_align_empty(N, 16, 'complex128')
b = fft.n_byte_align_empty(N, 16, 'complex128')
c = fft.n_byte_align_empty(N, 16, 'complex128')


fft_object = fft.FFTW(a, b)
r[0:1] = 1
a[:] = r + 1j*0


ifft_object = fft.FFTW(b, c, direction='FFTW_BACKWARD')

fft_a = fft_object()
b = b*np.conj(b)
ifft_b = ifft_object()

print np.real(ifft_b[0:64]), np.imag(c)
