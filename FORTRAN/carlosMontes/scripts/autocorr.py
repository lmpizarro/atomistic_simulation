#http://stackoverflow.com/questions/16044491/statistical-scaling-of-autocorrelation-using-numpy-fft

import numpy as np
fft = np.fft

def autocorrelation(x):
    """
    Compute autocorrelation using FFT
    The idea comes from 
    http://dsp.stackexchange.com/a/1923/4363 (Hilmar)
    """
    x = np.asarray(x)
    N = len(x)
    x = x-x.mean()
    s = fft.fft(x, N*2-1)
    result = np.real(fft.ifft(s * np.conjugate(s), N*2-1))
    result = result[:N]
    result /= result[0]
    return result

def AutoCorrelation(x):
    x = np.asarray(x)
    y = x-x.mean()
    result = np.correlate(y, y, mode='full')
    result = result[len(result)//2:]
    result /= result[0]
    return result 

def autocorrelate(x):
    fftx = fft.fft(x)
    fftx_mean = np.mean(fftx)
    fftx_std = np.std(fftx)

    ffty = np.conjugate(fftx)
    ffty_mean = np.mean(ffty)
    ffty_std = np.std(ffty)

    result = fft.ifft((fftx - fftx_mean) * (ffty - ffty_mean))
    result = fft.fftshift(result)
    return [i / (fftx_std * ffty_std) for i in result.real]


np.set_printoptions(precision=3, suppress=True)

def autocorrelacion_w(a):
    #la ventana
    aver = np.average(a)
    sigma=np.real(np.size(a)) /2.5
    idx=np.arange(0,np.size(a),1,dtype="int")
    window=np.exp(-idx*idx/(2*sigma*sigma)) / (sigma*np.sqrt(2*np.pi))
    #el padding
    b=np.zeros(2*np.size(a))
    aw = (a - aver)*window
    b[:np.size(a)] = aw


    c=fft.fft(b)
    c=c*np.conj(c)
    d=np.fft.ifft(c)
    return d

def psd_autocorr(a):
    d = autocorrelacion_w(a)	
    h=fft.fft(d)
    h=np.real(h*np.conj(h))
    return h


def main():
    """
    These tests come from
    http://www.maplesoft.com/support/help/Maple/view.aspx?path=Statistics/AutoCorrelation
    http://www.maplesoft.com/support/help/Maple/view.aspx?path=updates%2fMaple15%2fcomputation
    """
    tests = [
        ([1,2,1,2,1,2,1,2], [1,-0.875,0.75,-0.625,0.5,-0.375,0.25,-0.125]),
        ([1,-1,1,-1], [1, -0.75, 0.5, -0.25]),
        ]

    for x, answer in tests:
        x = np.array(x)
        answer = np.array(answer)
        # print(autocorrelate(x)) 
        print(autocorrelation(x))
        print(AutoCorrelation(x))
        assert np.allclose(AutoCorrelation(x), answer)
        print

    """
    Test that autocorrelation() agrees with AutoCorrelation()
    """
    for i in range(10):
        x = np.random.random(np.random.randint(2,10))
        assert np.allclose(autocorrelation(x), AutoCorrelation(x))


    import matplotlib.pyplot as plt

    N = 1000
    x1=np.linspace(0,6*np.pi,N)
    x2=np.linspace(0,12*np.pi,N)
    #la seal
    a=np.sin(x1) + np.sin(x2) + np.sin(np.linspace(0,24*np.pi,N)) + + np.sin(np.linspace(0,4*np.pi,N))
    a += 3 * np.random.rand(N)

    plt.plot(a)

    plt.show()


    x= AutoCorrelation(a)

    d = autocorrelacion_w(a)

    p = psd_autocorr(a)

    plt.plot(np.real(p[:np.size(a)/10]))
    #plt.plot(np.real(d))

    #plt.plot(a)
    plt.show()

    plt.plot(x)
    plt.show()


import matplotlib.pyplot as plt 

def graficar_autocorr(arch, Nl):
    datos = np.loadtxt(arch,skiprows=1)
    averx = np.average(datos[:,0])
    avery = np.average(datos[:,1])
    averz = np.average(datos[:,2])
    print averx, avery, averz
    px = psd_autocorr(datos[:,0]) 
    py = psd_autocorr(datos[:,1]) 
    pz = psd_autocorr(datos[:,2]) 
    plt.plot(px[:Nl])
    plt.plot(py[:Nl])
    plt.plot(pz[:Nl])
    plt.show()

def grafica_archivos():
    arch = "veloc_fac_1.dat" 	
    graficar_autocorr(arch)
    arch = "veloc_fac_2.dat" 	
    graficar_autocorr(arch)
    arch = "veloc_fac_3.dat" 	
    graficar_autocorr(arch)
    arch = "veloc_ver_1.dat" 	
    graficar_autocorr(arch)


def grafica_psd_componentes(Nl):
    archivos = [ "veloc_fac_1.dat",  "veloc_fac_2.dat",  "veloc_fac_3.dat",  "veloc_ver_1.dat"] 
    
    vxx = np.empty(0)
    vyy = np.empty(0)
    vzz = np.empty(0)
    for arch in archivos:
        datos = np.loadtxt(arch,skiprows=1)
	vx = datos[:,0]
	vy = datos[:,1]
	vz = datos[:,2]
	vxx =np.hstack((vxx,vx))
	vyy =np.hstack((vyy,vy))
	vzz =np.hstack((vzz,vz))
        averx = np.average(vx)
        avery = np.average(vy)
        averz = np.average(vz)
        print averx, avery, averz, averz + averx + avery
    averx = np.average(vxx)

    px = psd_autocorr(vxx)
    plt.plot(px[:Nl])
    px = psd_autocorr(vyy)
    plt.plot(px[:Nl])
    px = psd_autocorr(vzz)
    plt.plot(px[:Nl])

    plt.show()

def grafica_auto_corr_componentes():
    archivos = [ "veloc_fac_1.dat",  "veloc_fac_2.dat",  "veloc_fac_3.dat",  "veloc_ver_1.dat"] 
    
    vxx = np.empty(0)
    vyy = np.empty(0)
    vzz = np.empty(0)
    for arch in archivos:
        datos = np.loadtxt(arch,skiprows=1)
	vx = datos[:,0]
	vy = datos[:,1]
	vz = datos[:,2]
	vxx =np.hstack((vxx,vx))
	vyy =np.hstack((vyy,vy))
	vzz =np.hstack((vzz,vz))
        averx = np.average(vx)
        avery = np.average(vy)
        averz = np.average(vz)
        print averx, avery, averz, averz + averx + avery
    averx = np.average(vxx)
    
    Nl = np.size(vxx) / 2
    px = autocorrelacion_w(vxx)
    plt.plot(px[:Nl])
    px = autocorrelacion_w(vyy)
    plt.plot(px[:Nl])
    px = autocorrelacion_w(vzz)
    plt.plot(px[:Nl])

    plt.show()




if __name__ == "__main__":
    archivos = [ "veloc_fac_1.dat",  "veloc_fac_2.dat",  "veloc_fac_3.dat",  "veloc_ver_1.dat"] 
    datos = np.loadtxt(archivos[0],skiprows=1)
    N = 2*np.size(datos[:,0])
   
    vxx = np.zeros(N)
    vyy = np.zeros(N)
    vzz = np.zeros(N)
    for arch in archivos:
        datos = np.loadtxt(arch,skiprows=1)
	vx = datos[:,0]
	vy = datos[:,1]
	vz = datos[:,2]
    
        vxx = vxx + autocorrelacion_w(vx)
        vyy = vyy + autocorrelacion_w(vy)
        vzz = vzz + autocorrelacion_w(vz)

    h=fft.fft(vxx)
    h=np.real(h*np.conj(h))

    N = 300
    plt.plot(np.real(h[:N]))

    h=fft.fft(vyy)
    h=np.real(h*np.conj(h))
    plt.plot(np.real(h[:N]))

    h=fft.fft(vzz)
    h=np.real(h*np.conj(h))
    plt.plot(np.real(h[:N]))

    plt.show()


