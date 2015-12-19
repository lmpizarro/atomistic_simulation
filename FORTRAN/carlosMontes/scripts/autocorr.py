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

if __name__ == "__main__":
    main()
