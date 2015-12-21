import autocorr as ac
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    dirs=["masa_0.6", "masa_1.6",  "masa_2.6", "masa_3.6", "masa_4.6", "masa_5.6", "masa_6.6", "masa_7.6",
		    "masa_1.0",  "masa_2.0", "masa_3.0", "masa_4.0", "masa_5.0", "masa_6.0", "masa_7.0", "masa_8.0" ]
    archivos = [ "veloc_fac_1.dat",  "veloc_fac_2.dat",  "veloc_fac_3.dat",  "veloc_ver_1.dat"]

    name = dirs[0] + "/" + archivos[0]
    Ntran = 10000
    datos = np.loadtxt(name,skiprows=1+Ntran)
    N = np.size(datos[:,0])

    print N

    vxx = np.zeros(2*N)
    vyy = np.zeros(2*N)
    vzz = np.zeros(2*N)
    i=1
    for d in dirs:
        for arch in archivos:
            name = d+"/"+arch
            datos = np.loadtxt(name,skiprows=1 + Ntran)

	    vx = datos[:,0]
	    vy = datos[:,1]
	    vz = datos[:,2]
    
            vxx = vxx + ac.autocorrelacion_w(vx)
            vyy = vyy + ac.autocorrelacion_w(vy)
            vzz = vzz + ac.autocorrelacion_w(vz)
	    figname =("%s%03d")%("vxx",i)
	    print np.max(vxx)/np.real(i), np.min(vyy)/np.real(i), np.average(vzz)/np.real(i) 
	    #plt.ylim(ymax = 0.000003 , ymin = -0.000003 )
	    #0.000203522824707+0j) (-0.000150
	    #plt.plot(vxx[:2000]/np.real(i))
	    #plt.savefig(figname)
	    #plt.close()
	    i = i + 1

    N = 1500
    h=np.fft.fft(vxx / np.real(i))
    h=np.real(h*np.conj(h))
    print np.sum(h)
    plt.plot(np.real(h[:N]))
    figname = "psd_vxx.pdf"
    np.savetxt("psd_vxx.txt",h)
    plt.savefig(figname)
    plt.close()


    h=np.fft.fft(vyy / np.real(i))
    h=np.real(h*np.conj(h))
    print np.sum(h)
    plt.plot(np.real(h[:N]))
    figname = "psd_vyy.pdf"
    np.savetxt("psd_vyy.txt",h)
    plt.savefig(figname)
    plt.close()


    h=np.fft.fft(vzz / np.real(i))
    h=np.real(h*np.conj(h))
    print np.sum(h)
    plt.plot(np.real(h[:N]))
    figname = "psd_vzz.pdf"
    np.savetxt("psd_vzz.txt",h)
    plt.savefig(figname)
    plt.close()




