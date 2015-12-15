#!/usr/bin/python
# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt

axis_font = { 'size':'20'}
plt.figure(1,figsize=(8,6))


x = [5.0 ,    
     5.3 ,
     5.5 ,
     5.6 ,
     5.7 ,
     5.8 ,
     5.9 ,
     5.95,
     6.0 ,
     6.1 ,
     6.2 ,
     6.3 ,
     6.4 ,
     6.5 ,
     6.6 ,
     6.7 ]


y = [ 16.8   ,
       0.388 ,
       -4.03 ,
       -5.2  ,
       -5.92 ,
       -6.32 ,
       -6.488,
       -6.5  ,
       -6.483,
       -6.37 ,
       -6.17 ,
       -5.92 ,
       -5.54 ,
       -5.34 ,
       -5.72 ,
       -5.59 ]

x = np.asarray(x) / 4.0
y = np.asarray(y)


plt.plot(x,y,'k.',markersize=10)
plt.title(u'Energía potencial a $T^*=0$')
plt.xlabel('$a^*$',**axis_font)
plt.ylabel(u'$U^*$',**axis_font)
plt.ylim(-6.6,-5.1)
plt.xlim(1.35,1.7)
plt.grid(True)

plt.show()
