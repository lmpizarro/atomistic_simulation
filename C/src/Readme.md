
[Molecular Simulations](http://www.pages.drexel.edu/~cfa22/msim/msim.html)

compilar

gcc -c *.c

linkear

gcc -o mclj_npt mclj_npt.o -lm -lgsl -lblas

gcc -o mdlj_berp mdlj_berp.o -lm -lgsl -lblas
