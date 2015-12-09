# autoCorrelation

## Purpose

Computes the autocorrelation function of a *d*-dimensional array signal for multiple particles.

## Background and theory

Please refer to manual.pdf for detailed instructions and a discussion of the mathematical background.

## Author

Written by [Maximilien Levesque](http://www.chimie.ens.fr/?q=pasteur/pct/Maximilien_Levesque/bio), [École Normale Supérieure](http://www.ens.fr), Paris, France.

## Important contributors

* [Vincent Dahirel](http://www.phenix.cnrs.fr/spip.php?rubrique91), UPMC, UMR 8234 PHENIX, Paris, France, for rewritting the brute force algorithm in a faster and more readable way, Oct. 2013.
* [Marie Jardat](http://www.phenix.cnrs.fr/spip.php?rubrique44), UPMC, UMR 8234 PHENIX, Paris, France, for discussions, extensive testing and bug reports. Sept. 2013.
* Xudong Zhao, UPMC, UMR 8234 PHENIX, Paris, France, for providing test cases, for discussions and Apple specific bug reports. Sept. 2013

## Changelog

* Sept. 2013:   Basic version of the code, with bruteforce algorithm only.
* Oct. 2013:    Various improvements of the code and of the bruteforce algorithm.
* Feb. 2014:    Introduction of the Fourier space algorithm, inducing runtimes divided by 100.

## How to make it work

### Download it

The best and easiest way is to use git:
```
    git clone https://github.com/maxlevesque/autoCorrelation  
```
will create a folder called autoCorrelation and download the latest version of the sources automatically.

You may also go to https://github.com/maxlevesque/autoCorrelation/archive/master.zip from your web brower.

Or use `get`, `wget`, ... or whatever seems good to you. I repeat: the best and easiest way to download and stay up-to-date is to use `git`.

### Compile it

The Fourier space algorithm requires the `FFTW3` library. Please visite their [website](www.fftw.org) for more informations.

```bash
$ make
```

## How to use `autoCorrelation`

### Inputs
All you need is a file, whatever its name, containing all your data, velocities (or forces, or ...) in the following format, e.g., for 3dim signal:
``` 
x(1,t) y(1,t) z(1,t)  
x(2,t) y(2,t) z(2,t)  
x(i,t) y(i,t) z(i,t)  
...    
x(Nat,t) y(Nat,t) z(Nat,t)    
x(1,t+1) y(1,t+1) z(1,t+1)  
x(2,t+1) y(2,t+1) z(2,t+1)  
x(i,t+1) y(i,t+1) z(i,t+1)  
...  
x(Nat,t+1) y(Nat,t+1) z(Nat,t+1)  
...  
...  
x(Nat,t+Nstep) y(Nat,t+Nstep) z(Nat,t+Nstep)  
```  

where `Nat` is the total number of atoms you have in your supercell,  
and `Nstep` is the number of timesteps in your trajectory. The signal has three columns: it is 3-dimensional.

In other words, you print the coordinates of all sites for a given timestep, then for the next one, etc.  
You do not print blank lines anywhere. You can add as many spaces you want between the columns.

## Execution

The executable is waiting for arguments:  
1. `Nat`, defined above  
2. `a`, an integer with value `1` or `2`, for bruteforce or fourierspace algorithm.  
3. `d`, an integer that declares how many columns in the file, i.e., if you're computing autocorrelation of a signal (1), a force (3), etc.  
3. `filename` of the trajectory in the format discussed above.  
  
So, you have to execute:  
`$ autoCorrelation 800 2 3 ./analysis/velocities.out`

## Outputs

Execution of `autoCorrelation` will result in a single file: `acf.out`.  
This file has a very simple ASCII format:
``` 
1 1.12931233  
2 1.24034134  
3 1.5908O123  
...  
Nstep 123987.12383  
``` 
where the first column is the timestep, and the second column is the autocorrelation function.  

### Units
Both timesteps and autocorrelation functions are in units of your simulation.  


## Disclaimer

This program has been thoroughly tested and validated.
I would not share a tool that is not production ready. I am very confident in this program.
*Nevertheless*, this program may contain bugs or restrictions that would lead to unexpected results.
Please, read carefully this readme file, and test it thoroughly on data for which you already know the results.

Comments, bug reports or even thanks ;) are welcome!

===================================================
Some equations (just to show how to type equations)
===================================================

Example ...

.. math::

    H = \sum_i h (i) + \frac{1}{2} \sum_{i \ne j} g (i, j) + V_{NN}; \quad V_{NN} = \frac{1}{2} \sum_{A \ne B} \frac{Z_A Z_B}{R_{AB}}

Here is some inline example: :math:`V_{eN}`.
And another example:

.. math::

    h_D = \beta mc^2 + c (\alpha \cdot \mathbf{p} ) + V_{eN}

And yet another example:

.. math::

   g^{Coulomb} (1, 2) = \frac{1}{r_{12}}
