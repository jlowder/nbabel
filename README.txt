Filename:           nbabel.lisp
Date:               August 2016
Author:             Jason Lowdermilk <jlowdermilk@gmail.com>
Integration scheme: Predictor-corrector leapfrog
Compiler:           SBCL 1.2.13
Operating           system: Ubuntu 14.04 x86_64 3.13.0-76-generic
Hardware:           Intel(R) Core(TM) i7-4600M CPU @ 2.90GHz

Performance:

Time step: 1.e-3 N-body unit
Integration from t=0 to t=1.0
Initial Conditions: input`N` (Plummer distribution of `N` equal mass particles)

N    tCPU     |dE/E|
===  =======  =========
 64   1.551s  3.6361e-6
128   8.280s  2.1276e-6
256  52.169s  6.6877e-4

Building:

This code does not rely on any external packages or libraries. Therefore it should build and run on
practically any standard lisp implementation. I have tested with sbcl, ecl, and clisp.

To compile using sbcl:

sbcl --noinform --no-sysinit --no-userinit --disable-debugger --load nbabel.lisp
* (save-lisp-and-die "nbabel" :toplevel #'main :executable t)

Running:

./nbabel <input128
Time                    Kinetic Energy          Potential Energy        Total Energy            Energy Error
======================= ======================= ======================= ======================= =======================
0.0                     0.24996077365592212     -0.4999607732812352     -0.2499999996253131     -1.4987469043603607e-9
0.010000000000000002    0.2495915276127595      -0.49959152231817877    -0.24999999470541928    -2.1178322207049418e-8
0.02000000000000001     0.24926710693327656     -0.499267095986324      -0.24999998905304746    -4.3787809489437306e-8
0.03000000000000002     0.24899248355448364     -0.4989924671799096     -0.24999998362542594    -6.549829556501411e-8
0.04000000000000003     0.24877280416254166     -0.4987727826567151     -0.24999997849417344    -8.602330558993293e-8
.
<snip>
.
0.9700000000000008      0.2392809015858227      -0.48928045423486705    -0.24999955264904436    -1.789403821894632e-6
0.9800000000000008      0.23916524637005485     -0.4891648592304759     -0.24999961286042108    -1.5485583150232085e-6
0.9900000000000008      0.23961000669864183     -0.48960987828873453    -0.2499998715900927     -5.136396284965986e-7
1.0000000000000007      0.2402892846248178      -0.49028981652335946    -0.25000053189854166    2.1275941672893607e-6


License:

GNU GPL v3
