# nbabel
CL code submission for nbabel (http://nbabel.org), n-body code using leapfrog integration.

## Overview

    Filename:           nbabel.lisp
    Date:               July 2016
    Author:             Jason Lowdermilk <jlowdermilk@gmail.com>
    Integration scheme: Predictor-corrector leapfrog
    Compiler:           SBCL 1.2.13
    Operating           system: Ubuntu 14.04 x86_64 3.13.0-76-generic
    Hardware:           Intel(R) Core(TM) i7-4600M CPU @ 2.90GHz

## Performance

    Time step: 1.e-3 N-body unit
    Integration from t=0 to t=1.0
    Initial Conditions: input`N` (Plummer distribution of `N` equal mass particles)


    N      tCPU     dE/E
    ====== ======== ======================
    16       0.050s  6.388468209946375e-7
    32       0.145s -1.2211915069442633e-7
    64       0.451s  5.41408800014409e-7
    128      2.405s  4.6473401360595493e-7
    256      9.879s  3.3276633353171325e-7
    512     54.293s -4.270974619816691e-5
    1024   344.144s  6.046204287957143e-7

## Building

This code relies on cl-ppcre for regular expression support, which it loads with quicklisp.  It also
uses uiop for command line parsing.  Therefore, your lisp installation needs to have quicklisp and
asdf3 enabled (most modern lisps already do).

To compile using sbcl:

    sbcl --noinform --load nbabel.lisp
    * (save-lisp-and-die "nbabel" :toplevel #'main :executable t)
    
## Running

    ./nbabel input16
    Time                    Kinetic Energy          Potential Energy        Total Energy            Energy Error
    ======================= ======================= ======================= ======================= =======================
    0                       0.25000003              0.5                     -0.24999997             -0.0
    0.010000000000000002d0  0.2501344339915326d0    0.5001344418222696d0    -0.250000007830737d0    1.5053225550873724d-7
    0.02000000000000001d0   0.2502690298421953d0    0.5002690383206694d0    -0.2500000084784741d0   1.5312320416236297d-7
    0.03000000000000002d0   0.25040375922390457d0   0.5004037683904836d0    -0.25000000916657905d0  1.5587562433807292d-7
    0.04000000000000003d0   0.25053859153464725d0   0.5005386014297196d0    -0.2500000098950724d0   1.5878959797033604d-7
    .
    <snip>
    .
    0.9600000000000007d0    0.2598396053865419d0    0.5098396928180761d0    -0.2500000874315342d0   4.689354821691571d-7
    0.9700000000000008d0    0.2599005894241833d0    0.5099006765062735d0    -0.25000008708209015d0  4.6753770588545563d-7
    0.9800000000000008d0    0.2599615339033171d0    0.5099616206345376d0    -0.25000008673122054d0  4.6613422727357527d-7
    0.9900000000000008d0    0.2600224880671878d0    0.510022574448355d0     -0.25000008638116716d0  4.6473401360595493d-7

## Implementation Notes

This code was intended more for flexibility and precision than raw
speed. For example, the vector functions are written using the
recursive-folding style that will work with any number of vectors, of
any number of elements. The leapfrog algorithm itself is implemented
as a closure which can easily be substituted for different algorithms
for the sake of comparison. Double-floats are used for improved
precision.

## License

GNU GPL v3
