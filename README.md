This repository contains code accompanying the paper ``Sharp isoperimetric inequalities on the Hamming cube near the critical exponent'' by P. Durcik, P. Ivanisvili, J. Roos.

The code produces provably correct bounds using interval arithmetic. It and all libraries it depends on are open source.


Installation
=============

1. Install a Python 3 distribution, Anaconda is recommended: https://www.anaconda.com/

2. Install Python-FLINT. The current latest release (0.6.0) lacks the inverse error function. Thus the library must currently be built from source:

    * Install git, see https://git-scm.com/

    * Clone the python-flint repository

          git clone https://github.com/flintlib/python-flint.git
    
    * Follow build instructions here: https://fredrikj.net/python-flint/setup.html

(It is strongly recommended, but not strictly necessary, to do this on Linux.)

3. Clone this repository

Usage
=======

To verify all computer-assisted claims in the paper use

    python run.py

