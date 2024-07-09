This repository contains code accompanying the paper *Sharp isoperimetric inequalities on the Hamming cube near the critical exponent* by Polona Durcik, Paata Ivanisvili and Joris Roos.

The code produces provably correct bounds using interval/ball arithmetic relying on [FLINT/Arb](https://flintlib.org/doc/arb.html).

Installation
=============

1. Install a [Python 3](https://www.python.org/downloads/) distribution (at least 3.9)

2. Install the latest version of [python-flint](https://github.com/flintlib/python-flint). As of this writing, the latest stable release (0.6.0) lacks the inverse error function, but the function is available in recent pre-releases >=0.7.0a1. To install use:

        pip install python-flint>=0.7.0a1

3. Clone this repository

        git clone https://github.com/roos-j/dir24-isoperim.git

Usage
=======

To verify all computer-assisted claims in the paper use

    python run.py --beta 0.50057 --c 0.997

To save partition data to file use

    python run.py --filename partitions.py

To view all command line options run

    python run.py -h
