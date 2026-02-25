This repository contains code verifying numerical claims in the papers [arXiv:2407.12674](https://arxiv.org/abs/2407.12674) and [arXiv:2602.20462](https://arxiv.org/abs/2602.20462).

The code produces provably correct bounds using interval/ball arithmetic relying on [FLINT/Arb](https://flintlib.org/doc/arb.html).

Installation
=============

1. Install a [Python 3](https://www.python.org/downloads/) distribution (at least 3.9).

2. Install the latest version of [python-flint](https://github.com/flintlib/python-flint) (at least 0.7.0).

        pip install python-flint

3. Install [Git](https://git-scm.com/downloads).

4. Clone this repository:

        git clone https://github.com/roos-j/dir24-isoperim.git

Usage
=======

To verify all computer-assisted claims in the papers use

    python run.py

To save partition data to file use

    python run.py --filename partitions.py

To view all command line options run

    python run.py -h
