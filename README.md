This repository contains code accompanying the paper ``Sharp isoperimetric inequalities on the Hamming cube near the critical exponent'' by P. Durcik, P. Ivanisvili, J. Roos.

The code produces provably correct bounds using interval arithmetic. It and the libraries it depends on are open source.


Installation
=============

Instructions for generic platform
----------------------
This will be most straightforward in a Unix-like environment.


1. Install a Python 3 distribution (at least 3.9)

2. Build and install python-flint. The current latest release (0.6.0) lacks the inverse error function. Thus the library must currently be built from source:

    * Install git: https://git-scm.com/

    * Clone the python-flint repository

          git clone https://github.com/flintlib/python-flint.git
    
    * Follow build instructions here: https://fredrikj.net/python-flint/setup.html

(It is strongly recommended, but not strictly necessary, to do this on Linux.)

3. Clone this repository

        git clone https://github.com/roos-j/dir24-isoperim.git


Specific instructions for Windows (x86-64)
--------------------

1. Install MSYS2 (https://www.msys2.org); it is recommended to use the default directory: `C:\msys64`

2. (Optional) Add the following to PATH environment variable:
        
        C:\msys64\; C:\msys64\usr\bin; C:\msys64\ucrt64\bin

3. Install required packages

    * Open UCRT64 shell: type UCRT64 in the Windows search bar or in a command prompt
    
    * Shell will open in the default home directory at `C:\msys64\home\[USERNAME]`
    
    * To install required packages enter:
    
            pacman -S mingw-w64-ucrt-x86_64-gcc
            pacman -S mingw-w64-ucrt-x86_64-cmake
            pacman -S mingw-w64-ucrt-x86_64-flint
            pacman -S mingw-w64-ucrt-x86_64-python
            pacman -S mingw-w64-ucrt-x86_64-python-pip
            pacman -S git

    (Packages will be installed to C:\msys64\ucrt64 or C:\msys64\usr)
    
4. Build and install python-flint (must be built from source):
        
        git clone https://github.com/flintlib/python-flint.git
        
        cd python-flint

        pip install .
        
        cd ..

5. Clone this repository

        git clone https://github.com/roos-j/dir24-isoperim.git

Usage
=======

To verify all computer-assisted claims in the paper use

    python run.py --beta 0.50057 --c 0.997

To view all command line options run

    python run.py -h
