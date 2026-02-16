# (c) 2024 Joris Roos <jroos.math@gmail.com>
# pylint: disable=invalid-name,no-name-in-module

''' Some general definitions '''

import sys
from flint import arb

from .util import Log, log, FMT_FAIL, FMT_PASS, err, warn, Output

from .labels import lbl_dict


def left(i):
    '''Left half of interval.'''
    return (i[0], (.5*(i[0]+i[1])).upper())

def right(i):
    '''Right half of interval.'''
    return ((.5*(i[0]+i[1])).lower(), i[1])

def intvl_exact(x):
    '''Return True if interval endpoints are exact'''
    return x[0].is_exact() and x[1].is_exact()

def contains_root(f, x):
    '''Return True if f has different sign at interval endpoints'''
    return f(x[0])*f(x[1])<0

def _find_root_rec(f, x):
    '''See `find_root`'''
    assert intvl_exact(x)
    assert contains_root(f, x)
    l = left(x)
    r = right(x)
    if contains_root(f, l):
        x = _find_root_rec(f, l)
    elif contains_root(f, r):
        x = _find_root_rec(f, r)
    # Recursion will stop when root can no longer be guaranteed to lie in one of the two
    # i.e. when current precision limit is reached
    return x

def find_root(f, x):
    '''Naive root finding by interval bisection of an initial guess; 
    roughly up to current precision.
    
    f -- Function taking an `arb` 
    x -- Initial interval, assumed to be a 2-tuple of exact `arb`s

    Assumptions:
    - f is an enclosure of a continuous function
    - f has at least one root on x and x is not 
    
    Output is an `arb` that is guaranteed to contain the left-most root 
    '''
    x = _find_root_rec(f, x)
    assert intvl_exact(x)
    ans = arb.union(x[0], x[1])
    assert 0 in f(ans) # Sanity check
    return ans

def part_rect(g, x, y, depth=0, maxDepth=12):
    r'''Recursive dyadic partitioning on a given rectangle to prove positivity of given function.
    
    Return empty list on failure.
    If successful, return admissible partition as list of rectangles, 
    each rectangle given by a pair of exact intervals.

    Parameters:
    g --- Lower bound function that takes rectangle parameter *xm, xM, ym, yM*
    x -- Interval in x coordinate; must be exact
    y -- Interval in y coordinate; must be exact
    depth -- Initial depth, used for recursion (default: 0)
    maxDepth -- Maximum depth (default: 12)

    This implementation is written for simplicity and readability,
    not for best possible performance.
    '''
    assert intvl_exact(x) and intvl_exact(y)
    if g(*x, *y) > 0:
        return True, [(x, y)]
    if depth >= maxDepth:
        return False, [(x, y)]
    rv = []
    for (cx, cy) in [(left(x), left(y)), (right(x), left(y)),
                     (left(x), right(y)), (right(x), right(y))]:
        s, t = part_rect(g, cx, cy, depth+1, maxDepth)
        if not s:
            return s, t
        rv += t
    return True, rv

def part_intvl(g, x, depth=0, maxDepth=12):
    '''
    Partition interval to show positivity of g.
    Same as `part_rect` but in one dimension.

    On success, return partition of given interval.
    '''
    assert intvl_exact(x)
    rv = [x[0]] if depth == 0 else []
    if g(*x) > 0:
        return True, rv + [x[1]]
    if depth >= maxDepth:
        return False, [x[0], x[1]]
    s, t = part_intvl(g, left(x), depth+1, maxDepth)
    if not s:
        return s, t
    rv += t
    s, t = part_intvl(g, right(x), depth+1, maxDepth)
    if not s:
        return s, t
    rv += t
    return True, rv

def min_val_rect(g, rects):
    '''Return minimum value of g on given partition of rectangles.'''
    return min(g(*x, *y) for (x,y) in rects).lower()

def min_val_intvl(g, intvls):
    '''Return minimum value of g on given partition of intervals.'''
    return min(g(intvls[i], intvls[i+1]) for i in range(len(intvls)-1)).lower()

def tuple_to_arb(x):
    '''Convert interval in tuple format to an arb.'''
    return arb((x[0]+x[1])/2, (x[1]-x[0])/2)

def arb_to_tuple(x):
    '''Convert an arb to an interval as a tuple.'''
    return (x.lower(), x.upper())

# Default values of constants
b0 = arb(.50057) # 0.5+arb("37/65536")
b1 = 0.5+arb("31/1024")
c0 = arb("0.997") # 1.-arb("3/1024")

def L(x: arb, b: arb) -> arb:
    '''Logarithmic function :math:`L_b(x)`'''
    if x == arb(0):
        return arb(0)
    return x*(arb.log(1/x)/arb.log(arb(2)))**b

def Q(x: arb, b: arb) -> arb:
    '''Cubic function :math:`Q_b(x)`'''
    return 2*x/3*(1-x)*(2**(2+b)-3 + (12-2**(3+b))*x)

def alpha0(b: arb) -> arb:
    '''Constant alpha0'''
    return 2**(2+b)-5

def alpha1(b: arb) -> arb:
    '''Constant alpha1'''
    return 3-2**(1+b)

def DQ(x: arb, b: arb) -> arb:
    '''Derivative of cubic function'''
    return (-3+2**(2+b))*2/3 - 4*alpha0(b)*x - 8*alpha1(b)*x**2

def phi(t: arb) -> arb:
    '''Gaussian distribution function'''
    return (2*arb.pi())**(-.5)*arb.exp(-t*t/2)

def PhiInv(t: arb) -> arb:
    '''Inverse Gaussian cdf'''
    return arb(2)**.5*arb.erfinv(2*t-1)

def bobkovI(x: arb) -> arb:
    '''Gaussian isoperimetric profile'''
    return phi(PhiInv(x))

def wtox(w: arb) -> arb:
    return 1-w/2

class Jconst: # pylint: disable=too-few-public-methods
    '''Container for w0, x0.'''


# Verification

def verify(v):
    '''
    Call a generic verification routine (expected to return a boolean) and log result.
    Return result of verification.
    '''
    log(v.__name__ + ": ", end="")
    rv = v()
    log((FMT_PASS%"ok" if rv else FMT_FAIL%"fail"), indent=0)
    return rv

def verify_positive(g, x, y=None, maxDepth=12, verbose=1, tag="", **args):
    '''
    Verify that given lower bound function is positive using partitioning 
    on a given rectangle or interval and output result.
    '''
    msg = g.__name__
    if g.__name__ in lbl_dict:
        msg += f" [display ({lbl_dict[g.__name__]})] "
    if len(args) > 0:
        msg += " with "+", ".join([f"{k}={float(v)}" for k,v in args.items()])
    if verbose:
        log(msg + ": ", end="")

    def G(*p):
        return g(*p, **args)

    if y is None:
        success, part = part_intvl(G, x, maxDepth=maxDepth)
        if success:
            Output.get_instance().write_comment(msg)
            cmt = f"{len(part)-1:d} intervals, min. val = {min_val_intvl(G, part)}"
            Output.get_instance().write_part(g.__name__+tag, part, cmt)
            if verbose:
                log(FMT_PASS%"ok", indent=0)
                log(f"   {cmt}")
    else:
        success, part = part_rect(G, x, y, maxDepth=maxDepth)
        if success:
            Output.get_instance().write_comment(msg)
            cmt = f"{len(part):d} rectangles, min. val = {min_val_rect(G, part)}"
            Output.get_instance().write_part(g.__name__+tag, part, cmt)
            if verbose:
                log(FMT_PASS%"ok", indent=0)
                log(f"   {cmt}")
    if not success and verbose:
        log(FMT_FAIL%"fail", indent=0)
        log(f"   at {part}")
    return success, part

def batch_verify(label, methods, verbose=1):
    '''Run a list of verification methods.'''
    if verbose:
        log(label)
    Log.lvl = 1
    for m in methods:
        m()
    Log.lvl = 0

if not hasattr(arb, "erfinv"):
    err("version of python-flint too old: erfinv missing")
    sys.exit()
