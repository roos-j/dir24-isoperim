from flint import arb

def G1(x, y, B, b):
    '''Function :math:`G^1_b[B](x,y)`

    Meant to take `arb` or `float` parameters.
    '''
    return ((y-x)**(1/b)+B(y)**(1/b))**b + B(x) - 2*B((x+y)/2)

def G2(x, y, B, b):
    '''Function :math:`G^2_b[B](x,y)`
    
    Meant to take `arb` or `float` parameters.
    '''
    return y-x + (2**b-1)*B(y) + B(x) - 2*B((x+y)/2)

def G(x, y, B, b):
    '''Function `G_b[B](x,y)`'''
    return max(G1(x, y, B, b), G2(x, y, B, b))

def left(i): 
    '''Left half of interval.'''
    return (i[0], (i[0]+i[1])/2)

def right(i): 
    '''Right half of interval.'''
    return ((i[0]+i[1])/2, i[1])

def intvl_exact(x):
    return x[0].is_exact() and x[1].is_exact()
    
def contains_root(f, x):
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
    '''Naive root finding by interval bisection of an initial guess; roughly up to current precision.
    
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
    If successful, return admissible partition as list of rectangles, each rectangle given by a pair of exact intervals.

    Parameters:
    g --- Lower bound function that takes rectangle parameter *xm, xM, ym, yM*
    x -- Interval in x coordinate; must be exact
    y -- Interval in y coordinate; must be exact
    depth -- Initial depth, used for recursion (default: 0)
    maxDepth -- Maximum depth (default: 12)

    This implementation is written for simplicity and readability, not for best possible performance.
    '''
    assert intvl_exact(x) and intvl_exact(y)
    if g(*x, *y) > 0:
        return True, [(x, y)]
    elif depth >= maxDepth:
        return False, [(x, y)]
    else:
        rv = []
        s, t = part_rect(g, left(x), left(y), depth+1, maxDepth)
        if not s: return s, t
        else: rv += t
        s, t = part_rect(g, right(x), left(y), depth+1, maxDepth)
        if not s: return s, t
        else: rv += t
        s, t = part_rect(g, left(x), right(y), depth+1, maxDepth)
        if not s: return s, t
        else: rv += t
        s, t = part_rect(g, right(x), right(y), depth+1, maxDepth)
        if not s: return s, t
        else: rv += t
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
    elif depth >= maxDepth:
        return False, [x[0], x[1]]
    else:
        s, t = part_intvl(g, left(x), depth+1, maxDepth)
        if not s: return s, t
        else: rv += t
        s, t = part_intvl(g, right(x), depth+1, maxDepth)
        if not s: return s, t
        else: rv += t
        return True, rv

def min_val_rect(g, rects):
    '''Return minimum value of g on given partition of rectangles.'''
    return min([g(*x, *y) for (x,y) in rects])

def min_val_intvl(g, intvls):
    '''Return minimum value of g on given partition of intervals.'''
    return min([g(intvls[i], intvls[i+1]) for i in range(len(intvls)-1)])

# def fake_erfinv(x):
#     '''Shouldn't be used.'''
#     if x < -1 or x > 1: return arb("nan")
#     elif x == 1: return arb("inf")
#     elif x == -1: return arb("-inf")
#     else:
#         return arb(float(sympy.erfinv(x)), rad=1E-15) # Should be fine but need to check

b0 = 0.5+arb("19/32768")
b1 = 0.5+arb("31/1024")
c0 = 1-arb("3/1024")

def L(x: arb, b: arb) -> arb:
    '''Logarithmic function :math:`L_b(x)`'''
    if x == arb(0): return arb(0)
    else: return x*(arb.log(1/x)/arb.log(arb(2)))**b

def Q(x: arb, b: arb) -> arb:
    '''Cubic function :math:`Q_b(x)`'''
    return 2*x/3*(1-x)*(2**(2+b)-3 + (12-2**(3+b))*x)

def alpha0(b: arb) -> arb:
    return 2**(2+b)-5

def alpha1(b: arb) -> arb:
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

def Jw(x: arb, w: arb) -> arb:
    '''Rescaled Gaussian isoperimetric profile'''
    return arb(2)**.5*arb(w)*bobkovI((1-arb(x))/arb(w))

w0 = find_root(lambda w: Jw(arb(.5), w)-.5, (arb(.75), arb(1)))
x0 = 1-w0/2

def J(x: arb) -> arb:
    '''Specific rescaling that we use'''
    return Jw(x, w0)

def DJ(x: arb) -> arb:
    '''Derivative of J'''
    return arb(2)**.5*PhiInv((1-x)/w0)
