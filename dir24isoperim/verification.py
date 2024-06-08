from flint import arb

from .general import b0, b1, c0, Jconst, L, Q, DQ, J, DJ, G, alpha0, alpha1, \
                    part_rect, part_intvl, min_val_intvl, min_val_rect

from .util import Log, log, FMT_FAIL, FMT_PASS, err, warn, Output

#
#  Define functions for verification
# 

# Case J

# Manually constructed enclosures for J, |J'|
def Jm(xm: arb, xM: arb) -> arb: 
    return arb.min(J(xm), J(xM))

def JM(xm: arb, xM: arb) -> arb:
    if xM < Jconst.x0: return J(xM)
    elif xm > Jconst.x0: return J(xm)
    else: return J(Jconst.x0)

def absDJm(xm: arb, xM: arb) -> arb:
    if xM < Jconst.x0: return DJ(xM)
    elif xm > Jconst.x0: return -DJ(xm)
    else: return arb(0)

def absDJM(xm: arb, xM: arb) -> arb:
    return arb.max(abs(DJ(xm)), abs(DJ(xM)))

def g_J_1(xm: arb, xM: arb, hm: arb, hM: arb, b: arb, c: arb) -> arb:
    rv = b*c**(1-1/b)*JM(xm+hm, xM+hM)**(1-1/b)
    rv -= .5*b*(1-b)*c**(1-2/b)*Jm(xm+hm, xM+hM)**(1-2/b)*hM**(1/b)
    rv -= c/2*Jm(xm, xM)**(-1)*hM**(2-1/b)
    if xM < Jconst.x0:
        rv += c/4*DJ(xM)*J(xM)**(-2)*hm**(3-1/b)
        rv += c/32*DJ(xM)*(7+3*DJ(xM)**2)*J(xM)**(-4)*hm**(5-1/b)
    else: # *not* equivalent to xM >= x0 since x0 is not exact
        rv -= c/4*absDJM(xm, xM)*Jm(xm, xM)**(-2)*hM**(3-1/b)
        rv -= c/32*absDJM(xm, xM)*(7+3*absDJM(xm, xM)**2)*Jm(xm, xM)**(-4)*hM**(5-1/b)
    rv -= 7*c/48*(1+absDJM(xm, xM)**2)*Jm(xm, xM)**(-3)*hM**(4-1/b)
    rv -= c/90*(7+23*absDJM(xm, xM+hM)**2+6*absDJM(xm, xM+hM)**4)*Jm(xm, xM+hM)**(-5)*hM**(6-1/b)
    rv += c/2880*(7+23*absDJm(xm, xM+hM/2)**2+6*absDJm(xm, xM+hM/2)**4)*JM(xm, xM+hM/2)**(-5)*hm**(6-1/b)
    return rv

def g_J_2(xm: arb, xM: arb, ym: arb, yM: arb):
    return (ym-xM)**2 + J(yM)**2 - (2*J((xm+ym)/2)-J(xm))**2

# Case Q

def g_Q_1(hm: arb, hM: arb, ym: arb, yM: arb, b: arb) -> arb:
    rv = b*Q(yM,b)**(1-1/b)
    rv -= (alpha0(b)- 2*alpha1(b)*hm+4*alpha1(b)*yM)*hM**(2-1/b)
    rv -= b/2*(1-b)*Q(ym,b)**(1-2/b)*hM**(1/b)
    return rv    

def g_Q_1_y1_4(*p, b): return g_Q_1(arb(.25), arb(.25), *p, b)

def g_Q_2(hm: arb, hM: arb, ym: arb, yM: arb, b: arb) -> arb:
    rv = -12*alpha1(b)*hM**2
    rv += (2*alpha0(b)+8*alpha1(b)*ym)*hm 
    rv -= 6*(3-1/b)*alpha1(b)*Q(yM, b)**(1/b)*hM**(2-1/b)
    rv += (2-1/b)*(2*alpha0(b)+8*alpha1(b)*ym)*Q(ym, b)**(1/b)*hM**(1-1/b)
    return rv

# Case LJQ

def g_LJQ_1(xm: arb, xM: arb, ym: arb, yM: arb, b: arb) -> arb:
    l = ((ym-xM)**(1/b) + Jm(ym, yM)**(1/b))**b
    r = ym - xM + (2**b-1)*Jm(ym, yM)
    rv = arb.max(l, r)
    rv += L(xm, b)
    rv -= 2*Q((xM+yM)/2, b)
    return rv

def g_LJQ_2(ym: arb, yM: arb, bm: arb, bM: arb) -> arb:
    return ym - arb(1/16) + (2**bm-1)*Jm(ym, yM)+L(arb(1/16), bm) - 2*Q(yM/2+arb(1/32), bM)

# Case QJQ

def g_QJQ(xm: arb, xM: arb, ym: arb, yM: arb, b: arb) -> arb:
    return ym - xM + J(ym)*DJ(ym) - (2*Q((xm+yM)/2, b) - Q(xm, b))*DQ((xm+ym)/2, b)

# Case QJ

def g_QJ_1(xm: arb, xM: arb, ym: arb, yM: arb, b: arb, c: arb) -> arb:
    rv = (ym-xM)**(1/b-1)
    if xM < Jconst.x0:
        rv += c**(1/b)*(2*Jm((xm+ym)/2, (xM+yM)/2)-Q(xM,b))**(1/b-1)*DJ((xM+yM)/2)
    else:
        rv -= c**(1/b)*(2*JM((xm+ym)/2, (xM+yM)/2)-Q(xm,b))**(1/b-1)*absDJM((xm+ym)/2, (xM+yM)/2)
    rv -= c**(1/b)*(2*JM((xm+ym)/2, (xM+yM)/2)-Q(xm,b))**(1/b-1)*DQ(xm, b)
    return rv

def g_QJ_2(xm: arb, xM: arb, ym: arb, yM: arb) -> arb:
    return ((ym - xM)**2 + J(yM)**2)**.5 + Q(xm,arb(.5)) - 2*JM((xm+ym)/2, (xM+yM)/2)

# Poincare

def g_P_1(x: arb, b: arb=b1) -> arb:
    return 2**(-2*b)*(arb.log(1/x)/arb.log(arb(2)))**.5+2**(-2*b)*2**arb(.5)*arb(.77)*arb.log(Jconst.w0/x)**.5-2

def g_P_1_at1_32() -> arb:
    return g_P_1(arb(1/32)) > 0

b0P = .5+3*2**(-12)
def g_P_2(xm: arb, xM: arb, b: arb=b0P) -> arb:
    return 2**(-2*b)*(L(xm, arb(.5))+J(1-xm))-2*xM*(1-xm)

def g_P_3(xm: arb, xM: arb, b: arb=b0P) -> arb:
    rv = -.5*DQ(xm, arb(.5))
    if 1-xm < Jconst.x0:
        rv += 2**(-2*b)*DJ(1-xm)
    else:
        rv -= .5*absDJM(1-xM, 1-xm)
    rv += xm**(2*b-1)*(1-xM)
    rv -= xM
    rv += (1-xM)**(2*b)
    rv -= 2*b*xM
    return rv

# Auxiliary

def g_JL(xm: arb, xM: arb, b: arb) -> arb:
    return 2/JM(xm, xM)-b*arb.log(arb(2))**(-b)*(1-xM)**(-1)*(1-b+arb.log(1/(1-xM)))*arb.log(1/(1-xm))**(-2+b)

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
    Verify that given lower bound function is positive using partitioning on a given rectangle or interval and output result.
    '''
    if verbose: 
        msg = g.__name__
        if len(args) > 0:
            msg += " with "+", ".join(["%s=%f"%(k,float(v)) for k,v in args.items()])
        log(msg + ": ", end="")
    G = lambda *p: g(*p, **args)
    if y == None: 
        success, part = part_intvl(G, x, maxDepth=maxDepth) 
        if success:
            cmt = "%d intervals, min. val = %s"%(len(part)-1, min_val_intvl(G, part))
            Output.get_instance().write_part(g.__name__+tag, part, cmt)
            if verbose: 
                log(FMT_PASS%"ok", indent=0)
                log("   %s"%cmt)
    else: 
        success, part = part_rect(G, x, y, maxDepth=maxDepth)
        if success:
            cmt = "%d rectangles, min. val = %s"%(len(part), min_val_rect(G, part))
            Output.get_instance().write_part(g.__name__+tag, part, cmt)
            if verbose:
                log(FMT_PASS%"ok", indent=0)
                log("   %s"%cmt)
    if not success and verbose:
        log(FMT_FAIL%"fail", indent=0)
        log("   at %s"%part)
    return success, part

def batch_verify(label, methods, verbose=1):
    if verbose:log(label)
    Log.lvl = 1
    for i, m in enumerate(methods):
        m()
    Log.lvl = 0

def verify_all(b0=b0, c0=c0):
    Output.get_instance().write("# Partition data for beta0=%s, c0=%s\n\n"%(repr(b0), repr(c0)))
    if not (arb(.5) <= b0 <= 1):
        err("b0 must lie in [0.5, 1]")
        return
    if not c0 <= 1 and c0 > 0:
        err("c0 must lie in (0,1]")
        return
    if b0 > b1:
        warn("b0>%f: running only case J"%float(b1))
    batch_verify("case J", [
        lambda: verify_positive(g_J_1, (arb(1/2), arb(5/8)), (arb(0), arb(3/16)), b=b0, c=arb(1)),
        lambda: verify_positive(g_J_1, (arb(1/2), arb(5/8)), (arb(0), arb(3/16)), tag="h", b=arb(.5), c=c0),
        lambda: verify_positive(g_J_2, (arb(1/2), arb(9/16)), (arb(11/16), arb(1)))
    ])
    if b0 > b1: return
    batch_verify("case Q", [
        lambda: verify_positive(g_Q_1, (arb(0), arb(1/4)), (arb(1/4), arb(1/2)), b=b0),
        lambda: verify_positive(g_Q_1_y1_4, (arb(1/4), arb(1/2)), b=arb(.5)),
        lambda: verify_positive(g_Q_2, (arb(1/4), arb(1/2)), (arb(1/4), arb(1/2)), b=b0),
        lambda: verify_positive(g_Q_2, (arb(1/4), arb(1/2)), (arb(1/4), arb(1/2)), tag="h", b=arb(.5))
    ])
    batch_verify("case LJQ", [
        lambda: verify_positive(g_LJQ_1, (arb(1/16), arb(1/4)), (arb(1/2), arb(3/4)), b=b0),
        lambda: verify_positive(g_LJQ_1, (arb(1/16), arb(1/4)), (arb(1/2), arb(3/4)), tag="h", b=arb(.5)),
        lambda: verify_positive(g_LJQ_2, (arb(1/2), arb(3/4)), (arb(1/2), arb(1)))
    ])
    batch_verify("case QJQ", [
        lambda: verify_positive(g_QJQ, (arb(1/4), arb(1/2)), (arb(1/2), arb(3/4)), b=b1),
        lambda: verify_positive(g_QJQ, (arb(1/4), arb(1/2)), (arb(1/2), arb(3/4)), tag="h", b=arb(.5))
    ])
    batch_verify("case QJ", [
        lambda: verify_positive(g_QJ_1, (arb(1/4), arb(1/2)), (arb(1/2), arb(5/8)), b=b0, c=arb(1)),
        lambda: verify_positive(g_QJ_1, (arb(1/4), arb(1/2)), (arb(1/2), arb(5/8)), tag="h", b=arb(.5), c=c0),
        lambda: verify_positive(g_QJ_2, (arb(1/4), arb(1/2)), (arb(5/8), arb(1)))
    ])
    if (b0 > b0P):
        warn("b0>%f: skipping Poincare"%float(b0P))
        return
    batch_verify("Poincare", [
        lambda: verify(g_P_1_at1_32),
        lambda: verify_positive(g_P_2, (arb(1/32), arb(1/4))),
        lambda: verify_positive(g_P_3, (arb(1/4), arb(1/2)))
    ])
    batch_verify("Auxiliary", [
        lambda: verify_positive(g_JL, (arb(1/2), arb(7/8)), b=b0)
    ])


if __name__ == '__main__':
    verify_all()
