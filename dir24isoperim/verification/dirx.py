# (c) 2026 Joris Roos <jroos.math@gmail.com>
# pylint: disable=too-many-arguments,invalid-name,no-name-in-module

''' Verification of numerical inequalities in DIRX26 '''

from flint import arb

from ..general import Jconst, L, Q, DQ, bobkovI, PhiInv, \
                    wtox, batch_verify, verify_positive

from ..util import Output

#
#  Define functions for verification
#

def Jw(x: arb, w: arb) -> arb:
    '''Rescaled Gaussian isoperimetric profile'''
    return .5*bobkovI((1-arb(x))/arb(w))/bobkovI(1/(2*w))

def J(x: arb) -> arb:
    '''Specific rescaling that we use'''
    return Jw(x, Jconst.w1)

def DJ(x: arb) -> arb:
    '''Derivative of J'''
    return .5*1/Jconst.w1*1/bobkovI(1/(2*Jconst.w1))*PhiInv((1-x)/Jconst.w1)

def gamma() -> arb:
    '''Value of $\gamma_w$ for $w=w_1$'''
    return 1/(4*Jconst.w1**2*bobkovI(1/(2*Jconst.w1))**2)

def Jm(xm: arb, xM: arb) -> arb:
    '''Lower bound for J'''
    return arb.min(J(xm), J(xM))

def JM(xm: arb, xM: arb) -> arb:
    '''Upper bound for J'''
    if xM < Jconst.x1:
        return J(xM)
    if xm > Jconst.x1:
        return J(xm)
    return J(Jconst.x1)

def absDJm(xm: arb, xM: arb) -> arb:
    '''Lower bound for |J'|'''
    if xM < Jconst.x1:
        return DJ(xM)
    if xm > Jconst.x1:
        return -DJ(xm)
    return arb(0)

def absDJM(xm: arb, xM: arb) -> arb:
    '''Upper bound for |J'|'''
    return arb.max(abs(DJ(xm)), abs(DJ(xM)))

def h_LJ_1(cm: arb, cM: arb) -> arb:
    '''Case LJ'''
    return 2*cm - .5 + (2**arb(.5)-1)*Jm(2*cm-.25,2*cM-.25) + L(arb(.25)) - 2*JM(cm, cM)

def h_LJQ_1(xm: arb, xM: arb, ym: arb, yM: arb) -> arb:
    '''Case LJQ'''
    l = ((ym-xM)**2 + Jm(ym, yM)**2)**.5
    r = ym - xM + (2**arb(.5)-1)*Jm(ym, yM)
    rv = arb.max(l, r)
    rv += L(xm, arb(.5))
    rv -= 2*Q((xM+yM)/2, arb(.5))
    return rv

def h_LJQ_2(ym: arb, yM: arb) -> arb:
    '''Case LJQ'''
    return ym - arb(1/16) + (2**arb(.5)-1)*Jm(ym, yM)+L(arb(1/16), arb(.5)) \
            - 2*Q(yM/2+arb(1/32), arb(.5))

def h_QJQ_1(xm: arb, xM: arb, ym: arb, yM: arb) -> arb:
    '''Case QJQ.1'''
    return ym - xM + J(yM)*DJ(yM) - (2*Q((xm+yM)/2) - Q(xm))*DQ((xm+ym)/2)

def h_QJQ_2(xm: arb, xM: arb, ym: arb, yM: arb) -> arb:
    '''Case QJQ.2'''
    return (ym - xM)**2 + Jm(ym, yM)**2 - (2*Q((xm+yM)/2)-Q(xm))**2

def h_QJ_1(xm: arb, xM: arb, ym: arb, yM: arb) -> arb:
    '''Case QJ'''
    rv = ym-xM
    if (xM+yM)/2 < Jconst.x1:
        rv += (2*Jm((xm+ym)/2, (xM+yM)/2)-Q(xM))*DJ((xM+yM)/2)
    else:
        rv -= (2*JM((xm+ym)/2, (xM+yM)/2)-Q(xm))*absDJM((xm+ym)/2, (xM+yM)/2)
    rv -= (2*JM((xm+ym)/2, (xM+yM)/2)-Q(xm))*DQ(xm)
    return rv

def h_QJ_2(xm: arb, xM: arb, ym: arb, yM: arb) -> arb:
    '''Case QJ'''
    return ((ym - xM)**2 + J(yM)**2)**.5 + Q(xm) - 2*JM((xm+ym)/2, (xM+yM)/2)

def h_P_1(xm: arb, xM: arb) -> arb:
    '''Poincare'''
    return .5*(L(xm, arb(.5))+J(1-xm)) - 2*xM*(1-xm)

def h_P_2(xm: arb, xM: arb) -> arb:
    '''Poincare'''
    rv = -.5*DQ(xm)
    rv += .5*DJ(1-xm)
    rv += -4*xM + 2
    return rv

def verify_all():
    '''Verify all claims in the paper.'''
    Output.get_instance().write(f"# Partition data for DIRX26\n\n")
    batch_verify("case LJQ", [
        lambda: verify_positive(h_LJQ_1, (arb(1/16), arb(1/4)), (arb(1/2), arb(3/4))),
        lambda: verify_positive(h_LJQ_2, (arb(1/2), arb(3/4)))
    ])
    batch_verify("case LJ", [
        lambda: verify_positive(h_LJ_1, (arb(1/2), arb(5/8)))
    ])
    batch_verify("case QJQ", [
        lambda: verify_positive(h_QJQ_1, (arb(1/4), arb(1/2)), (arb(1/2), arb(33/64))),
        lambda: verify_positive(h_QJQ_2, (arb(1/4), arb(1/2)), (arb(33/64), arb(3/4)))
    ])
    batch_verify("case QJ", [
        lambda: verify_positive(h_QJ_1, (arb(1/4), arb(1/2)), (arb(1/2), arb(5/8))),
        lambda: verify_positive(h_QJ_2, (arb(1/4), arb(1/2)), (arb(5/8), arb(1)))
    ])
    batch_verify("Poincare", [
        lambda: verify_positive(h_P_1, (arb(1/64), arb(1/4))),
        lambda: verify_positive(h_P_2, (arb(1/4), arb(1/2)))
    ])

Jconst.w1 = arb(29/32)
Jconst.x1 = wtox(Jconst.w1)
