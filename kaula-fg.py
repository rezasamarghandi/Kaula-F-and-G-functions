from sympy import *
import math as math


def ffun(l, m, p):
    c, s, t = symbols('c s t', integer=True)
    inc = symbols('inc')
    k = math.trunc((l - m) / 2)
    y = Sum(((factorial(2 * l - 2 * t)) / (
            factorial(t) * factorial(l - t) * factorial(l - m - 2 * t) * 2 ** (2 * l - 2 * t)) * (
                 (sin(inc) ** (l - m - 2 * t)))) * Sum(binomial(m, s) * ((cos(inc)) ** s) * Sum(
        binomial(l - m - 2 * t + s, c) * binomial(m - s, p - t - c) * (-1) ** (c - k), (c, 0, 5)), (s, 0, m)),
            (t, 0, min(p, k))).doit()
    return y


def gfun(l, p, q):
    e, d, k, s, t, = symbols('e d k s t', integer=True)
    if p > l / 2:
        pp = l - p
        qp = -q
    else:
        pp = p
        qp = q

    beta = e / (1 + sqrt(1 - e ** 2))

    qqq = 0
    for k in range(0, 10 + 1):
        if qp < 0:
            hp = k
            hq = k - qp
        else:
            hp = k + qp
            hq = k
        plpq = 0
        qlpq = 0
        for rp in range(0, hp + 1):
            plpq = plpq + binomial(2 * pp - 2 * l, hp - rp) * (((-1) ** rp) / factorial(rp)) * (
                    (l - 2 * pp + qp) * e / (2 * beta)) ** rp
        for rq in range(0, hq + 1):
            qlpq = qlpq + binomial(-2 * pp, hq - rq) / factorial(rq) * (
                    (l - 2 * pp + qp) * e / (2 * beta)) ** rq
        qqq = qqq + plpq * qlpq * beta ** (2 * k)

    z = (-1) ** abs(q) * (1 + beta ** 2) ** l * beta ** abs(q) * qqq

    return z



l = int(input("l = "))
m = int(input("m = "))
p = int(input("p = "))
q = int(input("q = "))

f=ffun(l,m,p)
g=gfun(l,p,q)
print("F = ",f)
print("G = ",g)


