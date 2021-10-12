#include "polynomial.h"

quintuple_poly::quintuple_poly() {}
quintuple_poly::quintuple_poly(const F &aa, const F &bb, const F &cc, const F &dd, const F &ee, const F &ff)
{
    a = aa;
    b = bb;
    c = cc;
    d = dd;
    e = ee;
    f = ff;
}

quintuple_poly quintuple_poly::operator + (const quintuple_poly &x) const
{
    return quintuple_poly(a + x.a, b + x.b, c + x.c, d + x.d, e + x.e, f + x.f);
}

F quintuple_poly::eval(const F &x) const
{
    return (((((a * x) + b) * x + c) * x + d) * x + e) * x + f;
}

quadruple_poly::quadruple_poly() {}
quadruple_poly::quadruple_poly(const F &aa, const F &bb, const F &cc, const F &dd, const F &ee)
{
    a = aa;
    b = bb;
    c = cc;
    d = dd;
    e = ee;
}

quadruple_poly quadruple_poly::operator + (const quadruple_poly &x) const
{
    return quadruple_poly(a + x.a, b + x.b, c + x.c, d + x.d, e + x.e);
}

F quadruple_poly::eval(const F &x) const
{
    return ((((a * x) + b) * x + c) * x + d) * x + e;
}

cubic_poly::cubic_poly() {}
cubic_poly::cubic_poly(const F &aa, const F &bb, const F &cc, const F &dd)
{
    a = aa;
    b = bb;
    c = cc;
    d = dd;
}

cubic_poly cubic_poly::operator + (const cubic_poly &x) const
{
    return cubic_poly(a + x.a, b + x.b, c + x.c, d + x.d);
}

F cubic_poly::eval(const F &x) const
{
    return (((a * x) + b) * x + c) * x + d;
}


quadratic_poly::quadratic_poly() {}
quadratic_poly::quadratic_poly(const F &aa, const F &bb, const F &cc)
{
    a = aa;
    b = bb;
    c = cc;
}

quadratic_poly quadratic_poly::operator + (const quadratic_poly &x) const
{
    return quadratic_poly(a + x.a, b + x.b, c + x.c);
}

quadratic_poly quadratic_poly::operator+(const linear_poly &x) const
{
    return quadratic_poly(a, b + x.a, c + x.b);
}

cubic_poly quadratic_poly::operator * (const linear_poly &x) const
{
    return cubic_poly(a * x.a, a * x.b + b * x.a, b * x.b + c * x.a, c * x.b);
}

quadratic_poly quadratic_poly::operator*(const F &x) const
{
    return quadratic_poly(a * x, b * x, c * x);
}

F quadratic_poly::eval(const F &x) const
{
    return ((a * x) + b) * x + c;
}





linear_poly::linear_poly() { a = F_ZERO; b = F_ZERO; }
linear_poly::linear_poly(const F &aa, const F &bb)
{
    a = aa;
    b = bb;
}
linear_poly::linear_poly(const F &x)
{
    a = F_ZERO;
    b = x;
}

linear_poly linear_poly::operator + (const linear_poly &x) const
{
    return linear_poly(a + x.a, b + x.b);
}

quadratic_poly linear_poly::operator * (const linear_poly &x) const
{
    return quadratic_poly(a * x.a, a * x.b + b * x.a, b * x.b);
}

linear_poly linear_poly::operator*(const F &x) const
{
    return linear_poly(a * x, b * x);
}

F linear_poly::eval(const F &x) const
{
    return a * x + b;
}
