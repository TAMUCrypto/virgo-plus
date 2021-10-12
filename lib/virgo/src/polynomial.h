#pragma once

#include <vector>
#include "fieldElement.hpp"

namespace virgo {
    class linear_poly;

//ax^3 + bx^2 + cx + d
    class cubic_poly {
    public:
        virgo::fieldElement a, b, c, d;

        cubic_poly();

        cubic_poly(const virgo::fieldElement &, const virgo::fieldElement &, const virgo::fieldElement &,
                   const virgo::fieldElement &);

        cubic_poly operator+(const cubic_poly &) const;

        virgo::fieldElement eval(const virgo::fieldElement &) const;
    };

//ax^2 + bx + c
    class quadratic_poly {
    public:
        virgo::fieldElement a, b, c;

        quadratic_poly();

        quadratic_poly(const virgo::fieldElement &, const virgo::fieldElement &, const virgo::fieldElement &);

        quadratic_poly operator+(const quadratic_poly &) const;

        quadratic_poly operator+(const linear_poly &) const;

        cubic_poly operator*(const linear_poly &) const;

        quadratic_poly operator*(const virgo::fieldElement &) const;

        virgo::fieldElement eval(const virgo::fieldElement &) const;
    };


//ax + b
    class linear_poly {
    public:
        virgo::fieldElement a, b;

        linear_poly();

        linear_poly(const virgo::fieldElement &, const virgo::fieldElement &);

        linear_poly(const virgo::fieldElement &);

        linear_poly operator+(const linear_poly &) const;

        quadratic_poly operator*(const linear_poly &) const;

        linear_poly operator*(const virgo::fieldElement &) const;

        virgo::fieldElement eval(const virgo::fieldElement &) const;
    };


//ax^4 + bx^3 + cx^2 + dx + e
    class quadruple_poly {
    public:
        virgo::fieldElement a, b, c, d, e;

        quadruple_poly();

        quadruple_poly(const virgo::fieldElement &, const virgo::fieldElement &, const virgo::fieldElement &,
                       const virgo::fieldElement &, const virgo::fieldElement &);

        quadruple_poly operator+(const quadruple_poly &) const;

        virgo::fieldElement eval(const virgo::fieldElement &) const;
    };

//ax^5 + bx^4 + cx^3 + dx^2 + ex + f
    class quintuple_poly {
    public:
        virgo::fieldElement a, b, c, d, e, f;

        quintuple_poly();

        quintuple_poly(const virgo::fieldElement &, const virgo::fieldElement &, const virgo::fieldElement &,
                       const virgo::fieldElement &, const virgo::fieldElement &, const virgo::fieldElement &);

        quintuple_poly operator+(const quintuple_poly &) const;

        virgo::fieldElement eval(const virgo::fieldElement &) const;
    };
}