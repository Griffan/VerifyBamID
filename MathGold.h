#ifndef _MATHGOLD_H_
#define _MATHGOLD_H_

#include "MathConstant.h"
#include "MathVector.h"

// Minimizes functions of one variable in one dimension
class ScalarMinimizer
{
public:
    double(*func)(double);          // function to be minimized
    double a, b, c, min;
    double fa, fb, fc, fmin;

    ScalarMinimizer()
    {
        func = NULL;
    };
    virtual ~ScalarMinimizer() { }

    virtual double f(double x);

    void   Bracket(double a, double b);      // bracket a minimum near a and b
    virtual double Brent(double tol = TOL);  // return minimum, to precision TOL
    // result stored in min
};

class LineMinimizer : public ScalarMinimizer
// Minimizes f(P) along the line define by P = point + x * line
// Stores the best point (in point) along the line
// and the displacement from the original (in line)
{
private:
    bool         garbage;
public:
    VectorFunc * func;      // function to be minimized
    Vector       line, point, temp;

    LineMinimizer();
    LineMinimizer(VectorFunc & vfunc);
    LineMinimizer(double(*vfunc)(Vector & v));

    virtual ~LineMinimizer()
    {
        if (garbage) delete func;
    }

    virtual double f(double x);

    virtual double Brent(double tol = TOL);
};

#endif


