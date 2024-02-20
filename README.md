[//]: <> (README.md)
[//]: <> (fnt: Numerical Toolbox)
[//]: <> ()
[//]: <> (Copyright [c] 2024 Bryan Franklin. All rights reserved.)

# FNT
FNT: Numerical Toolbox

## Background
Many numerical method libraries (e.g., [Dlib](http://dlib.net)) require the caller to pass the objective
function as a parameter and the library will call that objective function as it
sees fit.  In many cases the additional overhead needed to allow such
libraries to make the call can be quite cumbersome or require esoteric
knowledge of the language (e.g., function pointers or closures) to get right.

The goal of this library is to decouple the objective function from the
library.  The resulting API allows the caller to ask what the next input to
the objective function should be, call the objective function as it normally
would be called, then update the library with the value returned.

## Getting Started

### Linux or macOS

Build from source with `cmake`:
```text
$ git clone https://github.com/fra99le/fnt.git
$ cd fnt/src
$ cmake .
$ make
```

## Example

To help illustrate how this library differs from other libraries, below is an
example.  As with any library of this type there is some initialization and
teardown that needs to happen, but the actual numerical method (bisection in
this case) happens in the while loop.

Excluding setup and teardown, the execution of any method in this library works as follows:
1. Check if complete.
2. Get next input to evaluate.
3. Evaluate the objective function for provided input.
4. Pass objective function value back to method.
5. Loop.

```c
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "fnt.h"

double polynomial(double x) {
    // 3x^3 - 5x^2 - 6x + 5
    return 3*pow(x, 3.0) - 5*pow(x,2.0) - 6*x + 5;
}

int main() {

    void *fnt = NULL;
    fnt_init(&fnt, "./methods");

    /* load bisection method for a one dimensional input */
    fnt_set_method(fnt, "bisection", 1);

    /* display info about the method */
    fnt_info(fnt);

    /* set threshold for completion */
    double f_tol = 1e-5;
    double x_tol = 1e-5;
    fnt_hparam_set(fnt, "f_tol", &f_tol);
    fnt_hparam_set(fnt, "x_tol", &x_tol);

    /* place initial bounds for search */
    double x_0 = 2.0;
    double x_1 = 3.0;
    fnt_hparam_set(fnt, "upper", &x_0);
    fnt_hparam_set(fnt, "lower", &x_1);

    /* allocate one dimensional input vector for objective function */
    fnt_vect_t x;
    fnt_vect_calloc(&x, 1);

    /* MARK: End of setup */

    /* loop as long as method is not complete */
    while( fnt_done(fnt) == FNT_CONTINUE ) {

        /* get vector to try */
        if( fnt_next(fnt, &x) != FNT_SUCCESS ) { break; }

        /* call objective function */
        double fx = polynomial(FNT_VECT_ELEM(x, 0));

        /* update method */
        if( fnt_set_value(fnt, &x, fx) != FNT_SUCCESS ) { break; }
    }

    /* Get the root. */
    double x_root;
    if( fnt_root(fnt, &x, &x_root) == FNT_SUCCESS ) {
        fnt_vect_println(&x, "Root found at: x=", "%.3f");
    }

    /* MARK: Begin teardown */

    /* free input vector */
    fnt_vect_free(&x);

    /* free the method */
    fnt_free(&fnt);

    return 0;
}
```

## Why use this library?

By splitting iterations of numerical methods into two phases where the
caller is tasked with evaluating the objective function every iteration, the
caller maintains much more control over how the overall method is run.
This strategy also eliminates the need to create a wrapper function around
your objective function that conforms to any restrictions that other libraries
might place on its inputs.

For example, the function `polynomial` could have been written as:
```c
double polynomial(double x, int n, double *coefficients) {
    double sum = 0.0;
    for(int i=0; i<=n; ++i) {
        sum += coefficients[i] * pow(x, (double)i);
    }
    return sum;
}
```

If that were the case, then it wouldn't be compatible with a library that
expects a function that takes only a single `double` as its parameter, and would
require a wrapper function.
However, with this library one need only update the call
to the objective function to use this alternate polynomial function.

Thus the objective function call in the while loop:
```c
/* call objective function */
double fx = polynomial(FNT_VECT_ELEM(x, 0));
```

Becomes:
```c
/* call objective function */
double[] coeffs = [5, -6.0, -5.0, 3.0]; // 3x^3 - 5x^2 - 6x + 5
double fx = polynomial(FNT_VECT_ELEM(x, 0), 3, coeffs);
```

No other changes are required.

This is of course a contrived example, but in the real world it is often the
case that objective functions have many more inputs than just the parameters
being optimized.

[//]: <> (TODO: Add a section on methods available.)
[//]: <> (TODO: Add information on how to get a list of methods.)
