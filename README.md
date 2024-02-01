# fnt
Numerical Toolbox

## Background
Many numberical method libraries require the caller to pass the objective
function as a parameter and the library will call the objective function as it
sees fit.  In many cases the additional overhead needed to allow such
libraries to make the call can be quite cumbersome.

The goal of this library is to decouple the objective function from the
library.  The resulting API allows the caller to ask what the next input to
the objective function should be, call the objective function as it normally
would be called, then update the library with the value returned.
