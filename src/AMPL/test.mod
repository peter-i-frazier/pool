set I;
set J;
set K;

param a {i in I};
param ita {I,K,J};

var x {K,J} binary;
var y {i in I};

maximize Total: sum {i in I} y[i] * a[i];

subject to star {i in I}: y[i] <= (prod {j in J} (sum {k in K} x[k,j] * ita[i,k,j]))/(1+prod {j in J} (sum {k in K} x[k,j] * ita[i,k,j]));
subject to xk {j in J}: sum {k in K} x[k,j] = 1;