/* Comparator only. This is a convex interpolation, not a rank square root. */
mix(th,z)=(1-th)*exp(log(sexp(1))*z)+th*sexp(z);
mixjetset(th,x,r,N)={Pmix=pentcauchy(t->mix(th,t),x,r,N);0;};
