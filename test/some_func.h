
template < typename Real, int DIM >
void rosenbrock(Real *x, Real *y, Real *c) {
  y[0] = 0;
  for(int i=0;i!=DIM-1;++i){
    y[0] += 100 * pow((x[i+1] - x[i] * x[i]), 2) +
          pow((1 - x[i]), 2);
  }
}
template < typename Real>
void G9(Real *x, Real *y, Real *c) {
  c[0] = -(2 * pow(x[0], 2) + 3 * pow(x[1], 4) + x[2] + 4 * pow(x[3], 2) +
           5 * x[4] - 127);
  c[1] = -(7 * x[0] + 3 * x[1] + 10 * pow(x[2], 2) + x[3] - x[4] - 282);
  c[2] = -(23 * x[0] + pow(x[1], 2) + 6 * pow(x[5], 2) - 8 * x[6] - 196);
  c[3] = -(4 * pow(x[0], 2) + pow(x[1], 2) - 3 * x[0] * x[1] +
           2 * pow(x[2], 2) + 5 * x[5] - 11 * x[6]);
  y[0] = pow((x[0] - 10), 2) + 5 * pow((x[1] - 12), 2) + pow(x[2], 4) +
         3 * pow((x[3] - 11), 2) + 10 * pow(x[4], 6) + 7 * pow(x[5], 2) +
         pow(x[6], 4) - 4 * x[5] * x[6] - 10 * x[5] - 8 * x[6];
}
template < typename Real>
void test_fun(Real *x, Real *y, Real *c){
  y[0] = exp(-pow(x[0] - 2, 2)) + exp(-pow(x[0] - 6, 2)/10) + 1 / (x[0]*x[0] + 1);
}
