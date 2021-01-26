/* Brazil Salt Policy model: a decision support tool for primary prevention of NCDs
   Copyright (C) 2019 Chris Kypridemos

   Brazil Sodium Policy model is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have a copy of the GNU General Public License
   along with this program; if not, see <http://www.gnu.org/licenses/>
   or write to the Free Software Foundation, Inc., 51 Franklin Street,
   Fifth Floor, Boston, MA 02110-1301  USA.*/

#include <omp.h>
#include <Rcpp.h>
#include <string>
#include <math.h>
#include <Rmath.h>
#include "pcg_random.h" // from http://www.pcg-random.org
#include "randutils.h" //from https://gist.githubusercontent.com/imneme/540829265469e673d045/raw/8486a610a954a8248c12485fb4cfc390a5f5f854/randutils.hpp and for explanation http://www.pcg-random.org/posts/simple-portable-cpp-seed-entropy.html
#include <algorithm>    // std::max

using namespace Rcpp;
using namespace std;

// [[Rcpp::plugins(openmp)]]


//quantile implementation for default R method (type 7)
//' @export
// [[Rcpp::export]]
NumericVector fquantile(NumericVector x,
                        NumericVector probs,
                        bool na_rm = true) {
  if (all(is_na(x)))
  {
    NumericVector out(probs.size(), NA_REAL);
    return(out);
  }
  if (na_rm) x = na_omit(x);
  const int n = x.size();
  NumericVector out(probs.size());
  IntegerVector ii(probs.size());
  NumericVector h(probs.size());
  NumericVector index = 1 + (n - 1) * probs;
  NumericVector lo = floor(index); //floor
  //ceiling
  NumericVector hi = ceiling(index);
  for(int i = 0; i < probs.size(); i++)
  {//catch corner case when index element is int and ceiling = floor
    h[i] = index[i] - lo[i];
  }
  sort(x.begin(), x.end());
  out = x[as<IntegerVector>(lo) - 1];
  if (all(is_na(ii)))
  {
    return(out);
  } else
  {
    x = x[as<IntegerVector>(hi) - 1];
    for(int i = 0; i < probs.size(); i++)
    {
      out[i] = (1 - h[i]) * out[i] + h[i] * x[i];
    }
    return(out);
  }
}

//' @export
// [[Rcpp::export]]
List fquantile_byid(NumericVector x,
                    NumericVector q,
                    StringVector id,
                    bool rounding = false,
                    bool na_rm = true) {
  // Need to be sorted by id
  const int n = x.size();
  const int m = unique(id).size();
  NumericMatrix z(m, q.size());
  StringVector id_nam(m);
  int start = 0;
  int counter = 0;
  int end = 0;
  int counter_row = 0;

  for (int i = 1; i < n; i++) { // start from 2nd element
    if (id[i] == id[i-1]) counter++;
    else
    {
      start = i - counter - 1;
      end = i - 1;
      counter = 0;
      if (rounding) z.row(counter_row) = round(fquantile(x[seq(start, end)], q, na_rm), 0);
      else z.row(counter_row) = fquantile(x[seq(start, end)], q, na_rm);
      id_nam[counter_row] = id[end];
      counter_row++;
    }
  }
  // take care the last group
    if (rounding) z.row(counter_row) = round(fquantile(x[seq(n - 1 - counter, n - 1)], q, na_rm), 0);
    else z.row(counter_row) = fquantile(x[seq(n - 1 - counter, n - 1)], q, na_rm);
    id_nam[counter_row] = id[n - 1];

  // return(z);
  const int tt = 1 + q.size();
  List outputList(tt);
  outputList[0] = id_nam;
  for (int i = 1; i < tt; i++) {
    outputList[i] = z(_, i - 1);
  }

  return outputList;
}

//' @export
// [[Rcpp::export]]
NumericVector fbound(const NumericVector &x, NumericVector &a, NumericVector &b) {
  const int n = x.size();
  NumericVector out(n);
  for(int i = 0; i < n; i++) {
    if (NumericVector::is_na(x[i])) out[i] = NA_REAL;
    else
    {
      if (a[i] > b[i]) {double c = a[i]; a[i] = b[i]; b[i] = c;}; // ensure a < b

      if (x[i] < a[i]) out[i] = a[i];
      else if (x[i] > b[i]) out[i] = b[i];
      else out[i] = x[i];
    }
  }
  return out;
}

//' @export
// [[Rcpp::export]]
IntegerVector fbound_int(const IntegerVector &x, IntegerVector &a, IntegerVector &b) {
  const int n = x.size();
  IntegerVector out(n);
  for(int i = 0; i < n; i++) {
    if (IntegerVector::is_na(x[i])) out[i] = NA_INTEGER;
    else
    {
      if (a[i] > b[i]) {double c = a[i]; a[i] = b[i]; b[i] = c;}; // ensure a < b

      if (x[i] < a[i]) out[i] = a[i];
      else if (x[i] > b[i]) out[i] = b[i];
      else out[i] = x[i];
    }
  }
  return out;
}

//' @export
// [[Rcpp::export]]
NumericVector fbound_inplace(NumericVector x, NumericVector &a, NumericVector &b) {
  const int n = x.size();
  for(int i = 0; i < n; i++) {
    if (!NumericVector::is_na(x[i]))
    {
      if (a[i] > b[i]) {double c = a[i]; a[i] = b[i]; b[i] = c;}; // ensure a < b

      if (x[i] < a[i]) x[i] = a[i];
      else if (x[i] > b[i]) x[i] = b[i];
    }
  }

  return x;
}

//' @export
// [[Rcpp::export]]
IntegerVector fbound_inplace_int(IntegerVector x, IntegerVector &a, IntegerVector &b) {
  const int n = x.size();
  for(int i = 0; i < n; i++) {
    if (!IntegerVector::is_na(x[i]))
    {
       if (a[i] > b[i]) {double c = a[i]; a[i] = b[i]; b[i] = c;}; // ensure a < b

       if (x[i] < a[i]) x[i] = a[i];
       else if (x[i] > b[i]) x[i] = b[i];
    }
  }

  return x;
}

//' @export
// [[Rcpp::export]]
LogicalVector fequal(NumericVector x, double tol) {
  NumericVector y = na_omit(x);
  const int n = y.size();
  for (int i = 0; i < n; ++i) {
    if (y[i] - y[0] > tol || y[0] - y[i] > tol)
      return wrap(false);
  }

  return wrap(true);
}

//' @export
// [[Rcpp::export]]
NumericVector fnormalise(const NumericVector& x) { // between 0, 1
  const int n = x.size();
  const double minx = min(na_omit(x));
  const double maxx = max(na_omit(x));
  NumericVector out(n);
  for (int i = 0; i < n; i++) {
    out[i] = (x[i] - minx) / (maxx - minx);
  }

  return out;
}

//' @export
// [[Rcpp::export]]
NumericVector fcompress(const NumericVector& x, float limit) {
  const int n = x.size();
  NumericVector out = clone(x);
  float sum_out = sum(na_omit(out));
  int counter = 0;
  while(pow(sum_out - limit, 2) > 0.1)
  {
    if (counter == 1000) break;
    sum_out = sum(na_omit(out));
    for(int i = 0; i < n; i++)
    {
      out[i] *= limit/sum_out;
      if (out[i] < 0.0) out[i] = 0.0;
      else if (out[i] > 1.0) out[i] = 1.0;
    }
    counter++;
  }
  return out; // transforms input vector inplace
}

//' @export
// [[Rcpp::export]]
NumericVector roll_mean_left(const NumericVector& dat,  int window) {
  const int n = dat.size();
  window =  ((window > n ) ? n : window);

  NumericVector out(n);
  double summed = 0.0;
  for (int i=0; i < window; i++) {
    summed += dat[i];
    out[i] = summed/(i+1);
  }
  for (int i=window; i < n; i++) {
    summed += (dat[i] - dat[i-window]);
    out[i] = summed / window;
  }
  return out;
}

//' @export
// [[Rcpp::export]]
NumericVector shift_byidNum(const NumericVector& x, int lag,
                            float replace, IntegerVector id) {
  // id should be sorted and same length as x
  const int n = x.size();
  NumericVector out(n);
  if (id[0] == id[lag]) {
    for (int i = 0; i < lag; i++) out[i] = replace;
    for (int i = lag; i < n; i++)
    {
      if (id[i] == id[i-lag]) {
        out[i] = x[i-lag];
      } else {
        out[i] = replace;
      }
    }
    return out;
  } else {
    stop("This function does not work because number of ids is too small!");
  }
}

//' @export
// [[Rcpp::export]]
IntegerVector shift_byidInt(const IntegerVector& x, int lag,
                            int replace, IntegerVector id) {
  // id should be sorted and same length as x
  const int n = x.size();
  IntegerVector out(n);
  if (id[0] == id[lag]) {
    for (int i = 0; i < lag; i++) out[i] = replace;
    for (int i = lag; i < n; i++)
    {
      if (id[i] == id[i-lag]) {
        out[i] = x[i-lag];
      } else {
        out[i] = replace;
      }
    }
    return out;
  } else {
    stop("This function does not work because number of ids is too small!");
  }
}

//' @export
// [[Rcpp::export]]
StringVector  shift_byidStr(const CharacterVector& x, int lag,
                            std::string replace, IntegerVector id) {
  // id should be sorted and same length as x
  const int n = x.size();
  StringVector out(n);
  if (id[0] == id[lag]) {
    for (int i = 0; i < lag; i++) out[i] = replace;
    for (int i = lag; i < n; i++)
    {
      if (id[i] == id[i-lag]) {
        out[i] = x[i-lag];
      } else {
        out[i] = replace;
      }
    }
    return out;
  } else {
    stop("This function does not work because number of ids too small!");
  }
}

//' @export
// [[Rcpp::export]]
IntegerVector sort_byidInt(const IntegerVector& x, IntegerVector id) {
  // Need to be sorted by id
  const int n = x.size();
  IntegerVector out(n);
  int start = 0;
  int counter = 0;
  int end = 0;
  out(0) = x(0); // because ext loop starts from 2nd element
  for (int i = 1; i < n; i++) { // start from 2nd element
    if (id[i] == id[i-1]) counter++;
    else
    {
      start = i - counter - 1;
      end = i - 1;
      counter = 0;
      IntegerVector tmp = x[seq(start, end)];
      out[seq(start, end)] = tmp.sort();
    }
  }
  // take care the last group
  IntegerVector tmp =  x[seq(n-counter-1, n-1)];
  out[seq(n-counter-1, n-1)] = tmp.sort();
  return(out);
}

//' @export
// [[Rcpp::export]]
NumericVector sort_byidNum(const NumericVector& x, IntegerVector id) {
  // Need to be sorted by id
  const int n = x.size();
  NumericVector out(n);
  int start = 0;
  int counter = 0;
  int end = 0;
  out(0) = x(0); // because ext loop starts from 2nd element
  for (int i = 1; i < n; i++) { // start from 2nd element
    if (id[i] == id[i-1]) counter++;
    else
    {
      start = i - counter - 1;
      end = i - 1;
      counter = 0;
      NumericVector tmp = x[seq(start, end)];
      out[seq(start, end)] = tmp.sort();
    }
  }
  // take care the last group
  NumericVector tmp =  x[seq(n-counter-1, n-1)];
  out[seq(n-counter-1, n-1)] = tmp.sort();
  return(out);
}

//' @export
// [[Rcpp::export]]
NumericMatrix df2mat_numeric(const DataFrame& x) {
  int xsize=x.size();
  NumericMatrix y(x.nrows(),xsize);
  for(int i=0; i<xsize; i++){
    y(_,i) = NumericVector(x[i]);
  }
  y.attr("dimnames") = List::create( R_NilValue, x.attr("names") ) ;
  return(y);
}

//' @export
// [[Rcpp::export]]
IntegerMatrix df2mat_integer(const DataFrame& x) {
  int xsize=x.size();
  IntegerMatrix y(x.nrows(),xsize);
  for(int i=0; i<xsize; i++){
    y(_,i) = IntegerVector(x[i]);
  }
  y.attr("dimnames") = List::create( R_NilValue, x.attr("names") ) ;
  return(y);
}


//' @export
// [[Rcpp::export]]
IntegerVector collect_output(
    const IntegerVector&  m1,// other mort
    const IntegerVector&  m2,// chd mort
    const IntegerVector&  m3,// stroke mort
    const NumericVector&  d1,// dice to resolve discrete time mortality bias
    const IntegerVector&  p1,// htn incid
    const IntegerVector&  p2,// chd incid
    const IntegerVector&  p3,// stroke incid
    const IntegerVector&  id
) {
  const int n = id.size();
  IntegerVector lifecourse(n);

  int aux[] = {0, 0, 0};

  // Resolve 1st row
  aux[1]   = m1[0] + m2[0] + m3[0]; // if dead
  aux[2]   = p1[0] + p2[0] + p3[0]; //if not dead estimate incidence

  // resolve discrete time mortality bias
  if (aux[1] == -123)
  {
  	if (d1[0] < 1/3) aux[1] = -3;
  	else if (d1[0] >= 2/3) aux[1] = -20;
  	else aux[1] = -100;
  };
  if (aux[1] == -120)
  {
  	if (d1[0] < 1/2) aux[1] = -100;
  	else aux[1] = -20;
  };
  if (aux[1] == -103)
  {
  	if (d1[0] < 1/2) aux[1] = -100;
  	else aux[1] = -3;
  };
  if (aux[1] == -23)
  {
  	if (d1[0] < 1/2) aux[1] = -3;
  	else aux[1] = -20;
  };

  lifecourse[0]         = ((aux[1] < 0) ? aux[1] : aux[2]);

  for (int i = 1; i < n; i++) //start loop over id (from 2nd row)
  {
    aux[0] = lifecourse[i-1]; // not accurate id id-1 != id
    aux[1] = m1[i] + m2[i] + m3[i]; // if dead
    aux[2] = p1[i] + p2[i] + p3[i]; //if not dead estimate incidence

    // resolve discrete time mortality bias
    if (aux[1] == -123)
    {
    	if (d1[i] < 1/3) aux[1] = -3;
    	else if (d1[i] >= 2/3) aux[1] = -20;
    	else aux[1] = -100;
    };
    if (aux[1] == -120)
    {
    	if (d1[i] < 1/2) aux[1] = -100;
    	else aux[1] = -20;
    };
    if (aux[1] == -103)
    {
    	if (d1[i] < 1/2) aux[1] = -100;
    	else aux[1] = -3;
    };
    if (aux[1] == -23)
    {
    	if (d1[i] < 1/2) aux[1] = -3;
    	else aux[1] = -20;
    };

    if (id[i] == id[i-1]) // if same person
    {
      if (aux[0] < 0) lifecourse[i] = NA_INTEGER; ///those dead remain dead NA after death. NA_integer < x always true
      else lifecourse[i] = ((aux[1] < 0) ? aux[1] : aux[2]);
    }
    else // if different person
    {
      lifecourse[i] = ((aux[1] < 0) ? aux[1] : aux[2]);
    }
  }

  return lifecourse;
}


// optimise BCPEo (and BCPEo trunc)
int sign(const double& x) {
  return (x > 0) - (x < 0);
}

double F_T(const double& t, const double& tau) {
  double log_c = 0.5 * (-(2.0/tau) * log(2.0) + lgamma(1.0/tau) - lgamma(3.0/tau));
  double c = exp(log_c);
  double s = 0.5 * (pow((abs(t/c)), tau));
  double F_s = R::pgamma(s, 1/tau, 1.0, true, false);
  double cdf = 0.5 * (1.0 + F_s * sign(t));
  return cdf;
}


double my_pBCPEo_scalar(const double& q,
                        const double& mu,
                        const double& sigma,
                        const double& nu,
                        const double& tau)
{
  double z = 0.0;
  if (nu == 0.0) z = log(q/mu)/sigma; else z = ((pow((q/mu),nu) - 1.0)/(nu * sigma));


  double FYy1 = F_T(z, tau);

  double FYy2 = 0.0;
  if (nu > 0.0)
  {
    FYy2 = F_T(-1.0/(sigma * abs(nu)), tau);
  }
  double FYy3 = F_T(1/(sigma * abs(nu)), tau);
  return (FYy1 - FYy2)/FYy3;
}

double q_T(const double& p, const double& tau) {
  double log_c = 0.5 * (-(2.0/tau) * log(2.0) + lgamma(1.0/tau) - lgamma(3.0/tau));
  double c = exp(log_c);
  double s = R::qgamma((2.0 * p - 1.0) * sign(p - 0.5), (1.0/tau), 1.0, true, false);
  return sign(p - 0.5) * pow((2.0 * s),(1.0/tau)) * c;
}

double my_qBCPEo_scalar(const double& p,
                        const double& mu,
                        const double& sigma,
                        const double& nu,
                        const double& tau)
{
  double za = 0.0;
  if (nu < 0)
  {
    za = q_T(p * F_T(1.0/(sigma * abs(nu)), tau), tau);
  } else if (nu == 0)
  {
    za = q_T(p, tau);
  } else // if nu > 0
  {
    za = q_T((1.0 - (1.0 - p) * F_T(1.0/(sigma * abs(nu)), tau)), tau);
  }
  double ya = 0.0;
  if (nu == 0.0) ya = mu * exp(sigma * za); else ya = mu * pow((nu * sigma * za + 1.0),(1.0/nu));

  return ya;
}

double my_qBCPEo_trunc_scalar(const double& p,
                              const double& mu,
                              const double& sigma,
                              const double& nu,
                              const double& tau,
                              const double& lower_lim,
                              const double& upper_lim) {
  double pp1 = my_pBCPEo_scalar(lower_lim, mu, sigma, nu, tau);
  double pp2 = my_pBCPEo_scalar(upper_lim, mu, sigma, nu, tau);
  return my_qBCPEo_scalar((p * (pp2 - pp1) + pp1),  mu, sigma, nu, tau);
}

//' @export
// [[Rcpp::export]]
NumericVector my_pBCPEo(const NumericVector& q,
                              const NumericVector& mu,
                              const NumericVector& sigma,
                              const NumericVector& nu,
                              const NumericVector& tau,
                              const int& n_cpu) {
  // TODO check all have same length and are > 0 etc.
  const int n = q.length();
  NumericVector out(n);

  omp_set_num_threads(n_cpu); // Use n_cpu threads for all

  // consecutive parallel regions
#pragma omp parallel for default(shared)
  for (int i = 0; i < n; i++)
  {
    out[i] = my_pBCPEo_scalar(q[i], mu[i], sigma[i], nu[i], tau[i]);
  }

  return out;
}

//' @export
// [[Rcpp::export]]
NumericVector my_qBCPEo(const NumericVector& p,
                                    const NumericVector& mu,
                                    const NumericVector& sigma,
                                    const NumericVector& nu,
                                    const NumericVector& tau,
                                    const int& n_cpu) {
  // TODO check all have same length and are > 0 etc.
  const int n = p.length();
  NumericVector out(n);

  omp_set_num_threads(n_cpu); // Use n_cpu threads for all

  // consecutive parallel regions
#pragma omp parallel for default(shared)
  for (int i = 0; i < n; i++)
  {
    out[i] = my_qBCPEo_scalar(p[i], mu[i], sigma[i], nu[i], tau[i]);
  }

  return out;
}

//' @export
// [[Rcpp::export]]
NumericVector my_qBCPEo_trunc(const NumericVector& p,
                                    const NumericVector& mu,
                                    const NumericVector& sigma,
                                    const NumericVector& nu,
                                    const NumericVector& tau,
                                    const double& lower_lim,
                                    const double& upper_lim,
                                    const int& n_cpu) {
  // TODO check all have same length and are > 0 etc.
  const int n = p.length();
  NumericVector out(n);

  omp_set_num_threads(n_cpu); // Use n_cpu threads for all

  // consecutive parallel regions
#pragma omp parallel for default(shared)
  for (int i = 0; i < n; i++)
  {
    out[i] = my_qBCPEo_trunc_scalar(p[i], mu[i], sigma[i], nu[i], tau[i], lower_lim, upper_lim);
  }

  return out;
}


double my_pBCT_scalar(const double& q,
                      const double& mu,
                      const double& sigma,
                      const double& nu,
                      const double& tau)
  {
// TODO test input within accepted range
   double z;
   if (nu != 0) z = ((pow(q/mu, nu) - 1)/(nu * sigma));
   else z = log(q/mu)/sigma;
   double FYy1 = R::pt(z, tau, true, false);
   double FYy2 = 0;
   if (nu > 0) FYy2 = R::pt(-1/(sigma * abs(nu)), tau, true, false);
   double FYy3 = R::pt(1/(sigma * abs(nu)), tau, true, false);
   return (FYy1 - FYy2)/FYy3;
  }


double my_qBCT_scalar(const double& p,
                      const double& mu,
                      const double& sigma,
                      const double& nu,
                      const double& tau)
{
    // TODO test input within accepted range
  double z;
  if (nu <= 0) z = R::qt(p * R::pt(1/(sigma * abs(nu)), tau, true, false), tau, true, false);
  else z = R::qt(1 - (1 - p) * R::pt(1/(sigma * abs(nu)), tau, true, false), tau, true, false);
  double ya;
  if (nu != 0) ya = mu * pow((nu * sigma * z + 1), (1/nu));
  else ya = mu * exp(sigma * z);
  return ya;
}

double my_qBCT_trunc_scalar(const double& p,
                              const double& mu,
                              const double& sigma,
                              const double& nu,
                              const double& tau,
                              const double& lower_lim,
                              const double& upper_lim) {
  double pp1 = my_pBCT_scalar(lower_lim, mu, sigma, nu, tau);
  double pp2 = my_pBCT_scalar(upper_lim, mu, sigma, nu, tau);
  return my_qBCT_scalar((p * (pp2 - pp1) + pp1),  mu, sigma, nu, tau);
}

//' @export
// [[Rcpp::export]]
NumericVector my_pBCT(const NumericVector& q,
                              const NumericVector& mu,
                              const NumericVector& sigma,
                              const NumericVector& nu,
                              const NumericVector& tau,
                              const int& n_cpu) {
  // TODO check all have same length and are > 0 etc.
  const int n = q.length();
  NumericVector out(n);

  omp_set_num_threads(n_cpu); // Use n_cpu threads for all

  // consecutive parallel regions
#pragma omp parallel for default(shared)
  for (int i = 0; i < n; i++)
  {
    out[i] = my_pBCT_scalar(q[i], mu[i], sigma[i], nu[i], tau[i]);
  }

  return out;
}

//' @export
// [[Rcpp::export]]
NumericVector my_qBCT(const NumericVector& p,
                              const NumericVector& mu,
                              const NumericVector& sigma,
                              const NumericVector& nu,
                              const NumericVector& tau,
                              const int& n_cpu) {
  // TODO check all have same length and are > 0 etc.
  const int n = p.length();
  NumericVector out(n);

  omp_set_num_threads(n_cpu); // Use n_cpu threads for all

  // consecutive parallel regions
#pragma omp parallel for default(shared)
  for (int i = 0; i < n; i++)
  {
    out[i] = my_qBCT_scalar(p[i], mu[i], sigma[i], nu[i], tau[i]);
  }

  return out;
}

//' @export
// [[Rcpp::export]]
NumericVector my_qBCT_trunc(const NumericVector& p,
                                    const NumericVector& mu,
                                    const NumericVector& sigma,
                                    const NumericVector& nu,
                                    const NumericVector& tau,
                                    const double& lower_lim,
                                    const double& upper_lim,
                                    const int& n_cpu) {
  // TODO check all have same length and are > 0 etc.
  const int n = p.length();
  NumericVector out(n);

  omp_set_num_threads(n_cpu); // Use n_cpu threads for all

  // consecutive parallel regions
#pragma omp parallel for default(shared)
  for (int i = 0; i < n; i++)
  {
    out[i] = my_qBCT_trunc_scalar(p[i], mu[i], sigma[i], nu[i], tau[i], lower_lim, upper_lim);
  }

  return out;
}


//' @export
// [[Rcpp::export]]
NumericVector my_runif_seeded(const int& n,
                       const double& lower_bound,
                       const double& upper_bound,
                       const int& seed)
{
  // Make a random number engine
  pcg32 rng(seed);

  // Choose a random number between lower_bound and upper_bound
  uniform_real_distribution<double> uniform_dist(lower_bound, upper_bound);
  NumericVector out(n);


  for (int i = 0; i < n; i++)
    out[i] = uniform_dist(rng);

  return out;
}

//' @export
// [[Rcpp::export]]
NumericVector my_runif(const int& n,
                       const double& lower_bound,
                       const double& upper_bound)
{
  // Make a random number engine
    pcg32 rng(randutils::auto_seed_128{}.base());

  // Choose a random number between lower_bound and upper_bound
  uniform_real_distribution<double> uniform_dist(lower_bound, upper_bound);
  NumericVector out(n);


  for (int i = 0; i < n; i++)
    out[i] = uniform_dist(rng);

  return out;
}


//' @export
// [[Rcpp::export]]
LogicalVector mk_new_simulant_markers(const IntegerVector &pid)
{
  const int n = pid.size();
  LogicalVector new_simulant_markers(n);
  new_simulant_markers[0] = true;
  int previous_pid = pid[0];
  // Loop with no conditional branches in the body (therefore branch predictor should get it right almost every time) and minimal memory access by retaining previous_pid.
  for (int i = 1; i < n; i++)
  {
    new_simulant_markers[i] = pid[i] != previous_pid;
    previous_pid = pid[i];
  }
  return new_simulant_markers;
}

//' @export
// [[Rcpp::export]]
IntegerVector incidence_type_1(
    const IntegerVector &year,
    const LogicalVector &is_new_simulant,
    const NumericVector &prob,
    const IntegerVector &prevalence,
    const int init_year)
{
  // pid will be sorted by year. TODO force and check in R side
  const int n = is_new_simulant.size();
  IntegerVector out = clone(prevalence);

  for (int i = 0; i < n; i++)
    if (year[i] >= init_year)  // start logic after init year
    {
      if (is_new_simulant[i])
      {
        if (prob[i] == 1)
          out[i]++;
      }
      else
      {
        out[i] = out[i - 1];
        if (prob[i] == 1)
          out[i]++;
      }
    }
    return out;
}

//' @export
// [[Rcpp::export]]
IntegerVector incidence_type_2(
    const IntegerVector &year,
    const LogicalVector &is_new_simulant,
    const NumericVector &prob,
    const NumericVector &rn,
    const IntegerVector &prevalence,
    const int init_year)
{
  // pid will be sorted by year. TODO force and check in R side
  const int n = is_new_simulant.size();
  IntegerVector out = clone(prevalence);

  for (int i=0; i<n; i++)
    if (year[i] >= init_year)  // start logic after init year
    {
      if (is_new_simulant[i])
      {
        if (out[i] == 0 && prob[i] > rn[i])
          out[i]++;
      }
      else
      {
        if (out[i-1] > 0)
          out[i] = out[i-1] + 1;
        else if (out[i-1] == 0 && prob[i] > rn[i])
          out[i]++;
      }
    }
    return(out);
}

//' @export
// [[Rcpp::export]]
DataFrame incidence_type_3(
    const IntegerVector &year,
    const LogicalVector &is_new_simulant,
    const NumericVector &t2dm_prob,
    const NumericVector &t2dm_rn,
    const IntegerVector &t2dm_prvl,
    const NumericVector &t2dm_cvd_multiplier,
    const NumericVector &chd_prob,
    const NumericVector &chd_rn,
    const IntegerVector &chd_prvl,
    const NumericVector &stroke_prob,
    const NumericVector &stroke_rn,
    const IntegerVector &stroke_prvl,
    const NumericVector &cvd_t2dm_multiplier,
    const int t2dm_lag,
    const int cvd_lag,
    const int init_year)
{
  // pid will be sorted by year. TODO force and check in R side

  const int n = is_new_simulant.size();
  IntegerVector out_t2dm   = clone(t2dm_prvl);
  IntegerVector out_chd    = clone(chd_prvl);
  IntegerVector out_stroke = clone(stroke_prvl);
  IntegerVector out_cvd(n);

  for (int i=0; i<n; i++)
    if (year[i] >= init_year)  // start logic after init year
    {
      if (is_new_simulant[i])
      {
        out_cvd[i] = max(chd_prvl[i], stroke_prvl[i]);
      }
      else
      {
        if (out_cvd[i-1]    > 0) out_cvd[i]    = out_cvd[i-1] + 1;
        if (out_chd[i-1]    > 0) out_chd[i]    = out_chd[i-1] + 1;
        if (out_stroke[i-1] > 0) out_stroke[i] = out_stroke[i-1] + 1;
        if (out_t2dm[i-1]   > 0) out_t2dm[i]   = out_t2dm[i-1] + 1;
      }

      double used_cvd_t2dm_multiplier = out_t2dm[i] < t2dm_lag ? 1.0 : cvd_t2dm_multiplier[i];
      double used_t2dm_cvd_multiplier = out_cvd[i] < cvd_lag ? 1.0 : t2dm_cvd_multiplier[i];
      if (out_t2dm[i] == 0 && (t2dm_prob[i] * used_t2dm_cvd_multiplier) > t2dm_rn[i])
        out_t2dm[i]++;
      if (out_chd[i] == 0 && (chd_prob[i] * used_cvd_t2dm_multiplier) > chd_rn[i])
        out_chd[i]++;
      if (out_stroke[i] == 0 && (stroke_prob[i] * used_cvd_t2dm_multiplier) > stroke_rn[i])
        out_stroke[i]++;
      out_cvd[i] = max(out_chd[i], out_stroke[i]); // update cvd after all the calculations
    }

    return DataFrame::create(
      _["t2dm_prvl"]  = out_t2dm,
      _["chd_prvl"]   = out_chd,
      _["stroke_prvl"]= out_stroke,
      _["cvd_prvl"]   = out_cvd);
}

//' @export
// [[Rcpp::export]]
IntegerVector mortality_type_1(
    const IntegerVector &year,
    const LogicalVector &is_new_simulant,
    const NumericVector &prob,
    const NumericVector &rn,
    const IntegerVector &prevalence,
    const int init_year,
    const int cod)
{
  // pid will be sorted by year. TODO force and check in R side
  const int n = is_new_simulant.size();
  IntegerVector out(n);

  for (int i=0; i<n; i++)
    if (year[i] >= init_year)  // start logic after init year
    {
      if (is_new_simulant[i])
      {
        if (prevalence[i] > 0 && prob[i] > rn[i])
          out[i] = cod;
      }
      else
      {
        if (out[i - 1] == cod)
          out[i] = cod;
        else if (prevalence[i] > 0 && prob[i] > rn[i])
          out[i] = cod;
      }
    }
    return(out);
}

//' @export
// [[Rcpp::export]]
DataFrame mortality_type_2(
    const IntegerVector &year,
    const LogicalVector &is_new_simulant,
    const NumericVector &prob,
    const NumericVector &rn,
    const IntegerVector &diagnosed,
    const IntegerVector &prevalence,
    const int init_year,
    const int cod,
    const int cure
)
{
  // pid will be sorted by year. TODO force and check in R side
  const int n = is_new_simulant.size();
  IntegerVector out_prvl = clone(prevalence);
  IntegerVector out_dgn = clone(diagnosed);
  IntegerVector out_mrtl(n);

  for (int i = 0; i < n; i++)
    if (year[i] >= init_year)  // start logic after init year
    {
      if (is_new_simulant[i])
      {
        if (prevalence[i] > 0 && prevalence[i] <= cure && prob[i] > rn[i])
          out_mrtl[i] = cod;
        else
        { // Known: out_mrtl[i] != cod
          if (prevalence[i] > cure)
          {
            out_prvl[i] = 0;
            out_dgn[i] = 0;
          }
        }
      }
      else
      {
        if (out_mrtl[i - 1] == cod)
          out_mrtl[i] = cod;
        else
        { // Known: out_mrtl[i - 1] != cod
          if (prevalence[i] > 0 && prevalence[i] <= cure && prob[i] > rn[i])
            out_mrtl[i] = cod;
          else
          { // Known: out_mrtl[i] != cod
            if (prevalence[i] > cure)
            {
              out_prvl[i] = 0;
              out_dgn[i] = 0;
            }
          }
        }
      }
    }
    return DataFrame::create(
      _["disease_prvl"]  = out_prvl,
      _["disease_dgn"]   = out_dgn,
      _["disease_mrtl"]  = out_mrtl);
}

inline int calculate_cod_for_row(
    const vector<IntegerVector> &cods,
    const NumericVector &rn,
    const int row)
{
  const int conditionCount = cods.size();
  // Aggregate the possible causes.
  int possible_causes[conditionCount];
  int next_possible_cause = 0;
  for (int cod_index = 0; cod_index < conditionCount; cod_index++)
  {
    if (cods[cod_index][row] != 0)
      possible_causes[next_possible_cause++] = cods[cod_index][row];
  }
  if (next_possible_cause == 0)
  {
    // No death - the very common case!
    return 0;
  }
  else if (next_possible_cause == 1)
  {
    // Unique cause of death (the commonest case when dying by far) - just put the sole cause into the output and we're done.
    return possible_causes[0];
  }
  else
  {
    // More than one cause of death; select randomly (and equally) between them.
    // TODO: This assumes rn is in the range [0,1) i.e. does not include 1.0!
    return possible_causes[(int)(rn[row] * next_possible_cause)];
  }
}

//' @export
// [[Rcpp::export]]
IntegerVector mortality_resolve(
    const IntegerVector &year,
    const LogicalVector &is_new_simulant,
    const IntegerVector &cod1,
    const IntegerVector &cod2,
    const IntegerVector &cod3,
    const NumericVector &rn,
    int init_year)
{
  // pid will be sorted by year. TODO force and check in R side
  const int n = is_new_simulant.size();
  IntegerVector out(n);
  const vector<IntegerVector> cods({cod1, cod2, cod3});

  for (int i = 0; i < n; i++)
    if (year[i] >= init_year)  // start logic after init year
    {
      if (is_new_simulant[i]) // if the same simulant
        out[i] = calculate_cod_for_row(cods, rn, i);
      else // same simulant - which guarantees i > 0 and hence [i - 1] is legal
      {
        if (out[i - 1] != 0)
          out[i] = out[i - 1];
        else
          out[i] = calculate_cod_for_row(cods, rn, i);
      }
    }
    return out;
}

inline string calculate_morb_for_row(
    const vector<IntegerVector> &cods,
    const CharacterVector disease_names,
    const CharacterVector cod_names,
    const int row)
{
  const int conditionCount = cods.size();
  // TODO check cods.size() - 1 == disease_names.size()
  // Concatenate diseases.
  string possible_causes = "";
  for (int cod_index = 0; cod_index < conditionCount; cod_index++)
  {
    if (cods[cod_index][row] != 0)
    {
      if (!possible_causes.empty()) possible_causes += ", ";
      // handle diseases
      if (cod_index < (conditionCount - 1))
      {
        possible_causes += disease_names[cod_index];
        possible_causes += to_string(cods[cod_index][row]);
      }
      // handle deaths
      else
      {
        possible_causes += ("Died of " + cod_names[cods[cod_index][row] - 1]);
      }
    }
  }
  if (possible_causes.empty()) possible_causes = "Healthy";

  return possible_causes;
}

//' @export
// [[Rcpp::export]]
CharacterVector morbidity_resolve(
    const IntegerVector   &year,
    const IntegerVector   &disease1,
    const IntegerVector   &disease2,
    const IntegerVector   &death,
    const CharacterVector &disease_names,
    const CharacterVector &cod_names,
    int init_year)
{
  const int n = year.size();
  CharacterVector out(n);
  const vector<IntegerVector> cods({disease1, disease2, death});

  for (int i = 0; i < n; i++)
    if (year[i] >= init_year)  // start logic after init year
    {
      out[i] = calculate_morb_for_row(cods, disease_names, cod_names, i);
    }
    return out;
}

//' @export
// [[Rcpp::export]]
int count_if(LogicalVector x, bool na_rm = false) {
  if (na_rm) x = na_omit(x); // remove NA from denominator
  const int n = x.size();
  int counter = 0;
  for(int i = 0; i < n; i++) {
    if(x[i] == TRUE) {
      counter++;
    }
  }
  return counter;
}

//' @export
// [[Rcpp::export]]
double prop_if(LogicalVector x, bool na_rm = false) {
  if (na_rm) x = na_omit(x); // remove NA from denominator
  const int n = x.size();
  int counter = 0;
  for(int i = 0; i < n; i++) {
    if(x[i] == TRUE) {
      counter++;
    }
  }
  return counter/(double)n;
}


double my_runif_scalar(const double& lower_bound,
                              const double& upper_bound)
{
  // Make a random number engine
  pcg32 rng(randutils::auto_seed_128{}.base());

  // Choose a random number between lower_bound and upper_bound
  uniform_real_distribution<double> uniform_dist(lower_bound, upper_bound);

  return uniform_dist(rng);
}

// helper function that makes a double jump within x - jump and x + jump
// special care to avoid values outside (0, 1)
double fscramble_hlp (const double x, const double jump)
{
  double out = my_runif_scalar(x - jump, x + jump);
  if ((out <= 0.0) | (out >= 1.0)) out = x;

  return out;
}

//' @export
// [[Rcpp::export]]
NumericVector fscramble_trajectories(const NumericVector &x, const IntegerVector &pid, const double &jump) {
  // pid should be sorted and same length as x
  int n = x.size();
  NumericVector out(n);
  out[0] = x[0];
  for (int i = 1; i < n; ++i)
  {
    if (pid[i] == pid[i-1])
    {
      out[i] = fscramble_hlp(out[i-1], jump);
    }
    else
    {
      out[i] = x[i];
    }
  }
  return out;
}

//' @export
// [[Rcpp::export]]
LogicalVector identify_longdeads(const IntegerVector& x, const LogicalVector& pid) {
  const int n = x.size();
  LogicalVector out(n);
  for (int i = 0; i < n; i++)
  {
    if (!pid[i] && x[i-1] != 0) out[i] = true; else out[i] = false;
  }
  return out;
}

