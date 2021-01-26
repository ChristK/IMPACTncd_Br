#include <Rcpp.h>
using namespace Rcpp;

// This is borrowed from Kmisc package for convenience, since this is not in CRAN anymore
// https://github.com/kevinushey/Kmisc

inline bool IsNA( const Rcpp::internal::const_string_proxy<STRSXP>& x ) {
  return x == NA_STRING;
}

inline bool IsNA(int x) {
  return x == NA_INTEGER;
}

inline bool IsNA(SEXP x) {
  return x == NA_STRING;
}

inline bool IsNA(double x) {
  return Rcpp::internal::Rcpp_IsNA(x);
}

inline bool IsNaN(double x) {
  return Rcpp::internal::Rcpp_IsNaN(x);
}

// borrowed from data.table package;
// see https://github.com/arunsrinivasan/datatable/blob/master/pkg/src/countingcharacter.c
inline int StrCmp(SEXP x, SEXP y)
{
  if (x == NA_STRING) return (y == NA_STRING ? 0 : 1);
  else if (y == NA_STRING) return -1;
  else if (x == y) return 0;  // same string in cache
  else return strcmp(char_nocheck(x), char_nocheck(y));
}

template <typename T>
struct NACompare;

template <>
struct NACompare<int> {
  inline bool operator()(int left, int right) const {
    if (left == NA_INTEGER) return false;
    if (right == NA_INTEGER) return true;
    return left < right;
  }
};

template <>
struct NACompare<double> {
  inline bool operator()(double left, double right) const {

    bool leftNaN = (left != left);
    bool rightNaN = (right != right);

    // this branch inspired by data.table: see
    // https://github.com/arunsrinivasan/datatable/commit/1a3e476d3f746e18261662f484d2afa84ac7a146#commitcomment-4885242
    if (IsNaN(right) and IsNA(left)) return true;

    if (leftNaN != rightNaN) {
      return leftNaN < rightNaN;
    } else {
      return left < right;
    }

  }

};

template <>
struct NACompare<SEXP> {
  inline bool operator()(SEXP left, SEXP right) const {
    return StrCmp(left, right) < 0;
  }
};

template <typename T, typename U>
IntegerVector do_counts(const T&);

template <>
inline IntegerVector do_counts<IntegerVector, int>(const IntegerVector& x) {
  std::map< int, int, NACompare<int> > counts;
  int n = x.size();
  for (int i=0; i < n; ++i) {
    ++counts[ x[i] ];
  }
  return wrap(counts);
}

template <>
inline IntegerVector do_counts<LogicalVector, int>(const LogicalVector& x) {
  std::map< int, int, NACompare<int> > counts;
  int n = x.size();
  for (int i=0; i < n; ++i) {
    ++counts[ x[i] ];
  }
  IntegerVector output = wrap(counts);

  // yuck
  SEXP namesptr = Rf_getAttrib(output, R_NamesSymbol);
  for (int i=0; i < output.size(); ++i) {
    if (strcmp(CHAR(STRING_ELT(namesptr, i)), "0") == 0) {
      SET_STRING_ELT(namesptr, i, Rf_mkChar("FALSE"));
    }
    if (strcmp(CHAR(STRING_ELT(namesptr, i)), "1") == 0) {
      SET_STRING_ELT(namesptr, i, Rf_mkChar("TRUE"));
    }
  }

  return output;

}

template <>
inline IntegerVector do_counts<CharacterVector, SEXP>(const CharacterVector& x) {
  IntegerVector tmp = table(x);
  std::map< SEXP, int, NACompare<SEXP> > counts;
  for (int i=0; i < tmp.size(); ++i) {
    counts[ STRING_ELT( tmp.attr("names"), i ) ] = tmp[i];
  }
  IntegerVector output = wrap(counts);

  CharacterVector names = output.attr("names");
  CharacterVector::iterator it = std::find( names.begin(), names.end(), "NA" );
  if (it != names.end()) {
    *it = NA_STRING;
  }
  return output;

}

template <>
inline IntegerVector do_counts<NumericVector, double>(const NumericVector& x) {
  std::map< double, int, NACompare<double> > counts;
  int n = x.size();
  for (int i=0; i < n; ++i) {
    ++counts[ x[i] ];
  }
  IntegerVector output = wrap(counts);

  // explicitly use R's double-to-character coercion to get good names
  int m = counts.size();
  NumericVector keys = no_init(m);
  typedef std::map< double, int, NACompare<double> >::iterator MapItr;
  int i = 0;
  for (MapItr it = counts.begin(); it != counts.end(); ++it) {
    keys[i] = it->first;
    ++i;
  }
  CharacterVector names = Rf_coerceVector(keys, STRSXP);
  output.attr("names") = names;
  // fix names
  for (int i=0; i < output.size(); ++i) {
    bool samestr = strcmp(
      CHAR(STRING_ELT(output.attr("names"), i)),
      "-0"
    ) == 0;
    if (samestr) {
      SET_STRING_ELT(output.attr("names"), i, Rf_mkChar("0"));
    }
  }
  return output;
}

//' @export
// [[Rcpp::export]]
IntegerVector tableRcpp(SEXP x) {
  switch (TYPEOF(x)) {
  case INTSXP: return table(as<IntegerVector>(x));
  case REALSXP: return table(as<NumericVector>(x));
  case STRSXP: return table(as<CharacterVector>(x));
  case LGLSXP: return table(as<LogicalVector>(x));
  default: {
    stop("unrecognized SEXP type");
    return R_NilValue;
  }
  }
}

//' @export
// [[Rcpp::export]]
IntegerVector counts(SEXP x) {
  switch (TYPEOF(x)) {
  case REALSXP: return do_counts<NumericVector, double>(x);
  case STRSXP: return do_counts<CharacterVector, SEXP>(x);
  case INTSXP: return do_counts<IntegerVector, int>(x);
  case LGLSXP: return do_counts<LogicalVector, int>(x);
  default: {
    stop("unrecognized SEXP type");
    return R_NilValue;
  }
  }
}

