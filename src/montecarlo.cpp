#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <string>
#include <optional>
#include <vector>
#include <random>
#include "tinyexpr.h"

using namespace Rcpp;

// =============================================================
// UNIVERSAL FUNCTION WRAPPER (R function / tinyexpr)
// =============================================================
struct FunctionWrapper {
  bool use_tinyexpr;
  std::string expr;
  te_expr *compiled_expr = nullptr;
  std::vector<std::string> var_names;
  std::vector<double> var_values;
  std::optional<Rcpp::Function> f;

  FunctionWrapper(Rcpp::Function f_) : use_tinyexpr(false), f(f_) {}
  FunctionWrapper(const std::string& expr_, const std::vector<std::string>& vars)
    : use_tinyexpr(true), expr(expr_), var_names(vars), var_values(vars.size()) {
    if (expr.empty()) Rcpp::stop("Empty expression passed to FunctionWrapper.");
    std::vector<te_variable> te_vars;
    for (size_t i = 0; i < vars.size(); i++)
      te_vars.push_back({vars[i].c_str(), &var_values[i]});
    compiled_expr = te_compile(expr.c_str(), te_vars.data(), te_vars.size(), nullptr);
    if (!compiled_expr) Rcpp::stop("Failed to compile expression with tinyexpr.");
  }

  ~FunctionWrapper() { if (compiled_expr) te_free(compiled_expr); }

  double eval(const NumericVector& x) {
    if (!use_tinyexpr) {
      if (!f.has_value()) Rcpp::stop("Function not initialized.");
      return Rcpp::as<double>((*f)(x));
    } else {
      if (x.size() != var_values.size()) Rcpp::stop("Input vector length mismatch.");
      for (size_t i = 0; i < var_values.size(); i++) var_values[i] = x[i];
      return te_eval(compiled_expr);
    }
  }
};

// =============================================================
// SUPPORT FUNCTIONS
// =============================================================
double sample_truncated_normal(double mean, double sd, double lower, double upper) {
  double alpha = R::pnorm(lower, mean, sd, 1, 0);
  double beta  = R::pnorm(upper, mean, sd, 1, 0);
  double u = R::runif(alpha, beta);
  return R::qnorm(u, mean, sd, 1, 0);
}

double truncated_normal_density(double x, double mean, double sd, double lower, double upper) {
  double Z = R::pnorm(upper, mean, sd, 1, 0) - R::pnorm(lower, mean, sd, 1, 0);
  if (Z < 1e-12) Z = 1e-12;
  return R::dnorm(x, mean, sd, 0) / Z;
}

double sample_truncated_beta(double alpha, double beta, double lower, double upper) {
  double u = R::runif(0,1);
  double q = R::qbeta(u, alpha, beta, 1, 0);
  return lower + q * (upper - lower);
}

double truncated_beta_density(double x, double alpha, double beta, double lower, double upper) {
  double z = (x - lower) / (upper - lower);
  double dens = R::dbeta(z, alpha, beta, 0) / (upper - lower);
  return dens;
}

double sample_truncated_exponential(double rate, double lower, double upper) {
  double u = R::runif(0,1);
  double cdf_lower = 1 - exp(-rate * lower);
  double cdf_upper = 1 - exp(-rate * upper);
  double z = cdf_lower + u * (cdf_upper - cdf_lower);
  return -log(1 - z)/rate;
}

double truncated_exponential_density(double x, double rate, double lower, double upper) {
  double Z = exp(-rate * lower) - exp(-rate * upper);
  if (Z < 1e-12) Z = 1e-12;
  return rate * exp(-rate * x) / Z;
}

// =============================================================
// ADAPTIVE PARTITIONING
// =============================================================
List auto_partition_adaptive(Function f, double lower, double upper, int n_grid = 1000, int pilot = 10) {
  NumericVector x(n_grid), y(n_grid);
  double dx = (upper - lower) / (n_grid - 1);
  for (int i = 0; i < n_grid; i++) {
    x[i] = lower + i * dx;
    y[i] = as<double>(f(x[i]));
  }

  NumericVector y_smooth = clone(y);
  for (int i = 1; i < n_grid - 1; i++) y_smooth[i] = (y[i-1]+y[i]+y[i+1])/3.0;

  std::vector<double> troughs = {lower};
  for (int i = 1; i < n_grid - 1; i++)
    if (y_smooth[i] < y_smooth[i-1] && y_smooth[i] < y_smooth[i+1])
      troughs.push_back(x[i]);
    troughs.push_back(upper);

    List partitions;
    for (size_t i = 0; i < troughs.size() - 1; i++) {
      double p_lower = troughs[i], p_upper = troughs[i+1];
      double max_val = R_NegInf, max_x = p_lower;
      for (int j = 0; j < pilot; j++) {
        double xi = p_lower + (p_upper - p_lower)*R::runif(0,1);
        double yi = as<double>(f(xi));
        if (yi > max_val) { max_val = yi; max_x = xi; }
      }
      partitions.push_back(List::create(
          _["lower"] = NumericVector::create(p_lower),
          _["upper"] = NumericVector::create(p_upper),
          _["center"] = max_x,
          _["max_value"] = max_val
      ));
    }
    return partitions;
}

List auto_partition_nd_adaptive(Function f, NumericVector lower, NumericVector upper, int n_per_dim = 4, int pilot = 10) {
  int dim = lower.size();
  if (dim == 1) return auto_partition_adaptive(f, lower[0], upper[0], n_per_dim*100, pilot);

  List partitions;
  std::vector<double> dx(dim);
  for (int d = 0; d < dim; d++) dx[d] = (upper[d]-lower[d])/n_per_dim;

  IntegerVector idx(dim, 0);
  NumericVector part_lower(dim), part_upper(dim);
  bool done = false;

  while(!done) {
    for (int d = 0; d < dim; d++) {
      part_lower[d] = lower[d]+idx[d]*dx[d];
      part_upper[d] = lower[d]+(idx[d]+1)*dx[d];
    }
    double max_val = R_NegInf;
    NumericVector center(dim);
    for (int i=0;i<pilot;i++){
      NumericVector x(dim);
      for (int d=0; d<dim; d++) x[d]=part_lower[d]+(part_upper[d]-part_lower[d])*R::runif(0,1);
      double val = as<double>(f(x));
      if (val>max_val) { max_val=val; center = clone(x); }
    }
    partitions.push_back(List::create(
        _["lower"] = clone(part_lower),
        _["upper"] = clone(part_upper),
        _["center"] = clone(center),
        _["max_value"] = max_val
    ));

    for (int d=dim-1; d>=0; d--){
      idx[d]++;
      if (idx[d]<n_per_dim) break;
      else if(d==0) {done=true; break;} else idx[d]=0;
    }
  }
  return partitions;
}

// =============================================================
// IMPORTANCE PARAMETER OPTIMIZATION
// =============================================================
List optimize_importance_params(FunctionWrapper &fw, const NumericVector &lower, const NumericVector &upper,
                                int dim, const std::string &distribution, int n_pilot=500) {
  NumericMatrix samples(n_pilot, dim);
  NumericVector values(n_pilot);
  for (int i=0;i<n_pilot;i++){
    NumericVector x(dim);
    for(int d=0;d<dim;d++){x[d]=lower[d]+(upper[d]-lower[d])*R::runif(0,1); samples(i,d)=x[d];}
    values[i] = std::abs(fw.eval(x));
  }
  double sum_values = std::max(double(sum(values)),1e-12);
  NumericVector weights = values/sum_values;
  List params;
  if(distribution=="normal"){
    NumericVector means(dim), sds(dim);
    for(int d=0;d<dim;d++){
      double mean_est=0,var_est=0;
      for(int i=0;i<n_pilot;i++) mean_est+=weights[i]*samples(i,d);
      for(int i=0;i<n_pilot;i++) var_est+=weights[i]*pow(samples(i,d)-mean_est,2);
      means[d]=mean_est;
      sds[d]=std::max(std::sqrt(var_est),0.01*(upper[d]-lower[d]));
    }
    params["mean"]=means; params["sd"]=sds;
  } else if(distribution=="beta"){
    NumericVector alpha(dim,2.0), beta(dim,2.0);
    params["alpha"]=alpha; params["beta"]=beta;
  } else if(distribution=="exponential"){
    NumericVector rate(dim,1.0);
    params["rate"]=rate;
  } else if(distribution=="mixture_normal"){
    NumericVector means(dim), sds(dim);
    for(int d=0;d<dim;d++){means[d]=(upper[d]+lower[d])/2; sds[d]=0.3*(upper[d]-lower[d]);}
    params["mean"]=means; params["sd"]=sds;
  }
  return params;
}

// =============================================================
// IMPORTANCE SAMPLER
// =============================================================
List sample_importance(const NumericVector &lower, const NumericVector &upper, int dim,
                       const std::string &distribution, const List &params){
  NumericVector x(dim); double dens=1.0;
  if(distribution=="normal"){
    NumericVector mean=params["mean"], sd=params["sd"];
    for(int d=0;d<dim;d++){
      x[d]=sample_truncated_normal(mean[d],sd[d],lower[d],upper[d]);
      dens*=truncated_normal_density(x[d],mean[d],sd[d],lower[d],upper[d]);
    }
  } else if(distribution=="beta"){
    NumericVector alpha=params["alpha"], beta=params["beta"];
    for(int d=0;d<dim;d++){
      x[d]=sample_truncated_beta(alpha[d],beta[d],lower[d],upper[d]);
      dens*=truncated_beta_density(x[d],alpha[d],beta[d],lower[d],upper[d]);
    }
  } else if(distribution=="exponential"){
    NumericVector rate=params["rate"];
    for(int d=0;d<dim;d++){
      x[d]=sample_truncated_exponential(rate[d],lower[d],upper[d]);
      dens*=truncated_exponential_density(x[d],rate[d],lower[d],upper[d]);
    }
  } else if(distribution=="mixture_normal"){
    NumericVector mean=params["mean"], sd=params["sd"];
    for(int d=0;d<dim;d++){
      if(R::runif(0,1)<0.5) x[d]=sample_truncated_normal(mean[d]-sd[d],sd[d],lower[d],upper[d]);
      else x[d]=sample_truncated_normal(mean[d]+sd[d],sd[d],lower[d],upper[d]);
      dens*=0.5*truncated_normal_density(x[d],mean[d]-sd[d],sd[d],lower[d],upper[d])
        +0.5*truncated_normal_density(x[d],mean[d]+sd[d],sd[d],lower[d],upper[d]);
    }
  }
  return List::create(_["x"]=x,_["density"]=dens);
}

// =============================================================
// MONTE CARLO INTEGRATOR
// [[Rcpp::export]]
List montecarlo_integrate_cpp(SEXP f_input, NumericVector lower, NumericVector upper,
                              int n_samples, bool partition=false, int dim=1,
                              std::string expr="", Nullable<CharacterVector> vars=R_NilValue,
                              bool importance_sampling=false, std::string is_distribution="normal",
                              Nullable<List> is_params=R_NilValue, int step=100){
  std::vector<std::string> varnames;
  if(!Rf_isNull(vars)){
    CharacterVector v(vars);
    for(int i=0;i<v.size();i++) varnames.push_back(as<std::string>(v[i]));
  } else for(int d=0;d<dim;d++) varnames.push_back("x"+std::to_string(d+1));

  FunctionWrapper fw = Rf_isFunction(f_input)? FunctionWrapper(as<Function>(f_input)) : FunctionWrapper(expr,varnames);

  List params;
  if(importance_sampling){
    if(is_params.isNull()) params = optimize_importance_params(fw, lower, upper, dim, is_distribution);
    else params = as<List>(is_params);
  }

  List partitions;
  if(partition){
    if(dim==1 && Rf_isFunction(f_input)) partitions = auto_partition_adaptive(as<Function>(f_input), lower[0], upper[0]);
    else if(dim>1 && Rf_isFunction(f_input)) partitions = auto_partition_nd_adaptive(as<Function>(f_input), lower, upper);
    else partitions = List::create(List::create(_["lower"]=lower,_["upper"]=upper));
  } else partitions = List::create(List::create(_["lower"]=lower,_["upper"]=upper));

  int n_partitions = partitions.size();
  int samples_per_partition = std::max(1,n_samples/n_partitions);

  NumericVector estimates; double sum=0.0; int total_samples=0;
  for(int p=0;p<n_partitions;p++){
    List part=partitions[p];
    NumericVector part_lower=part["lower"], part_upper=part["upper"];
    for(int i=0;i<samples_per_partition;i++){
      NumericVector x(dim);
      if(!importance_sampling){
        for(int d=0;d<dim;d++) x[d]=part_lower[d]+(part_upper[d]-part_lower[d])*R::runif(0,1);
        sum+=fw.eval(x);
      } else {
        List s=sample_importance(part_lower, part_upper, dim, is_distribution, params);
        NumericVector xs = s["x"];
        double qx = s["density"];
        double val = fw.eval(xs);
        sum+=val/qx;
      }
      total_samples++;
      if(total_samples%step==0){
        double volume=1.0; for(int d=0;d<dim;d++) volume*=(upper[d]-lower[d]);
        double current_estimate = importance_sampling ? sum/total_samples : volume*sum/total_samples;
        estimates.push_back(current_estimate);
      }
    }
  }

  double volume=1.0; for(int d=0;d<dim;d++) volume*=(upper[d]-lower[d]);
  double final_estimate = importance_sampling ? sum/total_samples : volume*sum/total_samples;

  return List::create(
    _["estimate"]=final_estimate,
    _["estimates"]=estimates,
    _["importance_sampling"]=importance_sampling,
    _["is_distribution"]=is_distribution,
    _["params_used"]=importance_sampling?Rcpp::wrap(params):R_NilValue,
    _["expr"]=fw.use_tinyexpr?Rcpp::wrap(fw.expr):R_NilValue,
    _["n_partitions"]=n_partitions
  );
}
