#include <Rcpp.h>
using namespace Rcpp;

// -------------------- 1D Adaptive Signed Partition --------------------
// [[Rcpp::export]]
List adaptive_partition_signed_1d(Function f, double lower, double upper, int n_grid = 1000) {
  NumericVector x(n_grid);
  NumericVector y(n_grid);
  double dx = (upper - lower)/(n_grid-1);

  for(int i=0;i<n_grid;i++){
    x[i] = lower + i*dx;
    y[i] = as<double>(f(x[i]));
  }

  std::vector<double> breaks;
  breaks.push_back(lower);

  for(int i=1;i<n_grid-1;i++){
    // split at zero crossings or local extrema
    if((y[i] < y[i-1] && y[i] < y[i+1]) || (y[i] > y[i-1] && y[i] > y[i+1]) || (y[i]*y[i-1]<0)){
      breaks.push_back(x[i]);
    }
  }
  breaks.push_back(upper);

  List partitions;
  for(size_t i=0;i<breaks.size()-1;i++){
    partitions.push_back(
      List::create(
        _["lower"] = NumericVector::create(breaks[i]),
        _["upper"] = NumericVector::create(breaks[i+1])
      )
    );
  }
  return partitions;
}

// -------------------- ND Grid Generator --------------------
void generate_grid_nd(const NumericVector& lower, const NumericVector& upper,
                      const IntegerVector& n_grid, int dim,
                      NumericVector& current, List& all_points) {
  if(dim == lower.size()){
    all_points.push_back(clone(current));
    return;
  }

  double step = (upper[dim] - lower[dim]) / (n_grid[dim]-1);
  for(int i=0;i<n_grid[dim];i++){
    current[dim] = lower[dim] + i*step;
    generate_grid_nd(lower, upper, n_grid, dim+1, current, all_points);
  }
}

// -------------------- ND Signed Partition --------------------
List adaptive_partition_signed_nd(Function f, NumericVector lower, NumericVector upper, int n_grid_per_dim = 5){
  int dim = lower.size();
  IntegerVector n_grid(dim, n_grid_per_dim);

  List all_points;
  NumericVector current(dim);
  generate_grid_nd(lower, upper, n_grid, 0, current, all_points);

  List partitions;
  for(int i=0;i<all_points.size();i++){
    NumericVector pt = all_points[i];
    NumericVector delta = (upper - lower) / (2.0*n_grid_per_dim);

    // define lower and upper for this partition
    NumericVector p_lower = pt - delta;
    NumericVector p_upper = pt + delta;

    // make sure bounds don't exceed original bounds
    for(int d=0;d<dim;d++){
      if(p_lower[d] < lower[d]) p_lower[d] = lower[d];
      if(p_upper[d] > upper[d]) p_upper[d] = upper[d];
    }

    partitions.push_back(
      List::create(
        _["lower"] = p_lower,
        _["upper"] = p_upper
      )
    );
  }
  return partitions;
}

// -------------------- Function Evaluator --------------------
double eval_f(Function f, const NumericVector& x, int dim){
  if(dim==1) return as<double>(f(x[0]));
  else if(dim==2) return as<double>(f(x[0], x[1]));
  else if(dim==3) return as<double>(f(x[0], x[1], x[2]));
  else {
    List args(dim);
    for(int d=0;d<dim;d++) args[d] = x[d];
    return as<double>(Rcpp::Language(f, args).eval());
  }
}

// -------------------- Sample Allocation --------------------
IntegerVector allocate_samples(const List& partitions, int n_samples){
  int n = partitions.size();
  IntegerVector alloc(n);
  int base = n_samples / n;
  int rem = n_samples % n;
  for(int i=0;i<n;i++){
    alloc[i] = base + (i<rem?1:0);
  }
  return alloc;
}

// -------------------- Integrate Partition --------------------
double integrate_partition(Function f, const List& partition, int samples, int dim){
  NumericVector lower = partition["lower"];
  NumericVector upper = partition["upper"];
  NumericVector x(dim);

  double volume = 1.0;
  for(int d=0;d<dim;d++) volume *= (upper[d]-lower[d]);

  double sum = 0.0;
  for(int i=0;i<samples;i++){
    for(int d=0;d<dim;d++){
      x[d] = lower[d] + (upper[d]-lower[d])*R::runif(0,1);
    }
    sum += eval_f(f, x, dim);
  }

  return volume * sum / samples;
}

// -------------------- Main Monte Carlo Function --------------------
// [[Rcpp::export(name="montecarlo_integrate_cpp")]]
double montecarlo_integrate(Function f, NumericVector lower, NumericVector upper,
                            int n_samples, bool partition=false, int dim=1,
                            int n_grid=1000, int n_grid_per_dim=5){

  if(!partition){
    NumericVector x(dim);
    double sum=0.0, volume=1.0;
    for(int d=0;d<dim;d++) volume *= (upper[d]-lower[d]);
    for(int i=0;i<n_samples;i++){
      for(int d=0;d<dim;d++) x[d] = lower[d] + (upper[d]-lower[d])*R::runif(0,1);
      sum += eval_f(f, x, dim);
    }
    return volume * sum / n_samples;
  } else {
    List partitions;
    if(dim==1){
      partitions = adaptive_partition_signed_1d(f, lower[0], upper[0], n_grid);
    } else {
      partitions = adaptive_partition_signed_nd(f, lower, upper, n_grid_per_dim);
    }

    if(partitions.size()==0) return montecarlo_integrate(f, lower, upper, n_samples,false,dim);

    IntegerVector alloc = allocate_samples(partitions, n_samples);
    double total=0.0;

    for(int i=0;i<partitions.size();i++){
      if(alloc[i]>0) total += integrate_partition(f, partitions[i], alloc[i], dim);
    }

    return total;
  }
}
