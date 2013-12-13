# Generate an order-1 linear ODE system
#   x'(t) = A x(t) + b
# and its curve.

generation.lode1 <- function (
  dimension
  , timepoint
  , orthogonal_transformation = list()
  , row_column_permutation = TRUE
  , constant = list(0,0)
)

# INPUT:
# dimension: Dimension of the system.
#   Currently the system can contian complex eigen-values,
#   which appear in pairs,
#   and at most one real eigen-value,
#   since currently we have not found a good way to control
#   the correlation between curves by real eigen-values.
#   Hence if "dimension" is even, all eigen-values are complex,
#   otherwise there exists one real eigen-value in the spectrum.
# timepoint: Time points for observation data.
# orthogonal_transformation: A list which applies transformation
#   to coefficient matrix to adjust the sparsity and structure.
#   Each component of the list is a 2-tuple (M,N)
#   A random orthogonal matrix of dimension (N-M+1) is left multipled
#   to Mth-Nth row of the coefficient matrix,
#   and its transpose right multipled to Mth-Nth column of the
#   coefficient matrix.
# row_column_permutation: Make the sparsity structure less obvious
#   by permuting rows and columns of the coefficient matrix.
# constant: A list of two vectors of length "dimension"
#   indicating lower and upper bounds of constant term.
#   Each component can also be a scalar,
#   which will be automatically expanded to a vector.
#   The default NULL value disables the constant term.
#   (Notice that NULL is not mathematically equivalent to "list(0,0)".)

# OUTPUT:
# time: Input argument "timepoint".
# parameter$initial: Initial condition of the system.
# parameter$linear: Linear term of the system.
#   It will be identical to the first column of "curve".
# parameter$constant: Constant term of the system.
#   Value NULL indicates the system does not contain a constant term.
# curve: Data matrix for value at given time points.
#   Each row is a curve.
#   Each column is for one time point.

{

ret <- list()
class(ret) <- 'lode1'

# Sanity Check#{{{

dimension <- as.integer(dimension)
if ( length(dimension)!=1 || dimension<=0 )
{
  stop('Argument "dimension" must be a positive integer.')
}

if ( length(timepoint)<2 )
{
  stop('Argument "timepoint" must be a vector ' ,
    'longer than 1 with strictly ascending elements.')
}
temp <- as.numeric(unique(sort(timepoint)))
if ( length(timepoint)!=length(temp) || !all(timepoint==temp) )
{
  stop('Argument "timepoint" must be a vector ' ,
    'longer than 1 with strictly ascending elements.')
}

row_column_permutation <- as.logical(row_column_permutation)
if ( length(row_column_permutation)!=1 )
{
  stop('Argument "row_column_permutation" must be a logical scalar.')
}

constant_exist <- !is.null(constant)
if ( constant_exist )
{
  constant <- as.list(constant)
  if ( length(constant)!=2 )
  {
    stop('Argument "constant" must be '
      ,'a list of two vectors or scalars.')
  }
  constant[[1]] <- as.numeric(constant[[1]])
  if ( length(constant[[1]])==1 )
  {
    constant[[1]] <- rep ( constant[[1]] , dimension )
  }
  if ( length(constant[[1]])!=dimension )
  {
    stop('Each element of argument "constant" must be a scalar or ' ,
      'a vector of length "dimension".')
  }
  constant[[2]] <- as.numeric(constant[[2]])
  if ( length(constant[[2]])==1 )
  {
    constant[[2]] <- rep ( constant[[2]] , dimension )
  }
  if ( length(constant[[2]])!=dimension )
  {
    stop('Each element of argument "constant" must be a scalar or ' ,
      'a vector of length "dimension".')
  }
}
#}}}

# Generate eigenvalues#{{{

# Generate a block-diagonal matrix.
# Each block is a scalar "a",
# indicating a real eigen-value,
# or of form:
#     a b
#    -b a
# which indicates a pair of complex eigen-values: a+-bi
# The analytic solution is:
#     exp(at)
# for real eigen-value, and
#     exp(at)sin(bt)
#     exp(at)cos(bt)
# for a pair of complex eigen-values.

# "a"s are stored in "eigen_real",
# and "b"s in "eigen_imaginary".

# If there exists one real eigen-value,
# it is the last in "a".

# Generally "a" should be negative, or at least not too positive,
# to make the system stable.
# Each "a" is generated independently and uniform-randomly in range:
#     [real_min,real_max].
# The "real_min" value satisfies that the magnitude at the last
# time point is approximately half of that at the first.
# "real_max" is half of the lower-bound.

# Each "b" should be bounded away from each other
# to decrease correlation of the system.
# Currently they are chosen as an arithmetic sequence,
# plus a small random white noise.
# The smallest "b" is chosen to contain approximately
# one period in the time span when no real eigen-value exists,
# and double that value when there exists one real eigen-value.

num_real_eigen <- dimension%%2
num_complex_eigen <- floor(dimension/2)

time_span <- tail(timepoint,1) - timepoint[1]

real_min = -0.7/time_span
real_max = real_min/2

eigen_real <- runif (
  num_complex_eigen + num_real_eigen
  , real_min
  , real_max
)

eigen_imaginary_random_level <- 0.1
eigen_imaginary <- 1:num_complex_eigen
if ( num_real_eigen==1 )
{
  eigen_imaginary <- eigen_imaginary + 1
}
eigen_imaginary <- (
  eigen_imaginary
  + rnorm ( num_complex_eigen , 0 , eigen_imaginary_random_level )
)
eigen_imaginary <- eigen_imaginary * 2 * pi / time_span

require('permute')
permute_index <- permute::shuffle(dimension/2)
eigen_imaginary <- eigen_imaginary[permute_index]
#}}}

# Coefficient Matrix and Observation#{{{

temp <- (1:num_complex_eigen) * 2

ret$parameter$linear <- matrix ( 0 , dimension , dimension )

ret$parameter$linear [ cbind(temp,temp) ] <-
  eigen_real[1:num_complex_eigen]
ret$parameter$linear [ cbind(temp-1,temp-1) ] <-
  eigen_real[1:num_complex_eigen]
ret$parameter$linear [ cbind(temp-1,temp) ] <-
  eigen_imaginary
ret$parameter$linear [ cbind(temp,temp-1) ] <-
  -eigen_imaginary

ret$curve <- matrix ( 0 , length(timepoint) , dimension )

lapply ( 1 : num_complex_eigen ,
  function(index)
  {
    temp <- exp ( eigen_real[index] * timepoint )
    ret$curve[,2*index-1] <<-
      temp * sin ( eigen_imaginary[index] * timepoint )
    ret$curve[,2*index] <<-
      temp * cos ( eigen_imaginary[index] * timepoint )
  }
)

if ( num_real_eigen == 1 )
{
  ret$parameter$linear[dimension,dimension] <- tail(eigen_real,1)
  ret$curve[,dimension] <- exp ( tail(eigen_real,1) * timepoint )
}
#}}}

# Orthogonal Transformation#{{{

orthogonal_transformation <- as.list(orthogonal_transformation)
for ( item in orthogonal_transformation )
{
  item <- as.integer(item)
  if ( length(item)!=2 )
  {
    stop('Each element of argument "orthogonal_transformation" '
      ,'must be length 2.')
  }
  if ( item[1]>=item[2] )
  {
    stop('In each element of argument "orthogonal_transformation", ' ,
      'The second component must be larger than first.'
    )
  }
  if ( item[1]<1 || item[2]>dimension )
  {
    stop('Out-of-bound index in elements of '
      ,'argument "orthogonal_transformation".')
  }

  require('pracma')
  temp <- diag(dimension)
  temp [ item[1]:item[2] , item[1]:item[2] ] <-
    pracma::rortho ( item[2]-item[1]+1 )
  ret$curve <- ret$curve %*% t(temp)
  ret$parameter$linear <- temp %*% ret$parameter$linear %*% t(temp)
}
#}}}

# Row-Column Permutation#{{{

# Permute rows and columns
# by left multiplying a permutation matrix
# and right multiplying its transpose
# to the coefficient matrix.
# This will make the sparsity structure less obvious,
# however it does not change the property,
# which means the system is still unconnected,
# if it is unconnected before this permutation.

if ( row_column_permutation == TRUE )
{
  require('permute')
  permute_index <- permute::shuffle ( dimension )
  ret$parameter$linear <-
    ret$parameter$linear [ permute_index , permute_index ]
  ret$curve <- ret$curve [ , permute_index ]
}
#}}}

# Constant#{{{

ret$constant <- NULL

ret$curve <- t(ret$curve)

if ( constant_exist )
{
  ret$parameter$constant <-
    runif ( dimension , constant[[1]] , constant[[2]] )
  temp <- solve ( ret$parameter$linear , ret$parameter$constant )
  lapply ( 1 : length(time_point) , function(index)
  {
    ret$curve[,index] <<- ret$curve[,index] - temp
    return()
  } )
}

#}}}

# Return#{{{

ret$time <- timepoint
ret$parameter$initial <- ret$curve[,1]

return(ret)
#}}}

}
