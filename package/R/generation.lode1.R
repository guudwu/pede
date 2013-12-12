generation.lode1 <- function (
  dimension
  , timepoint
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

# OUTPUT:
# time: Input argument "timepoint".
# parameter$initial: Initial condition of the system.
# parameter$linear: Linear term of the system.
#   It will be identical to the first column of "curve".
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

# Return#{{{

ret$time <- timepoint
ret$curve <- t(ret$curve)
ret$parameter$initial <- ret$curve[,1]

return(ret)
#}}}

}
