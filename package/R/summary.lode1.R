# Check property of generated ODE system.

# Currently only two properties are checked:
# Pair-wise correlation of curves.
# Variance Inflation Factor.

summary.lode1 <- function
(
  lode1
  , correlation = TRUE
  , vif = TRUE
  , verbose = 1
)

# INPUT:
# lode1: An order-1 linear ODE system.
#   Return value of "generation.lode1" can be directly used.
#   It should contain element "curve".
# correlation: Whether to check pair-wise correlation.
# vif: Whether to check VIF.
# verbose: If positive, output the property results to STDOUT.

{

ret <- lode1

# Sanity check#{{{

if ( is.null(ret$curve) )
{
  stop('Argument "lode1" should contain element "curve".')
}
if ( is.null(dim(ret$curve)) )
{
  stop('Element "curve" of argument "lode1" '
    ,'should be a numeric matrix.')
}
ret$curve <- as.numeric(ret$curve)
dim(ret$curve) <- dim(lode1$curve)

correlation <- as.logical(correlation)
if ( length(correlation)!=1 )
{
  stop('Argument "correlation" should be a logical scalar.')
}

vif <- as.logical(vif)
if ( length(vif)!=1 )
{
  stop('Argument "vif" should be a logical scalar.')
}

verbose <- as.integer(verbose)
if ( length(verbose)!=1 )
{
  stop('Argument "verbose" should be an integer scalar.')
}
#}}}

# Correlation#{{{

if ( correlation )
{
  if ( !is.null(ret$property$correlation) )
  {
    warning('Correlation property exists.  Skip repeated computing.')
  }
  else
  {
    ret$property$correlation <-
      cor ( t(lode1$curve) )
    ret$property$max_correlation <-
      max ( abs (
        ret$property$correlation - diag(nrow(lode1$curve))
      ) )
  }

  if ( verbose >= 1 )
  {
    cat('Pair-wise correlation matrix: \n')
    print(ret$property$correlation)
    cat('Maximal correlation: ')
    cat(ret$property$max_correlation)
    cat('\n')
  }
}
#}}}

# VIF#{{{

if ( vif )
{
  if ( !is.null(ret$property$vif) )
  {
    warning('VIF property exists.  Skip repeated computing.')
  }
  else
  {
    require('HH')
    temp <- lapply ( 1 : nrow(lode1$curve) , function(index)
    {
      return ( lode1$curve[index,] )
    } )
    ret$property$vif <-
      as.vector (
        HH::vif ( data.frame(temp) )
      )
  }

  if ( verbose >= 1 )
  {
    cat('Variance Inflation Factor: \n')
    cat(ret$property$vif)
    cat('\n')
  }
}
#}}}

# Return#{{{

return(ret)
#}}}

}
