library('pede')

set.seed(0)
dimension <- 5
time_point <- 0:10

# Generate model
truth <-
  generation.lode1 (
    dimension
    , time_point
    , list(c(2,3),c(4,5))
  )

# Check curve with internal "ode" function.
linODE <- function ( time , state , pars )
{
  res <- pars[[1]] %*% state
  if ( ! is.null ( pars[[2]] ) )
    res <- res + pars[[2]]
  return ( list(res) )
}

library('deSolve')
ode_res <- deSolve::ode (
  truth$parameter$initial
  , truth$time
  , linODE
  , list ( truth$parameter$linear , NULL )
)

rss <- norm ( t(ode_res[,-1]) - truth$curve , 'F' )
cat('rss:')
cat(rss)
cat('\n')
print(truth)
