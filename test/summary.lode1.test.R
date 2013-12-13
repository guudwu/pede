library('pede')

set.seed(0)
dimension <- 5
time_point <- 0:10

# Generate model
truth <-
  generation.lode1 (
    dimension
    , time_point
    , list()
    , TRUE
    , list(0,1e-2)
  )

summary.lode1(truth)
