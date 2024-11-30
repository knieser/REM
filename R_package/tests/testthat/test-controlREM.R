
# steps must be a positive integer greater than 1
expect_error(controlREM(steps = 0))
expect_error(controlREM(steps = 3.5))
expect_error(controlREM(steps = -3))

# tol values should be between 0 and 1
expect_error(controlREM(tol = 0))
expect_error(controlREM(tol = 10))

# maxiter must be greater than 2
expect_error(controlREM(maxiter = 1))
expect_error(controlREM(maxiter = 0))
expect_error(controlREM(maxiter = -1))

# min weights must be between 0 and 1
expect_error(controlREM(min_weights = 0))
expect_error(controlREM(min_weights = 1))

# max ueps must be between 0 and 1
expect_error(controlREM(max_ueps = -1))
expect_error(controlREM(max_ueps = 1.1))

# chk_gamma must be between 0 and 1
expect_error(controlREM(chk_gamma = 1))
expect_error(controlREM(chk_gamma = 0))
expect_error(controlREM(chk_gamma = -1))

