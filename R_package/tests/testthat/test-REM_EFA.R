# Holzinger data
library(lavaan)
HS.full <- HolzingerSwineford1939
HS.df <- HolzingerSwineford1939[,c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9")]

# test case when data contains non-numeric entries
expect_error(REM_EFA(X = HS.full, delta = 0, k_range = 2))

# test when delta < 0
expect_error(REM_EFA(X = HS.df, delta = -1, k_range = 2))

# test when delta > 1
expect_error(REM_EFA(X = HS.df, delta = 10, k_range = 2))

# test when k_range < 1
expect_error(REM_EFA(X = HS.df, delta = 0.05, k_range = 0))

# test when k_range < 1 AND delta is < 0
expect_error(REM_EFA(X = HS.df, delta = -10, k_range = 0))

# test when k_range < 1 AND delta is > 1
expect_error(REM_EFA(X = HS.df, delta = 10, k_range = 0))

# test when k_range is non-integer
expect_error(REM_EFA(X = HS.df, delta = 10, k_range = 0.3))
