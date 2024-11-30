# Holzinger data
library(lavaan)
HS.full <- HolzingerSwineford1939
HS.df <- HolzingerSwineford1939[,c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9")]
HS.model <-  'visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9'

# test case when data contains non-numeric columns
expect_error(REM_CFA(X = HS.full, model = HS.model))

# test case when dimension of model does not match dimension of data
expect_error(REM_CFA(X = HS.df[,-1], model = HS.model))

# test case when delta < 0
expect_error(REM_CFA(X = HS.df, delta = -0.05, model = HS.model))

# test case when delta > 1
expect_error(REM_CFA(X = HS.df, delta = 1.05, model = HS.model))

# test case with unexpected model form
HS.model2 <-  'visual  =~ x1 + 2 x2 + x3
              speed   =~ x7 + x8 + x9
              textual =~ x4 + x5+ x6'
expect_error(REM_CFA(X = HS.df, model = HS.model2))

