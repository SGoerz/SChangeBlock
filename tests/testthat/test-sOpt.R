test_that("sOpt returns the correct exponent", {
  n <- c(10, 20, 50, 100, 200, 400)
  s.opt <- sOpt(n, 0.6)
  l <- round(n^s.opt)
  b <- floor(n / l)
  
  expect_equal(l * b, n)
})

test_that("build_blocks output has the correct format", {
  s <- build_blocks(c(20, 20), c(4, 4))
  expect_equal(dim(s), c(20, 20))
  
  s <- build_blocks(c(20, 20), c(3, 3))
  expect_equal(dim(s), c(20, 20))
  
  s <- build_blocks(20, 2)
  expect_equal(length(s), 20)
})

test_that("build_blocks produces the correct output", {
  s <- build_blocks(c(20, 20), c(2, 2))
  expect_equal(s, kronecker(matrix(1:4, ncol = 2), matrix(1, ncol = 10, nrow = 10)))
  
  s <- build_blocks(20, 3)
  expect_equal(s, matrix(c(rep(1, 7), rep(2, 7), rep(3, 6))))
})