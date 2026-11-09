test_that("Graphics returns the correct output", {
  X <- genField(c(20, 20))
  p <- plot(X)
  
  expect_true(is(p, "ggplot"))
  expect_type(p, "object")
})