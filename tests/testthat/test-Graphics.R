test_that("Graphics returns the correct output", {
  X <- genField(c(40, 20))
  p <- plot(X, facet_id = factor(rep(1:2, each = 20^2)))
  
  expect_true(is(p, "ggplot"))
  expect_type(p, "object")

  x <- arima.sim(list(ar = 0.3), 100)
  class(x) <- "RandomField"
  
  p <- plot(x, facet_id = factor(rep(1:2, 50)), main = "XYZ", 
       p = theme(panel.grid.major = element_line(color = "white"))) +
    theme(panel.background = element_rect(fill = "grey", color = NA))
  
  expect_true(is(p, "ggplot"))
  expect_type(p, "object")
})