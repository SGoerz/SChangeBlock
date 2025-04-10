test_that("alternatives returns the correct format", {
  expect_equal(class(alternatives(c(50, 50), type = "1a")), "integer")
  expect_equal(class(alternatives(c(50, 50), type = "1b")), "integer")
  expect_equal(class(alternatives(c(50, 50), type = "1c")), "integer")
  expect_equal(class(alternatives(c(50, 50), type = 2)), "integer")
  expect_equal(class(alternatives(c(50, 50), type = 4)), "integer")
  expect_equal(class(alternatives(c(50, 50), type = 5)), "integer")
  
  expect_equal(class(alternatives(c(50, 50), type = 3)), "numeric")
  expect_length(alternatives(c(50, 50), type = 3), 50^2)
})
