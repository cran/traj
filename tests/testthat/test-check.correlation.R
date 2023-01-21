test_that("test correlation", {
  output = matrix(c(.5, .6, .2, .7), nr = 2)
  expect_equal(check.correlation(output), NULL)
})