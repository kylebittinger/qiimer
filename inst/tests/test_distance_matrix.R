context('Paired distances')

sample_names <- c("a", "b", "c", "d")
distmat <- matrix(
  c(0, 1, 2, 3,
    1, 0, 4, 5,
    2, 4, 0, 6,
    3, 5, 6, 0), 
  ncol=4, 
  dimnames=list(sample_names, sample_names))
p_levels <- c("Within Pairs", "Between Pairs")

test_that("paired_distances labels pairs correctly", {
  pairs <- paired_distances(distmat, c("a", "b"), c("c", "d"))

  category_ac <- with(pairs, Category[(SampleA == "a") & (SampleB == "c")])
  expect_that(category_ac, equals(factor("Within Pairs", levels=p_levels)))

  category_ad <- with(pairs, Category[(SampleA == "a") & (SampleB == "d")])
  expect_that(category_ad, equals(factor("Between Pairs", levels=p_levels)))
})

test_that("paired_distances works with length-1 sample lists", {
  pairs <- paired_distances(distmat, "a", "b")

  expect_that(pairs$SampleA, equals(factor("a")))
  expect_that(pairs$SampleB, equals(factor("b")))
  expect_that(pairs$Category, equals(factor("Within Pairs", levels=p_levels)))
  expect_that(pairs$Distance, equals(1))
})

context("Grouped distances")

g_levels <- c("Within1", "Between", "Within2")

test_that("grouped_distances labels groups correctly", {
  groups <- grouped_distances(distmat, c("a", "b"), c("c", "d"))
  
  category_ab <- with(groups, Category[(SampleA == "a") & (SampleB == "b")])
  expect_that(category_ab, equals(factor("Within1", levels=g_levels)))

  category_ac <- with(groups, Category[(SampleA == "a") & (SampleB == "c")])
  expect_that(category_ac, equals(factor("Between", levels=g_levels)))
  
  category_cd <- with(groups, Category[(SampleA == "c") & (SampleB == "d")])
  expect_that(category_cd, equals(factor("Within2", levels=g_levels)))
})

test_that("grouped_distances works with length-1 sample lists", {
  groups <- grouped_distances(distmat, "a", "b")
  
  expect_that(groups$SampleA, equals(factor("a")))
  expect_that(groups$SampleB, equals(factor("b")))
  expect_that(groups$Category, equals(factor("Between", levels=g_levels)))
  expect_that(groups$Distance, equals(1))
})
