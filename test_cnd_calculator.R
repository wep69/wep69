# Load testthat library
library(testthat)

# Source the R CND calculator functions
# Assuming test_cnd_calculator.R is in the same directory as cnd_calculator.R
# In a real package structure, you'd use devtools::load_all() or similar.
source("cnd_calculator.R")

# Define a small tolerance for floating point comparisons
TOL <- 1e-7

context("Test calculate_Rd_R Function")

test_that("calculate_Rd_R works with typical list input", {
  concentrations <- list(N=3, P=0.5, K=2)
  expect_equal(calculate_Rd_R(concentrations), 100 - (3+0.5+2), tolerance = TOL)
})

test_that("calculate_Rd_R works with typical vector input", {
  concentrations <- c(N=3, P=0.5, K=2)
  expect_equal(calculate_Rd_R(concentrations), 100 - (3+0.5+2), tolerance = TOL)
  expect_named(calculate_Rd_R(concentrations), NULL) # Rd is a single numeric value
})

test_that("calculate_Rd_R handles empty input", {
  expect_equal(calculate_Rd_R(list()), 100, tolerance = TOL)
  expect_equal(calculate_Rd_R(c()), 100, tolerance = TOL)
  expect_equal(calculate_Rd_R(NULL), 100, tolerance = TOL)
})

test_that("calculate_Rd_R handles single nutrient input", {
  expect_equal(calculate_Rd_R(list(N=5)), 95, tolerance = TOL)
  expect_equal(calculate_Rd_R(c(N=5)), 95, tolerance = TOL)
})

test_that("calculate_Rd_R handles input with NA", {
  concentrations_list <- list(N=3, P=NA, K=2)
  expect_equal(calculate_Rd_R(concentrations_list), 100 - (3+2), tolerance = TOL) # NA is ignored
  concentrations_vec <- c(N=3, P=NA, K=2)
  expect_equal(calculate_Rd_R(concentrations_vec), 100 - (3+2), tolerance = TOL)
})

test_that("calculate_Rd_R handles all NA input", {
  expect_equal(calculate_Rd_R(list(N=NA, P=NA)), 100, tolerance = TOL)
})

test_that("calculate_Rd_R handles Rd = 0 and Rd < 0", {
  expect_equal(calculate_Rd_R(list(N=100)), 0, tolerance = TOL)
  expect_equal(calculate_Rd_R(list(N=110)), -10, tolerance = TOL)
})


context("Test calculate_geometric_mean_R Function")

test_that("calculate_geometric_mean_R with typical input", {
  data_list <- list(N=3.0, P=0.5, K=2.0, Rd=94.5)
  expected_list <- (3.0 * 0.5 * 2.0 * 94.5)^(1/4)
  expect_equal(calculate_geometric_mean_R(data_list), expected_list, tolerance = TOL)

  data_vec <- c(N=3.0, P=0.5, K=2.0, Rd=94.5)
  expected_vec <- (3.0 * 0.5 * 2.0 * 94.5)^(1/4)
  expect_equal(calculate_geometric_mean_R(data_vec), expected_vec, tolerance = TOL)
})

test_that("calculate_geometric_mean_R with single component", {
  expect_equal(calculate_geometric_mean_R(list(N=10.0)), 10.0, tolerance = TOL)
  expect_equal(calculate_geometric_mean_R(c(N=10.0)), 10.0, tolerance = TOL)
})

test_that("calculate_geometric_mean_R with zero concentration", {
  # If any concentration is 0, G should be 0
  expect_equal(calculate_geometric_mean_R(list(N=3.0, P=0, K=2.0, Rd=95.0)), 0, tolerance = TOL)
})

test_that("calculate_geometric_mean_R with empty input", {
  # Should return NaN as per function spec
  expect_true(is.nan(calculate_geometric_mean_R(list())))
  expect_true(is.nan(calculate_geometric_mean_R(c())))
})

test_that("calculate_geometric_mean_R with NA input", {
  # Should return NA_real_ as per function spec
  expect_true(is.na(calculate_geometric_mean_R(list(N=10, P=NA))))
  expect_true(is.na(calculate_geometric_mean_R(c(N=10, P=NA))))
})

test_that("calculate_geometric_mean_R with negative numbers", {
  # Produces NaN for even root of negative product
  expect_warning(result_neg1 <- calculate_geometric_mean_R(list(N=-2, P=4, K=1, Rd=1)), "Negative concentrations found") # (-8)^(1/4) -> NaN
  expect_true(is.nan(result_neg1))

  # Produces real for odd root of negative product
  expect_warning(result_neg2 <- calculate_geometric_mean_R(list(N=-8, P=1, K=1)), "Negative concentrations found") # (-8)^(1/3) -> -2
  expect_equal(result_neg2, -2, tolerance = TOL)
  
  # Produces real for even root of positive product (from two negatives)
  expect_warning(result_neg3 <- calculate_geometric_mean_R(list(N=-4, P=-2, K=1, Rd=1)), "Negative concentrations found") # (8)^(1/4)
  expect_equal(result_neg3, 8^(1/4), tolerance = TOL)
})


context("Test calculate_V_values_R Function")

test_that("calculate_V_values_R with typical input (list and vector)", {
  concs_list <- list(N=3.0, P=0.5, K=2.0, Rd=94.5)
  g_val <- calculate_geometric_mean_R(concs_list)
  v_list <- calculate_V_values_R(concs_list, g_val)
  expect_type(v_list, "list")
  expect_named(v_list, names(concs_list))
  expect_equal(v_list$N, log(3.0/g_val), tolerance = TOL)
  expect_equal(v_list$P, log(0.5/g_val), tolerance = TOL)

  concs_vec <- c(N=3.0, P=0.5, K=2.0, Rd=94.5)
  g_val_vec <- calculate_geometric_mean_R(concs_vec)
  v_vec <- calculate_V_values_R(concs_vec, g_val_vec)
  expect_type(v_vec, "double") # unlist returns a vector
  expect_named(v_vec, names(concs_vec))
  expect_equal(v_vec["N"], log(3.0/g_val_vec), tolerance = TOL)
})

test_that("calculate_V_values_R when G is 0, NA, or NaN", {
  concs <- list(N=1, P=2)
  v_g_zero <- calculate_V_values_R(concs, 0)
  expect_true(all(sapply(v_g_zero, is.nan)))
  expect_named(v_g_zero, names(concs))

  v_g_na <- calculate_V_values_R(concs, NA_real_)
  expect_true(all(sapply(v_g_na, is.nan)))
  
  v_g_nan <- calculate_V_values_R(concs, NaN)
  expect_true(all(sapply(v_g_nan, is.nan)))
})

test_that("calculate_V_values_R with concentration 0, negative, or NA", {
  g_val <- 10.0 # Assume a valid G
  
  concs_zero <- list(N=0, P=1, K=2)
  v_zero <- calculate_V_values_R(concs_zero, g_val)
  expect_true(is.nan(v_zero$N))
  expect_false(is.nan(v_zero$P))

  concs_neg <- list(N=-1, P=1, K=2)
  v_neg <- calculate_V_values_R(concs_neg, g_val)
  expect_true(is.nan(v_neg$N))
  expect_false(is.nan(v_neg$P))
  
  concs_na <- list(N=NA, P=1, K=2)
  v_na <- calculate_V_values_R(concs_na, g_val)
  expect_true(is.nan(v_na$N))
  expect_false(is.nan(v_na$P))
})

test_that("calculate_V_values_R when G = 1", {
  concs <- list(N=3.0, P=0.5)
  v_vals <- calculate_V_values_R(concs, 1.0)
  expect_equal(v_vals$N, log(3.0), tolerance = TOL)
  expect_equal(v_vals$P, log(0.5), tolerance = TOL)
})

test_that("calculate_V_values_R with empty input", {
  expect_equal(calculate_V_values_R(list(), 10.0), list())
  expect_equal(calculate_V_values_R(c(), 10.0), vector("double")) # unlist of empty c() is NULL, mode is NULL, vector("NULL") error
                                                                  # function returns vector(mode(concentrations_vec))
                                                                  # mode(unlist(c())) is "NULL", so vector("NULL")
                                                                  # The R function returns vector(mode(concentrations_vec)) which is logical() for c()
  expect_length(calculate_V_values_R(c(), 10.0),0)
  expect_type(calculate_V_values_R(c(), 10.0), "logical") # because mode(unlist(c())) is "logical" in R 4.x
})


context("Test calculate_CND_indices_R Function")

test_that("calculate_CND_indices_R with typical input", {
  v_vals <- list(N=1.0, P=-0.5, K=0.0)
  v_star <- list(N=0.8, P=-0.6, K=0.1)
  sd_star <- list(N=0.2, P=0.1, K=0.05)
  
  indices <- calculate_CND_indices_R(v_vals, v_star, sd_star)
  expect_type(indices, "list")
  expect_named(indices, names(v_vals))
  expect_equal(indices$N, (1.0 - 0.8) / 0.2, tolerance = TOL)
  expect_equal(indices$P, (-0.5 - (-0.6)) / 0.1, tolerance = TOL)
  expect_equal(indices$K, (0.0 - 0.1) / 0.05, tolerance = TOL)
})

test_that("calculate_CND_indices_R when SD_star is 0", {
  v_vals <- list(N=1.0)
  v_star <- list(N=0.8)
  sd_star_zero <- list(N=0)
  indices <- calculate_CND_indices_R(v_vals, v_star, sd_star_zero)
  expect_true(is.nan(indices$N))
})

test_that("calculate_CND_indices_R with missing reference data", {
  v_vals <- list(N=1.0, P=-0.5) # P is not in v_star_missing
  v_star_missing <- list(N=0.8)
  sd_star_missing <- list(N=0.2)
  indices <- calculate_CND_indices_R(v_vals, v_star_missing, sd_star_missing)
  expect_false(is.nan(indices$N))
  expect_true(is.nan(indices$P)) # P missing in ref
})

test_that("calculate_CND_indices_R when V_value is NA or NaN", {
  v_vals_na <- list(N=NA, P=-0.5)
  v_vals_nan <- list(N=NaN, P=-0.5)
  v_star <- list(N=0.8, P=-0.6)
  sd_star <- list(N=0.2, P=0.1)

  indices_na <- calculate_CND_indices_R(v_vals_na, v_star, sd_star)
  expect_true(is.nan(indices_na$N))
  expect_false(is.nan(indices_na$P))

  indices_nan <- calculate_CND_indices_R(v_vals_nan, v_star, sd_star)
  expect_true(is.nan(indices_nan$N))
})

test_that("calculate_CND_indices_R with empty input", {
  expect_equal(calculate_CND_indices_R(list(), list(), list()), list())
  expect_length(calculate_CND_indices_R(c(), c(), c()),0)
  expect_type(calculate_CND_indices_R(c(), c(), c()), "logical") # Based on typeof(V_values_vec)
})

test_that("calculate_CND_indices_R with reference data NA/NaN", {
  v_vals <- list(N=1.0)
  v_star_na <- list(N=NA)
  sd_star_na <- list(N=NA)
  v_star_ok <- list(N=0.8)
  sd_star_ok <- list(N=0.1)

  expect_true(is.nan(calculate_CND_indices_R(v_vals, v_star_na, sd_star_ok)$N))
  expect_true(is.nan(calculate_CND_indices_R(v_vals, v_star_ok, sd_star_na)$N))
})


context("Test calculate_CND_r_squared_R Function")

test_that("calculate_CND_r_squared_R with typical input", {
  indices <- list(N=1.0, P=-2.0, K=0.5)
  expected <- 1.0^2 + (-2.0)^2 + 0.5^2
  expect_equal(calculate_CND_r_squared_R(indices), expected, tolerance = TOL)
})

test_that("calculate_CND_r_squared_R with one index NA/NaN", {
  indices_na <- list(N=1.0, P=NA, K=0.5)
  expected_na <- 1.0^2 + 0.5^2
  expect_equal(calculate_CND_r_squared_R(indices_na), expected_na, tolerance = TOL)

  indices_nan <- list(N=1.0, P=NaN, K=0.5)
  expected_nan <- 1.0^2 + 0.5^2
  expect_equal(calculate_CND_r_squared_R(indices_nan), expected_nan, tolerance = TOL)
})

test_that("calculate_CND_r_squared_R with all indices NA/NaN", {
  expect_equal(calculate_CND_r_squared_R(list(N=NA, P=NA)), 0, tolerance = TOL)
  expect_equal(calculate_CND_r_squared_R(list(N=NaN, P=NaN)), 0, tolerance = TOL)
})

test_that("calculate_CND_r_squared_R with empty input", {
  expect_equal(calculate_CND_r_squared_R(list()), 0, tolerance = TOL)
  expect_equal(calculate_CND_r_squared_R(c()), 0, tolerance = TOL)
})


context("Test calculate_cnd_R Main Function")

test_that("calculate_cnd_R comprehensive test (list input)", {
  sample_conc <- list(N=2.8, P=0.25, K=2.5) 
  # Rd = 100 - (2.8 + 0.25 + 2.5) = 100 - 5.55 = 94.45
  # nutrient_concentrations_with_Rd = list(N=2.8, P=0.25, K=2.5, Rd=94.45)
  # G = (2.8 * 0.25 * 2.5 * 94.45)^(1/4) approx 3.585900199
  # V_N = log(2.8 / G)   approx -0.247511
  # V_P = log(0.25 / G)  approx -2.663697
  # V_K = log(2.5 / G)   approx -0.360787
  # V_Rd = log(94.45 / G) approx 3.271995

  high_yield_means <- list(N= -0.20, P= -2.70, K= -0.30, Rd= 3.25)
  high_yield_sds <- list(N= 0.1, P= 0.15, K= 0.12, Rd= 0.10)

  # Expected CND Indices:
  # I_N = (-0.247511 - (-0.20)) / 0.1 = -0.47511
  # I_P = (-2.663697 - (-2.70)) / 0.15 = 0.24202
  # I_K = (-0.360787 - (-0.30)) / 0.12 = -0.506558
  # I_Rd = (3.271995 - 3.25) / 0.10 = 0.21995
  
  # Expected r_squared = (-0.47511)^2 + (0.24202)^2 + (-0.506558)^2 + (0.21995)^2
  #                   = 0.2257295 + 0.0585737 + 0.256601 + 0.048378 approx 0.5892822

  result <- calculate_cnd_R(sample_conc, high_yield_means, high_yield_sds)
  
  expect_type(result, "list")
  expect_named(result, c("CND_indices", "r_squared"))
  expect_type(result$CND_indices, "list")
  expect_named(result$CND_indices, c("N", "P", "K", "Rd")) # Order from nutrient_concentrations_with_Rd

  expect_equal(result$CND_indices$N, -0.47511, tolerance = 1e-5)
  expect_equal(result$CND_indices$P, 0.24202, tolerance = 1e-5)
  expect_equal(result$CND_indices$K, -0.506558, tolerance = 1e-5)
  expect_equal(result$CND_indices$Rd, 0.21995, tolerance = 1e-5)
  expect_equal(result$r_squared, 0.5892822, tolerance = 1e-5)
})

test_that("calculate_cnd_R with empty sample concentrations (vector input)", {
  empty_sample_conc <- c()
  high_yield_means <- c(N= -0.20, P= -2.70, K= -0.30, Rd= 3.25) # Must have Rd
  high_yield_sds <- c(N= 0.1, P= 0.15, K= 0.12, Rd= 0.10)   # Must have Rd

  # Rd = 100
  # G = 100
  # V_Rd = log(100/100) = 0
  # I_Rd = (0 - 3.25) / 0.10 = -32.5
  # r_squared = (-32.5)^2 = 1056.25
  # Other nutrients N, P, K are not in sample, so they won't be in V_values, 
  # and thus not in CND_indices from the loop in calculate_CND_indices_R
  
  result <- calculate_cnd_R(empty_sample_conc, high_yield_means, high_yield_sds)
  
  expect_type(result$CND_indices, "double") # vector
  expect_named(result$CND_indices, "Rd") # Only Rd should be present
  expect_length(result$CND_indices, 1)
  expect_equal(result$CND_indices[["Rd"]], -32.5, tolerance = TOL)
  expect_equal(result$r_squared, (-32.5)^2, tolerance = TOL)
})

test_that("calculate_cnd_R propagation of NaN if G is NaN", {
  # If G is NaN (e.g., empty concentrations_with_Rd, which calculate_geometric_mean_R handles by returning NaN)
  # then V_values will be all NaN, CND_indices all NaN, r_squared will be 0.
  
  # To force G to be NaN, we can pass empty list to calculate_geometric_mean_R
  # but main function constructs nutrient_concentrations_with_Rd_vec.
  # If sample_nutrient_concentrations is empty, Rd is 100, G is 100.
  # If sample_nutrient_concentrations has an NA, G is NA.
  
  sample_na <- list(N=NA, P=1) # This will make G = NA_real_
  high_yield_means <- list(N=0.1, P=-1.5, Rd=3.0)
  high_yield_sds <- list(N=0.2, P=0.1, Rd=0.3)
  
  result <- calculate_cnd_R(sample_na, high_yield_means, high_yield_sds)
  
  # Rd = 100 - 1 = 99
  # nutrient_concentrations_with_Rd = list(N=NA, P=1, Rd=99)
  # G = NA (due to N=NA)
  # V_values = list(N=NaN, P=NaN, Rd=NaN) (because G is NA)
  # CND_indices = list(N=NaN, P=NaN, Rd=NaN)
  # r_squared = 0
  
  expect_true(is.nan(result$CND_indices$N))
  expect_true(is.nan(result$CND_indices$P))
  expect_true(is.nan(result$CND_indices$Rd)) # Rd also becomes NaN because G is NA
  expect_equal(result$r_squared, 0, tolerance = TOL)
})

test_that("calculate_cnd_R nutrient in sample but not in reference data", {
  sample_conc = list(N=2.8, P=0.25, XTRA=1.0) # XTRA not in reference
  high_yield_means <- list(N= -0.20, P= -2.70, Rd= 3.25)
  high_yield_sds <- list(N= 0.1, P= 0.15, Rd= 0.10)
  
  result <- calculate_cnd_R(sample_conc, high_yield_means, high_yield_sds)
  
  expect_true(is.nan(result$CND_indices$XTRA))
  expect_false(is.nan(result$CND_indices$N))
  expect_false(is.nan(result$CND_indices$P))
  expect_false(is.nan(result$CND_indices$Rd))
  
  # r_squared should be calculated from N, P, Rd only
  expected_r_sq <- result$CND_indices$N^2 + result$CND_indices$P^2 + result$CND_indices$Rd^2
  expect_equal(result$r_squared, expected_r_sq, tolerance = TOL)
})
