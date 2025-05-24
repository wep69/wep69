# Functions for Compositional Nutrient Diagnosis (CND) in R

#' Calculate Residual Value (Rd)
#'
#' Calculates Rd = 100 - sum of nutrient concentrations.
#' Handles empty or NULL input by returning 100.
#' Ignores NA values in the sum.
#'
#' @param nutrient_concentrations A named list or vector of nutrient concentrations
#'                                (e.g., `list(N=3.5, P=0.4)` or `c(N=3.5, P=0.4)`).
#' @return The calculated Rd value (numeric).
#' @export
#' @examples
#' calculate_Rd_R(list(N=3, P=0.5, K=2)) # Expected: 100 - (3+0.5+2) = 94.5
#' calculate_Rd_R(c(N=3, P=0.5, K=2))   # Expected: 94.5
#' calculate_Rd_R(list())                # Expected: 100
#' calculate_Rd_R(NULL)                  # Expected: 100
#' calculate_Rd_R(c())                   # Expected: 100
#' calculate_Rd_R(list(N=3, P=NA, K=2))  # Expected: 100 - (3+2) = 95
calculate_Rd_R <- function(nutrient_concentrations) {
  # Check for NULL or empty (length 0) input
  if (is.null(nutrient_concentrations) || length(nutrient_concentrations) == 0) {
    return(100)
  }
  
  # Calculate sum of concentrations, ignoring NA values
  # unlist is important if input is a list
  sum_of_concentrations <- sum(unlist(nutrient_concentrations), na.rm = TRUE)
  
  Rd <- 100 - sum_of_concentrations
  return(Rd)
}

#' Calculate Geometric Mean (G)
#'
#' Calculates G = (product of concentrations)^(1 / number of components).
#' Concentrations are typically expected to be positive.
#'
#' @param nutrient_concentrations_with_Rd A named list or vector of nutrient concentrations including Rd.
#' @return The calculated G value (numeric). Returns `NaN` if input is empty,
#'         or if calculation results in a non-real number (e.g., fractional power of a negative product).
#'         Returns `NA_real_` if any input concentration is `NA`. Returns `0` if any concentration is `0`.
#' @export
#' @examples
#' calculate_geometric_mean_R(list(N=30, P=20, K=40, Rd=10)) # Expected: (30*20*40*10)^(1/4)
#' calculate_geometric_mean_R(c(N=0, P=1, K=1, Rd=1))      # Expected: 0
#' calculate_geometric_mean_R(list())                       # Expected: NaN
#' calculate_geometric_mean_R(list(N=10, P=NA))             # Expected: NA_real_
#' calculate_geometric_mean_R(list(N=-2, P=4, K=1, Rd=1))   # Expected: NaN (sqrt of negative)
#' calculate_geometric_mean_R(list(N=-8, P=1, K=1))         # Expected: -2 ((-8)^(1/3)) (if simplified)
#'                                                          # R's default (-8)^(1/3) is -2
#' calculate_geometric_mean_R(list(N=-4, P=-2, K=1, Rd=1))  # Expected: (8)^(1/4) approx 1.681793
calculate_geometric_mean_R <- function(nutrient_concentrations_with_Rd) {
  concentrations <- unlist(nutrient_concentrations_with_Rd)
  
  if (length(concentrations) == 0) {
    return(NaN) 
  }
  
  if (any(is.na(concentrations))) {
    return(NA_real_)
  }

  if (any(concentrations < 0)) {
    warning("Negative concentrations found. Geometric mean calculation may result in NaN if product is negative and root is even.")
  }

  product_of_concentrations <- prod(concentrations) # na.rm = FALSE by default
  num_components <- length(concentrations)
  
  # If product is 0, G is 0, unless num_components is 0 (already handled)
  if (product_of_concentrations == 0) {
    return(0)
  }
  
  # Standard calculation for G
  # R handles fractional powers of negative numbers correctly if the result is real (e.g. (-8)^(1/3) = -2)
  # If the result is complex (e.g. (-4)^(1/2)), it returns NaN, which is desired.
  G <- product_of_concentrations^(1 / num_components)
  
  # Ensure that if product was negative and root was even, result is NaN, not potentially a complex number's magnitude or principal root.
  # R's default behavior for `^` already handles this by returning NaN for e.g. (-2)^(0.5)
  # However, an explicit check for clarity or to enforce specific behavior might be:
  # if (product_of_concentrations < 0 && num_components %% 2 == 0 && (1/num_components) %% 1 != 0.5) {
  # This check is tricky because (1/num_components) %% 1 != 0 means it's not an integer root.
  # And num_components %% 2 == 0 means an even root.
  # R's `^` operator is generally reliable for this. For example:
  # (-4)^(1/2) is NaN
  # (-8)^(1/3) is -2
  # (-16)^(1/4) is NaN (R gives NaN for even roots of negative numbers)
  
  return(G)
}

#' Calculate Centered Log-Ratio (CLR) Transformed Variables (V values)
#'
#' Calculates V_nutrient = ln(concentration / G) for each nutrient.
#' ln is the natural logarithm.
#'
#' @param nutrient_concentrations_with_Rd A named list or vector of nutrient concentrations including Rd.
#' @param geometric_mean The geometric mean (G) of all components, a single numeric value.
#' @return A named list or vector of V values, matching input type.
#'         V values will be `NaN` if G is invalid (0, NA, NaN) or if the
#'         corresponding concentration is not positive (<=0, NA).
#' @export
#' @examples
#' concs1 <- list(N=30, P=20, K=40, Rd=10) # G = (240000)^(1/4) approx 22.13364
#' G1 <- calculate_geometric_mean_R(concs1)
#' calculate_V_values_R(concs1, G1)
#' # Expected: N=log(30/G1), P=log(20/G1), K=log(40/G1), Rd=log(10/G1)
#'
#' concs2 <- c(N=1, P=0, Rd=99) # G = 0
#' G2 <- calculate_geometric_mean_R(concs2) # Should be 0
#' calculate_V_values_R(concs2, G2) # Expected: list(N=NaN, P=NaN, Rd=NaN)
#'
#' calculate_V_values_R(list(N=1, P=NA, Rd=99), 10) # Conc P is NA
#' # Expected: N=log(1/10), P=NaN, Rd=log(99/10)
#'
#' calculate_V_values_R(list(N=-1, P=10, Rd=89), 10) # Conc N is negative
#' # Expected: N=NaN, P=log(10/10), Rd=log(89/10)
#'
#' calculate_V_values_R(list(), 10) # Empty input
#' # Expected: list()
#'
#' calculate_V_values_R(c(N=1, P=1), NaN) # G is NaN
#' # Expected: c(N=NaN, P=NaN)
calculate_V_values_R <- function(nutrient_concentrations_with_Rd, geometric_mean) {
  # Store original type to return consistently
  is_input_list <- is.list(nutrient_concentrations_with_Rd)
  
  # Convert to simple named vector for processing, names are preserved
  concentrations_vec <- unlist(nutrient_concentrations_with_Rd)
  original_names <- names(concentrations_vec) # Get names after unlisting

  # Handle empty input early
  if (length(concentrations_vec) == 0) {
    if (is_input_list) return(list()) else return(vector(mode(concentrations_vec))) # Use mode for type
  }

  # If G is invalid, all V values are NaN
  if (is.na(geometric_mean) || is.nan(geometric_mean) || geometric_mean == 0) {
    V_values <- rep(NaN, length(concentrations_vec))
    names(V_values) <- original_names
    if (is_input_list) return(as.list(V_values)) else return(V_values)
  }

  # Perform calculation using vectorized operations
  # Note: R's log() is the natural logarithm
  # ratios <- concentrations_vec / geometric_mean # Can produce Inf if G is tiny, or if conc is Inf
  # V_values_vec <- log(ratios)                   # log(Inf) is Inf, log(-Inf) is NaN
                                                # log(0/G) is -Inf if G > 0
                                                # log(<0/G) is NaN if G > 0 (or complex if G < 0)
  
  # Direct calculation to avoid intermediate Inf from ratio if conc is 0.
  # log(0) is -Inf. log(negative) is NaN. These are handled by the explicit check below.
  V_values_vec <- log(concentrations_vec / geometric_mean)

  # Ensure conditions like conc <= 0 or NA conc result in NaN for V-value.
  # CND implies concentrations > 0.
  # is.na() check is important for original NA concentrations.
  # concentrations_vec <= 0 check handles original zero or negative concentrations.
  invalid_conc_indices <- is.na(concentrations_vec) | (concentrations_vec <= 0 & !is.na(concentrations_vec))
  V_values_vec[invalid_conc_indices] <- NaN
  
  names(V_values_vec) <- original_names # Ensure names are reassigned

  if (is_input_list) {
    return(as.list(V_values_vec))
  } else {
    return(V_values_vec)
  }
}

#' Calculate CND Indices (I)
#'
#' Calculates I_nutrient = (V_nutrient - V_star_nutrient) / SD_star_nutrient.
#' @param V_values A named list or vector of V values for the sample.
#' @param V_star_values A named list or vector of mean V values (V*) for the high-yield population.
#' @param SD_star_values A named list or vector of standard deviations of V values (SD*) for the high-yield population.
#' @return A named list or vector of CND indices, matching type of V_values.
#' @export
#' @examples
#' v_vals <- list(N=0.5, P=-0.2, K=0.1, Rd=-0.4)
#' v_star <- list(N=0.6, P=-0.1, K=0.2, Rd=-0.5, X=1) # X is extra
#' sd_star <- list(N=0.1, P=0.05, K=0.05, Rd=0.1, X=0.1)
#' calculate_CND_indices_R(v_vals, v_star, sd_star)
#' 
#' # Example with SD_star = 0 for one nutrient
#' sd_star_zero <- list(N=0.1, P=0, K=0.05, Rd=0.1)
#' calculate_CND_indices_R(v_vals, v_star, sd_star_zero)
#' 
#' # Example with V_value being NA
#' v_vals_na <- list(N=NA, P=-0.2, K=0.1, Rd=-0.4)
#' calculate_CND_indices_R(v_vals_na, v_star, sd_star)
#'
#' # Example with missing reference data for 'K'
#' v_star_missing <- list(N=0.6, P=-0.1, Rd=-0.5)
#' sd_star_missing <- list(N=0.1, P=0.05, Rd=0.1)
#' calculate_CND_indices_R(v_vals, v_star_missing, sd_star_missing)
calculate_CND_indices_R <- function(V_values, V_star_values, SD_star_values) {
  V_values_vec <- unlist(V_values)
  V_star_values_vec <- unlist(V_star_values)
  SD_star_values_vec <- unlist(SD_star_values)
  
  nutrient_names <- names(V_values_vec)
  if (length(nutrient_names) == 0 && length(V_values_vec) == 0) { # Check length of vector too for empty unnamed vector
    if (is.list(V_values)) return(list()) else return(vector(typeof(V_values_vec)))
  }
  
  CND_indices_vec <- numeric(length(nutrient_names))
  names(CND_indices_vec) <- nutrient_names
  
  for (nutrient in nutrient_names) {
    v_nutrient <- V_values_vec[[nutrient]]
    
    # Check if nutrient exists in reference data, and V_nutrient is valid
    # The condition `!nutrient %in% names(V_star_values_vec)` checks if the key exists.
    # Accessing with `V_star_values_vec[[nutrient]]` would return NULL if key doesn't exist, 
    # which `is.na()` or `is.nan()` might not catch as intended for "missing key".
    # So, checking name presence first is more robust for "missing key".
    if (is.na(v_nutrient) || is.nan(v_nutrient) || 
        !(nutrient %in% names(V_star_values_vec)) || 
        !(nutrient %in% names(SD_star_values_vec))) {
      CND_indices_vec[[nutrient]] <- NaN
      next
    }
    
    v_star_nutrient <- V_star_values_vec[[nutrient]]
    sd_star_nutrient <- SD_star_values_vec[[nutrient]]
    
    # Check for NA/NaN in reference values or SD_star being zero
    # This check is now for the *values* after confirming keys exist.
    if (is.na(v_star_nutrient) || is.nan(v_star_nutrient) ||
        is.na(sd_star_nutrient) || is.nan(sd_star_nutrient) || 
        (!is.na(sd_star_nutrient) && sd_star_nutrient == 0)) { # Ensure sd_star_nutrient is not NA before comparing to 0
      CND_indices_vec[[nutrient]] <- NaN
    } else {
      CND_indices_vec[[nutrient]] <- (v_nutrient - v_star_nutrient) / sd_star_nutrient
    }
  }
  
  if (is.list(V_values)) {
    return(as.list(CND_indices_vec))
  } else {
    return(CND_indices_vec)
  }
}

#' Calculate CND Nutritional Unbalance Index (r^2)
#'
#' Calculates r^2 = sum of (CND_index^2), ignoring NA/NaN values.
#' @param CND_indices A named list or vector of CND indices.
#' @return The calculated r^2 value. Returns 0 if all indices are NA/NaN or input is empty.
#' @export
#' @examples
#' calculate_CND_r_squared_R(list(N=1, P=-2, K=0.5))
#' calculate_CND_r_squared_R(list(N=1, P=NA, K=0.5)) # NA should be ignored
#' calculate_CND_r_squared_R(list(N=NaN, P=NA))     # All NA/NaN, should be 0
#' calculate_CND_r_squared_R(list())               # Empty, should be 0
calculate_CND_r_squared_R <- function(CND_indices) {
  indices_vec <- unlist(CND_indices)
  
  if (length(indices_vec) == 0) {
    # This check is technically redundant if all NA/NaN also leads to sum(...na.rm=TRUE) = 0
    # but it's explicit for an empty input.
    return(0)
  }
  
  # Square the indices, then sum, removing NAs from the sum.
  # If all values after squaring are NA (e.g. original were all NA/NaN), sum(..., na.rm=TRUE) will be 0.
  r_squared <- sum(indices_vec^2, na.rm = TRUE)
  
  return(r_squared)
}

#' Calculate Compositional Nutrient Diagnosis (CND)
#'
#' Main function to calculate CND indices and r^2 value.
#' @param sample_nutrient_concentrations A named list or vector of sample nutrient concentrations.
#' @param high_yield_v_means A named list or vector of mean V values (V*) from high-yield population.
#' @param high_yield_v_std_devs A named list or vector of standard deviations of V values (SD*) from high-yield population.
#' @return A list containing `CND_indices` (named list/vector matching type of sample_nutrient_concentrations) and `r_squared` (numeric).
#' @export
#' @examples
#' sample_conc <- list(N=2.5, P=0.25, K=1.5, Ca=0.8, Mg=0.3)
#' # V_star and SD_star would typically come from a reference dataset (e.g. csv, or pre-calculated)
#' # For this example, let's assume some placeholder V_star and SD_star values
#' # These must include 'Rd' as well.
#' V_star_example <- list(N=0.1, P=-1.5, K=-0.5, Ca=-0.8, Mg=-1.2, Rd=3.0)
#' SD_star_example <- list(N=0.2, P=0.1, K=0.2, Ca=0.15, Mg=0.1, Rd=0.3)
#' calculate_cnd_R(sample_conc, V_star_example, SD_star_example)
#' 
#' # Example with empty sample concentrations
#' calculate_cnd_R(list(), V_star_example, SD_star_example)
#' 
#' # Example with vector input
#' sample_conc_vec <- c(N=2.5, P=0.25, K=1.5, Ca=0.8, Mg=0.3)
#' calculate_cnd_R(sample_conc_vec, V_star_example, SD_star_example)
calculate_cnd_R <- function(sample_nutrient_concentrations, 
                            high_yield_v_means, 
                            high_yield_v_std_devs) {

  # Store original type of sample_nutrient_concentrations to match for CND_indices output
  is_input_list <- is.list(sample_nutrient_concentrations)

  # Ensure inputs are of a consistent type (vectors) for easier processing internally
  sample_conc_vec <- unlist(sample_nutrient_concentrations)
  v_means_vec <- unlist(high_yield_v_means)
  v_std_devs_vec <- unlist(high_yield_v_std_devs)

  # 1. Calculate Rd
  # calculate_Rd_R handles list or vector input for sample_conc_vec correctly
  Rd <- calculate_Rd_R(sample_conc_vec) # Pass the vector form
  
  # 2. Create nutrient_concentrations_with_Rd
  # Ensure Rd is added with the name "Rd"
  nutrient_concentrations_with_Rd_vec <- c(sample_conc_vec, Rd = Rd)
  
  # Handle case where sample_nutrient_concentrations might be empty
  # If sample_conc_vec was empty, nutrient_concentrations_with_Rd_vec is just c(Rd=100)
  
  # 3. Calculate Geometric Mean (G)
  # calculate_geometric_mean_R handles vector input
  G <- calculate_geometric_mean_R(nutrient_concentrations_with_Rd_vec)
  
  # 4. Calculate V values
  # calculate_V_values_R handles vector input for concentrations and returns vector
  V_values_vec <- calculate_V_values_R(nutrient_concentrations_with_Rd_vec, G)
  
  # 5. Calculate CND indices
  # calculate_CND_indices_R handles vector inputs and returns a vector
  CND_indices_vec <- calculate_CND_indices_R(V_values_vec, v_means_vec, v_std_devs_vec)
  
  # 6. Calculate CND r_squared
  # calculate_CND_r_squared_R handles vector input
  r_squared <- calculate_CND_r_squared_R(CND_indices_vec)
  
  # Return results as a list
  # Convert CND_indices_vec to a list if the original sample_nutrient_concentrations was a list
  final_cnd_indices <- if (is_input_list) {
    as.list(CND_indices_vec) 
  } else {
    CND_indices_vec # it's already a vector
  }

  return(list(
    CND_indices = final_cnd_indices,
    r_squared = r_squared
  ))
}
