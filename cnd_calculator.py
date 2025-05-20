import math

def calculate_Rd(nutrient_concentrations: dict) -> float:
    """
    Calculates the residual dry matter (Rd) from nutrient concentrations.

    Args:
        nutrient_concentrations: A dictionary where keys are nutrient names
                                 (e.g., 'N', 'P', 'K') and values are their
                                 concentrations (e.g., float).

    Returns:
        The calculated Rd value (float).
    """
    total_nutrient_concentration = sum(nutrient_concentrations.values())
    rd = 100.0 - total_nutrient_concentration
    return rd

def calculate_geometric_mean(nutrient_concentrations_with_Rd: dict) -> float:
    """
    Calculates the geometric mean (G) of nutrient concentrations including Rd.

    Args:
        nutrient_concentrations_with_Rd: A dictionary where keys are nutrient names
                                         (e.g., 'N', 'P', 'K', 'Rd') and values are
                                         their concentrations (float).

    Returns:
        The calculated geometric mean (G) value (float).
    """
    product_of_concentrations = 1.0
    for concentration in nutrient_concentrations_with_Rd.values():
        product_of_concentrations *= concentration
    
    num_components = len(nutrient_concentrations_with_Rd)
    if num_components == 0:
        return 0.0  # Or raise an error, depending on desired behavior for empty input

    geometric_mean = product_of_concentrations ** (1 / num_components)
    return geometric_mean

def calculate_V_values(nutrient_concentrations_with_Rd: dict, geometric_mean: float) -> dict:
    """
    Calculates the V values for each nutrient.

    V_nutrient = ln(concentration_nutrient / geometric_mean)

    Args:
        nutrient_concentrations_with_Rd: A dictionary where keys are nutrient names
                                         (e.g., 'N', 'P', 'K', 'Rd') and values are
                                         their concentrations (float).
        geometric_mean: The geometric mean (G) of the nutrient concentrations (float).

    Returns:
        A dictionary where keys are nutrient names and values are the
        calculated V values (float).
    """
    v_values = {}
    if geometric_mean == 0: # Avoid division by zero if G is 0
        # Handle this case as appropriate for the application
        # For now, return empty or raise an error, or return V values as 0 or NaN
        for nutrient_name in nutrient_concentrations_with_Rd.keys():
            v_values[nutrient_name] = math.nan # Or some other indicator of an issue
        return v_values

    for nutrient_name, concentration in nutrient_concentrations_with_Rd.items():
        if concentration <= 0: # Avoid log of non-positive number
            # Handle this case as appropriate
            v_values[nutrient_name] = math.nan # Or some other indicator
            continue
        v_values[nutrient_name] = math.log(concentration / geometric_mean)
    return v_values

def calculate_CND_indices(V_values: dict, V_star_values: dict, SD_star_values: dict) -> dict:
    """
    Calculates the CND indices for each nutrient.

    CND_index_nutrient = (V_nutrient - V_star_nutrient) / SD_star_nutrient

    Args:
        V_values: A dictionary of nutrient names and their calculated V values (float).
        V_star_values: A dictionary of nutrient names and their mean V values (V*)
                       for the high-yield reference population (float).
        SD_star_values: A dictionary of nutrient names and their standard deviation
                        of V values (SD*) for the high-yield reference
                        population (float).

    Returns:
        A dictionary where keys are nutrient names and values are the
        calculated CND indices (float).
    """
    cnd_indices = {}
    for nutrient_name, v_nutrient in V_values.items():
        if nutrient_name not in V_star_values or nutrient_name not in SD_star_values:
            cnd_indices[nutrient_name] = math.nan
            # Optionally, log a warning here:
            # print(f"Warning: Missing V* or SD* for nutrient {nutrient_name}")
            continue

        v_star_nutrient = V_star_values[nutrient_name]
        sd_star_nutrient = SD_star_values[nutrient_name]

        if sd_star_nutrient == 0:
            cnd_indices[nutrient_name] = math.nan # Avoid division by zero
            # Optionally, log a warning here:
            # print(f"Warning: SD* is zero for nutrient {nutrient_name}")
            continue
        
        # Ensure v_nutrient is not nan before calculation
        if math.isnan(v_nutrient):
            cnd_indices[nutrient_name] = math.nan
            continue

        cnd_indices[nutrient_name] = (v_nutrient - v_star_nutrient) / sd_star_nutrient
    return cnd_indices

def calculate_CND_r_squared(CND_indices: dict) -> float:
    """
    Calculates the nutritional unbalance index (r^2) from CND indices.

    r^2 = IN^2 + IP^2 + ... + IR^2

    Args:
        CND_indices: A dictionary where keys are nutrient names (e.g., 'N', 'P', 'K', 'Rd')
                     and values are their calculated CND indices (float).

    Returns:
        The calculated r^2 value (float). If a CND index is math.nan,
        it's excluded from the sum.
    """
    r_squared = 0.0
    for cnd_index in CND_indices.values():
        if not math.isnan(cnd_index):
            r_squared += cnd_index ** 2
    return r_squared

def calculate_cnd(sample_nutrient_concentrations: dict, high_yield_v_means: dict, high_yield_v_std_devs: dict):
    """
    Calculates the Compositional Nutrient Diagnosis (CND) indices for a given sample.

    Args:
        sample_nutrient_concentrations: A dictionary where keys are nutrient names
                                       (e.g., 'N', 'P', 'K') and values are their
                                       concentrations (e.g., float).
        high_yield_v_means: A dictionary where keys are nutrient names (including 'Rd')
                            and values are the mean V values (V*) for the
                            high-yield reference population.
        high_yield_v_std_devs: A dictionary where keys are nutrient names (including 'Rd')
                               and values are the standard deviation of V values (SD*)
                               for the high-yield reference population.

    Returns:
        A tuple containing:
            - CND_indices (dict): A dictionary where keys are nutrient names
                                  (e.g., 'N', 'P', 'K', 'Rd') and values are their
                                  calculated CND indices (float).
            - r_squared (float): The calculated nutritional unbalance index (r^2).
    """
    # 1. Calculate Rd
    rd_value = calculate_Rd(sample_nutrient_concentrations)

    # 2. Create nutrient_concentrations_with_Rd
    # Make a copy to avoid modifying the input dictionary
    nutrient_concentrations_with_Rd = sample_nutrient_concentrations.copy()
    nutrient_concentrations_with_Rd['Rd'] = rd_value

    # Handle cases where a nutrient concentration might be zero or negative,
    # which could lead to issues in log operations later.
    # calculate_geometric_mean will produce 0 if any concentration is 0.
    # calculate_V_values will produce nan if concentration is <= 0 or G is 0.

    # 3. Calculate geometric mean (G)
    geometric_mean_g = calculate_geometric_mean(nutrient_concentrations_with_Rd)

    # 4. Calculate V values
    v_values = calculate_V_values(nutrient_concentrations_with_Rd, geometric_mean_g)

    # 5. Calculate CND indices
    # These means and std_devs must contain 'Rd' as a key if it's in V_values.
    cnd_indices = calculate_CND_indices(v_values, high_yield_v_means, high_yield_v_std_devs)

    # 6. Calculate CND r_squared
    r_squared = calculate_CND_r_squared(cnd_indices)

    return cnd_indices, r_squared
