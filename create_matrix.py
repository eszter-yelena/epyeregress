import pandas as pd
import numpy as np

def normalize_data(accumulator_matrix):
    normalised_matrix = accumulator_matrix.copy()
    for col in normalised_matrix.columns:
        col_variance = accumulator_matrix[col].var()
        if col_variance > 0:
            normalised_matrix[col] = (accumulator_matrix[col] - accumulator_matrix[col].mean()) / np.sqrt(col_variance)
        else:
            normalised_matrix[col] = 0
    return normalised_matrix

# Read the mobility data from the CSV file
mobility_data = pd.read_csv("data/mobility_report.csv")
vaccination_data = pd.read_csv("data/vaccination_scores.csv")
delta_data = pd.read_csv("data/delta_count_NZ.csv")
oxford_data = pd.read_csv("data/oxford_scores.csv")

# Select the last six columns
accumulator_matrix = mobility_data.iloc[:, -6:].copy()

# accumulator_matrix = pd.DataFrame()

# Add new columns with zeros
accumulator_matrix['vaccination'] = vaccination_data['per 100 people']
accumulator_matrix['delta'] = delta_data['percentage']

# Fill the matrix with the oxford scores
for col in oxford_data.columns:
    if col.startswith('V'):
        accumulator_matrix[col] = oxford_data[col]

# Create a new matrix based on the current matrix
accumulator_normalized = normalize_data(accumulator_matrix)

# Calculate standard deviations
factor = np.sqrt(accumulator_normalized.var(axis=0))

accumulator_normalized_np = accumulator_normalized.to_numpy()

# Calculate variances for each column
variances = np.var(accumulator_normalized_np, axis=0, ddof=1)
factor = np.sqrt(variances)

# Write the results to a CSV file
factor_df = pd.DataFrame(factor)
factor_df.to_csv(f'data/factor.csv', index=False)

# Save standard deviations to CSV
np.savetxt("data/factor.csv", factor, delimiter=",", fmt="%f", header="")

# Select columns with non-zero variance
accumulator_normalized = accumulator_normalized.loc[:, factor != 0]  # Select columns where corresponding factor value is not zero

# Print the result
# print(accumulator_normalized)
