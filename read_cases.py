import csv

# Define the file path
file_path = "data/cases_perday.csv"

# Initialize empty lists to store local and imported cases
local_cases = []
imported_cases = []

def read_cases_data(file_path):
  # Read the CSV file
  with open(file_path, newline='') as csvfile:
    reader = csv.DictReader(csvfile)

    # Iterate through each row in the CSV file
    for row in reader:
      # Extract local and imported cases (assuming they are integers)
      try:
        local_cases.append(int(row["local"]))
        imported_cases.append(int(row["imported"]))
      except ValueError:
        # Handle potential errors if the values are not numerical
        print(f"Error: Could not convert values in row {row} to integers")

# Print statements to verify (optional)

read_cases_data(file_path)
# print("Local cases:", local_cases)
# print("Imported cases:", imported_cases)
