import pandas as pd

# Load the Excel file
input_file = 'marinum_pav.xlsx'  # Replace with your actual file name
output_file = 'output_file.xlsx'  # Replace with your desired output file name

# Read the Excel file
data = pd.read_excel(input_file)

# Process the data
# Skip the first row (headers) and first column (Cluster ID)
processed_data = data.iloc[:, 1:].applymap(lambda x: 0 if x == '-' else 1)

# Combine the first column (Cluster ID) with the processed data
processed_data.insert(0, data.columns[0], data.iloc[:, 0])

# Save the processed data to a new Excel file
processed_data.to_excel(output_file, index=False)

print(f"Data has been processed and saved to {output_file}")