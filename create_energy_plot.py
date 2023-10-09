import re
import matplotlib.pyplot as plt


# Define the range of file numbers you want to import
num_files = 117  # Ending file number (adjust as needed)

# List to store extracted decimal values
decimal_values = []
# Loop through the file numbers and import data from each file
for i in range(1, num_files + 1):
    # Construct the file name based on the file number
    file_number = str(i).zfill(5)
    file_name = f'houdini/dam_break_small_APIC/energy{file_number}.obj'
    
    try:
        with open(file_name, 'r') as file:
            # Read the content of the file
            file_content = file.read()
            
            # Use regular expression to find decimal values following 'e'
            matches = file_content[2]
            file_content = file_content.replace('e ', '')
            file_content = file_content.replace(' "\\"n', '')
            
            # Add the extracted decimal values to the list
            decimal_values.append(float(file_content))
            
    except FileNotFoundError:
        print(f"File {file_name} not found.")

# Create the plot with customization
plt.plot(decimal_values, linestyle='-', color='b')
plt.xlabel('Time')
plt.xticks([])  # Pass an empty list to hide the ticks
plt.ylabel('Energy')
plt.title('Energy over time')
plt.grid(True)
plt.legend()

# Display the plot
plt.show()