#!/bin/bash

# Check if a prefix is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <prefix>"
    exit 1
fi

# Assign the user-provided prefix
prefix=$1
output_file="${prefix}.cfg"

# Define the template file (ensure this is in the same directory or adjust path)
template_file="config.txt"

# Check if the template file exists
if [ ! -f "$template_file" ]; then
    echo "Error: Template file '$template_file' not found!"
    exit 1
fi

# Replace {prefix} with the user-provided prefix and save as a new config file
sed "s/{prefix}/$prefix/g" "$template_file" > "$output_file"

echo "Configuration file '$output_file' has been generated successfully."
