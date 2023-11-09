#!/bin/bash

# Define arrays of values for each variable
ny_values=(64)
selected_case_values=(0 1)
bf_values=(5.0)
fact_values=(0.2)

# Loop through each combination of values and run the program
for ny in "${ny_values[@]}"; do
    for selected_case in "${selected_case_values[@]}"; do
        for bf in "${bf_values[@]}"; do
            for fact in "${fact_values[@]}"; do
                # Build the command
                command="./2d_analysis ny=${ny} selected_case=${selected_case} bf=${bf} fact=${fact}"
                echo "Running command: $command"
                # Run the command
                $command
                # Check if the command was successful
                if [ $? -ne 0 ]; then
                    echo "Command failed: $command"
                    # Uncomment the next line if you want the script to exit on failure
                    # exit 1
                fi
            done
        done
    done
done

echo "All combinations have been run."
