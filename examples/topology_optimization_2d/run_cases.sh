#!/bin/bash

# Define arrays of values for each variable
ny_values=(64)
selected_case_values=(2)
bf_values=(5.0)
fact_values=(0.2)
vol_values=(0.2 0.4 0.8) 
trac_values=(-100.0 -500.0 -1000.0) 

# Loop through each combination of values and run the program
for ny in "${ny_values[@]}"; do
    for selected_case in "${selected_case_values[@]}"; do
        for bf in "${bf_values[@]}"; do
            for fact in "${fact_values[@]}"; do
                for vol in "${vol_values[@]}"; do
                    for trac in "${trac_values[@]}"; do
                        # Build the command
                        command="./2d_analysis ny=${ny} selected_case=${selected_case} bf=${bf} fact=${fact} vol=${vol} trac=${trac}"
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
    done
done

echo "All combinations have been run."
