import subprocess
import pandas as pd
import numpy as np

# Define the command to run as a list
command = ["nextflow","run","main.nf","-c","test/test.config","-profile","ci"]

# Run the command and capture the output
result = subprocess.run(command, capture_output=True, text=True)

# Check if the command was successful
if result.returncode == 0:
    # Print the standard output from the command
    print("Command output:\n", result.stdout)
else:
    # Print the error if the command failed
    print("Command failed with error:\n", result.stderr)
