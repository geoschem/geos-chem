# Verify_KPP_Standalone.py
# This script is used to verify the KPP Standalone ability to replicate chemistry of
# grid cells from a 3D model using the GEOS-Chem chemical mechanism.

import os
import subprocess
import pandas as pd
import matplotlib.pyplot as plt

# 1) Read all the files in a samples local directory
sample_dir = 'samples'  # replace with your directory path
standalone_dir = '/Users/psturm/Desktop/Twilight_KPP/KPP-Standalone'
files = os.listdir(sample_dir)

# Create a pandas dataframe to store the warning files
df = None


# 2) Run the KPP Standalone on all of them
for file in files:
    process = subprocess.Popen([os.path.join(standalone_dir, 'kpp_standalone.exe'), 
                                os.path.join(sample_dir, file)], 
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()

    # 3) If the KPP Standalone writes any output with "Warning" in it, print a message
    if 'Warning' in stdout.decode() or 'Warning' in stderr.decode():
        print(f'Warning found in the output of file: {file}')

        # Show the output
        print(stdout.decode())
        print(stderr.decode())
        
        # Extract the filename and the other fields from the stdout
        lines = stdout.decode().split('\n')
        data = {'filename': lines[0].split(': ')[1]}
        for line in lines[1:]:
            if ': ' in line:
                key, value = line.split(': ')
                key = key.strip()
                data[key] = value
        
        # Write the data to the pandas DataFrame
        if df is None:
            df = pd.DataFrame(data, index=[0])
        else:
            df = pd.concat([df, pd.DataFrame(data, index=[0])], ignore_index=True)

# Print some diagnostics
# print the length of the data frame and the number of files checked
print(f'Number of files checked: {len(files)}')
if df is None:
    print('No files with warnings')
else:
    print(f'Number of files with warnings: {len(df)}')
    df['Number of internal timesteps ( standalone)'] = pd.to_numeric(df['Number of internal timesteps ( standalone)'])
    df['Number of internal timesteps (from 3D run)'] = pd.to_numeric(df['Number of internal timesteps (from 3D run)'])
    df['Difference'] = abs(df['Number of internal timesteps ( standalone)'] - df['Number of internal timesteps (from 3D run)'])
    # find number of non-zero differences
    print(f'Number of files differing in internal timesteps: {len(df[df["Difference"] > 0])}')
    print(f'Maximum difference in the number of internal timesteps: {df["Difference"].max()}')
    print(f'File with the maximum difference: {df.loc[df["Difference"].idxmax()]["filename"]}')

    # Save the data frame to a csv file
    # df.to_csv('warnings.csv', index=False)

    # plt.hist((df['Number of internal timesteps ( standalone)'] - df['Number of internal timesteps (from 3D run)']))
    # plt.xlabel('Difference in internal timesteps, Standalone to 3D')

# write the list of validated files to a new text file
with open('filelist_validated.txt', 'w') as f:
    for file in files:
        if "samples/"+file not in df['filename'].values:
            f.write(file + '\n')