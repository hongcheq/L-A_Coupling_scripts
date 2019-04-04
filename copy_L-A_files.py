"""
Function: based on files names in L-A_Coupling_scripts, copy files here. Then
            Upload to Github.
Date: 20190404
"""

import subprocess


with open('L-A_Coupling_scripts.txt') as f:
    data = f.readlines()

file_names = [x.strip() for x in data]
print(file_names)

pre_path = "/Users/hongcheq/Documents/Data/Data/Greenplanet/Programming/NCL/scripts/"
for i_file in file_names:
    print(i_file)
    bash_command = "cp " + pre_path+ i_file + " ."
    process = subprocess.Popen(bash_command.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
