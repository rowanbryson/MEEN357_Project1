# I need to get the data from the logs and put in into a json file

import json
import os
import sys
from rich import print

# Get the path to the logs
input = 'opt_logs_12-5-2022.txt'
output = 'opt_logs_12-5-2022.json'

# the file is a loguru log file
# I need to get the data from the message field and put it into a list of dictionaries and then into a json file

# I need to get the data from the logs and put in into a json file

data = []
with open(input, 'r') as f:
    for line in f:
        timestamp, loglevel, message = [piece.strip() for piece in line.split('|')]
        if loglevel == 'INFO' and not message == 'Starting optimization session...':
            data.append(message)

data = [line.replace("'", '"') for line in data]
data = [line.replace('array([', '[') for line in data]
data = [line.replace('])', ']') for line in data]

# load the data into python dictionaries
data = [json.loads(line) for line in data]

# write the data to a json file
with open(output, 'w') as f:
    json.dump(data, f, indent=4)
