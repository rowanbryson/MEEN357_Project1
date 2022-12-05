import json

def submit_design_candidate(design_info, filepath):
    """
    opens the json file at filepath and submits the design candidate

    if it is one of the top 10 designs, it adds to the list
    if there were already 10 designs, it also removes the worst one
    """
    design = design_info['design']
    design_nums = design_info['design_nums']
    time = design_info['time']
    percent_success = design_info['percent_success']

    with open(filepath, 'r') as f:
        # check if the design qualifies as top 10
        data = json.load(f)
        # sort the data by time where the highest time comes first
        data = sorted(data, key=lambda x: x['time'] * -1)
        if len(data) < 10:
            data.append({
                'design_nums': list(design_nums),
                'design': design,
                'time': time,
                'percent_success': percent_success
            })
            # sort the data by time where the lowest time comes first
            data = sorted(data, key=lambda x: x['time'])
        else:
            # check if the design is better than the worst design
            if time < data[0]['time']:
                # replace the worst design with the new design
                data[0] = {
                    'design_nums': list(design_nums),
                    'design': design,
                    'time': time,
                    'percent_success': percent_success
                }
        # sort the data by time where the lowest time comes first
        data = sorted(data, key=lambda x: x['time'])
        # write the data to the json file
        with open(filepath, 'w') as f:
            json.dump(data, f, indent=4)


if __name__ == '__main__':
    input_path = 'opt_logs_12-5-2022.json'
    output_path = 'best_designs_multi_terrain.json'

    with open(input_path, 'r') as f:
        # the json file starts with [ and ends with ]
        # handle this and load the data into a list of dictionaries
        data = json.load(f)['designs']

    for design_info in data:
        submit_design_candidate(design_info, output_path)


