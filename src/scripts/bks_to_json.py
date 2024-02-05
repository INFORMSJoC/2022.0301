
import json
import pandas as pd


def create_json_file(data, num_jobs, instance):
    json_file_name = f'wt{num_jobs:0>3}_{int(instance):0>3}.json'
    with open(json_file_name, 'w') as json_file:
        json.dump(data, json_file)
    print(f"JSON file '{json_file_name}' created successfully.")


def read_csv_and_create_json(csv_file):
    # Extracting number of jobs and instance from the CSV file name
    num_machine = []
    num_jobs = []
    for i in csv_file:
        _, _, num_machine_local, num_jobs_local = i.split('_')
        num_jobs_local = num_jobs_local.split('.')[0]
        num_machine.append(num_machine_local)
        num_jobs.append(num_jobs_local)

    pandas_data_frames = []
    for i in csv_file:
        pandas_data_frames.append(pd.read_csv(i))

    # zip the lines of the data frames together
    readers = [df.iterrows() for df in pandas_data_frames]
    # iterate over the rows of the zipped data frames

    for rows in zip(pandas_data_frames[0].iterrows(), pandas_data_frames[1].iterrows()):
        instance = rows[0][1]['Instance']
        bks0 = rows[0][1]['UILS']
        bks1 = rows[1][1]['UILS']
        data = {num_machine[0]: {'bks': int(bks0)}, num_machine[1]: {'bks': int(bks1)}}
        create_json_file(data, num_jobs[0],instance)

    print("JSON files created successfully.")


# Provide the path to your CSV file
csv_file_path = ['./instances/wt100/sol/sol_machine_2_100.txt',
                 './instances/wt100/sol/sol_machine_4_100.txt']
read_csv_and_create_json(csv_file_path)
