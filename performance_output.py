"""
This program takes the performance data from the phyCD algorithm and
consolidates it into a readable form
"""

import os
import pprint
import pandas as pd

def main():

    # finds directories that have performanceOutput in the name
    directories = []
    for name in os.listdir():
        if "performanceOutput" in name:
            directories.append(name)

    # makes output directory if it does not exist
    if not os.path.isdir('Combined Performance Files'):
        os.mkdir('Combined Performance Files')

    # runs through all directories that have performanceOutput in the name
    for dir in directories:

        # runs through each file in the output directory
        rowList = []
        for filename in os.listdir(dir):

            # gets the testing parameters
            parameters = filename[:-4].split(",")
            num_p = int(parameters[0])
            num_t = int(parameters[1])
            if len(parameters) == 3:
                str_len = int(parameters[2])
            else:
                str_len = None

            # gets performance value
            p_file = open(f"{dir}/{filename}", "r")
            time = float(p_file.read())

            # stores parameters and performance value in row list to be stored in dataframe
            rowList.append([num_p, num_t, str_len, time])

        # takes data from row list and makes dataframe
        performanceDf = pd.DataFrame(rowList, columns=["num_p", "num_t", "str_len", "time"])
        performanceDf.sort_values(by=["num_p", "num_t", "str_len"], inplace=True)

        # outputs dataframe to csv
        performanceDf.to_csv(f"Combined Performance Files/{dir}.csv", index=False)


if __name__ == '__main__':
    main()
