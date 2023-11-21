# Modified from: https://github.com/henriasv/lammps-logfile/tree/master/lammps_logfile
import pandas as pd
import numpy as np
from io import BytesIO, StringIO

def parse_logfile(logfile, return_df=True):
    start_thermo_strings = ['Memory usage per processor', 'Per MPI rank memory allocation']
    stop_thermo_strings = ['Loop time', 'ERROR']
    partial_logs = []
    
    with open(logfile, 'r') as log:
        contents = log.readlines()
    keyword_flag = False
    i = 0
    while i < len(contents):
        line = contents[i]

        if keyword_flag:
            # keywords = line.split()
            tmpString = ''
            # Check wheter any of the thermo stop strigs are in the present line
            while not sum([string in line for string in stop_thermo_strings]) >= 1:
                if "\n" in line:
                    tmpString += line
                i += 1
                if i < len(contents):
                    line = contents[i]
                else:
                    break
            partialLogContents = pd.read_table(StringIO(tmpString), sep=r'\s+')

            # if (self.keywords != keywords):
            #     # If the log keyword changes, i.e. the thermo data to be outputted chages,
            #     # we flush all prevous log data. This is a limitation of this implementation. 
            #     self.flush_dict_and_set_new_keyword(keywords)

            partial_dict = {}
            # for name in keywords:
            #     # data_dict[name] = np.append(data_dict[name],partialLogContents[name])
            #     partial_dict[name] = list(partialLogContents[name].values)#np.append(np.asarray([]), partialLogContents[name])
            if return_df:
                partial_logs.append(partialLogContents)
            else:
                partial_logs.append(partialLogContents.to_dict(orient = 'list'))
            keyword_flag = False

        # Check whether the string matches any of the start string identifiers
        if sum([line.startswith(string) for string in start_thermo_strings]) >= 1:
            keyword_flag = True
        i += 1
    return partial_logs
