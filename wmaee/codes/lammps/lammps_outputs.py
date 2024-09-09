# Modified from: https://github.com/henriasv/lammps-logfile/tree/master/lammps_logfile
import pandas as pd
from io import StringIO
from typing import Union, List, Dict, Optional

def parse_logfile(logfile: str, return_df: bool = True) -> Union[pd.DataFrame, List[Dict[str, List[Optional[float]]]]]:
    """
    Parse a LAMMPS logfile and extract thermo data.

    Parameters
    ----------
    logfile : str
        Path to the LAMMPS logfile.
    return_df : bool, optional
        If True, returns a pandas DataFrame; if False, returns a list of dictionaries,
        by default True.

    Returns
    -------
    Union[pd.DataFrame, List[Dict[str, List[Optional[float]]]]]
        Parsed thermo data either as a pandas DataFrame or a list of dictionaries.
    """
    # Strings that indicate the start and stop of thermo data in the logfile
    start_thermo_strings = ['Memory usage per processor', 'Per MPI rank memory allocation']
    stop_thermo_strings = ['Loop time', 'ERROR']
    
    # List to store parsed thermo data
    partial_logs = []

    with open(logfile, 'r') as log:
        contents = log.readlines()

    keyword_flag = False
    i = 0

    while i < len(contents):
        line = contents[i]

        if keyword_flag:
            tmpString = ''
            # Check whether any of the thermo stop strings are in the present line
            while not sum([string in line for string in stop_thermo_strings]) >= 1:
                if "\n" in line:
                    tmpString += line
                i += 1
                if i < len(contents):
                    line = contents[i]
                else:
                    break

            # Parse the thermo data using pandas
            partialLogContents = pd.read_table(StringIO(tmpString), sep=r'\s+')

            if return_df:
                partial_logs.append(partialLogContents)
            else:
                # Convert DataFrame to a list of dictionaries if return_df is False
                partial_logs.append(partialLogContents.to_dict(orient='list'))

            keyword_flag = False

        # Check whether the string matches any of the start string identifiers
        if sum([line.startswith(string) for string in start_thermo_strings]) >= 1:
            keyword_flag = True

        i += 1

    # Return either a single DataFrame or a list of dictionaries
    if len(partial_logs) == 1:
        return partial_logs[0]
    return partial_logs