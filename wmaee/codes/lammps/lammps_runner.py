from wmaee.core.io import working_directory
from wmaee.core.config import Config
import subprocess
from typing import Optional, Dict, Union

def run_lammps(script_name: Optional[str] = 'lammps.in',
               lmp_in: Optional[str] = None,
               command: Optional[str] = None,
               args: Optional[Dict[str, Union[str, int]]] = None,
               directory: Optional[str] = None,
               log: Union[bool, str] = True) -> None:
    """
    Run LAMMPS (Large-scale Atomic/Molecular Massively Parallel Simulator) using the specified command,
    arguments, working directory, and LAMMPS input script.

    Parameters
    ----------
    script_name : str, optional
        The name of the LAMMPS input script. Defaults to 'lammps.in'.
    lmp_in : str, optional
        The entire LAMMPS input script as a string. If provided, it will be written
        into the file <script_name> before executing lammps. If this file exists, it will be
        overwritten.
    command : str, optional
        The LAMMPS command. If not provided, it will be obtained from the
        configuration file.
    args : dict, optional
        Additional arguments to be passed to the LAMMPS command.
    directory : str, optional
        The working directory for running LAMMPS.
    log : bool or str, optional
        If True, capture the output to the screen. If a string is provided,
        capture the output to the specified log file. If False, run silently
        without capturing output.

    Returns
    -------
    None
    """
    
    if command == None:
        # get the lammps command from wmaee.conf.yaml
        cfg = Config()
        command = cfg.get('applications').get('lammps').get('runner').get('script')
        args_template = cfg.get('applications').get('lammps').get('runner').get('args')
        # replace defaults with whatever was provided
        if args == None:
            args = dict()
        for a in args:
            if a in args_template:
                args_template[a] = args[a]
        for a in args_template:
            command = command.replace('{{ ' + str(a) + ' }}', str(args_template[a]))
    command += ' < '+script_name
    
    logfile = False
    if not log:
        # silent
        stdout = subprocess.DEVNULL
        stderr = subprocess.DEVNULL
    elif isinstance(log, str):
        # capture to file
        stdout = subprocess.PIPE
        stderr = subprocess.STDOUT
    else:
        # don't capture -> screen
        stdout = None
        stderr = None
    
    if directory == None:
        directory = '.'
    with working_directory(directory):        
        # write input script
        if lmp_in != None:
            with open(script_name, 'w') as infile:
                infile.write(lmp_in)
        
        # run lammps
        run = subprocess.run(command, shell=True, stdout = stdout,
                             stderr = stderr, text = True)
        if isinstance(log, str):
            with open(log, 'w') as output:
                output.write(run.stdout)