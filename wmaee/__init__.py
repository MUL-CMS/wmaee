import sys

from ase import Atoms
from wmaee.core import GpawCalculation, VaspCalculation, MDCalculation, UnmetRequirement, show, show_traj, \
    available_models, available_lammps_models, AbinitCalculation, md as md


def monkey_patch_abinit_calculator():
    import io
    import os
    import time
    import subprocess
    from ase.calculators.calculator import Calculator, CalculatorSetupError, CalculationFailed
    from ase.calculators.calculator import all_changes

    def _calculate(self, atoms=None, properties=['energy'],
                  system_changes=all_changes, interval: float = 0.1):
        Calculator.calculate(self, atoms, properties, system_changes)
        self.write_input(self.atoms, properties, system_changes)
        if self.command is None:
            raise CalculatorSetupError(
                'Please set ${} environment variable '
                .format('ASE_' + self.name.upper() + '_COMMAND') +
                'or supply the command keyword')
        command = self.command
        if 'PREFIX' in command:
            command = command.replace('PREFIX', self.prefix)

        if self.output == '-':
            stdout, stderr = (sys.stdout, sys.stderr)
            popen_kwargs = dict(stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            do_forward = True
            have_opened_file = False
        elif self.output is None:
            stdout, stderr = (os.devnull, os.devnull)
            popen_kwargs = dict(stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            do_forward = False
            have_opened_file = False
        elif isinstance(self.output, str):
            stdout, stderr = (open(f"{self.output}.out", 'w'), open(f"{self.output}.err", 'w'))
            popen_kwargs = dict(stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            do_forward = True
            have_opened_file = True
        else:
            stdout, stderr = (self.output, self.output)
            popen_kwargs = dict(stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            do_forward = True
            have_opened_file = False

        try:
            proc = subprocess.Popen(command, shell=True, cwd=self.directory, **popen_kwargs)
        except OSError as err:
            # Actually this may never happen with shell=True, since
            # probably the shell launches successfully.  But we soon want
            # to allow calling the subprocess directly, and then this
            # distinction (failed to launch vs failed to run) is useful.
            msg = 'Failed to execute "{}"'.format(command)
            raise EnvironmentError(msg) from err

        if do_forward:
            stdout_src = io.TextIOWrapper(proc.stdout)
            stderr_src = io.TextIOWrapper(proc.stderr)
            while proc.poll() is None:
                stdout.write(stdout_src.read())
                stderr.write(stderr_src.read())
                time.sleep(interval)
            if have_opened_file:
                stdout.close()
                stderr.close()
        errorcode = proc.wait()

        if errorcode:
            path = os.path.abspath(self.directory)
            msg = ('Calculator "{}" failed with command "{}" failed in '
                   '{} with error code {}'.format(self.name, command,
                                                  path, errorcode))
            raise CalculationFailed(msg)

        self.read_results()

    from ase.calculators.abinit import Abinit
    setattr(Abinit, "calculate", _calculate)


monkey_patch_abinit_calculator()

__all__ = ["Atoms", "GpawCalculation", "VaspCalculation", "MDCalculation", "UnmetRequirement", "show", "show_traj",
           "available_models", "available_lammps_models", "md", "AbinitCalculation"]
