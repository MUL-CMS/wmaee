import sys
import logging
import json
from wmaee.extensions.unsupported.worker.filelock import FileLock
from tinydb import Query, TinyDB
from threading import Thread
from enum import Enum
from glob import glob
from os.path import exists, basename, dirname, join
from os import getcwd
from tqdm import tqdm
from wmaee import working_directory, Poscar, vasp, parse_output

logger = logging.getLogger(__file__)


class LockedTinyDB(TinyDB):

    def __init__(self, filename, *args, **kwargs):
        super(LockedTinyDB, self).__init__(filename, *args, **kwargs)
        self._filelock = FileLock(filename)

    def __enter__(self):
        result = super(LockedTinyDB, self).__enter__()
        self._filelock.acquire()
        return result

    def __exit__(self, exc_type, exc_val, exc_tb):
        result = super(LockedTinyDB, self).__exit__(exc_type, exc_val, exc_tb)
        self._filelock.release()
        return result


def database(fname):
    return LockedTinyDB(fname)


Calculations = Query()


class Status(Enum):
    Queue = 0
    Finished = 1
    Crashed = 2
    Running = 3


def initialize(directory, dbfname, prefix):
    with working_directory(directory):
        if not exists(dbfname):
            # collect all the structures
            folder_prefix = prefix
            calcs = glob(folder_prefix)
            if len(calcs) > 0:
                bulk = []
                for c_id, fd in tqdm(enumerate(calcs), total=len(calcs)):
                    data = {
                        'id': c_id,
                        'folder': fd,
                        'structure': Poscar.from_file(join(fd, 'POSCAR')).structure.as_dict(),
                        'status': Status.Queue.value,
                        'results': {}
                    }
                    bulk.append(data)
                with database(fname=dbfname) as db:
                    db.insert_multiple(bulk)
        else:
            logger.info('Skipped initialization')
            pass


def fetch_next(directory, dbfname):
    with working_directory(directory):
        with database(fname=dbfname) as db:
            remaining = db.search(Calculations.status == Status.Queue.value)
            if len(remaining) < 1:
                logger.info('Worker exiting, since no calculations are available any more')
                result = None
            else:
                first_element = next(iter(remaining))
                first_id = first_element['id']
                # tell the database that it is running now
                db.update(dict(status=Status.Running.value), Calculations.id == first_id)
                first_element['status'] = Status.Running.value
                result = first_element

        return result


def calculate(calculation, directory):
    try:
        with working_directory(directory):
            if exists('call.json'):
                with open('call.json', 'r') as h:
                    kwargs = json.loads(h.read())
            else:
                kwargs = {}
            with working_directory(calculation['folder']):
                from time import sleep
                sleep(0.1)
                # vasp(**kwargs)
                output = parse_output()
                results = {
                    'E': output.final_energy,
                    'forces': list(output.final_forces)
                }
                calculation['results'] = results
                calculation['status'] = Status.Finished.value
                success = True
    except Exception as e:
        logger.exception('An error occured while processing folder "%s"' % calculation['folder'], exc_info=e)
        calculation['status'] = Status.Crashed.value
        success = False
    return success, calculation


def update_calculation(calculation, directory, dbfname):
    cid = calculation['id']
    with working_directory(directory):
        with database(fname=dbfname) as db:
            db.update(dict(status=calculation['status'], results=calculation['results']), Calculations.id == cid)


class CalculationWorker(Thread):

    def __init__(self, directory, dbfname=None, prefix=None, *args, **kwargs):
        super(CalculationWorker, self).__init__(*args, **kwargs)
        self.prefix, self.name = dirname(directory), basename(directory)
        if not self.prefix:
            self.prefix = getcwd()
        if prefix is None:
            self.folder_prefix = '%s-*' % self.name
        else:
            self.folder_prefix = prefix
        if dbfname is None:
            self._dbfname = '%s.tinydb.json' % self.name
        else:
            self._dbfname = dbfname
        self._directory = directory

    def run(self) -> None:
        initialize(self._directory, self._dbfname, self.folder_prefix)
        calculation = fetch_next(self._directory, self._dbfname)
        while calculation is not None:
            success, result = calculate(calculation, self._directory)
            logger.info('Finished calculation "%i" in folder "%s"' % (result['id'], result['folder']))
            update_calculation(result, self._directory, self._dbfname)
            calculation = fetch_next(self._directory, self._dbfname)


if __name__ == '__main__':
    if len(sys.argv) == 2:
        _, directory = sys.argv
        dbfname = None
    elif len(sys.argv) == 3:
        _, directory, dbfname = sys.argv
    else:
        print(sys.argv)
        raise ValueError('Expects either one or two arguments')
    runner = CalculationWorker(directory, dbfname=dbfname)
    runner.setDaemon(True)
    runner.start()
    runner.join()
