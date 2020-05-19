import sys
import logging
import json
from threading import Thread
from enum import Enum
from glob import glob
from os.path import exists, basename, dirname, join
from os import getcwd
from tqdm import tqdm
from wmaee import working_directory, Poscar, vasp, parse_output
from wmaee.core.common import fullname
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, String, Integer, DateTime, create_engine
from sqlalchemy.orm import sessionmaker
from contextlib import contextmanager
from datetime import datetime
from collections import namedtuple
from time import sleep
logger = logging.getLogger(__file__)
# One day of timeout limit
LOCK_TIMEOUT = (24*2600)

__DATABASE_SESSION = None
__DATABASE_ENGINE = None

Base = declarative_base()



def get_connection_string(fname: str) -> str:
    return 'sqlite:///%s' % fname


class Calculation(Base):

    __tablename__ = 'calculations'
    calculation_id = Column(Integer, primary_key=True)
    folder = Column(String(512), nullable=False)
    status = Column(Integer, nullable=False)
    started = Column(DateTime, nullable=True)
    finished = Column(DateTime, nullable=True)

    def __repr__(self):
        return '%s(id=%i, folder=%s, status=%i)' % (fullname(self), self.calculation_id, self.folder, self.status)


def make_engine(fname=None):
    global __DATABASE_ENGINE
    if __DATABASE_ENGINE is None:
        __DATABASE_ENGINE = create_engine(get_connection_string(fname))
        # bind the metadata to the engine
        Base.metadata.bind = __DATABASE_ENGINE
    return __DATABASE_ENGINE


def get_session(engine=None):
    global __DATABASE_SESSION
    if __DATABASE_SESSION is None:
        __DATABASE_SESSION = sessionmaker(bind=engine)
    return __DATABASE_SESSION


@contextmanager
def database_session():
    """Provide a transactional scope around a series of operations."""
    session = get_session()()
    try:
        yield session
        session.commit()
    except Exception as exc:
        session.rollback()
        logger.exception('Transaction failed!', exc_info=exc)
        raise
    finally:
        session.close()


class Status(Enum):
    Queue = 0
    Finished = 1
    Crashed = 2
    Running = 3


def show_calculations(fname):
    engine = create_engine(get_connection_string(fname))
    import pandas as pd
    with engine.connect() as connection:
        return pd.read_sql_table(Calculation.__tablename__, connection)


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
            self._dbfname = '%s.db' % self.name
        else:
            self._dbfname = dbfname
        self._directory = directory
        self._engine = make_engine(fname=join(self._directory, self._dbfname))
        self._session_meta = get_session(self._engine)

    def initialize_database(self):
        Base.metadata.create_all(self._engine)
        # After that we can bind it

    def initialize(self):
        with working_directory(self._directory):
            if not exists(self._dbfname):
                # collect all the structures
                self.initialize_database()
                # Now the database does not yet exist we have to create the datatable
                folder_prefix = self.folder_prefix
                calcs = glob(folder_prefix)
                if len(calcs) > 0:
                    bulk = []
                    for c_id, fd in tqdm(enumerate(calcs), total=len(calcs)):
                        xml_path = join(fd, 'vasprun.xml')
                        if exists(xml_path):
                            try:
                                with working_directory(fd):
                                    parse_output()
                            except Exception:
                                status = Status.Crashed.value
                            else:
                                status = Status.Finished.value
                        else:
                            status = Status.Queue.value
                        finished = datetime.now() if status == Status.Finished.value else None
                        bulk.append(Calculation(calculation_id=c_id, folder=fd, status=status, started=None, finished=finished))
                    with database_session() as session:
                        for calculation in bulk:
                            session.add(calculation)
                    initialized = True
            else:
                logger.info('Skipped initialization')
                initialized = False
        return initialized

    def fetch_next(self):
        with database_session() as session:
            remaining = session.query(Calculation).filter(Calculation.status == Status.Queue.value).count()
            if remaining < 1:
                logger.info('Worker exiting, since no calculations are available any more')
                result = None
            else:
                first_calculation = session.query(Calculation).filter(Calculation.status == Status.Queue.value).first()
                first_calculation.started = datetime.now()
                first_calculation.status = Status.Running.value
                result = (first_calculation.calculation_id, first_calculation.folder, first_calculation.status)
        return result

    def calculate(self, calculation):
        calculation_id, folder, _ = calculation
        print('Running calculation directory "%s"' % join(self._directory, folder))
        try:
            with working_directory(directory):
                if exists('call.json'):
                    with open('call.json', 'r') as h:
                        kwargs = json.loads(h.read())
                else:
                    kwargs = {}
                with working_directory(folder):
                    vasp(**kwargs)
                    success = True
                status = Status.Finished.value
        except Exception as e:
            logger.exception('An error occurred while processing folder "%s"' % folder, exc_info=e)
            status = Status.Crashed.value
            success = False
        # logger
        with database_session() as session:
            calc = session.query(Calculation).filter(Calculation.calculation_id == calculation_id).one()
            calc.status = status
            calc.finished = datetime.now()
        return success




class SimpleCalculationWorker(Thread):

    def __init__(self, directory, dbfname=None, prefix=None, *args, **kwargs):
        super(SimpleCalculationWorker, self).__init__(*args, **kwargs)
        self.prefix, self.name = dirname(directory), basename(directory)
        if not self.prefix:
            self.prefix = getcwd()
        if prefix is None:
            self.folder_prefix = '%s-*' % self.name
        else:
            self.folder_prefix = prefix
        if dbfname is None:
            self._dbfname = '%s.db' % self.name
        else:
            self._dbfname = dbfname
        self._directory = directory
        self._done = []
        self._remaining = []
        self._crashed = []
        self._lock_file_name = None
        self._calc_meta = namedtuple('CalculationTuple', ['calculation_id', 'folder', 'status'])

    def _get_status(self, pth):
        xml_path = join(pth, 'vasprun.xml')
        if exists(xml_path):
            try:
                with working_directory(pth):
                    parse_output()
            except Exception:
                status = Status.Crashed.value
            else:
                status = Status.Finished.value
        else:
            status = Status.Queue.value
        return status

    def initialize(self):
        with working_directory(self._directory):
            # Now the database does not yet exist we have to create the datatable
            folder_prefix = self.folder_prefix
            calcs = glob(folder_prefix)
            if len(calcs) > 0:
                for c_id, fd in tqdm(enumerate(calcs), total=len(calcs)):
                    status = self._get_status(fd)
                    calc = self._calc_meta(calculation_id=c_id, folder=fd, status=status)
                    if status == Status.Queue.value:
                        self._remaining.append(calc)
                    elif status == Status.Finished.value:
                        self._done.append(calc)
                    elif status == Status.Crashed.value:
                        self._crashed.append(calc)

    def _update_info_remaining(self):
        with working_directory(self._directory):
            for r in self._remaining:
                r.status = self._get_status(r.folder)
            finished = [s for s in self._remaining if  s.status == Status.Finished.value]
            crashed = [s for s in self._remaining if  s.status == Status.Crashed.value]
            for fin in finished:
                self._remaining.remove(fin)
                self._done.append(fin)
            for cr in crashed:
                self._remaining.remove(cr)
                self._crashed.append(cr)

    def fetch_next(self):
        if len(self._remaining) < 1:
            logger.info('Worker exiting, since no calculations are available any more')
            result = None
        else:
            candidate = self._remaining.pop(0)
            def is_valid(direct, cand):
                with working_directory(join(direct, cand.folder)):
                    locks = glob('*.lock')
                    if len(locks) != 0:
                        return False
                    # Lets wait if probably another worker want's to have this directory

    return result


    def run(self) -> None:
        self.initialize()
        calculation = self.fetch_next()
        while calculation is not None:
            success = self.calculate(calculation)
            logger.info('Finished calculation "%i" in folder "%s" = %s' % (calculation[0], calculation[1], success))
            calculation = self.fetch_next()

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
