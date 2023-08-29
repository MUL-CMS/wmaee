#!/usr/bin/python
import re
import os
import pwd
import shlex
import fnmatch
import datetime
import tempfile
import operator
import functools
from typing import *

BACKUP_DIR="/media/backup"
KEEP = 7
EXCLUDE_USERS = frozenset({"nobody"})
EXCLUDE_PATTERNS = {"*.git*", "*.tar.gz", "*.local/*", "*/.bash*", "*.kim-api/*", "*.cache/*", "*.npm/*", "*.jupyter/*", "*.ipython/*", "*.db", "*.profile"}

users = {pw.pw_name: pw.pw_dir for pw in pwd.getpwall() if pw.pw_uid >= 1000 and pw.pw_name not in EXCLUDE_USERS}
snd = operator.itemgetter(1)
fst = operator.itemgetter(0)
date = datetime.datetime.now()


def build_archive_name_regex() -> re.Pattern:
    data_part = r"(?P<year>\d+)-(?P<month>\d+)-(?P<day>\d+)-(?P<hour>\d+)-(?P<minute>\d+)\.tar\.gz"
    return re.compile(r"(?P<user>" + '|'.join(users) + ")-" + data_part)


def get_archive_name(username: str, suffix: str = "tar.gz") -> str:
    return f"{username}-{date.year}-{date.month}-{date.day}-{date.hour}-{date.minute}.{suffix}"


def files_in_directory(user: str) ->  Iterator[str]:
    user_dir = users.get(user)
    for root, dirs, files in os.walk(users.get(user)):
        yield from (os.path.relpath(os.path.join(root, f), start=os.path.dirname(user_dir)) for f in files)


def discard_excluded(fs: Iterator[str]) -> bool:
    for f in fs:
        if any(fnmatch.fnmatch(f, pat) for pat in EXCLUDE_PATTERNS):
            continue
        else:
            yield f


def run_command(*command: str):
    result = os.system(shlex.join(command))
    if result != 0:
        raise RuntimeError(f"Return code {result}")


def pack_user_archive(user: str, filename: str) -> str:
    user_dir = users.get(user)
    curr_dir = os.getcwd()
    user_parent_dir = os.path.dirname(user_dir)
    os.chdir(user_parent_dir)
    archive_name = os.path.join(BACKUP_DIR, get_archive_name(user))
    try:
        run_command("tar", "-czvf", archive_name, "-T", filename)
    except RuntimeError:
        raise
    finally:
        os.chdir(curr_dir)

    return archive_name


def parse_archive_name(user: str, year: str, month: str, day: str, hour: str, minute: str) -> Tuple[str, datetime.datetime]:
    return user, datetime.datetime(year=int(year), month=int(month), day=int(day), hour=int(hour), minute=int(minute))


@functools.lru_cache(maxsize=1)
def get_archives(root: str) -> Dict[str, Dict[str, datetime.datetime]]:
    pattern = build_archive_name_regex()
    result = dict()
    for f in os.listdir(root):
        m = pattern.match(os.path.basename(f))
        if m:
            user, t = parse_archive_name(**m.groupdict())
            files = result.get(user, {}) 
            files[f] = t
            result[user] = files
    return result


def classify_archives(u: str) -> Tuple[List[str], List[str]]:
    files = list(map(fst, sorted(get_archives(BACKUP_DIR).get(u, {}).items(), key=snd)))
    return files[:-KEEP], files[-KEEP:]


if __name__ == "__main__":
    for user in users:
        with tempfile.NamedTemporaryFile('w') as file_list:
            with file_list.file as handle:
                for f in discard_excluded(files_in_directory(user)):
                    print(f, file=handle)
            pack_user_archive(user, file_list.name)

            delete, _ = classify_archives(user)
            for outdated_archive in delete:
                os.remove(os.path.join(BACKUP_DIR, outdated_archive))
