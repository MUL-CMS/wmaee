from pyiron_base import state

def assign_partition_mulhpc():
    """
    Assigns a free partition from available partitions `p11` or `p12`.

    This function queries the remote system for idle partitions and assigns one of the available 
    partitions (`p11` or `p12`). If neither partition is available, it defaults to returning `p11`.

    Returns
    -------
    str
        The name of the assigned partition (`p11` or `p12`).
        If no free partition is found, returns `p11` by default.

    Notes
    -----
    - The function relies on the `sinfo` command to retrieve the state of partitions.
    - Partitions are considered free if they are marked as "idle" in the output of `sinfo`.
    - If the partition name ends with an asterisk (`*`), it is stripped out before returning.

    Examples
    --------
    >>> assign_partition_mulhpc()
    'p11'  # or 'p12', depending on availability
    """
    free_partitions = state.queue_adapter._adapter._execute_remote_command("sinfo | grep idle").split("\n")
    for partition in free_partitions:
        if partition:
            parts = partition.split()
            if "p11" in parts[0] or "p12" in parts[0]:
                use_part = parts[0].split("*")[0]  # p11 is always written like p11*
            else:
                use_part = "p11"
            return use_part
