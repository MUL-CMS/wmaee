from pyiron_base import state

def assign_partition_mulhpc():
    # returns a free partition, either p11 or p12
    # if there is no free partition the p11 partition is returned automatically
    free_partitions = state.queue_adapter._adapter._execute_remote_command("sinfo | grep idle").split("\n")
    for partition in free_partitions:
        if partition:
            parts = partition.split()
            if "p11" in parts[0] or "p12" in parts[0]:
                use_part = parts[0].split("*")[0] # p11 is always written like p11*
            else:
                use_part = "p11"
            return use_part
