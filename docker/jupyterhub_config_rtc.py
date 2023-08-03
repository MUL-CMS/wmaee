

# https://discourse.jupyter.org/t/plans-on-bringing-rtc-to-jupyterhub/9813/7
NON_ADMIN_USERS = ALL_USERS - ADMIN_USERS


def make_access_group_name(user: str) -> str:
    return f"rtc-access-{user}"


def make_manage_role_name(user: str) -> str:
    return f"rtc-manage-{user}"


# allow the admins to access non admins servers
load_groups = { make_access_group_name(user): list(ADMIN_USERS) for user in NON_ADMIN_USERS }

load_roles = [dict(
    name=make_access_group_name(user),
    description=f"RTC access role to {user}",
    scopes=[f'access:servers!user={user}'],
    groups=[make_access_group_name(user)]
) for user in NON_ADMIN_USERS]

# create the manage roles
managing_roles = [dict(
    name=make_manage_role_name(user),
    description=f"Manage users in group {make_access_group_name(user)}",
    scopes=[f'groups!group={make_access_group_name(user)}'],
    users=[user]
) for user in NON_ADMIN_USERS]

load_roles.extend(managing_roles)

c.JupyterHub.load_groups = load_groups
c.JupyterHub.load_roles = load_roles
c.Spawner.args = ['--LabApp.collaborative=True']