# apply home-made patches
# see folder patches for details
from wmaee.core.config import is_pmg_avail
if is_pmg_avail:
    import wmaee.codes.patches.pmg_vasprun