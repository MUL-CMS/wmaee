# apply home-made patches
# see folder patches for details
from wmaee.core.requirements import test_pmg
if test_pmg():
    import wmaee.codes.patches.pmg_vasprun