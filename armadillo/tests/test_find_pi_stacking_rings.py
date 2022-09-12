# Import package, test suite, and other packages as needed
import armadillo
import pytest
import sys
import MDAnalysis as mda


def test_armadillo_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "armadillo" in sys.modules

u = mda.Universe('g3.pdb', guess_bonds=True)

ring_list = armadillo.find_all_rings(u)

thio_stacking, not_stacking = armadillo.find_pi_stacking_rings(ring_list, u, distmax=5.0, distmin=2.0,
                           angle_threshold_parallel=45, angle_threshold_tshaped=80)

def test_find_pi_stacking_rings():
    assert len(thio_stacking) == 135 and len(not_stacking) == 62
