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
ring_a = ring_list[0]
ring_b = ring_list[1]

def test_ring_normal_angle():
    assert armadillo.ring_normal_angle(ring_a, ring_b) == pytest.approx(67.65, 0.05)
