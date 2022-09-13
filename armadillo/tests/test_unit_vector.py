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
ag = ring_list[0]

def test_unit_vector():
    assert armadillo.unit_vector(ag) == "ERROR!"
