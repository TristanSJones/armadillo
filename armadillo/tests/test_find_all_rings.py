# Import package, test suite, and other packages as needed
import armadillo
import pytest
import sys
import MDAnalysis as mda


def test_armadillo_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "armadillo" in sys.modules

u = mda.Universe('g3.pdb', guess_bonds=True)
#test file with known amount of rings

def test_find_all_rings():
    assert len(armadillo.find_all_rings(u)) == 767
    """"Tests that the find_all_rings function finds the correct number of rings
    in the example universe"""
