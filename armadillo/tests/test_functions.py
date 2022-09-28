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
    """"Tests that the find_all_rings function finds the correct number of rings
    in the example universe"""
    assert len(armadillo.find_all_rings(u)) == 767


def test_unit_vector():
    """Test that unit_vector functions find a error when a atom group of more
    than two atoms is passed through it"""
    ring_list = armadillo.find_all_rings(u)
    ag = ring_list[0]
    assert armadillo.unit_vector(ag) == 'ERROR'


def test_normal_vector():
    """Test that the normal_vector function finds a three dimensional vector"""
    ring_list = armadillo.find_all_rings(u)
    ring = ring_list[0]
    assert len(armadillo.normal_vector(ring)) == 3


def test_ring_normal_angle():
    """Tests that ring_normal_angle returns the correct angle between two rings
    that have a known angle between them"""
    ring_list = armadillo.find_all_rings(u)
    ring_a = ring_list[0]
    ring_b = ring_list[1]
    assert armadillo.ring_normal_angle(ring_a, ring_b) == pytest.approx(67.65, 0.05)

def test_ring_distance():
    """Tests that ring_distance returns the correct distance between the centers
    of two rings that have a known distance between them"""
    ring_list = armadillo.find_all_rings(u)
    ring_a = ring_list[0]
    ring_b = ring_list[1]
    assert armadillo.ring_distance(ring_a, ring_b, u) == pytest.approx(3.948, 0.01)

def test_pi_stacking_individual():
    """Tests that when two rings that are known not to stack due to alpha angle
    being too large a false value is returned"""
    ring_list = armadillo.find_all_rings(u)
    ring_a = ring_list[0]
    ring_b = ring_list[1]
    assert armadillo.pi_stacking(ring_a, ring_b, u, dist_max=5.0, alpha_max=45) == False

def test_find_pi_stacking_rings():
    """Tests that when a list of rings is passed through the function two lists
    are created that have a known length for stacking and not stacking"""
    ring_list = armadillo.find_all_rings(u)
    thio_stacking, t_shaped = armadillo.find_pi_stacking_rings(ring_list, u, distmax=5.0, distmin=2.0,
                               angle_threshold_parallel=45, angle_threshold_tshaped=80)
    assert len(thio_stacking) == 118 and len(t_shaped) == 17
