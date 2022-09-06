"""
armadillo.py
pi stacking project

Handles the primary functions
"""
import numpy as np
import pandas as pd
import MDAnalysis as mda
import networkx as nx
import itertools
import nglview
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
from MDAnalysis.lib import distances

def unit_vector(ag):
    # for a 2 member atomgroup, return the unit vector
    if len(ag) != 2: #changed to a 5 memebered group?
        raise ValueError("ERROR!")
    vec = ag[1].position - ag[0].position
    return vec / np.sqrt((vec * vec).sum())
#Calculates the unit vector which is required for the normal vector used to calculate alpha

def find_all_rings(self):
    g = nx.Graph()
    g.add_edges_from(u.bonds.to_indices())
    # cycle_basis gives all rings
    # but a fused ring
    cycles = nx.cycle_basis(g)
    # all atoms that are part of a ring
    #all 5 and 6 member rings plus indole (9) rings
    ring_members = set(itertools.chain(*cycles))

    rings=[u.atoms[c] for c in cycles] # this line creates a list of atomgroups, one per ring
    return rings

def ring_center(self):
    return self.center_of_geometry() #Gives coordinates of the center of the geometry of the ring

def normal_vector(ring):
    ag1 = ring[[0, 2]] #Plane between atom 1 and atom 3
    ag2 = ring[[1, 3]] #Plane between atom 2 and atom 4

    v1 = unit_vector(ag1) #Calculates vector of the ag1 plane
    v2 = unit_vector(ag2) #Calculates vector of the ag2 plane

    normal = np.cross(v1, v2) #Calculates normal via cross product method
    return normal / np.linalg.norm(normal)

def ring_distance(ring_a, ring_b):
    Distance = distances.calc_bonds(center_a, center_b, box=u.dimensions)
    #Takes two center geometry coords and caluclates distance between
    return Distance

def ring_normal_angle(ring_a, ring_b):
    cos_normal = np.dot(vector_a, vector_b/(np.linalg.norm(vector_a)*np.linalg.norm(vector_b)))
    #calcs the cos of between the two rings
    cos_squar_normal = cos_normal*cos_normal
    #calcs the cos squar so that normalised to 90 deg
    alpha_angle = np.rad2deg(np.arccos(cos_squar_normal))
    #calcs the alpha angle
    return alpha_angle
#Calcs dot product of two  vectors divided by normal vectors of both multiplied together

def pi_stacking(ring_a, ring_b, dist_max=5.0, alpha_max=np.deg2rad(45)):
    return ring_distance(ring_a, ring_b) < dist_max and ring_normal_angle(ring_a, ring_b) < alpha_max
    #distance between them is smaller than max distance and angle between them is smaller than max angle then pi stacking is seen

def find_pi_stacking_rings(list_of_rings, distmax, distmin,
                           angle_threshold_parallel, angle_threshold_tshaped):
    '''
    loops over an array-like list of ring objects. check if they are pi-stacking based on
    distance, where the max=5.0 and min=2.0
    Rings are then checked to see if they are in the same fragment, if they are then discarded
    Rings that pass this check point then have their alpha angle checked to create two lists.
    One list is stacking and one in which they aren't. Stacking includes t-shaped configuration

    Parameters
    -----------
    list_of_rings :
        List of atom groups which are rings to be stacked
    distmax : float
        Maximum cutoff distance between a ring pair
    distmin : float
        Minimum cutoff distance betweem a ring pair
    angle_theshold_parallel : float
        Maximum cutoff for angle between rings for parallel and parallel displaced
    angle_theshold_tshaped : float
        Minumum cutoff for t-shaped stacking angle

    Returns
    --------
    Returns a lists with n = [0] being ix_of_actually_pi_stacking and n= [1] ix_of_none_pi_stacking
    '''
    # use capped distance to select rings within 5 A, return indeces of pairs of rings
    ix_rings = mda.lib.distances.self_capped_distance(ring_cogs_array,
                                                    max_cutoff=distmax,
                                                    min_cutoff=distmin,
                                                    box=u.dimensions, return_distances=False)

    # checks if ring pairs are in same fragment, if so, we don't want to calculate pi stacking!
    ring_couples_to_stack = []
    for a, b in ix_rings: #Seperates ringA and ringB
        if ring_list[a].atoms.fragments!=ring_list[b].atoms.fragments: # check if same fragment
            ring_couples_to_stack.append([a, b]) #Adds the falses to ring_couples_to_stack


    # now check angle between normal of ring a and distance vector
    ix_of_actually_pi_stacking=[]
    ix_of_none_pi_stacking=[]
    #ix_t_shaped=[]
    for a, b in ring_couples_to_stack:

        vector_a = normal_vector(ring_list[a])
        vector_b = normal_vector(ring_list[b])
        #calculate angle between 2 vectors
        cos = (np.dot(vector_a, vector_b/(np.linalg.norm(vector_a)*np.linalg.norm(vector_b))))
        #calcs cos for the two rings
        cos_squar = cos*cos
        #cosine squared so that the angle is within 0-90 deg for alpha
        alpha = np.rad2deg(np.arccos(cos_squar)) #calculates the angle


        if alpha < angle_threshold_parallel or alpha > angle_threshold_tshaped:
            #angle_theshold for stacking chosen by user
            ix_of_actually_pi_stacking.append([a, b])
        else:
            ix_of_none_pi_stacking.append([a, b])

    return ix_of_actually_pi_stacking, ix_of_none_pi_stacking
