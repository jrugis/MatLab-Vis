# -*- coding: utf-8 -*-
#
# mat2vtk.py

import numpy as np
import matplotlib.pyplot as plt
import os
import scipy.io as sc
import struct
import subprocess
from evtk.hl import unstructuredGridToVTK
from evtk.vtk import VtkTriangle

##################################################################
# functions
##################################################################

def write_tris(fname, verts, tris, data=None):
  nverts = verts.shape[0]
  xyz = np.empty([3, nverts]) # because it needs re-ordering
  for i in range(nverts):
    for j in range(3):
      xyz[j,i] = verts[i,j]  
  ntris = tris.shape[0]
  conn = np.empty(3*ntris)
  for i in range(ntris):
    conn[3*i] = tris[i,0]-1 # index from zero
    conn[3*i+1] = tris[i,1]-1
    conn[3*i+2] = tris[i,2]-1
  offset = np.zeros(ntris, dtype=int)
  for i in range(ntris):
    offset[i] = 3*(i+1) 
  ctype = np.zeros(ntris)
  for i in range(ntris):
    ctype[i] = VtkTriangle.tid 
  unstructuredGridToVTK(fname, \
    xyz[0,:], xyz[1,:], xyz[2,:], \
    connectivity=conn, offsets=offset, cell_types=ctype, \
    cellData=data, pointData=None)  # write out vtu file
  return

##################################################################
# main program
##################################################################

#---------------------------------------------------
# read Nathan's apical data file 
apical_name = 'conductancesFc0Fip0t200nc7adec59140448e185c5f2d9b2376a32ff.mat'
apical_key = 'w_IPR'

print('reading apical data file: ' + apical_name) 
dist = sc.loadmat(apical_name)
#print('keys:', dist.keys())

apical_data = dist[apical_key]
#print('apical data shape:')
#for i in range(7):
#  print(apical_data[i,0].shape)

#---------------------------------------------------
# read Nathan's mesh data file
mesh_name = 'mod_basal1data_smoothed_mesh.mat'
verts_key = 'p'
tris_key = 'triangles'

print('reading mesh data file: ' + mesh_name) 
dist = sc.loadmat(mesh_name)
#print('keys:', dist.keys())

verts = dist[verts_key]
#print('verts data shape:')
#for i in range(7):
#  print(verts[i,0].shape)

tris = dist[tris_key]
#print('tris data shape:')
#for i in range(7):
#  print(tris[i,0].shape)
#  print(np.amin(tris[i,0]))

#---------------------------------------------------
# write out an apical VTU file for each cell
aname = 'apical_'

print('writing VTU files:')
for i in range(7):
  tname = aname + str(i+1)
  print(tname)
  cell1_verts = verts[i,0]
  t1 = np.nonzero(apical_data[i,0]) # indices of all marked with 45
  t2 = np.array(t1[0])              # make it numpy array type
  cell1_tris = tris[i,0][t2,:]      # apical subset of tris (zero based indexing!)

  write_tris(tname, cell1_verts, cell1_tris)

#---------------------------------------------------







