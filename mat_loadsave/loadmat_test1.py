################################################################################
# loadmat_test1.py
# Matlab .mat file input & output test using scipy package
# Read vel_struct.mat
# Access contents
# Modify and save in .mat format
# - save1: save as is
# - save2: modify using python array
# - save3: modify using numpy array
################################################################################

import scipy.io as sio
import numpy as np


#path_mag = "C:/d/project/Case_test/001A00505620161214/3dpc/mrstruct_subject_20161214_user011/mag_struct.mat"
path_vel = "C:/d/project/Case_test/001A00505620161214/3dpc/mrstruct_subject_20161214_user011/vel_struct.mat"
#path_mask = "C:/d/project/Case_test/001A00505620161214/3dpc/mrstruct_subject_20161214_user011/LA_LV_volume.mat"

mrStructRaw = sio.loadmat(path_vel, squeeze_me=False)
print("[[[[[[[[[[[[[[[[[[[[RAW]]]]]]]]]]]]]]]]]]]]")
print(mrStructRaw)

mrStruct = mrStructRaw['mrStruct'][0,0]
print("[[[[[[[[[[[[[[[[[[[[mrStruct]]]]]]]]]]]]]]]]]]]]")
print(mrStruct)
print(mrStruct.dtype)
print("dim1", mrStruct['dim1'])
print("dim1", mrStruct['dim1'][0])
print("vox", mrStruct['vox'])
print("vox", mrStruct['vox'][0])
print("edges", mrStruct['edges'])

mrStructUser = mrStruct['user'][0,0]
print("[[[[[[[[[[[[[[[[[[[[user 0,0]]]]]]]]]]]]]]]]]]]]")
print(mrStructUser)
print(mrStructUser.dtype)

venc_in_plane = mrStructUser['venc_in_plane'][0,0]
print("[[[[[[[[[[[[[[[[[[[[venc_in_plane]]]]]]]]]]]]]]]]]]]]")
print(venc_in_plane)
print(venc_in_plane.dtype)
print(mrStructUser['venc_in_plane'])
print(mrStructUser['venc_in_plane'].dtype)

mrStructUser = np.squeeze(mrStruct['user'])
print("[[[[[[[[[[[[[[[[[[[[user 0,0]]]]]]]]]]]]]]]]]]]]")
print(mrStructUser)
print(mrStructUser.dtype)

venc_in_plane = mrStructUser['venc_in_plane']
print("[[[[[[[[[[[[[[[[[[[[venc_in_plane]]]]]]]]]]]]]]]]]]]]")
print(venc_in_plane)
print(venc_in_plane.dtype)
# print(venc_in_plane[0])
# print(venc_in_plane[0,0])


# dataAy = mrStruct['dataAy'][0, 0]
# dim1 = mrStruct['dim1'][0, 0][0]
# dim2 = mrStruct['dim2'][0, 0][0]
# dim3 = mrStruct['dim3'][0, 0][0]
# dim4 = mrStruct['dim4'][0, 0][0]
# dim5 = mrStruct['dim5'][0, 0][0]
# dim6 = mrStruct['dim6'][0, 0][0]
# vox = mrStruct['vox'][0, 0][0]
# edges = mrStruct['edges'][0, 0]
# patient = mrStruct['patient'][0, 0]
#
# mrStructRaw = sio.loadmat(path_vel, squeeze_me=False, struct_as_record=False)
print(mrStructRaw.keys())
# print(mrStructRaw)
print(mrStructRaw['mrStruct'])
print(mrStructRaw['mrStruct'].dtype)
print(mrStructRaw['mrStruct'][0,0]['dataAy'].size)
print(mrStructRaw['mrStruct'][0,0]['dataAy'].shape)
print(mrStructRaw['mrStruct'][0,0]['dataAy'][0,0].shape)

print("dim1", mrStruct['dim1'])
print("dim1", mrStruct['dim1'][0])
print("dim1", np.squeeze(mrStruct['dim1'][0]))
print("vox", np.squeeze(mrStruct['vox'][0]))
print("edges", mrStruct['edges'])
print("edges", np.squeeze(mrStruct['edges']))

mrStructUser = mrStruct['user'][0,0]
print(mrStructUser)
print(mrStructUser.dtype)
venc_in_plane = mrStructUser['venc_in_plane']
print(venc_in_plane)
venc_in_plane = mrStructUser['venc_in_plane'][0,0]
print(venc_in_plane)
venc_in_plane = np.squeeze(mrStructUser['venc_in_plane'])
print(venc_in_plane)


print("vox", mrStruct['vox'])
print("vox", mrStruct['vox'][0])

sio.savemat("save1_original.mat", mrStructRaw)

mrStruct['vox'] = [[1.2,2.2,3.2,4.2]]
print("vox", mrStruct['vox'])
print("vox", mrStruct['vox'][0])
mrStructUser['venc_in_plane'] = 32
mrStructUser['encoding_type'] = "nose"
sio.savemat("save2_pythonArray.mat", mrStructRaw)

mrStruct['vox'] = np.array([[1.3,2.3,3.3,4.3]])
mrStruct['edges'] = np.array([
    [1.3,2.3,3.3,4.3],
    [1.2,2.2,3.2,4.2],
    [1.4,2.4,3.4,4.4],
    [0,0,0,1]
    ])
print("vox", mrStruct['vox'])
print("vox", mrStruct['vox'][0])
print("edges", mrStruct['edges'])
print("edges", np.squeeze(mrStruct['edges']))
sio.savemat("save3_npArray.mat", mrStructRaw)
