import numpy as np
import numpy.ma as ma
import scipy.io as sio
import copy
import smoothn
import matplotlib.pyplot as plt
import os.path


class MrStruct:
    # Attributes
    #     mrStructRaw: loadmat output dictionary
    #     dataAy:     3-D, 4-D, or 5-D single array
    #     dim1:       'size_x'
    #     dim2:       'size_y'
    #     dim3:       'size_z'
    #     dim4:       'size_t' | 'echos'
    #     dim5:       'unused' | 'size_t'
    #     vox:        1x4 double
    #     edges:      4x4 double
    #
    # Matlab structure [mrStruct]
    # dataAy:     3-D, 4-D, or 5-D single array
    # memoryType: ''
    # dim1:       'size_x'
    # dim2:       'size_y'
    # dim3:       'size_z'
    # dim4:       'size_t' | 'echos'
    # dim5:       'unused' | 'size_t'
    # dim6:       'unused'
    # dim7:       'unused'
    # dim8:       'unused'
    # dim9:       'unused'
    # dim10:      'unused'
    # dim11:      'unused'
    # vox:        1x4 double
    # edges:      4x4 double
    # orient:     ''
    # method:     ''
    # te:         double
    # tr:         double
    # ti:         []
    # patient:    char array (case ID)
    # user:       1x1 struct
    #     size_t_step:        double
    #     encoding_type:      'velocity'
    #     venc_in_plane:      double
    #     venc_thorugh_plane: double
    def __init__(self, filepath=""):
        if (filepath != ""):
            self.importMat(filepath)

    def importMat(self, filepath):
        self.mrStructRaw = sio.loadmat(filepath, squeeze_me=False)

        mrStruct = self.mrStructRaw['mrStruct']
        self.dataAy = mrStruct['dataAy'][0, 0]
        self.dim1 = mrStruct['dim1'][0, 0][0]
        self.dim2 = mrStruct['dim2'][0, 0][0]
        self.dim3 = mrStruct['dim3'][0, 0][0]
        self.dim4 = mrStruct['dim4'][0, 0][0]
        self.dim5 = mrStruct['dim5'][0, 0][0]
        self.dim6 = mrStruct['dim6'][0, 0][0]
        self.vox = mrStruct['vox'][0, 0][0]
        self.edges = mrStruct['edges'][0, 0]
        self.patient = mrStruct['patient'][0, 0]

        # print (f"Raw: {mrStruct[0,0]}")
        # print (f"Len: {len(mrStruct[0,0])}")
        # print (f"1: {mrStruct[0,0][2]}")
        # print (f"10: {mrStruct[0,0][10]}")
        # print (f"dtype: {mrStruct[0,0].dtype}")
        # print (f"names: {mrStruct[0,0].dtype.names}")
        mrStruct_names = mrStruct[0,0].dtype.names
        if "user" in mrStruct_names:
            self.user = dict()
            user = mrStruct['user'][0,0][0,0]
            user_names = user.dtype.names
            for i in range(len(user)):
                if(np.issubdtype(user[i].dtype, np.number)):
                    self.user[user_names[i]] = user[i][0,0]
                else:
                    self.user[user_names[i]] = user[i][0]
            print(self.user)

        self.name = os.path.basename(filepath)

    def Smooth_Velocity(self):
        shape = self.dataAy.shape
        e = shape[3]
        t = shape[4]
        print(f"Applying Velocity Filter")
        for t_i in range(t):
            for e_i in range(e):
                # s = 0.15 in source matlab code
                # Experimentally, s = 0.3 here gives similar result.
                # self.dataAy[:, :, :, e_i, t_i], _, _, _ = smoothn.smoothn(self.dataAy[:, :, :, e_i, t_i], s=0.3)
                self.dataAy[:, :, :, e_i, t_i], _, _, _ = smoothn.smoothn(self.dataAy[:, :, :, e_i, t_i], s=0.3)
            print(f"Applying Velocity Filter {(t_i + 1) / t * 100}%")

    def exportMat(self, filepath):
        mrStruct = np.array( \
            [(self.dataAy, self.dim1, self.dim2, self.dim3, self.dim4, self.dim5, self.dim6, self.vox, self.edges, \
              self.patient)], \
            dtype=[('dataAy', 'O'), ('dim1', 'O'), ('dim2', 'O'), ('dim3', 'O'), ('dim4', 'O'), ('dim5', 'O'), \
                   ('dim6', 'O'), ('vox', 'O'), ('edges', 'O'), ('patient', 'O')])
        sio.savemat(filepath, {'mrStruct': mrStruct}, appendmat=False, do_compression=True)

    def generateSpeedArray(self, mrStructMag, flag):
        if not ((self.dataAy.ndim == 5) and (self.dim4 == 'echos') and (self.dim5 == 'size_t')):
            return None

        imaSPEED1 = np.sqrt(np.sum(self.dataAy ** 2, axis=3))
        imaMAG = mrStructMag.dataAy
        imgMask = np.ones(imaMAG.shape[0:3])
        t = imaMAG.shape[3]

        # [Omitted] Mask selection
        # calculate MIP of time-averaged magnitude data
        # ImgMag = np.max(np.mean(imaMAG, 3), 2)
        # plt.figure()
        # plt.gray()
        # plt.imshow(ImgMag)
        # plt.show()

        minVWhole = np.min(imaMAG)
        maxVWhole = np.max(imaMAG)
        maxV = maxVWhole - (maxVWhole - minVWhole) * 0.3
        minV = minVWhole + (maxVWhole - minVWhole) * 0.0

        imaMAG = np.clip(imaMAG, minV, maxV)
        imaMAG = (imaMAG - minV) / (maxV - minV)

        imaSPEED = np.mean(imaSPEED1 ** 2, axis=3)
        imaMAG1 = np.mean(imaMAG, axis=3)

        if flag == "correlation":
            imaSPEED = (np.sum(imaSPEED1 ** 2, axis=3)) ** 0.2
            print("Calculating Mask")
            pcmraData1 = np.zeros(imaSPEED1.shape)
            for t_ind in range(t):
                # settings for waitbar
                print(f"Calculating Mask {t_ind / t * 100}%")
                # pcmraData1[:, :, :, t_ind] = (
                #         ((imaSPEED1[:, :, :, t_ind] ** 2) * imaMAG[:, :, :, t_ind]) * imgMask).astype(float)
                temp_smooth, _, _, _ = smoothn.smoothn((imaSPEED1[:, :, :, t_ind] ** 2) * imaMAG[:, :, :, t_ind],
                                                       s=0.10)
                pcmraData1[:, :, :, t_ind] = (temp_smooth * imgMask).astype(float)

            # Correlate across time
            # Pad the mag array
            padmag = np.pad(pcmraData1, ((1, 1), (1, 1), (1, 1), (0, 0)), 'constant', constant_values=0)

            # Shift the mag matrix
            # X-1
            Xminus = pearsonCC(pcmraData1, padmag[:-2, 1:-1, 1:-1, :], 3);
            # X+1
            Xplus = pearsonCC(pcmraData1, padmag[2:, 1: -1, 1: -1, :], 3);
            # Y-1
            Yminus = pearsonCC(pcmraData1, padmag[1:-1, :-2, 1: -1, :], 3);
            # Y+1
            Yplus = pearsonCC(pcmraData1, padmag[1:-1, 2:, 1: -1, :], 3);
            # Z-1
            Zminus = pearsonCC(pcmraData1, padmag[1:-1, 1: -1, :-2, :], 3);
            # Z+1
            Zplus = pearsonCC(pcmraData1, padmag[1:-1, 1: -1, 2:, :], 3);

            # Create a matrix to normalize the data
            normalization = 6 * np.ones(pcmraData1.shape[0:3])

            # Compensate for the edges
            normalization[0, :, :] = normalization[0, :, :] - 1;
            normalization[-1, :, :] = normalization[-1, :, :] - 1;
            normalization[:, 0, :] = normalization[:, 0, :] - 1;
            normalization[:, -1, :] = normalization[:, -1, :] - 1;
            normalization[:, :, 0] = normalization[:, :, 0] - 1;
            normalization[:, :, -1] = normalization[:, :, -1] - 1;

            # Store the data after summation and normalization
            pcmraData = (Xplus + Xminus + Yplus + Yminus + Zplus + Zminus) / normalization;
            pcmraData[np.isnan(pcmraData)] = 0;
            pcmraData[np.isinf(pcmraData)] = 0;

        else:
            pcmraData = np.multiply(np.multiply(imaSPEED, imaMAG1), imgMask)

        mrStructSpeed = copy.deepcopy(mrStructMag)
        mrStructSpeed.dataAy = pcmraData
        mrStructSpeed.name = "speed"
        return mrStructSpeed


# Helper function
def pearsonCC(a, b, dim):
    # Computes the Pearson correlation coefficient between elements in a and b
    # along dim. For example, in a dataset of X,Y,Z,T: dim 4 will compute the
    # correlation coefficient for each x,y,z voxel along the time dimension.
    # Remove means
    az = a - np.mean(a, dim, keepdims=True)
    bz = b - np.mean(b, dim, keepdims=True)
    # Standard Pearson correlation coefficient formula
    a2 = az ** 2;
    b2 = bz ** 2;
    ab = az * bz;
    r = np.sum(ab, dim) / np.sqrt(np.sum(a2, dim) * np.sum(b2, dim));

    return r
