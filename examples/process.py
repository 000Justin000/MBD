import numpy as np
import sys

ori_dat  = np.loadtxt("graphene")
fal_dat  = np.zeros((0,4))
lat_vec  = [12.29750, 12.78000, 15.00000]
lat_vecs = np.array([[12.29750, 0.00000, 0.00000], [0.00000, 12.78000, 0.00000], [0.00000, 0.00000, 15.00000]])
lat_vecs[0,:] = lat_vecs[0,:] * int(sys.argv[1])
lat_vecs[1,:] = lat_vecs[1,:] * int(sys.argv[2])

for i in range(0, int(sys.argv[1])):
    for j in range(0, int(sys.argv[2])):
        cpy_dat = np.zeros((len(ori_dat), 4))
        cpy_dat[:,0] = ori_dat[:,0] + lat_vec[0] * i
        cpy_dat[:,1] = ori_dat[:,1] + lat_vec[1] * j
        cpy_dat[:,2] = ori_dat[:,2]
        cpy_dat[:,3] = ori_dat[:,3]
        fal_dat = np.concatenate((fal_dat, cpy_dat), axis=0)

outfilename = "graphene_" + sys.argv[1] + sys.argv[2] + "1.in"

np.savetxt(outfilename, fal_dat, delimiter="    ")

with open(outfilename, "a") as outfile:
        outfile.write("lattice_vector" + "    " + "{:.15e}".format(lat_vecs[0][0]) + "    " + "{:.15e}".format(lat_vecs[0][1]) + "    " + "{:.15e}".format(lat_vecs[0][2]) + "\n")
        outfile.write("lattice_vector" + "    " + "{:.15e}".format(lat_vecs[1][0]) + "    " + "{:.15e}".format(lat_vecs[1][1]) + "    " + "{:.15e}".format(lat_vecs[1][2]) + "\n")
        outfile.write("lattice_vector" + "    " + "{:.15e}".format(lat_vecs[2][0]) + "    " + "{:.15e}".format(lat_vecs[2][1]) + "    " + "{:.15e}".format(lat_vecs[2][2])       )
