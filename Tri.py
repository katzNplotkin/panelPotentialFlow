#! /usr/bin/python3

import numpy as np
# from mpl_toolkits.mplot3d import axes3d
# import matplotlib.pyplot as plt
import meshio


eps = 1.0E-6
inv4pi = 0.07957747155

class Tri:
    """ Class for handling triangular vortex elements """
    def __init__(self, p0, p1, p2, ncap, rc):
        self.p = np.zeros((3, 3))
        self.p[:, 0] = p0
        self.p[:, 1] = p1
        self.p[:, 2] = p2
        self.ncap = ncap
        self.rc = rc
        self.cp = (self.p[:, 0] + self.p[:, 1] + self.p[:, 2])/3.0

    def vind(self, p):
        vel = self.vindFilament(0, p) + self.vindFilament(1, p) + \
                self.vindFilament(2, p)
        return vel

    def vindFilament(self, filamentNum, p):
        if filamentNum == 0:
            p1 = self.p[:, 0]
            p2 = self.p[:, 1]
        elif filamentNum == 1:
            p1 = self.p[:, 1]
            p2 = self.p[:, 2]
        else:
            p1 = self.p[:, 2]
            p2 = self.p[:, 0]

        r1 = p-p1
        r2 = p-p2
        r0 = r1-r2

        r1xr2 = np.cross(r1, r2)
        r1xr2abs2 = np.dot(r1xr2, r1xr2)

        r1Unit = r1/np.linalg.norm(r1)
        r2Unit = r2/np.linalg.norm(r2)

        vel = 0.0
        if r1xr2abs2 > eps*eps:
            vel = (r1xr2*inv4pi*np.dot(r0, r1Unit-r2Unit))/ \
                    np.sqrt((self.rc*np.linalg.norm(r0))**4.0+r1xr2abs2**2.0)

        return vel

class Surface:
    def __init__(self, rc, STLfilename, isClosed):
        stlData = meshio.read(STLfilename)
        self.mesh = stlData
        self.nElem = stlData.cells[0].data.shape[0]
        self.gamma = np.zeros(self.nElem)
        self.ele = []
        self.isClosed = isClosed
        for i in range(self.nElem):
            vtx = stlData.cells_dict['triangle'][i, :]
            p0 = stlData.points[vtx[0], :]
            p1 = stlData.points[vtx[1], :]
            p2 = stlData.points[vtx[2], :]

            try:
                ncap = stlData.cell_data['facet_normals'][0][i, :]
            except KeyError:
                ncap = np.cross(p1-p0, p2-p1)

            ncap = ncap/np.linalg.norm(ncap)
            self.ele.append(Tri(p0, p1, p2, ncap, rc))
        print('Computing AIC ...')
        self.aic = self.getAIC()

    def getAIC(self):
        aic = np.zeros((self.nElem, self.nElem))
        for i in range(self.nElem):
            for j in range(self.nElem):
                aic[i, j] = np.dot(self.ele[j].vind(self.ele[i].cp), \
                                   self.ele[i].ncap)
        if self.isClosed:
            # Assumes first element has gamma = 0.0
            return aic[1:, 1:]
        else:
            return aic

    def vind(self, p):
        vel = 0.0
        for i in range(self.nElem):
            vel += self.ele[i].vind(p)*self.gamma[i]
        return vel

    def getRHS(self, vinf):
        RHS = np.zeros((self.nElem))

        for i in range(self.nElem):
            RHS[i] = np.dot(np.array(vinf),  self.ele[i].ncap)

        if self.isClosed:
            return RHS[1:]
        else:
            return RHS

    def assignGamma(self, gamma):
        if len(gamma) == self.nElem-1:
            self.gamma[0] = 0.0
            for i in range(1, self.nElem):
                self.gamma[i] = gamma[i-1]
        else:
            self.gamma = gamma
        self.mesh.cell_data['gamma'] = self.gamma

    def writeVTK(self, gamma, VTKfilename):
        self.mesh.write(VTKfilename)

    def writeGrid(self, vinf, xlim, ylim, zlim, nvec, TECfilename):
        x = np.linspace(xlim[0], xlim[1], nvec[0])
        y = np.linspace(ylim[0], ylim[1], nvec[1])
        z = np.linspace(zlim[0], zlim[1], nvec[2])

        xg, yg, zg = np.meshgrid(x, y, z)
        u = np.zeros((nvec[0], nvec[1], nvec[2]))
        v = np.zeros((nvec[0], nvec[1], nvec[2]))
        w = np.zeros((nvec[0], nvec[1], nvec[2]))

        for k in range(nvec[2]):
            for j in range(nvec[1]):
                for i in range(nvec[0]):
                    u[i,j,k], v[i,j,k], w[i,j,k] = \
                            self.vind([xg[i, j, k], \
                                       yg[i, j, k], \
                                       zg[i, j, k]])

        u += vinf[0]
        v += vinf[1]
        w += vinf[2]

        # Write to tec file
        with open(TECfilename, "w") as fh:
            fh.write('TITLE = "Grid"\n')
            fh.write('VARIABLES = "X" "Y" "Z" "U" "V" "W"\n')
            fh.write('ZONE I='+ str(nvec[0]) + ' J=' + str(nvec[1])+ \
                    ' K=' + str(nvec[2]) + ' DATAPACKING=BLOCK')

            for dt in [xg, yg, zg, u, v, w]:
                fh.write("\n")
                for k in range(nvec[2]):
                    for j in range(nvec[1]):
                        for i in range(nvec[0]):
                            fh.write("%e\n" % dt[i, j, k])


body = Surface(1E-6, "sphere.stl", isClosed = True)
vinf = [1,0,0]
RHS = body.getRHS(vinf)
print('Computing solution ...')
gamma = np.linalg.solve(body.aic, RHS)
body.assignGamma(gamma)
body.writeVTK(gamma, "body.vtk")
# body.writeGrid(vinf, [0,5], [0,5], [0,5], [30, 30, 30], 'grid.tec')
body.writeGrid(vinf, [-250, 260], [-560,-70], [-85,420], [30, 30, 30], 'grid.tec')
