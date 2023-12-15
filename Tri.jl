using WriteVTK
using VSPGeom
using LinearAlgebra
using Infiltrator

inv4pi = 0.07957747155
eps2 = eps()^2

""" Class for handling triangular vortex elements """
struct Tri
    p
    ncap
    rc
    cp
end

""" Constructor """
function Tri(p, ncap, rc)
    cp = zeros(3)
    cp .= sum(p; dims=2) / 3
    return Tri(p, ncap, rc, cp)
end

function Tri(p1, p2, p3, ncap, rc)
    p = zeros(3, 3)
    p[:, 1] = p1
    p[:, 2] = p2
    p[:, 3] = p3
    return Tri(p, ncap, rc)
end

""" Induced velocity """
function vind(self::Tri, p)
    return vindFil(self, 1, p) + vindFil(self, 2, p) + vindFil(self, 3, p)
end


""" Induced velocity by filament """
function vindFil(self::Tri, ifil::Int, p::Vector{Float64})
    if ifil < 3
        p1 = self.p[:, ifil]
        p2 = self.p[:, ifil+1]
    else
        p1 = self.p[:, 3]
        p2 = self.p[:, 1]
    end

    r1 = p-p1
    r2 = p-p2
    r0 = r1-r2

    r1xr2 = cross(r1, r2)
    r1xr2abs2 = dot(r1xr2, r1xr2)

    r1Unit = r1/norm(r1)
    r2Unit = r2/norm(r2)

    vel = 0.0
    if r1xr2abs2 > eps2
        vel = (r1xr2*inv4pi*dot(r0, r1Unit-r2Unit)) / sqrt((self.rc*norm(r0))^4.0+r1xr2abs2^2.0)
    end
    return vel
end

struct Surface{TI, TF}
    mesh::TriMesh
    ele::Vector{Tri}
    nElem::TI
    gamma::Vector{TF}
    isClosed::Bool
    aic::Matrix{TF}
end

""" Constructor """
function Surface(filename::String; rc=1e-6, isClosed=true, tol=1e-6)
    geom = readSTL(filename)
    mesh = geom[1]
    nElem = length(mesh.cells)
    p = zeros(3, 3)
    gamma = zeros(nElem)

    ele = Vector{Tri}(undef, nElem)
    ncap = zeros(3)
    for i in 1:nElem
        p[:, 1] = mesh.points[mesh.cells[i]][1]
        p[:, 2] = mesh.points[mesh.cells[i]][2]
        p[:, 3] = mesh.points[mesh.cells[i]][3]

        ncap .= cross(p[:, 2]-p[:, 1], p[:, 3]-p[:, 2])
        ele[i] = Tri(p, ncap, rc)
    end

    # Compute aic matrix
    aic = zeros(nElem, nElem)

    idxStart = isClosed ? 2 : 1
    for i in idxStart:nElem
        for j in idxStart:nElem
            aic[i, j] = dot(vind(ele[j], ele[i].cp), ele[i].ncap)
        end
    end

    return Surface(mesh, ele, nElem, gamma, isClosed, aic)
end

function vind(self::Surface, p)
    vel = 0.0
    for i in range(1:self.nElem)
        vel += vind(self.ele[i], p) .* self.gamma[i]
    end
    return vel
end

function getRHS(self::Surface, vinf)
    RHS = zeros(self.nElem)

    for i in range(1:self.nElem)
        RHS[i] = -1 .* dot(vinf, self.ele[i].ncap)
    end

    return self.isClosed ? RHS[2:end] : RHS
end

function writeMesh(s::Surface, filename::String)
    points, cells = getVTKElements(s.mesh)
    vtk_grid(filename, points, cells) do vtk
    end
end

#     def assignGamma(self, gamma):
#         if len(gamma) == self.nElem-1:
#             self.gamma[0] = 0.0
#             for i in range(1, self.nElem):
#                 self.gamma[i] = gamma[i-1]
#         else:
#             self.gamma = gamma
#         self.mesh.cell_data['gamma'] = self.gamma
#
#     def writeVTK(self, gamma, VTKfilename):
#         self.mesh.write(VTKfilename)
#
#     def writeGrid(self, vinf, xlim, ylim, zlim, nvec, TECfilename):
#         x = np.linspace(xlim[0], xlim[1], nvec[0])
#         y = np.linspace(ylim[0], ylim[1], nvec[1])
#         z = np.linspace(zlim[0], zlim[1], nvec[2])
#
#         xg, yg, zg = np.meshgrid(x, y, z)
#         u = np.zeros((nvec[0], nvec[1], nvec[2]))
    #         v = np.zeros((nvec[0], nvec[1], nvec[2]))
#         w = np.zeros((nvec[0], nvec[1], nvec[2]))
#
#         for k in range(nvec[2]):
#             for j in range(nvec[1]):
#                 for i in range(nvec[0]):
#                     u[i,j,k], v[i,j,k], w[i,j,k] = \
#                             self.vind([xg[i, j, k], \
#                                        yg[i, j, k], \
#                                        zg[i, j, k]])
#
#         u += vinf[0]
#         v += vinf[1]
#         w += vinf[2]
#
#         # Write to tec file
#         with open(TECfilename, "w") as fh:
#             fh.write('TITLE = "Grid"\n')
#             fh.write('VARIABLES = "X" "Y" "Z" "U" "V" "W"\n')
#             fh.write('ZONE I='+ str(nvec[0]) + ' J=' + str(nvec[1])+ \
#                     ' K=' + str(nvec[2]) + ' DATAPACKING=BLOCK')
#
#             for dt in [xg, yg, zg, u, v, w]:
#                 fh.write("\n")
#                 for k in range(nvec[2]):
#                     for j in range(nvec[1]):
#                         for i in range(nvec[0]):
#                             fh.write("%e\n" % dt[i, j, k])
#
#

body = Surface("sphere.stl")
writeMesh(body, "out")
# vinf = [1,0,0]
# RHS = body.getRHS(vinf)
# print('Computing solution ...')
# gamma = np.linalg.solve(body.aic, RHS)
# body.assignGamma(gamma)
# body.writeVTK(gamma, "body.vtk")
# # body.writeGrid(vinf, [0,5], [0,5], [0,5], [30, 30, 30], 'grid.tec')
# body.writeGrid(vinf, [-250, 260], [-560,-70], [-85,420], [30, 30, 30], 'grid.tec')
