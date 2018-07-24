import matplotlib.pyplot as plt
from input_file import vortex_length

def plot2d_subst(ax,biofilm):
    # Plots the concentration of substance at z=0 (plate)
    Lx, Ly, Lz = biofilm.length
    z_plane = np.zeros([int(Lx/vortex_length), int(Ly/vortex_length)])
    for vortex in biofilm.vortex_arr:
        if vortex.z == 1:
            z_plane[vortex.x, vortex.y] = vortex.conc_subst
    print(vortex.cqsm1)
    return ax.imshow(z_plane)


def plot2d_cellmass(ax,biofilm):
    # Plots the concentration of substance at z=0 (plate)
    Lx, Ly, Lz = biofilm.length
    Nx, Ny = int(Lx/vortex_length), int(Ly/vortex_length)
    z_plane = np.zeros([Nx, Ny])
    for vortex in biofilm.vortex_arr:
        if vortex.z == 0:
            z_plane[vortex.x, vortex.y] += vortex.get_mass()
    return ax.imshow(z_plane)

 
def plot3d(ax, biofilm):
    min_val = 1e8
    max_val = 1e-8
    Nx = int(biofilm.length[0] / vortex_length)
    Ny = int(biofilm.length[1] / vortex_length)
    Nz = int(biofilm.length[2] / vortex_length)
    array = [[[ None for z in range(Nz)] for y in range(Ny)] for x in range(Nz)]
    for vortex in biofilm.vortex_arr:
        x,y,z = vortex.x, vortex.y, vortex.z
        array[x][y][z] = vortex.conc_subst
        if max_val < vortex.conc_subst:
            max_val = vortex.conc_subst
        if min_val > vortex.conc_subst:
            min_val = vortex.conc_subst
    print("Minimum substrate concentration", min_val)
    print("Maximum substrate concentration", max_val)
    for x in range(Nx):
        for y in range(Ny):
            for z in range(Nz):
                if array[x][y][z] < 0:
                    c = '0.0'
                c = '%.1f' % ( 1 - (array[x][y][z]-min_val) / (max_val-min_val))            
                ax.scatter(x,y,z,color=c) 
