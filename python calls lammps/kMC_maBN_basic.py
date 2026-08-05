import numpy as np
from pymatgen.core import Structure, Lattice
from pymatgen.io.vasp import Poscar
from scipy.spatial import cKDTree
import os
import faulthandler
faulthandler.enable()

def writeFile(fn, s):
    with open(fn, 'w') as fout:
        fout.write(s)
    return

def writeLammpsData(structure, fn):
    out = [
        "maBN\n",
        "%d atoms\n"    % structure.num_sites,
        "%d atom types\n" % structure.ntypesp,
        "0.000000 %.6f  xlo xhi"%structure.lattice.a,
        "0.000000 %.6f  ylo yhi"%structure.lattice.b,
        "0.000000 %.6f  zlo zhi"%structure.lattice.c,
        "\nMasses\n",
        "1  14.0067", 
        "2  10.8110", 
        "\nAtoms\n"
    ]
    for i,site in enumerate(structure.sites):
        out.append("%d %d %.6f %.6f %.6f"%((i+1, 1 if str(site.specie)=="N" else 2)  + tuple(list(site.coords))))
    writeFile(fn, "\n".join(out))

def parselammps(fn):
    with open(fn,'r') as fin2:
        lines=fin2.readlines()
        num=int(lines[2].strip().split()[0])
        npairs=num/2
        [xmin,xmax]=lines[5].strip().split()[0:2]
        [ymin,ymax]=lines[6].strip().split()[0:2]
        [zmin,zmax]=lines[7].strip().split()[0:2]
        cord=np.array([i.strip().split() for i in lines[16:16+num]])
        Coords=cord[:,2:5].astype(np.float64)
        specs=['N' if i=='1' else "B" for i in cord[:,1] ]
        lat=[[float(xmax)-float(xmin),0,0],[0,float(ymax)-float(ymin),0],[0,0,float(zmax)-float(zmin)]]
    a=Structure(lattice=lat,species=specs,coords=Coords,coords_are_cartesian=True)
    a.sort()
    with open("./log.lammps", 'r') as fin:
        for line in fin:
            if "Energy initial, next-to-last, final =" in line: break
        line = fin.readline().strip().split()
        E_relaxed = float(line[2])
    return a,E_relaxed

def findAtomPair(structure, cutoff=2.0):
    import itertools
    lattice = structure.lattice
    coords = structure.cart_coords
    pos = []
    for i, j in itertools.product([0,1],[0,1]):
        pos+=list(coords+[i*lattice.a, j*lattice.b, 0])
    kd = cKDTree(pos)
    pairs = np.array(list(kd.query_pairs(cutoff)))
    pairs = np.array(list(filter(lambda x: np.any(x<structure.num_sites), pairs)))
    pairs = pairs%structure.num_sites 
    return pairs[np.random.randint(len(pairs))]

def rotateSingleBond(structure):
    gindex1, gindex2 = findAtomPair(structure)
    pos1 = structure.cart_coords[gindex1]
    pos2 = structure.cart_coords[gindex2]
    center = (pos1+pos2)/2
    pos1 -= center
    pos2 -= center
    pos1 = [pos1[1], -pos1[0], 0] + center
    pos2 = [pos2[1], -pos2[0], 0] + center        
    structure.sites[gindex1].coords = pos1
    structure.sites[gindex2].coords = pos2
    return gindex1, gindex2

def exchangeBond(structure):
    gindex2=gindex1=0
    while structure.species[gindex1]==structure.species[gindex2]:
        gindex1, gindex2 = findAtomPair(structure)
    pos1, pos2 = structure.cart_coords[[gindex1, gindex2]]
    structure.sites[gindex2].coords=pos1
    structure.sites[gindex1].coords=pos2

def runlammps():
    lammps_exe="lmp"
    os.system("%s -in maBN-mlp.in" %(lammps_exe))

def run(nstep):
    runlammps()
    structure, E_new = parselammps("relaxed.data")
    E_old = E_new
    kBT = 1.0
    log = open("maBN_MC.log", "a+")
    r = np.random.rand(nstep,2)
    for step in range(nstep):
        oldStructure = structure.copy()
        if r[step][0] >= 0.5:
            rotateSingleBond(structure)
        else:
            exchangeBond(structure)
        writeLammpsData(structure,"maBN.data")
        runlammps()
        structure, E_new = parselammps("relaxed.data")
        p = np.exp(-(E_new-E_old)/kBT)
        if r[step][1]>p:
            structure = oldStructure
        else:
            E_old = E_new
            writeLammpsData(structure,"./accepted/step-%d.data"%(step+1+i))
        s = ", ".join([
                "Step: %7d"%(step+1+i),
                "Eold: %12.5f"%E_old,
                "Enew: %12.5f"%E_new,
                "%s" %("Denied" if r[step][1]>p else "Accepted")
            ])
        log.write(s+"\n")
        log.flush()
    log.close()
    return

i = 0
run(30000)
