#! /usr/bin/env python
from math import pi
outside_lattitude = 10. * pi/180.

meshfile = open('data/mesh.msh')
lines = meshfile.readlines()
meshfile.close()
nb_nodes = int(lines[4].strip())
nb_elems = int(lines[4+nb_nodes+3].strip())

print 'Mesh contains:'
print 'nodes: ',nb_nodes
print 'elems: ',nb_elems

nodes = []
coords = []
max_lat = 0.
for l in range(5,5+nb_nodes):
    split = lines[l].strip().split()
    coords.append([float(split[1]),float(split[2])])
    #print str(float(split[2]))
    if float(split[2]) >= outside_lattitude:
        nodes.append(int(split[0]))
        max_lat = max(max_lat, float(split[2]))
print 'selected '+str(len(nodes))+' nodes out of '+str(nb_nodes)
print 'max_lat = '+str(max_lat)+ ' greater than '+str(outside_lattitude)
if (len(nodes) == 0):
    raise RuntimeError, 'No nodes were selected'

elems = []
connectivity = []
for l in range(8+nb_nodes,8+nb_nodes+nb_elems):
    split = lines[l].strip().split()
    idx = int(split[0])
    elem_nodes = [int(n) for n in split[5:]]
    keep = True
    for n in elem_nodes:
        if not n in nodes:
            keep = False
            break
    #if keep and (coords[elem_nodes[0]-1][1]-coords[elem_nodes[1]-1][1] == 0. and coords[elem_nodes[0]-1][1]-coords[elem_nodes[2]-1][1] == 0.):
    #    keep = False
    if keep:
        elems.append(idx)
        connectivity.append( elem_nodes )
        #print idx, elem_nodes

print 'selected '+str(len(elems))+' elems'

def latlon_to_stereo(coords):
    from math import pi,sin,cos
    theta = coords[0]
    phi   = (coords[1]+pi/2.)
    R     = sin(phi)/(1.-cos(phi))
    Theta = theta
    X = R*cos(Theta)
    Y = R*sin(Theta)
    return [X,Y]

meshfile = open('data/mesh_stereo.msh','w')

meshfile.write("""$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
""")
meshfile.write(str(len(nodes))+'\n')
for n in nodes:
    [X,Y] = latlon_to_stereo(coords[n-1])
    meshfile.write(str(n)+' '+str(X)+' '+str(Y)+' 0\n')

meshfile.write('$EndNodes\n$Elements\n')
meshfile.write(str(len(elems))+'\n')
for e, conn in zip(elems,connectivity):
    meshfile.write(str(e)+' 1 3 1 1 1 '+' '.join([str(n) for n in conn])+'\n')
meshfile.write('$EndElements\n')
meshfile.close()