import os
import math
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from Bio.PDB import PDBList, PDBParser, is_aa

def torsion_angles (p1,p2,p3,p4):
    """Calculate the torsion angle between plane p1-p2-p3 and plane p2-p3-p4. 
    P1,p2,p3,p4 are three dimensional coordinates"""
    b1 = p2 - p1
    b2 = p3 - p2
    b3 = p4 - p3
    if (np.cross(b1,b2)== (0.0, 0.0, 0.0)).all() or (np.cross(b2,b3)== (0.0, 0.0, 0.0)).all() :
        return 180.0
    else:
        n1 = np.cross(b1,b2) / np.linalg.norm(np.cross(b1,b2))
        n2 = np.cross(b2,b3) / np.linalg.norm(np.cross(b2,b3))
        m1 = np.cross(n1,b2 /np.linalg.norm(b2) )
        x = np.dot(n1,n2)
        y = np.dot(m1,n2)
        return -math.atan2(y,x) / math.pi * 180.0

def ramachandran_plot(phi,psi):
    fig = plt.figure(figsize=(8, 8))
    gs = GridSpec(4, 4, fig)
    # 主图：Ramachandran plot
    main = fig.add_subplot(gs[1:4, 0:3])
    sns.kdeplot(x=phi, y=psi, cmap='viridis')
    main.scatter(phi, psi, alpha=0.5, color='black', s=10)
    main.set_xlim(-180, 180)
    main.set_ylim(-180, 180)
    main.set_xlabel('φ (radians)')
    main.set_ylabel('ψ (radians)')
    # 上方的直方图
    top = fig.add_subplot(gs[0, 0:3], sharex=main)
    top.hist(phi, bins=30, color='gray', edgecolor='black')
    top.axis('off')
    # 右侧的直方图
    right = fig.add_subplot(gs[1:4, 3], sharey=main)
    right.hist(psi, bins=30, orientation='horizontal', color='gray', edgecolor='black')
    right.axis('off')
    plt.tight_layout()
    plt.show()

os.chdir(os.path.dirname(__file__))
# 研究对象的PDB id
pdb_ids = ['3j5p','6kxs','6l2t','6lmk','9flx','6lth','8yl4','8z3p','9b4p','9dio']
proteins = []
for pdb_id in pdb_ids:
    pdbl = PDBList()
    pdbl.retrieve_pdb_file(pdb_id, pdir='.', file_format='pdb')
    parser = PDBParser()
    structure = parser.get_structure(pdb_id, 'pdb'+pdb_id+'.ent')
    proteins.append(structure)

psi_list = []
phi_list = []
omega_list = []
for protein in proteins:
    model = protein[0]
    for chain in model:
        prev_residue = None
        for residue in chain:
            if not is_aa(residue):
                continue        
            if prev_residue:
                try:
                    next_residue = chain[residue.id[1] + 1]
                    ca_prev = prev_residue['CA'].coord
                    c_prev = prev_residue['C'].coord
                    n_prev = prev_residue['N'].coord
                    ca = residue['CA'].coord
                    c = residue['C'].coord
                    n = residue['N'].coord
                    ca_next = next_residue['CA'].coord
                    c_next = next_residue['C'].coord
                    n_next = next_residue['N'].coord
                    psi = torsion_angles(n,ca,c,n_next)
                    phi = torsion_angles(c_prev,n,ca,c)
                    omega = torsion_angles(ca_prev,c_prev,n,ca)
                    if residue.get_resname() == 'GLY':
                    # if residue.get_resname() == 'PRO':
                    # if prev_residue.get_resname() == 'PRO':
                        psi_list.append(psi)
                        phi_list.append(phi)
                        omega_list.append(omega)
                except (KeyError, IndexError):
                    pass
            prev_residue = residue
    # ramachandran_plot(np.array(phi_list),np.array(psi_list))
ramachandran_plot(np.array(phi_list),np.array(psi_list))
