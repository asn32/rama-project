#!/usr/bin/env python3

import urllib
import math
import os
import argparse
from matplotlib import pyplot as plt

class Atom:
    def __init__(self,atom_id,atom_name,res_name,chain_id,res_id,x,y,z):
        self.id = int(atom_id)
        self.name = atom_name
        self.res = res_name
        self.chain_id = chain_id
        self.res_id = int(res_id)
        self.coord = Vector3D([float(x),float(y),float(z)])

    def __repr__(self):
        message=f"Atom: {self.name} in Residue: {self.res} in Chain: {self.chain_id}"
        return(message)
    
class Residue: 
    def __init__(self):
        self.name = None
        self.id = 0
        self.chain_id = None
        self.N=None
        self.C_alpha = None
        self.C_beta = None
        self.sidechain=[]
    
    def __repr__(self):
        backbone = ''.join([x for i,x in enumerate(['N ','Ca ','Cb ']) if [self.N,self.C_alpha,self.C_beta][i] is not None])
        message= f"{self.name} {self.id} {backbone}"
        return message
        
    
    def add_atom(self,atom):
        if atom.name == "N" and self.C_alpha == None and self.C_beta == None and len(self.sidechain) == 0:
            self.N = atom
            self.update_name(atom)

        elif atom.name == "CA" and self.N != None and self.C_beta == None and len(self.sidechain) == 0:
            self.C_alpha = atom
            
        elif atom.name == "C" and self.N != None and self.C_alpha != None and len(self.sidechain) == 0:
            self.C_beta = atom
            
        else:
            self.sidechain.append(atom)
        
        
    def update_name(self,atom):
        if atom.res != self.name and self.name is not None:
            raise Exception("Atom added not part of current residue")
        else:
            self.name = atom.res
            self.id = int(atom.res_id)
            self.chain_id = atom.chain_id
    
class Chain:
    def __init__(self):
        self.name=None
        self.residues=[]
        
    def __repr__(self):
        message=f"Chain {self.name}, {len(self.residues)} residues"
        return(message)
    
    def add_residue(self,residue):
        self.residues.append(residue)
        self.update_name(residue)
        
    def update_name(self,residue):
        if residue.chain_id !=self.name  and self.name is not None:
            raise Exception("Residue added not part of chain")
        else:
            self.name = residue.chain_id
            
class Protein: 
    def __init__(self,name):
        self.name = name
        self.chains = {}
    
    def __repr__(self):
        message=''.join(x.__repr__()+"\n" for x in self.chains.values())
        return(message)
    
    def add_chain(self,chain):
        self.chains.update({chain.name : chain})
    
    
class Vector3D:
    def __init__(self,entries):
        self.entries = entries
    
    def __repr__(self):
        return(self.entries.__repr__())
        
    def __add__(self,other):
        return(Vector3D([x+y for x,y in zip(self.entries,other.entries)]))
    
    def __sub__(self,other):
        return(Vector3D([x-y for x,y in zip(self.entries,other.entries)]))
    
    def norm(self):
        norm = math.sqrt(sum([x**2 for x in self.entries]))
        return(Vector3D([x/norm for x in self.entries]))
    
    def dot(self,other):
        return(sum([x*y for x,y in zip(self.entries,other.entries)]))
    
    def cross(self,other):
        return(Vector3D(
            [
            self.entries[1]*other.entries[2] - self.entries[2]*other.entries[1],
            self.entries[2]*other.entries[0] - self.entries[0]*other.entries[2],
            self.entries[0]*other.entries[1] - self.entries[1]*other.entries[0]
            ]
        ))
    
    @staticmethod
    def calculate_dihedral_angle(v1,v2,v3,v4):
        b1 = v2-v1
        b2 = v3-v2
        b3 = v4-v3

        n1 = b1.cross(b2).norm()
        n2 = b2.cross(b3).norm()
        m1 = n1.cross(b2.norm())

        x = n1.dot(n2)
        y = m1.dot(n2)

        return(-180*math.atan2(y,x)/math.pi)


class PdbParser:
    def load_local_file(self,path_to_file):
        with open(path_to_file) as file:
            f = file.read().splitlines()
        return(f)
    
    def load_rcsb(self,pdb_id): ## double check this 
        get_request = "https://files.rcsb.org/download/"
        
        try:
            request = urllib.request.urlopen(get_request+pdb_id+".pdb")
        except HTTPError:
            raise Exception("PDB ID could not be found in RCSB database.")

        file = request.read().decode().splitlines()
        return(file)
    
    def raw_to_atoms(self,list_of_strings):
        list_of_atoms = []
        for string in list_of_strings:
            if string[0:6] == 'ATOM  ':
                list_of_atoms.append(
                                    Atom(atom_id=string[6:11].strip(),
                                         atom_name=string[12:16].strip(),
                                         res_name=string[17:20].strip(),
                                         chain_id=string[21].strip(),
                                         res_id=string[22:26].strip(),
                                         x=string[30:38].strip(),
                                         y=string[38:46].strip(),
                                         z=string[46:54].strip()
                                        )
                                    )

        return(list_of_atoms)
        
    def atoms_to_protein(self,list_of_atoms):
        prot = Protein('name')
        
        chains=list(set([x.chain_id for x in list_of_atoms]))
        chains.sort()

        for chain in chains:
            current_chain = Chain()
            atoms_in_chain = [x for x in list_of_atoms if x.chain_id == chain]

            residues = list(set([x.res_id for x in atoms_in_chain]))
            residues.sort()

            for residue in residues:
                current_residue = Residue()
                atoms_in_residue = [x for x in atoms_in_chain if x.res_id == residue]

                for atom in atoms_in_residue:
                    current_residue.add_atom(atom)

                current_chain.add_residue(current_residue)

            prot.add_chain(current_chain)
        
        return(prot)
    
    def file_to_protein(self,path):
        
        if path[-3:]=="pdb" or path[-3:]=="PDB":
            list_of_strings = self.load_local_file(path)
            name = os.path.basename(path)
            
        else:
            list_of_strings = self.load_rcsb(path)
            name=path

        list_of_atoms = self.raw_to_atoms(list_of_strings)
        protein = self.atoms_to_protein(list_of_atoms)
        protein.name = name
        return(protein)
        

def ramachandran_plot(x,y,title='Ramachandran Plot',save=None):
    plt.figure(figsize=(6,6))
    plt.axis([-190,190,-190,190])
    plt.xlabel('$\phi$')
    plt.ylabel('$\psi$')
    plt.grid(linestyle='--')
    plt.xticks([-180,-90,0,90,180])
    plt.yticks([-180,-90,0,90,180])
    plt.suptitle(title)
    plt.plot(x,y,color='blue',marker='.',ls='',markersize=5)
    
    if save != None:
        plt.savefig(save)
    return(plt)



if __name__=="__main__" :
    
    argument_parser=argparse.ArgumentParser()

    ## required argument, the parser file 
    argument_parser.add_argument("file",type=str,nargs="+",metavar='file.pdb or PDBID',
                                 help="The PDB(s) to be analyzed. Either .pdb files, or PDB IDs searchable on the RCSB Protein Data Bank (http://rcsb.org)")

    argument_parser.add_argument("-i","--interesting",dest='plot_interesting',action='store_true',
                                 help="Whether to plot interesting residues (pre-prolines, glycines) in separate Ramachandran plots")

    argument_parser.add_argument("-o","--output",
                             default=os.getcwd()+'/',type=str,dest='output_path',
                                 help="Specify the output file or directory for the generated Ramachandran plot(s)")

    
    args = argument_parser.parse_args()
    
    ## Convert first argument to files, if it is a directory.
    if os.path.isdir(args.file[0])==True:
        args.file=[f for f in os.listdir(args.file[0]) if f[-3:]=='.pdb']

        
    base_outfile_name = 'ramachandran_plot.png'
    outfile_names=['glycine',
                   'preproline',
                   'others']
    
    
    ## Determine the output file name
    if os.path.isdir(args.output_path)==True:
        outfile_names=[args.output_path+x+'_'+base_outfile_name for x in outfile_names]
        base_outfile_name= args.output_path + base_outfile_name
        
    else: ## assume the specified path is to a file 
        outfile_names=[os.path.dirname(args.output_path)+'/'+x+'_'+os.path.basename(args.output_path) for x in outfile_names]
        base_outfile_name= args.output_path
                
       
        
    ## initialize relevant variables
    phi=[]
    psi=[]
    glycine=[]
    preproline=[]
    
    ## Initialize the parser object
    pdb_parser = PdbParser()
    ## Parse the files into protein structure objects
    proteins =[pdb_parser.file_to_protein(x) for x in args.file]

    ## Iterate over the individual proteins and calculate Psi and Phi angles. 
    for protein in proteins:
        for chain in protein.chains.values():
            for i in range(1,len(chain.residues)-1):    
                #Ci–1 −Ni −Ciα −Ci
                phi.append(Vector3D.calculate_dihedral_angle(chain.residues[i-1].C_beta.coord,
                                                        chain.residues[i].N.coord,
                                                        chain.residues[i].C_alpha.coord,
                                                        chain.residues[i].C_beta.coord
                                                             ))

                #Ni −Ciα −Ci −Ni+1
                psi.append(Vector3D.calculate_dihedral_angle(chain.residues[i].N.coord,
                                                      chain.residues[i].C_alpha.coord,
                                                      chain.residues[i].C_beta.coord,
                                                    chain.residues[i+1].N.coord
                                                             ))

                if(args.plot_interesting==True):
                    if chain.residues[i].name == 'GLY':
                        glycine.append(True)
                    elif chain.residues[i+1].name == 'PRO':
                        preproline.append(True)
                    else:
                        glycine.append(False)
                        preproline.append(False)
                        
                        
    
    ## Plot Ramachandran Plot 
    
    ## add additionally identifying information if the number of proteins == 1
    if len(proteins) == 1:
        name_add = " of " + proteins[0].name
    else:
        name_add = " of " + str(len(proteins)) +" PDB structures"

    ## Generate the plots as required
    if(args.plot_interesting==True):
        plt_glycine = ramachandran_plot([x for x,y in zip(phi,glycine) if y == True],
                                        [x for x,y in zip(psi,glycine) if y == True],
                                        'Glycine Residues'+name_add,save=outfile_names[0])  

        plt_preproline = ramachandran_plot([x for x,y in zip(phi,preproline) if y == True],
                                           [x for x,y in zip(psi,preproline) if y == True],
                                          'Pre-proline Residues'+name_add,save=outfile_names[1])


        plt_normal = ramachandran_plot([x for x,y,z in zip(phi,preproline,glycine) 
                                        if y == False and z == False],
                                       [x for x,y,z in zip(psi,preproline,glycine) 
                                        if y == False and z == False],
                                        'Other residues'+name_add,save=outfile_names[2])

    ## If generating an 'interesting' plot is not specified, generate a normal plot
    else:    
        plt = ramachandran_plot(phi,psi,'Ramachandran Plot'+name_add,save=base_outfile_name)
    