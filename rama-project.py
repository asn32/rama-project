#!/usr/bin/env python3
"""Ramachandran Plot project

This script generates Ramachandran plots from a PBD structure file. Structure
files can be provided as paths to local .pdb files, or PDBIDs searchable on the RCSB
Protein Data Bank (http://www.rcsb.org).

If a directory is provided as the first argument to the script, then the directory is
searched for any .pdb files that can be used.

If multiple structure files are provided, in either .pdb or PDBID form, the generated
Ramachandran plot will contain all residues present in all structures.

If the "-i" flag is provided, separate Ramachandran plots for glycine, pre-proline
and all other residues will be generated. Plots for 'glycine', 'pre-proline', or 'others'
can be indicated to generate that specific plot.

This script requires matplotlib (>2.0) for plotting functionality, and
standard library modules urllib, math, os, argparse and f-strings (> Python 3.5) for functionality.

"""

import urllib
import math
import os
import argparse
from matplotlib import pyplot as plt


class Vector3D:
    """
    A  wrapper class for a list that implements basic 3D vector math.

    Attributes
    ----------
    entries : list of floats
        A list containing coordinate entries (3)

    Methods
    -------
    __repr__:
        calls the list.__repr__() function
    __add__:
        Implements vector addition as elementwise addition
    __sub__:
        Implements vector subtraction as the elementwise subtraction of the left vector from the right vector
    norm:
        Returns the Euclidean norm of a vector. Returns Vector3D object.
    dot:
        Implements the dot product between object instance and other Vector3D object.
        Returns Vector3D object.
    cross:
        Implements the 3D cross product between the object instance and other Vector3D object.
        Returns Vector3D object.
    calculate_dihedral_angle:
        Calculates the dihedral angle using four Vector3D objects, using the algorithim suggested
        in the project guidelines

    """
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



class Atom:
    """
    A class used to represent a single Atom in a PDB file.

    Attributes
    ----------
    atom_id : int
        The numeric ID of the atom within the PBD file
    atom_name : str
        The name of the atom ex. CA, or H
    res : str
        The name of the residue which the atom is a member of
    chain_id : str
        The chain of the protein that the atom is the member of
    res_id : int
        The residue number of the residue the atom is a member of
    coord: Vector3D
        A vector containing the Euclidean coordinates of the atom

    Methods
    -------
    __repr__:
        Returns an f-string formatted string indicating the atom and its residue and chain membership.

    """
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
    """
    A class used to represent a single Residue in a PDB file.

    Attributes
    ----------
    name : str
        The name of the residue ex. PHE, LEU etc.
    id : int
        The numeric ID of the residue based on the structure file.
    chain_id : str
        The chain of the protein that the atom is the member of
    N : Atom
        An Atom instance referring to the Nitrogen atom of the amino acid backbone of the residue
    C_alpha : Atom
        An Atom instance referring to the alpha carbon atom of the amino acid backbone of the residue
    C_beta : Atom
        An Atom instance referring to the beta carbon atom of the amino acid backbone of the residue
    C_alpha : List of Atom objects
        A list containing Atom instances referring to atoms part of the residue sidechain.

    Methods
    -------
    __repr__:
        Returns an f-string formatted string indicating the residue the atoms it contains.

    add_atom:
        Adds an atom to a residue, sorting the atom into the correct group - C_alpha, C_beta, N or sidechain

    update_name:
        Updates the name of the residue instance and throws an error if the atom added is derived
        from a different residue than the current one.

    """
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
        ## If the atom is a Nitrogen, and an Alpha or Beta carbon has not been added:
        if atom.name == "N" and self.C_alpha == None and self.C_beta == None and len(self.sidechain) == 0:
            self.N = atom
            self.update_name(atom)

        ## If the atom is an alpha carbon, and a Nitrogen has been added, but not a beta carbon:
        elif atom.name == "CA" and self.N != None and self.C_beta == None and len(self.sidechain) == 0:
            self.C_alpha = atom

        ## If the atom is an beta carbon, and a Nitrogen and a beta carbon have been added:
        elif atom.name == "C" and self.N != None and self.C_alpha != None and len(self.sidechain) == 0:
            self.C_beta = atom

        ## If the alpha carbon, beta carbon and nitrogen atoms have been added, atom must be part of sidechain
        else:
            self.sidechain.append(atom)


    def update_name(self,atom):
        ## If the added atom belongs to a different residue, and it is not the first update:
        if atom.res != self.name and self.name is not None:
            ## Then raise an exception
            raise Exception("Atom added not part of current residue")

        else: ## Otherwise, update the name, id, and chain_id of the residue.
            self.name = atom.res
            self.id = int(atom.res_id)
            self.chain_id = atom.chain_id


class Chain:
    """
    A class used to represent a single protein Chain in a PDB file.

    Attributes
    ----------
    name : str
        The name of the chain ex. A,B,C etc.
    residues : list of Residue objects
        A list of residues that belong to the chain, ordered accordig to their ID values.

    Methods
    -------
    __repr__:
        Returns an f-string formatted string indicating the chain and the residues it contains.

    add_residue:
       Adds a residue to the Chain object to the list of residues.

    update_name:
        Updates the name of the residue, throws an exception if residue does not share the same chain_id.

    """
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
        ## If the residue does not belong to the chain, and is not the first update
        if residue.chain_id !=self.name  and self.name is not None:
            ## Then raise an exception
            raise Exception("Residue added not part of chain")

        else: ## Otherwise, update the chain
            self.name = residue.chain_id



class Protein:
    """
    A class used to represent a Protein in a PDB file.

    Attributes
    ----------
    name : str
        The name of the protein.
    residues : dictionary of Residue objects
        A dictionary of residues that belong to the chain, ordered according to their ID values, with the chain ID as
        the key.

    Methods
    -------
    __repr__:
        Returns an f-string formatted string indicating the residue the atoms it contains.

    add_chain:
        Method to a chain to the Protein object. Stores chains in a dictionary, with the key as the chain ID.

    """
    def __init__(self,name):
        self.name = name
        self.chains = {}

    def __repr__(self):
        message=''.join(x.__repr__()+"\n" for x in self.chains.values())
        return(message)

    def add_chain(self,chain):
        self.chains.update({chain.name : chain})


class PdbParser:
    """
    A class implementing the methods to parse a PDB file to Protein object.

    Attributes
    ----------
    Does not use any class variables, so that the parsing function can be reused without initializing another
    class instance.

    Methods
    -------

    load_local_file:
        Opens a local file and reads it. Returns a list with each entry a string line of the PBD file

    load_rcsb:
        Opens a connection to the RCSB website and searches it for the specified PDB ID. If present, reads it.
        Returns a list with each entry a string line of the PBD file

    raw_to_atoms:
        Parses a list of string lines of a PDB file into a list of Atom objects, one for each ATOM line.
        Returns a list of Atom objects.

    atoms_to_proteins:
        Assembles a list of atoms into Residues, Chains, and Protein objects, preserving order. Returns a Protein object.

    file_to_protein:
        Takes either a local .pdb file or a PDBID, loads it and parses it into a Protein object. Returns a Protein object.

    """
    def load_local_file(self,path_to_file):
        with open(path_to_file) as file:
            f = file.read().splitlines()
        return(f)

    def load_rcsb(self,pdb_id): ## double check this
        get_request = "https://files.rcsb.org/download/"

        try: ## Try to load the pdb file from the RCSB website
            request = urllib.request.urlopen(get_request+pdb_id+".pdb")
            print(request)
        except urllib.error.HTTPError: ## If there's an error, raise an exception
            raise Exception("PDB ID could not be found in RCSB database.")

        ## File is recevied in binary, decode to default UTF-8.
        file = request.read().decode().splitlines()
        return(file)

    def raw_to_atoms(self,list_of_strings):
        list_of_atoms = []
        for string in list_of_strings: ## For every line in file
            if string[0:6] == 'ATOM  ': ## If the line is an Atom line
                list_of_atoms.append( ## Parse the line according to the PDB specification into an Atom object
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
        prot = Protein('name') ## Initialize the protein object with an empty name

        ## Get a unique list of the parsed chains and sort them alpha-numerically
        chains=list(set([x.chain_id for x in list_of_atoms]))
        chains.sort()

        for chain in chains: ## For each chain in that list
            current_chain = Chain() ## Initialize a new chain object

            ## Determine all the atoms that belong to that chain
            atoms_in_chain = [x for x in list_of_atoms if x.chain_id == chain]

            ## Get a unique list of residues in that chain and sort them alpha-numerically
            residues = list(set([x.res_id for x in atoms_in_chain]))
            residues.sort()

            for residue in residues: ## For each residue in that list
                current_residue = Residue() ## Initialize a new chain object

                ## Determine all the atoms that belong to that residue
                atoms_in_residue = [x for x in atoms_in_chain if x.res_id == residue]

                ## Add those atoms to the current residue
                for atom in atoms_in_residue:
                    current_residue.add_atom(atom)

                ## Add the current residue to the current chain
                current_chain.add_residue(current_residue)

            ## Add the current chain to the protein
            prot.add_chain(current_chain)

        ## Return the protein
        return(prot)

    def file_to_protein(self,path):
        ## If the input path terminates with a .pdb or .PDB format
        if path[-3:]=="pdb" or path[-3:]=="PDB":
            ## Load the file locally
            list_of_strings = self.load_local_file(path)
            name = os.path.basename(path)

        else: ## Else, treat it as a PDBID and load from RCSB
            list_of_strings = self.load_rcsb(path)
            name=path

        ## Transform parsed strings into list of Atom objects
        list_of_atoms = self.raw_to_atoms(list_of_strings)
        ## Assemble list of Atoms into a Protein object
        protein = self.atoms_to_protein(list_of_atoms)
        ## Rename the protein object
        protein.name = name
        return(protein)


def ramachandran_plot(x,y,title='Ramachandran Plot',save=None):
    """Generates a Ramachandran plot from x,y coordinates

    Parameters
    ----------
    x : list of floats
        The x coordinates, generally Phi
    y : list of floats
        The y coordiantes, generally Psi
    title : str
        String to title the figure. Defaults to 'Ramachandran Plot'
    save : str
        A string path indicating where to save the generated plot.

    Returns a matplotlib.pyplot object containing the Ramachandran plot
    """
    plt.figure(figsize=(6,6))
    plt.axis([-190,190,-190,190])
    plt.xlabel('$\phi$')
    plt.ylabel('$\psi$')
    plt.grid(linestyle='--')
    plt.xticks([-180,-90,0,90,180])
    plt.yticks([-180,-90,0,90,180])
    plt.suptitle(title)
    plt.plot(x,y,color='blue',marker='.',ls='',markersize=5)

    ## If a save-path is provided, save the figure at that path.
    if save != None:
        plt.savefig(save)
    return(plt)



if __name__=="__main__" :

    ## Initialize the argument parser
    argument_parser=argparse.ArgumentParser()

    ## Add argument for parsing pdb files
    argument_parser.add_argument("file",type=str,nargs="+",metavar='file.pdb or PDBID',
                                 help="The PDB(s) to be analyzed. Either local .pdb files, or PDB IDs searchable on the RCSB Protein Data Bank (http://rcsb.org)")

    ## Add argument for indicating whether to plot interesting residues
    argument_parser.add_argument("-i","--interesting",dest='plot_interesting',action='store_true',
                                 help="Controls whether to plot interesting residues (pre-prolines, glycines) in separate Ramachandran plots")

    ## Add argument for where to output file
    argument_parser.add_argument("-o","--output",
                             default=os.getcwd()+'/',type=str,dest='output_path',
                                 help="Specify the output file or directory for the generated Ramachandran plot(s)")

    ## Parse arguments
    args = argument_parser.parse_args()

    ## If first argument to "file" is a directory, then overwrite files with
    if os.path.isdir(args.file[0])==True:
        args.file=[f for f in os.listdir(args.file[0]) if f[-3:]=='.pdb']

    ## Initialize the filenames for the generated Ramachandran Plots
    base_outfile_name = 'ramachandran_plot.png'
    outfile_names=['glycine',
                   'preproline',
                   'others']

    ## Determine the output file name
    ## If the indicated output file is a directory:
    if os.path.isdir(args.output_path)==True:
        ## Construct full path using outfile_names.
        outfile_names=[args.output_path+x+'_'+base_outfile_name for x in outfile_names]
        base_outfile_name= args.output_path + base_outfile_name

    else: ## assume the specified path is to a file
        ## pre-append the outfile_names to the specified path
        outfile_names=[os.path.dirname(args.output_path)+'/'+x+'_'+os.path.basename(args.output_path) for x in outfile_names]
        base_outfile_name= args.output_path



    ## Initialize lists to store Phi,Psi values
    phi=[]
    psi=[]
    glycine=[]
    preproline=[]

    ## Initialize the parser object
    pdb_parser = PdbParser()
    ## Parse the files into protein structure objects
    proteins =[pdb_parser.file_to_protein(x) for x in args.file]

    ## Iterate over the individual proteins and calculate Psi and Phi angles.
    for protein in proteins: ## For each protein
        for chain in protein.chains.values(): ## For each chain in the protein
            for i in range(1,len(chain.residues)-1):  ## Starting with the 2nd residue, terminating at 2nd to last residue
                ## Calculate phi, using the dihedral angle of Ci–1 Ni Ciα Ci and store it in Phi
                phi.append(Vector3D.calculate_dihedral_angle(chain.residues[i-1].C_beta.coord,
                                                        chain.residues[i].N.coord,
                                                        chain.residues[i].C_alpha.coord,
                                                        chain.residues[i].C_beta.coord
                                                             ))

                ## Calculate psi, using the dihedral angle of Ni Ciα Ci Ni+1 and store it in Psi
                psi.append(Vector3D.calculate_dihedral_angle(chain.residues[i].N.coord,
                                                      chain.residues[i].C_alpha.coord,
                                                      chain.residues[i].C_beta.coord,
                                                    chain.residues[i+1].N.coord
                                                             ))
                ## If plotting 'interesting' residues
                if(args.plot_interesting==True):

                    ## Store a boolean indicating a residue is Gly or Pre-proline
                    if chain.residues[i].name == 'GLY':
                        glycine.append(True)
                    elif chain.residues[i+1].name == 'PRO':
                        preproline.append(True)
                    else:
                        glycine.append(False)
                        preproline.append(False)



    ## Plot Ramachandran Plot

    ## Add additionally identifying information if the number of proteins == 1, or a number of > 1
    if len(proteins) == 1:
        name_add = " of " + proteins[0].name
    else:
        name_add = " of " + str(len(proteins)) +" PDB structures"

    ## Generate the plots as required
    if(args.plot_interesting==True):
        ## For interesting residues, use boolean lists to index out Psi and Phi values belonging to:

        ## Glycine
        plt_glycine = ramachandran_plot([x for x,y in zip(phi,glycine) if y == True],
                                        [x for x,y in zip(psi,glycine) if y == True],
                                        'Glycine Residues'+name_add,save=outfile_names[0])
        ## Pre-proline
        plt_preproline = ramachandran_plot([x for x,y in zip(phi,preproline) if y == True],
                                           [x for x,y in zip(psi,preproline) if y == True],
                                          'Pre-proline Residues'+name_add,save=outfile_names[1])

        ## Other residues
        plt_normal = ramachandran_plot([x for x,y,z in zip(phi,preproline,glycine)
                                        if y == False and z == False],
                                       [x for x,y,z in zip(psi,preproline,glycine)
                                        if y == False and z == False],
                                        'Other residues'+name_add,save=outfile_names[2])

    ## If generating an 'interesting' plot is not specified, generate a normal plot
    else:
        plt = ramachandran_plot(phi,psi,'Ramachandran Plot'+name_add,save=base_outfile_name)
