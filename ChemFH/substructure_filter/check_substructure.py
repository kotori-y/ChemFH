'''
Description: the functions
Author: Kotori Y
Date: 2020-10-24 16:08:49
LastEditors: Kotori Y
LastEditTime: 2020-10-26 14:37:28
FilePath: \ChemFH\ChemFH\substructure_filter\check_substructure.py
AuthorMail: kotori@cbdd.me
'''

from collections import namedtuple
from functools import wraps
try:
    from load_pattern import loadpkl
except Exception:
    from .load_pattern import loadpkl


class InvalidInputError(Exception):
    pass


def withEndpoint(endpoint="PAINS"):
    def check(func):
        @wraps(func)
        def wrapper(mol, endpoint):
            pattl = loadpkl(endpoint)
            res = func(mol=mol, pattl=pattl)
            return res
        return wrapper
    return check


def checkValidMol(func):
    @wraps(func)
    def wrapper(mol, **kwgrs):
        if mol is not None:
            return func(mol, **kwgrs)
        else:
            return {}
    return wrapper

        


@withEndpoint()
@checkValidMol
def checkPattl(mol, pattl):
    
    checkRes = {}
    for name, patt in pattl:
        atoms = mol.GetSubstructMatches(patt)
        if atoms:
            checkRes[name] = atoms
    
    return checkRes



if '__main__' == __name__:
    from rdkit import Chem

    smi = 'O=c1cc(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12' 
    mol = Chem.MolFromSmiles(smi)

    print(checkPattl(mol=mol, endpoint="PAINS"))


        
        
