'''
Description: 
Author: Kotori Y
Date: 2020-10-24 16:02:15
LastEditors: Kotori Y
LastEditTime: 2020-10-26 14:41:40
FilePath: \ChemFH\ChemFH\__init__.py
AuthorMail: kotori@cbdd.me
'''

from substructure_filter import FrequentHitterFilter, highlightAtoms
from collections.abc import Iterable
import pandas as pd
from rdkit import Chem


class ChemFH(FrequentHitterFilter):

    def __init__(self, mols, endpoints, nJobs=1):

        self.mols = mols if isinstance(mols, Iterable) else (mols, )
        self.endpoints = endpoints if isinstance(
            endpoints, Iterable) else (endpoints, )
        self.nJobs = nJobs

    def _getMolInfo(self, mol):

        if mol:
            molSVG = highlightAtoms(mol, [])
            smiles = Chem.MolToSmiles(mol)
        else:
            molSVG, smiles = ['Invalid']*2
        return molSVG, smiles

    def _disposedRes(self, res, endpoint):
        num = []
        svg = []

        for mol, x in zip(mols, res):
            num.append(len(x))
            svg.append(
                ''.join(
                    [highlightAtoms(mol, atoms[0], [200, 100])
                     if x else 'Accepted' for atoms in x.values()]
                )
            )

            
        summary = pd.DataFrame(
            {
                f"{endpoint}": num,
            }
        )
        
        detail = pd.DataFrame(
            {
                f"{endpoint}": svg,
            }
        )
        yield summary, detail

    def getScreeningResult(self):

        smis = []
        molSVGs = []
        
        for mol in self.mols:
            molSVG, smiles = self._getMolInfo(mol)
            molSVGs.append(molSVG)
            smis.append(smiles)
            
        out1 = pd.DataFrame({"Molecule": molSVGs, "SMILES":smis})
        out2 = pd.DataFrame({"Molecule": molSVGs, "SMILES":smis})

        for endpoint in self.endpoints:
            res = super().screening(self.mols, endpoint=endpoint)
            for summary, detail in self._disposedRes(res, endpoint):
                out1 = pd.concat([out1, summary], axis=1)
                out2 = pd.concat([out2, detail], axis=1)

        pd.set_option('colheader_justify', 'center')   # FOR TABLE <th>
        html_string = """<html>
        <head><title>HTML Pandas Dataframe with CSS</title></head>
        <style>
        /* includes alternating gray and white with on-hover color */

        .mystyle {
            font-size: 11pt; 
            font-family: Arial;
            border-collapse: collapse; 
            border: 1px solid silver;

        }

        .mystyle td, th {
            padding: 5px;
        }

        .mystyle tr:nth-child(even) {
            background: #E0E0E0;
        }

        .mystyle tr:hover {
            background: silver;
            cursor: pointer;
        }
        </style>
        <body>
            %s
        </body>
        </html>"""

        out1 = html_string % out1.to_html(classes='mystyle', escape=False)
        out2 = html_string % out2.to_html(classes='mystyle', escape=False)

        return out1, out2


def readSmiles(smiles: list):
    
    mols = [Chem.MolFromSmiles(smi) for smi in smiles]
    return mols

def readSDFfile(sdfFile: str):

    suppl = Chem.FastSDMolSupplier(sdfFile)
    mols = [mol for mol in suppl]
    return mols

def main(mols, endpoints=["PAINS"]):

    fh = ChemFH(mols, endpoints=endpoints)
    out1, out2 = fh.getScreeningResult()
    return out1, out2



if '__main__' == __name__:
    from rdkit import Chem
    import pandas as pd

    smis = [
        'OC1=C[C-]2[OH+]C(c3ccc(O)c(O)c3)=C(O)C=C2C(O)=C1',
        'O=c1cc(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12',
        'OCC(O)C(O)C(O)C(O)CO',
        'test', # for testing invalid input
    ]
    mols = readSmiles(smis)
    out1, out2 = main(mols, endpoints=["PAINS", "ALARM_NMR"])
