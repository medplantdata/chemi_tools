from rdkit import Chem
import streamlit as st
import nispo
from pyopsin import PyOpsin
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Draw


opsin = PyOpsin()
smiles = opsin.to_smiles('ethanol')
print(smiles)

mol = Chem.MolFromSmiles(smiles[0])
iupac = nispo.mol_to_iupac(mol)
print(iupac)

def detect_input_type(input_string):
    input_string = input_string.strip()
    mol = Chem.MolFromSmiles(input_string)
    if mol:
        return mol, 'smiles'
    if not mol:
        mol = Chem.MolFromInchi(input_string)
        if mol:
            return mol, 'inchi'
    
    try:
        opsin = PyOpsin()
        smiles = opsin.to_smiles(input_string)
        mol = Chem.MolFromSmiles(smiles[0])
        if mol:
            return mol, 'iupac'
    except:
        return None, 'none'



st.title("Calitz Chemical Converter")

st.set_page_config(layout='wide')

left, right = st.columns(2)

with left:

    input_string = st.text_input('Input SMILES, InChi or iupac name', key='input_string_l', placeholder='CCCCCC1=CC(=C2[C@@H]3C=C(CC[C@H]3C(OC2=C1)(C)C)C)O')

    if input_string:
        mol, type = detect_input_type(input_string)
        if mol:
            c_smiles = Chem.MolToSmiles(mol,isomericSmiles=False,canonical=True)
            i_smiles = Chem.MolToSmiles(mol, isomericSmiles=True,canonical=True)
            InChi  = Chem.MolToInchi(mol)
            inchikey = Chem.MolToInchiKey(mol)
            iupac = nispo.mol_to_iupac(mol)
            mol_formula = Chem.rdMolDescriptors.CalcMolFormula(mol)

            st.success(f'We detected input as {type}')

            img = Chem.Draw.MolToImage(mol, size = (500,300))
            st.image(image=img, caption='Molecular Structure')

            st.caption('SMILES (RDkit canonical 2D))')
            st.code(c_smiles)

            st.caption('SMILES (RDkit isomeric canonical 3D))')
            st.code(i_smiles)

            st.caption('InChi')
            st.code(InChi)

            st.caption('InChiKey')
            st.code(inchikey)

            st.caption('Iupac Name')
            st.code(iupac)

            st.caption('Molecular Formula')
            st.code(mol_formula)
        else:
            st.error('Unfortunately we cannot parse this input')

    else:
        st.error('Please enter string')

with right:

    input_string = st.text_input('Input SMILES, InChi or iupac name', key='input_string_r', placeholder='CCO')

    if input_string:
        mol, type = detect_input_type(input_string)
        if mol:
            c_smiles = Chem.MolToSmiles(mol,canonical=True)
            i_smiles = Chem.MolToSmiles(mol, isomericSmiles=True,canonical=True)
            InChi  = Chem.MolToInchi(mol)
            inchikey = Chem.MolToInchiKey(mol)
            iupac = nispo.mol_to_iupac(mol)
            mol_formula = Chem.rdMolDescriptors.CalcMolFormula(mol)

            st.success(f'We detected input as {type}')

            img = Chem.Draw.MolToImage(mol, size = (500,300))
            st.image(image=img, caption='Molecular Structure')

            st.caption('SMILES (RDkit canonical))')
            st.code(c_smiles)

            st.caption('SMILES (RDkit isomeric canonical))')
            st.code(i_smiles)

            st.caption('InChi')
            st.code(InChi)

            st.caption('InChiKey')
            st.code(inchikey)

            st.caption('Iupac Name')
            st.code(iupac)

            st.caption('Molecular Formula')
            st.code(mol_formula)
        else:
            st.error('Unfortunately we cannot parse this input')

    else:
        st.error('Please enter string')



