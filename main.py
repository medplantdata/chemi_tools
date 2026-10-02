from rdkit import Chem
import streamlit as st
import nispo
from pyopsin import PyOpsin
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Draw
from rdkit import TautomerEnumerator


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

def compound_comparitor(input1,input2,type):
    if type == 'smiles':
        input1 = input1.strip()
        input2 = input2.strip()
        # Checks smiles strings are exactly identical
        if input1 == input2:
            return 'Both smiles are the same'
        m1 = Chem.MolFromSmiles(input1)
        m2 = Chem.MolFromSmiles(input2)
        # Checks mols made from smiles are substructs of each other ie are the same
        if m1.HasSubstructMatch(m2) and m2.HasSubstructMatch(m1):
            return 'Both smiles resolve to the same molecule'
        # Checks if mols are stereoisomers of each other (2D smiles are the same)
        non_isomeric1 = Chem.MolToSmiles(m1)
        non_isomeric2 = Chem.MolToSmiles(m2)
        if non_isomeric1 == non_isomeric2:
            return 'Smiles are stereoisomers of each other'
        #check tautomerism
        t1 = TautomerEnumerator().Canonicalize(m1)
        t2 = TautomerEnumerator().Canonicalize(m2)
        ts1 = Chem.MolFromSmiles(t1)
        ts2 = Chem.MolFromSmiles(t2)
        if ts1 == ts2:
            return 'Compounds are tautomers of each other'
        return 'Smiles refer to different compounds'        
    

st.title("Calitz Chemical Converter")

st.set_page_config(layout='wide')

left, right = st.columns(2)

with left:
     
    input_stringl = st.text_input('Input SMILES, InChi or iupac name', key='input_string_l', placeholder='CN1CCC[C@H]1C2=CN=CC=C2')

    if input_stringl:
        mol, typel = detect_input_type(input_stringl)
        if mol:
            c_smiles = Chem.MolToSmiles(mol,isomericSmiles=False,canonical=True)
            i_smiles = Chem.MolToSmiles(mol, isomericSmiles=True,canonical=True)
            InChi  = Chem.MolToInchi(mol)
            inchikey = Chem.MolToInchiKey(mol)
            iupac = nispo.mol_to_iupac(mol)
            mol_formula = Chem.rdMolDescriptors.CalcMolFormula(mol)

            st.success(f'We detected input as {typel}')

            img = Chem.Draw.MolToImage(mol, size = (500,300))
            st.image(image=img, caption='Molecular Structure')

            st.caption('SMILES (RDkit canonical 2D)')
            st.code(c_smiles)

            st.caption('SMILES (RDkit isomeric canonical 3D)')
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

    input_stringr = st.text_input('Input SMILES, InChi or iupac name', key='input_string_r', placeholder='CCO')

    if input_stringr:
        mol, typer = detect_input_type(input_stringr)
        if mol:
            c_smiles = Chem.MolToSmiles(mol, isomericSmiles=False,canonical=True)
            i_smiles = Chem.MolToSmiles(mol, isomericSmiles=True,canonical=True)
            InChi  = Chem.MolToInchi(mol)
            inchikey = Chem.MolToInchiKey(mol)
            iupac = nispo.mol_to_iupac(mol)
            mol_formula = Chem.rdMolDescriptors.CalcMolFormula(mol)

            st.success(f'We detected input as {typer}')

            img = Chem.Draw.MolToImage(mol, size = (500,300))
            st.image(image=img, caption='Molecular Structure')

            st.caption('SMILES (RDkit canonical 2D)')
            st.code(c_smiles)

            st.caption('SMILES (RDkit isomeric canonical)')
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



if input_stringl and input_stringr:
    if typel == 'smiles' and typer == 'smiles':
        st.sucess = compound_comparitor(input_stringl,input_stringr,'smiles')