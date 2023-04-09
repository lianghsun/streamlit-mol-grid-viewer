import pandas as pd
import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw, PandasTools

st.set_page_config(page_title='Molecule Grid Viewer', page_icon='👀')

def read_sdf(sdf):
    return PandasTools.LoadSDF(sdf)

# Function to display molecules in a grid format with 5 molecules per row
def display_molecules(frame):
    num_molecules = frame.__len__()
    row_size = st.session_state['row_size']
    rows = (num_molecules + row_size - 1) // row_size

    for row in range(rows):
        mols = frame[row * row_size:(row + 1) * row_size].copy()
        cols = st.columns(row_size)
        for col, (_, mol) in zip(cols, mols.iterrows()):
            img = Draw.MolToImage(mol['ROMol'], size=(200, 200))
            col.image(img, use_column_width=True)
            for prop in st.session_state['props']:
                col.write(f"__{prop}__: {mol[prop]}")

# Main app
def main():
    st.title("Molecule Grid Viewer")
    
    with st.expander('🤔 __What\'s this for?__'):
        st.markdown("Have you ever used RDKit's `MolsToGridImage()` function to draw molecular structures and tried to label relevant information underneath each molecule, like the example image below?")
        example_sms = pd.DataFrame(
            {
                'SMILES': ['O=C1NOC[C@@H]1/N=C/c1ccc(/C=N/[C@@H]2CONC2=O)cc1','CCO[C@H]1[C@@H](O)[C@H](CO)O[C@@H]1n1cnc2c(NC)nc(C)nc21','Cc1cc[n+]([O-])c(CS(=O)(=O)c2nnnn2-c2ccccc2)c1'],
                'MW': [302.29,323.353,331.357]
            }
        )
        PandasTools.AddMoleculeColumnToFrame(example_sms, 'SMILES')
        grid_image = Chem.Draw.MolsToGridImage(
            example_sms['ROMol'].to_list(),
            legends=[f"SMILES: {row['SMILES']}, MW: {row['MW']}" for idx, row in example_sms[['SMILES', 'MW']].iterrows()],
            molsPerRow=3)
        st.image(grid_image, use_column_width=True)
        st.markdown("_Doesn't this matter seem a bit off, right?_ Labels and such are all **crammed together 💥**, making it hard to read.")
        
    st.radio(
        label='See how this works',
        options=['Example molecules', 'Upload your SDF'],
        horizontal=True,
        key='method'
    )
    
    if st.session_state['method'] == 'Example molecules':         
        st.multiselect(
                label='Select which properties you would like to display.',
                options=example_sms.columns,
                default='SMILES',
                key='props'
        )
        st.number_input("Number of molecules per row", min_value=3, max_value=5, value=3, step=1, key='row_size')
        display_molecules(example_sms)
        
    elif st.session_state['method'] == 'Upload your SDF':
        uploaded_file = st.file_uploader("Upload an SDF file", type="sdf")

        if uploaded_file:
            frame = read_sdf(uploaded_file)

            st.multiselect(
                label='Select which properties you would like to display.',
                options=frame.columns,
                default='SMILES',
                key='props'
            )
            
            st.number_input("Number of molecules per row", min_value=3, max_value=5, value=5, step=1, key='row_size')

            if st.button(f'Display __{frame.__len__()}__ mols'):
                display_molecules(frame)

if __name__ == "__main__":
    main()
