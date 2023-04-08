import pandas as pd
import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw, PandasTools

st.set_page_config(page_title='Molecule Grid Viewer', page_icon='👀')

@st.cache_data
def read_sdf(sdf):
    return PandasTools.LoadSDF(sdf)

# Function to display molecules in a grid format with 5 molecules per row
def display_molecules(frame):
    num_molecules = frame.__len__()
    row_size = 5
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

    uploaded_file = st.file_uploader("Upload an SDF file", type="sdf")

    if uploaded_file:
        frame = read_sdf(uploaded_file)
        
        st.multiselect(
            label='Select which properties you would like to display.',
            options=frame.columns,
            default='SMILES',
            key='props'
        )
        
        if st.button(f'Display __{frame.__len__()}__ mols'):
            display_molecules(frame)

if __name__ == "__main__":
    main()
