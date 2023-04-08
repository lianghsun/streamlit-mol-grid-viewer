# Molecule Grid Viewer
[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://lianghsun-streamlit-mol-grid-viewer-main-hpf1jn.streamlit.app/)

*Say goodbye to the native RDKit's [MolsToGridImage](https://www.rdkit.org/docs/source/rdkit.Chem.Draw.html#rdkit.Chem.Draw.MolsToGridImage).* Have you ever used RDKit's built-in grid to display molecules, but when it comes to displaying label contents, it's very troublesome: _you can't segment them, or they may look unattractive._

Streamlit's built-in `st.columns` and `st.rows` allow you to easily create a grid viewer, and then you can assemble any content you want. In this repo, users can upload sdf files and choose the content they want to appear in the label, creating a clear and understandable grid image. Have fun!
# LICENSE
MIT