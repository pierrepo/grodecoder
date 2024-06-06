import os
from pathlib import Path
import tempfile

from loguru import logger
from PIL import Image
import streamlit as st

import grodecoder as gd


def last_line_log_file():
    with open("test.log", "rb") as file:
        try:
            file.seek(-2, os.SEEK_END)
            while file.read(1) != b'\n':
                file.seek(-2, os.SEEK_CUR)
        except OSError:
            file.seek(0)
        last_line = file.readline().decode()
    return last_line



if __name__ == "__main__":
    logger.remove()
    logger.add(st.write, format="{time:YYYY-MM-DD HH:mm:ss} | {message}", level="INFO")
    logger.add("test.log", format="{time:YYYY-MM-DD HH:mm:ss} | {level} | {message}")
    st.markdown(
        """
        <style>
        p {
            margin-bottom: 0px;
        }
        </style>
        """, unsafe_allow_html=True)
    logo = Image.open("assets/grodecoder_logo.png")
    st.sidebar.image(logo)
    st.sidebar.markdown("""
    **GroDecoder** extracts and identifies
    the molecular components of a structure file (PDB or GRO)
    issued from a molecular dynamics simulation.

    ---
    [Source code](https://github.com/pierrepo/grodecoder)                 
    """
    )

    st.markdown("# GroDecoder 📦")
    

    uploaded_file = st.file_uploader("Choose a structure file", type=["pdb", "gro"])

    st.markdown("""
    Examples:
    <a href="https://raw.githubusercontent.com/pierrepo/grodecoder/main/data/examples/barstar.gro" download target="_blank">Barstar</a>
    <br /><br />
    """, unsafe_allow_html=True)

    st.write("Options:")
    check_overlapping_residue_value = st.checkbox("Check overlapping residue", value=False)
    check_connectivity_value = st.checkbox("Check connectivity (add degree and number of edge in the fingerprint) ", value=False)
    bond_threshold_value = st.text_input("Threshold value ('auto' or any positive value):", value="auto")

    if (uploaded_file is not None) and st.button("Run analysis") : 
        temp_dir = tempfile.mkdtemp()
        path = Path(temp_dir) / uploaded_file.name
        with open(path, "wb") as f:
            f.write(uploaded_file.getvalue())
        st.write(f"Saved to {path}")
        with st.spinner("Running analysis..."):
            gd.main(
                path,
                draw_graph_option=False,
                check_overlapping_residue=check_overlapping_residue_value,
                check_connectivity=check_connectivity_value,
                bond_threshold=bond_threshold_value,
                query_pdb=True,
            )

        # Fetch the last line of the log file
        last_line = last_line_log_file()
        last_line_filename = last_line.split()[-1]

        # Test if the last line of the log file is a JSON
        # That also mean that the analysis terminated correctly
        if last_line_filename.split('.')[-1] == "json":
            with open(last_line_filename, 'r') as file:
                file_content = file.read()
            st.toast("Analysis completed successfully 🎉", icon="🥳")
            st.download_button(
                label="Download molecular inventory",
                data=file_content,
                file_name=last_line_filename,
                mime="text/plain"
            )
