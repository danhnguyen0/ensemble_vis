import os, json
import streamlit as st
import pandas as pd
import plotly.express as px
import py3Dmol

st.set_page_config(layout='wide')

# Initialize session state if it doesn't exist
if 'uploaded_file' not in st.session_state:
    st.session_state.uploaded_file = None
    st.session_state.df = None
    st.session_state.selected_chart = None
    st.session_state.selected_boxplot_column = None
    st.session_state.selected_scatterplot_x_column = None
    st.session_state.selected_scatterplot_y_column = None
    st.session_state.selected_barchart_x_column = None
    st.session_state.selected_barchart_y_column = None
    st.session_state.selected_model = None

def show_summary(df):
    st.subheader('Data Summary')
    summary = df.describe().transpose()
    st.write(summary)

def create_boxplot(df, column):
    fig = px.box(df, y=column, points='all')
    st.plotly_chart(fig, use_container_width=True, key='boxplot', config={'displayModeBar': True})

def create_scatterplot(df, x_column, y_column):
    fig = px.scatter(df, x=x_column, y=y_column)
    st.plotly_chart(fig, use_container_width=True)

def create_barchart(df, x_column, y_column):
    fig = px.bar(df, x=x_column, y=y_column)
    st.plotly_chart(fig, use_container_width=True)

def load_pdb_and_show_3d(pdb_file_path, style='stick', color='blue'):
    try:
        with open(pdb_file_path, 'r') as pdb_file:
            pdb_data = pdb_file.read()
    except Exception as e:
        st.error(f'Error loading PDB file: {e}')
        return

    viewer = py3Dmol.view(width=800, height=600)
    viewer.addModel(pdb_data, 'pdb')

    viewer.setStyle(
        {
            'stick': {},  # Represent atoms as sticks
            'cartoon': {},  # Represent secondary structures (alpha helix, beta sheets)
            'sphere': {'radius': 0.3, 'color': 'lightblue'},  # Show small spheres for atoms
            'surface': {'color': 'yellow'},  # Display molecular surfaces
            'ligand': {'stick': {'colorscheme': 'green'}}  # Display ligands in a green stick model
        }
    )
    viewer.zoomTo()
    viewer.setBackgroundColor('white')  # Set background color to white for better contrast

    viewer_html = viewer._make_html()
    st.components.v1.html(viewer_html, height=600)
    
def compare_model_parameters(df):
    st.subheader("Compare Model Parameters (Dmax, Gyration)")
    
    if not {'model', 'dmax', 'gyration'}.issubset(df.columns):
        st.error("Dataset must contain 'model', 'dmax', and 'gyration' columns.")
        return

    # Allow the user to select models for comparison
    selected_models = st.multiselect('Select models to compare:', df['model'].unique(), key='parameter_comparison_selector')
    
    if selected_models:
        filtered_df = df[df['model'].isin(selected_models)][['model', 'dmax', 'gyration']].reset_index(drop=True)
        
        st.write("### Parameter Table")
        st.dataframe(filtered_df)
        
        st.write("### Parameter Comparison Chart")
        fig = px.bar(
            filtered_df.melt(id_vars='model', value_vars=['dmax', 'gyration']),
            x='model', 
            y='value', 
            color='variable', 
            barmode='group',
            labels={'value': 'Parameter Value', 'variable': 'Parameter Type'},
            title="Dmax and Gyration Comparison"
        )
        st.plotly_chart(fig, use_container_width=True)
    
def show_3d_models_comparison(df):
    st.subheader('Compare Protein Structures')

    if not {'model', 'dmax', 'gyration'}.issubset(df.columns):
        st.error('Dataset must contain "model", "dmax", and "gyration" columns.')
        return

    models = df['model'].dropna().unique()
    selected_models = st.multiselect('Select protein models to compare:', models, key='model_comparison_selector')

    if selected_models:
        for i in range(0, len(selected_models), 3):
            cols = st.columns(3)  # Create 3 columns per row
            for j, selected_model in enumerate(selected_models[i:i + 3]):
                with cols[j]:
                    model_data = df[df['model'] == selected_model].iloc[0]
                    dmax = model_data['dmax']
                    gyration = model_data['gyration']

                    st.markdown(f"### Model: {selected_model}")
                    st.markdown(f"- **Dmax:** {dmax}")
                    st.markdown(f"- **Gyration:** {gyration}")

                    pdb_file_path = f'dataset/output_e00{int(st.session_state.uploaded_file.name.split(".")[0][-1])}/model_{selected_model}.pdb'

                    if os.path.exists(pdb_file_path):
                        load_pdb_and_show_3d(pdb_file_path)
                    else:
                        st.error(f'PDB file for {selected_model} not found at the expected location.')




def show_3d_model_from_selection(df):
    st.subheader('Individual Protein Structure')

    if 'model' not in df.columns:
        st.error('No model column in the dataset.')
        return

    models = df['model'].dropna().unique()
    selected_model = st.selectbox('Select a protein model to view:', models, key='model_selector')

    if selected_model:
        st.session_state.selected_model = selected_model

        pdb_file_path = f'dataset/output_e00{int(uploaded_file.name.split(".")[0][-1])}/model_{selected_model}.pdb'
        print(pdb_file_path)

        if os.path.exists(pdb_file_path):
            load_pdb_and_show_3d(pdb_file_path)
        else:
            st.error(f'PDB file for {selected_model} not found at the expected location.')

def show_3d_model_ensemble():
    st.subheader('Ensemble Protein Structure')

    pdb_file_path = f'dataset/ensembles/PED00020e00{int(uploaded_file.name.split(".")[0][-1])}.pdb'

    if os.path.exists(pdb_file_path):
        load_pdb_and_show_3d(pdb_file_path)
    else:
        st.error(f'PDB file for {pdb_file_path} not found at the expected location.')

st.title('Advanced Protein Structure Visualization')

# Upload CSV or JSON file
uploaded_file = st.file_uploader('Upload a CSV/JSON file', type=['csv', 'json'], accept_multiple_files=False)

if uploaded_file is not None:
    st.session_state.uploaded_file = uploaded_file

    if uploaded_file.name.endswith('.csv'):
        df = pd.read_csv(uploaded_file)
        st.session_state.df = df

    elif uploaded_file.name.endswith('.json'):
        try:
            file_content = uploaded_file.read().decode('utf-8')
            data = json.loads(file_content)

            if 'conformations' in data:
                df = pd.DataFrame(data['conformations'])
                st.session_state.df = df
            else:
                st.error('JSON does not contain the expected conformations key.')
        except Exception as e:
            st.error(f'Error processing file: {e}')
    else:
        st.error('Invalid file format! Please upload a CSV or JSON file.')

    # Show data summary
    show_summary(st.session_state.df)
    
    # Tabs for individual and comparison views
    tab1, tab2, tab3 = st.tabs(["Individual Protein View", "Compare Proteins", "Compare Parameters"])

    # Individual model view
    with tab1:
        show_3d_model_from_selection(st.session_state.df)

    # Comparison of protein structures
    with tab2:
        show_3d_models_comparison(st.session_state.df)
    
    with tab3:
        compare_model_parameters(st.session_state.df)
    

    # Select columns for plots
    st.subheader('Select Columns for Plots')
    x_column = st.selectbox(
        'Select X Column', 
        st.session_state.df.columns, 
        key='x_column',
        index=st.session_state.df.columns.get_loc(st.session_state.selected_scatterplot_x_column) 
            if st.session_state.selected_scatterplot_x_column in st.session_state.df.columns else 0
    )
    y_column = st.selectbox(
        'Select Y Column', 
        st.session_state.df.columns, 
        key='y_column',
        index=st.session_state.df.columns.get_loc(st.session_state.selected_scatterplot_y_column) 
            if st.session_state.selected_scatterplot_y_column in st.session_state.df.columns else 1
    )
    st.session_state.selected_scatterplot_x_column = x_column
    st.session_state.selected_scatterplot_y_column = y_column

    # Boxplot
    colss = st.columns(3)
    with colss[0]:
        st.subheader('Boxplot')
        create_boxplot(st.session_state.df, y_column)

    with colss[1]:
        # Scatterplot
        st.subheader('Scatterplot')
        create_scatterplot(st.session_state.df, x_column, y_column)
    with colss[2]:
        # Barchart
        st.subheader('Barchart')
        create_barchart(st.session_state.df, x_column, y_column)
        #show_3d_model_ensemble()
    
else:
    st.info('Please upload a CSV or JSON file to see the visualizations.')