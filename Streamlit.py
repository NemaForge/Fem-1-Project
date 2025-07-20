import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import os
import numpy as np
import glob # Needed for dynamic file discovery

# --- Page Configuration (Should be at the very top of your script) ---
st.set_page_config(
    page_title="Saurish and Xander's Biomart",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="collapsed"
)

# --- Global Data Loading and Constants ---
# Assuming AnalysisFile2.txt is in the root directory
ORIGINAL_DATA_FILE = "AnalysisFile2.txt"

# Original Germ Cell data files (assuming they are in the root directory)
GERM_CELL_FILES_MAP = {
    "Mature sperm": "finalMatureSperm.txt",
    "Meiotic germ cells": "finalMeiotic.txt",
    "Mitotic germ cells": "finalMitotic.txt",
    "Oocytes": "finalOocytes.txt",
    "Spermatids": "finalSpermatids.txt",
    "Spermatocytes": "finalSpermatocytes.txt",
    "Syncitial pachytene spermatocytes": "finalSyncitial.txt"
}

# Directory for the new processed somatic data
# IMPORTANT: This path is relative to your Streamlit app's root.
# Make sure your 'Processed Somatic Data Set_Grouped' folder is inside a 'data' folder
# within your GitHub repository (e.g., your_repo/data/Processed Somatic Data Set_Grouped/).
SOMATIC_CELL_DIRECTORY = "data/Processed Somatic Data Set_Grouped"

regular_dot_size = 5
fem1_dot_size = 10

@st.cache_data
def load_original_data(path):
    """Loads and preprocesses the original AnalysisFile2.txt data."""
    try:
        df_loaded = pd.read_csv(path, sep='\t')
        for col in ['Mean of Geneid Strains', 'Standard Deviation of Geneid Strains', 'Group']:
            df_loaded[col] = pd.to_numeric(df_loaded[col], errors='coerce')
        df_loaded.dropna(subset=['Mean of Geneid Strains', 'Standard Deviation of Geneid Strains', 'Group'], inplace=True)
        df_loaded['Group'] = df_loaded['Group'].astype(str)
        return df_loaded
    except FileNotFoundError:
        st.error(f"Error: Original data file not found at {path}. Please ensure the file exists in your repository.")
        st.stop()
    except Exception as e:
        st.error(f"Error loading original data: {e}")
        st.stop()

@st.cache_data
def load_all_single_cell_dataframes():
    """
    Loads both Germ and Somatic single-cell data, standardizing column names.
    Returns a nested dictionary: {'Germ': {df_name: df}, 'Somatic': {df_name: df}}.
    """
    all_single_cell_dfs = {'Germ': {}, 'Somatic': {}}
    
    # --- Load Germ Cell Data ---
    for display_name, filename in GERM_CELL_FILES_MAP.items():
        file_path = filename # Germ files assumed in root
        try:
            df_sc = pd.read_csv(file_path, sep='\t')
            
            # Standardize column names for Germ data
            # Assuming original germ files have 'gene name', 'Scaled_TPM', 'group number'
            if 'gene name' in df_sc.columns and 'Scaled_TPM' in df_sc.columns and 'group number' in df_sc.columns:
                df_sc = df_sc.rename(columns={
                    'gene name': 'gene_common',
                    'Scaled_TPM': 'expression_value',
                    'group number': 'group_identifier'
                })
                # Ensure numeric types
                df_sc['expression_value'] = pd.to_numeric(df_sc['expression_value'], errors='coerce')
                df_sc['group_identifier'] = pd.to_numeric(df_sc['group_identifier'], errors='coerce')
                df_sc.dropna(subset=['gene_common', 'expression_value', 'group_identifier'], inplace=True)
                df_sc['group_identifier'] = df_sc['group_identifier'].astype(int).astype(str) # Convert to string for grouping
                all_single_cell_dfs['Germ'][display_name] = df_sc
            else:
                st.warning(f"Germ cell file '{filename}' missing expected columns ('gene name', 'Scaled_TPM', 'group number'). Skipping.")

        except FileNotFoundError:
            st.warning(f"Germ cell file not found: {filename}. It will not be available. Please ensure the file exists in your repository.")
        except Exception as e:
            st.error(f"Error loading Germ cell data '{filename}': {e}")

    # --- Load Somatic Cell Data ---
    somatic_file_paths = glob.glob(os.path.join(SOMATIC_CELL_DIRECTORY, "*.txt"))
    if not somatic_file_paths:
        st.warning(f"No somatic cell data files found in '{SOMATIC_CELL_DIRECTORY}'. Please ensure the path is correct and files are .txt.")
    
    for file_path in somatic_file_paths:
        display_name = os.path.splitext(os.path.basename(file_path))[0] # Use filename without extension as display name
        try:
            df_sc = pd.read_csv(file_path, sep='\t')
            
            # Standardize column names for Somatic data
            # Assuming somatic files have 'gene_short_name', 'scaled_TPM', 'group_number'
            if 'gene_short_name' in df_sc.columns and 'scaled_TPM' in df_sc.columns and 'group_number' in df_sc.columns:
                df_sc = df_sc.rename(columns={
                    'gene_short_name': 'gene_common',
                    'scaled_TPM': 'expression_value',
                    'group_number': 'group_identifier'
                })
                # Ensure numeric types
                df_sc['expression_value'] = pd.to_numeric(df_sc['expression_value'], errors='coerce')
                df_sc['group_identifier'] = pd.to_numeric(df_sc['group_identifier'], errors='coerce')
                df_sc.dropna(subset=['gene_common', 'expression_value', 'group_identifier'], inplace=True)
                df_sc['group_identifier'] = df_sc['group_identifier'].astype(int).astype(str) # Convert to string for grouping
                all_single_cell_dfs['Somatic'][display_name] = df_sc
            else:
                st.warning(f"Somatic cell file '{display_name}' missing expected columns ('gene_short_name', 'scaled_TPM', 'group_number'). Skipping.")

        except FileNotFoundError:
            st.warning(f"Somatic cell file not found: {file_path}. It will not be available.")
        except Exception as e:
            st.error(f"Error loading Somatic cell data '{display_name}': {e}")
            
    return all_single_cell_dfs

# Load data globally
df_original = load_original_data(ORIGINAL_DATA_FILE)
# fem1_data_original is based on 'Gene Name' in df_original, no change here
fem1_data_original = df_original[df_original['Gene Name'] == 'fem-1']

# Load all single-cell data, now categorized
single_cell_dataframes_categorized = load_all_single_cell_dataframes()

# Helper function to get sorted unique groups for comparison legend
def get_sorted_groups(df, group_col):
    if df.empty or group_col not in df.columns:
        return []
    try:
        # Convert to numeric first for proper sorting, then back to string
        return sorted(df[group_col].dropna().astype(str).unique(), key=lambda x: int(x) if x.isdigit() else x)
    except:
        # Fallback for non-numeric group names
        return sorted(df[group_col].dropna().astype(str).unique())

# Populate these only if data exists for accurate sorting
early_embryo_groups_sorted = get_sorted_groups(df_original, 'Group')

# Find a representative single-cell dataframe to get sorted groups for legend
# Using the first available germ or somatic dataframe for common single cell groups
sample_sc_df = None
if single_cell_dataframes_categorized['Germ']:
    sample_sc_df = next(iter(single_cell_dataframes_categorized['Germ'].values()))
elif single_cell_dataframes_categorized['Somatic']:
    sample_sc_df = next(iter(single_cell_dataframes_categorized['Somatic'].values()))

if sample_sc_df is not None:
    single_cell_groups_sorted = get_sorted_groups(sample_sc_df, 'group_identifier')
else:
    single_cell_groups_sorted = []


# --- Helper Functions for Plotting ---
def create_aggregated_hover_data_flexible(df_to_process, gene_col, mean_col, std_dev_col, group_col, round_decimals=3):
    if df_to_process.empty:
        return pd.DataFrame(columns=[mean_col, std_dev_col, group_col, 'Aggregated Hover Text'])

    df_temp = df_to_process.copy()
    
    df_temp[f'{mean_col}_numeric_rounded'] = pd.to_numeric(df_temp[mean_col], errors='coerce').round(round_decimals)
    # std_dev_col might not exist for single-cell data, handle gracefully
    if std_dev_col and std_dev_col in df_temp.columns:
        df_temp[f'{std_dev_col}_numeric_rounded'] = pd.to_numeric(df_temp[std_dev_col], errors='coerce').round(round_decimals)
        subset_cols = [f'{mean_col}_numeric_rounded', f'{std_dev_col}_numeric_rounded', group_col]
        dropna_subset_cols = [f'{mean_col}_numeric_rounded', f'{std_dev_col}_numeric_rounded']
    else: # For single-cell data where std_dev is not provided
        df_temp[f'{std_dev_col}_numeric_rounded'] = np.nan # Placeholder if it doesn't exist
        subset_cols = [f'{mean_col}_numeric_rounded', group_col] # Group by mean and group only
        dropna_subset_cols = [f'{mean_col}_numeric_rounded']

    df_temp.dropna(subset=dropna_subset_cols, inplace=True)

    def generate_hover_text(x_group):
        lines = []
        for gene_name_val, mean_val, std_dev_val in zip(x_group[gene_col], x_group[mean_col], x_group[std_dev_col] if std_dev_col else [None]*len(x_group)):
            if std_dev_val is not None and pd.notna(std_dev_val):
                lines.append(f"<b>{gene_name_val}</b> (Expr: {float(mean_val):.2f}, Std Dev: {float(std_dev_val):.2f})")
            else:
                lines.append(f"<b>{gene_name_val}</b> (Expr: {float(mean_val):.2f})")
        return '<br>'.join(lines)

    grouped = df_temp.groupby(subset_cols).apply(generate_hover_text).reset_index(name='Aggregated Hover Text')


    grouped.rename(columns={
        f'{mean_col}_numeric_rounded': mean_col,
        f'{std_dev_col}_numeric_rounded': std_dev_col # This will exist but might be NaN for SC
    }, inplace=True)
    
    return grouped


def plot_gene_expression_set(df_data, fem1_data_subset, plot_title_prefix, gene_col, mean_col, std_dev_col, group_col, group_color_map):
    st.subheader(f"{plot_title_prefix}: All Genes Across All Groups")
    st.write(f"This plot shows the {mean_col} vs {std_dev_col} for all genes, categorized by their assigned group. Hover over a point to see all overlapping Gene Names and their data.")

    removed_gene_name_plot1 = None
    removed_gene_mean_value_plot1 = None
    df_filtered_for_plot1 = df_data.copy()

    if not df_data.empty and mean_col in df_data.columns:
        # Ensure 'fem-1' is checked against the correct gene column
        df_for_outlier_check = df_data[df_data[gene_col] != 'fem-1'].copy()
        if not df_for_outlier_check.empty and pd.to_numeric(df_for_outlier_check[mean_col], errors='coerce').notna().any():
            max_mean_gene_row = df_for_outlier_check.loc[pd.to_numeric(df_for_outlier_check[mean_col], errors='coerce').idxmax()]
            removed_gene_name_plot1 = max_mean_gene_row[gene_col]
            removed_gene_mean_value_plot1 = float(max_mean_gene_row[mean_col])
            df_filtered_for_plot1 = df_data[df_data[gene_col] != removed_gene_name_plot1].copy()

    # Pass the appropriate std_dev_col (can be None for single cell)
    plot1_data_for_hover = create_aggregated_hover_data_flexible(
        df_filtered_for_plot1[df_filtered_for_plot1[gene_col] != 'fem-1'],
        gene_col, mean_col, std_dev_col, group_col
    )

    fig1 = go.Figure()
    unique_groups = sorted(plot1_data_for_hover[group_col].unique(), key=lambda x: int(x))

    for group_name in unique_groups:
        group_df = plot1_data_for_hover[plot1_data_for_hover[group_col] == group_name]
        if not group_df.empty:
            # Check if std_dev_col is valid and exists in group_df before using it for y-axis
            y_axis_data = group_df[std_dev_col] if std_dev_col and std_dev_col in group_df.columns else [0] * len(group_df)
            mode = 'markers'
            if not std_dev_col or std_dev_col not in group_df.columns: # If no std dev, plot as a line of points at y=0
                y_axis_title_text = "" # Y-axis title will be empty for single cell
                fig1.update_layout(yaxis_showticklabels=False)
            else:
                y_axis_title_text = std_dev_col # Use original std dev name for early embryo
                fig1.update_layout(yaxis_showticklabels=True)


            fig1.add_trace(go.Scatter(
                x=group_df[mean_col],
                y=y_axis_data,
                mode=mode,
                name=f'{group_name}',
                marker=dict(size=regular_dot_size, color=group_color_map.get(str(group_name), 'lightgray')),
                hoverinfo='text',
                text=group_df['Aggregated Hover Text'],
                hovertemplate='%{text}<extra></extra>'
            ))

    if not fem1_data_subset.empty and fem1_data_subset[gene_col].iloc[0] != removed_gene_name_plot1:
        fem1_data_plot1 = fem1_data_subset[fem1_data_subset[gene_col] != removed_gene_name_plot1]
        if not fem1_data_plot1.empty:
            fem1_hover_text = (
                f"<b>{fem1_data_plot1[gene_col].iloc[0]}</b>"
                f"<br>{mean_col}: {float(fem1_data_plot1[mean_col].iloc[0]):.2f}"
                f"<br>{std_dev_col}: {float(fem1_data_plot1[std_dev_col].iloc[0]):.2f}"
                f"<br>Group: {fem1_data_plot1[group_col].iloc[0]}"
            )
            # Y-axis data for fem-1 (handle if std_dev_col is missing)
            fem1_y_data = fem1_data_plot1[std_dev_col] if std_dev_col and std_dev_col in fem1_data_plot1.columns else [0] * len(fem1_data_plot1)
            fig1.add_trace(go.Scatter(
                x=fem1_data_plot1[mean_col],
                y=fem1_y_data,
                mode='markers',
                marker=dict(size=fem1_dot_size, color='red', symbol='circle', line=dict(width=2, color='DarkRed')),
                name='fem-1',
                hoverinfo='text',
                text=[fem1_hover_text],
                hovertemplate='%{text}<extra></extra>'
            ))
            
    fig1.update_layout(
        title=f'{plot_title_prefix} Genes: {mean_col} vs {y_axis_title_text}',
        xaxis_title=mean_col,
        yaxis_title=y_axis_title_text, # Set dynamically
        font_family="Times New Roman",
        title_font_size=20,
        xaxis_title_font_size=14,
        yaxis_title_font_size=14,
        xaxis_type='log', # Set X-axis to logarithmic scale
        yaxis_type='log' if std_dev_col and std_dev_col in df_data.columns else 'linear', # Y-axis log only if std dev is present
        xaxis_exponentformat='power', # Display exponents
        xaxis_showexponent='all',    # Show exponent for all ticks
        xaxis_tickformat='e',        # Use scientific notation for ticks
        yaxis_exponentformat='power' if std_dev_col and std_dev_col in df_data.columns else 'none',
        yaxis_showexponent='all' if std_dev_col and std_dev_col in df_data.columns else 'none',
        yaxis_tickformat='e' if std_dev_col and std_dev_col in df_data.columns else '',
        xaxis_tickangle=90,
        yaxis_tickangle=0,
        width=900,
        height=600,
        legend_title_text='Group'
    )
    st.plotly_chart(fig1)
    if removed_gene_name_plot1:
        st.markdown(f'<p style="font-family:\'Times New Roman\', serif; font-size:11px; color:gray; text-align:center;">Note: The gene <b>{removed_gene_name_plot1}</b> (Mean: {removed_gene_mean_value_plot1:.2f}) was removed to improve plot clarity as it was the single highest outlier in "{mean_col}" across all genes.</p>', unsafe_allow_html=True)

    st.markdown("---")

    st.subheader(f"{plot_title_prefix}: Group 9 Genes")
    st.write(f"This plot focuses on genes within Group 9. Hover over a point to see all overlapping Gene Names and their data.")

    group9_data = df_data[df_data[group_col] == '9'].copy()
    plot2_data_for_hover = create_aggregated_hover_data_flexible(
        group9_data[group9_data[gene_col] != 'fem-1'],
        gene_col, mean_col, std_dev_col, group_col
    )

    fig2 = go.Figure()
    if not plot2_data_for_hover.empty:
        y_axis_data_2 = plot2_data_for_hover[std_dev_col] if std_dev_col and std_dev_col in plot2_data_for_hover.columns else [0] * len(plot2_data_for_hover)
        mode_2 = 'markers'
        if not std_dev_col or std_dev_col not in plot2_data_for_hover.columns:
            y_axis_title_text_2 = ""
            fig2.update_layout(yaxis_showticklabels=False)
        else:
            y_axis_title_text_2 = std_dev_col
            fig2.update_layout(yaxis_showticklabels=True)

        fig2.add_trace(go.Scatter(
            x=plot2_data_for_hover[mean_col],
            y=y_axis_data_2,
            mode=mode_2,
            marker=dict(size=regular_dot_size, color=group_color_map.get('9', 'blue')),
            hoverinfo='text',
            text=plot2_data_for_hover['Aggregated Hover Text'],
            hovertemplate='%{text}<extra></extra>'
        ))

    fem1_in_group9 = group9_data[group9_data[gene_col] == 'fem-1']
    if not fem1_in_group9.empty:
        fem1_hover_text_plot2 = (
            f"<b>{fem1_in_group9[gene_col].iloc[0]}</b>"
            f"<br>{mean_col}: {float(fem1_in_group9[mean_col].iloc[0]):.2f}"
            f"<br>{std_dev_col}: {float(fem1_in_group9[std_dev_col].iloc[0]):.2f}"
            f"<br>Group: {fem1_in_group9[group_col].iloc[0]}"
        )
        fem1_y_data_2 = fem1_in_group9[std_dev_col] if std_dev_col and std_dev_col in fem1_in_group9.columns else [0] * len(fem1_in_group9)
        fig2.add_trace(go.Scatter(
            x=fem1_in_group9[mean_col],
            y=fem1_y_data_2,
            mode='markers',
            marker=dict(size=fem1_dot_size, color='red', symbol='circle', line=dict(width=2, color='DarkRed')),
            name='fem-1',
            hoverinfo='text',
            text=[fem1_hover_text_plot2],
            hovertemplate='%{text}<extra></extra>'
        ))

    fig2.update_layout(
        title=f'{plot_title_prefix}: Group 9 Genes: {mean_col} vs {y_axis_title_text_2}',
        xaxis_title=mean_col,
        yaxis_title=y_axis_title_text_2,
        font_family="Times New Roman",
        title_font_size=20,
        xaxis_title_font_size=14,
        yaxis_title_font_size=14,
        xaxis_type='log', # Set X-axis to logarithmic scale
        yaxis_type='log' if std_dev_col and std_dev_col in df_data.columns else 'linear', # Y-axis log only if std dev is present
        xaxis_exponentformat='power',
        yaxis_showexponent='all' if std_dev_col and std_dev_col in df_data.columns else 'none',
        yaxis_tickformat='e' if std_dev_col and std_dev_col in df_data.columns else '',
        yaxis_exponentformat='power',
        xaxis_tickangle=90,
        yaxis_tickangle=0,
        width=900,
        height=600
    )
    st.plotly_chart(fig2)

    st.markdown("---")

    st.subheader(f"{plot_title_prefix}: Genes in Groups 8, 9, and 10")
    st.write(f"This plot shows genes from the top three groups. Hover over a point to see all overlapping Gene Names and their data.")

    selected_groups_data_raw = df_data[df_data[group_col].isin(['8', '9', '10'])].copy()

    removed_gene_name_plot3 = None
    removed_gene_mean_value_plot3 = None

    if not selected_groups_data_raw.empty and mean_col in selected_groups_data_raw.columns:
        df_for_outlier_check_plot3 = selected_groups_data_raw[selected_groups_data_raw[gene_col] != 'fem-1'].copy()
        if not df_for_outlier_check_plot3.empty and pd.to_numeric(df_for_outlier_check_plot3[mean_col], errors='coerce').notna().any():
            max_mean_gene_row_plot3 = df_for_outlier_check_plot3.loc[pd.to_numeric(df_for_outlier_check_plot3[mean_col], errors='coerce').idxmax()]
            removed_gene_name_plot3 = max_mean_gene_row_plot3[gene_col]
            removed_gene_mean_value_plot3 = float(max_mean_gene_row_plot3[mean_col])
            selected_groups_data_for_plot3 = selected_groups_data_raw[selected_groups_data_raw[gene_col] != removed_gene_name_plot3].copy()
        else:
            selected_groups_data_for_plot3 = selected_groups_data_raw.copy()
    else:
        selected_groups_data_for_plot3 = selected_groups_data_raw.copy()

    plot3_data_for_hover = create_aggregated_hover_data_flexible(
        selected_groups_data_for_plot3[selected_groups_data_for_plot3[gene_col] != 'fem-1'],
        gene_col, mean_col, std_dev_col, group_col
    )

    fig3 = go.Figure()
    unique_groups_plot3 = sorted(plot3_data_for_hover[group_col].unique(), key=lambda x: int(x))

    for group_name in unique_groups_plot3:
        group_df = plot3_data_for_hover[plot3_data_for_hover[group_col] == group_name]
        if not group_df.empty:
            y_axis_data_3 = group_df[std_dev_col] if std_dev_col and std_dev_col in group_df.columns else [0] * len(group_df)
            mode_3 = 'markers'
            if not std_dev_col or std_dev_col not in group_df.columns:
                y_axis_title_text_3 = ""
                fig3.update_layout(yaxis_showticklabels=False)
            else:
                y_axis_title_text_3 = std_dev_col
                fig3.update_layout(yaxis_showticklabels=True)

            fig3.add_trace(go.Scatter(
                x=group_df[mean_col],
                y=y_axis_data_3,
                mode=mode_3,
                name=f'{group_name}',
                marker=dict(size=regular_dot_size, color=group_color_map.get(str(group_name), 'gray')),
                hoverinfo='text',
                hovertext=group_df['Aggregated Hover Text'],
                hovertemplate='%{text}<extra></extra>'
            ))

    fem1_in_selected_groups = selected_groups_data_for_plot3[selected_groups_data_for_plot3[gene_col] == 'fem-1']
    if not fem1_in_selected_groups.empty:
        fem1_hover_text_plot3 = (
            f"<b>{fem1_in_selected_groups[gene_col].iloc[0]}</b>"
            f"<br>{mean_col}: {float(fem1_in_selected_groups[mean_col].iloc[0]):.2f}"
            f"<br>{std_dev_col}: {float(fem1_in_selected_groups[std_dev_col].iloc[0]):.2f}"
            f"<br>Group: {fem1_in_selected_groups[group_col].iloc[0]}"
        )
        fem1_y_data_3 = fem1_in_selected_groups[std_dev_col] if std_dev_col and std_dev_col in fem1_in_selected_groups.columns else [0] * len(fem1_in_selected_groups)
        fig3.add_trace(go.Scatter(
            x=fem1_in_selected_groups[mean_col],
            y=fem1_y_data_3,
            mode='markers',
            marker=dict(size=fem1_dot_size, color='red', symbol='circle', line=dict(width=2, color='DarkRed')),
            name='fem-1',
            hoverinfo='text',
            text=[fem1_hover_text_plot3],
            hovertemplate='%{text}<extra></extra>'
        ))

    fig3.update_layout(
        title=f'{plot_title_prefix}: Groups 8, 9, and 10 Genes: {mean_col} vs {y_axis_title_text_3}',
        xaxis_title=mean_col,
        yaxis_title=y_axis_title_text_3,
        font_family="Times New Roman",
        title_font_size=20,
        xaxis_title_font_size=14,
        yaxis_title_font_size=14,
        xaxis_type='log', # Log scale X-axis
        yaxis_type='log' if std_dev_col and std_dev_col in df_data.columns else 'linear', # Y-axis log only if std dev is present
        xaxis_exponentformat='power',
        yaxis_showexponent='all' if std_dev_col and std_dev_col in df_data.columns else 'none',
        yaxis_tickformat='e' if std_dev_col and std_dev_col in df_data.columns else '',
        yaxis_exponentformat='power',
        xaxis_tickangle=90,        # Rotate x-axis labels
        yaxis_tickangle=0,         # Keep y-axis labels horizontal
        width=900,
        height=600,
        legend_title_text='Group'
    )
    st.plotly_chart(fig3)
    if removed_gene_name_plot3:
        st.markdown(f'<p style="font-family:\'Times New Roman\', serif; font-size:11px; color:gray; text-align:center;">Note: The gene <b>{removed_gene_name_plot3}</b> (Mean: {removed_gene_mean_value_plot3:.2f}) was removed to improve plot clarity as it was the single highest outlier in "{mean_col}" within {group_col}s 8, 9, and 10.</p>', unsafe_allow_html=True)


def plot_single_cell_expression_set(df_data_sc, fem1_data_sc_subset, plot_title_prefix, gene_col, tpm_col, group_col, group_color_map):
    
    st.subheader(f"{plot_title_prefix}: All Genes (Sorted by {tpm_col})")
    st.write(f"This line graph shows {tpm_col} for all genes, sorted by expression, with points color-coded by their assigned group. Zoom in to see individual gene points.")
    
    df_sorted_by_tpm = df_data_sc.sort_values(by=tpm_col, ascending=True).reset_index(drop=True)
    
    removed_gene_name_plot_sc = None
    removed_gene_tpm_value_plot_sc = None

    if not df_sorted_by_tpm.empty and tpm_col in df_sorted_by_tpm.columns:
        df_for_outlier_check_sc = df_sorted_by_tpm[df_sorted_by_tpm[gene_col] != 'fem-1'].copy()
        if not df_for_outlier_check_sc.empty and pd.to_numeric(df_for_outlier_check_sc[tpm_col], errors='coerce').notna().any():
            max_tpm_gene_row_sc = df_for_outlier_check_sc.loc[pd.to_numeric(df_for_outlier_check_sc[tpm_col], errors='coerce').idxmax()]
            removed_gene_name_plot_sc = max_tpm_gene_row_sc[gene_col]
            removed_gene_tpm_value_plot_sc = float(max_tpm_gene_row_sc[tpm_col])
            df_sorted_by_tpm = df_sorted_by_tpm[df_sorted_by_tpm[gene_col] != removed_gene_name_plot_sc].copy()

    fig_sc = go.Figure()

    for group_name_sc in sorted(df_sorted_by_tpm[group_col].unique(), key=lambda x: int(x) if str(x).isdigit() else x):
        group_df_sc = df_sorted_by_tpm[df_sorted_by_tpm[group_col] == group_name_sc]
        if not group_df_sc.empty:
            fig_sc.add_trace(go.Scatter(
                x=group_df_sc[tpm_col],
                y=[0] * len(group_df_sc),
                mode='markers',
                name=f'Group {group_name_sc}', # Use "Group" prefix for clarity in legend
                marker=dict(size=regular_dot_size, color=group_color_map.get(str(group_name_sc), 'lightgray')),
                hoverinfo='text',
                text=[
                    f"<b>{row[gene_col]}</b><br>{tpm_col}: {float(row[tpm_col]):.2f}<br>Group: {row[group_col]}"
                    for idx, row in group_df_sc.iterrows()
                ],
                hovertemplate='%{text}<extra></extra>'
            ))

    fem1_in_sorted_sc = df_sorted_by_tpm[df_sorted_by_tpm[gene_col] == 'fem-1']
    if not fem1_in_sorted_sc.empty:
        fem1_hover_text_sc = (
            f"<b>{fem1_in_sorted_sc[gene_col].iloc[0]}</b>"
            f"<br>{tpm_col}: {float(fem1_in_sorted_sc[tpm_col].iloc[0]):.2f}"
            f"<br>Group: {fem1_in_sorted_sc[group_col].iloc[0]}"
        )
        fig_sc.add_trace(go.Scatter(
            x=fem1_in_sorted_sc[tpm_col],
            y=[0] * len(fem1_in_sorted_sc),
            mode='markers',
            marker=dict(size=fem1_dot_size, color='red', symbol='circle', line=dict(width=2, color='DarkRed')),
            name='fem-1',
            hoverinfo='text',
            text=[fem1_hover_text_sc],
            hovertemplate='%{text}<extra></extra>'
        ))
    
    fig_sc.update_layout(
        title=f'{plot_title_prefix} Genes: {tpm_col} Distribution',
        xaxis_title=tpm_col,
        yaxis_title='',
        yaxis_showticklabels=False,
        font_family="Times New Roman",
        title_font_size=20,
        xaxis_title_font_size=14,
        xaxis_type='log', # Set X-axis to logarithmic scale
        xaxis_exponentformat='power',
        xaxis_showexponent='all',
        xaxis_tickformat='e',
        xaxis_tickangle=90,
        width=1200,
        height=300,
        legend_title_text='Group'
    )
    st.plotly_chart(fig_sc)
    if removed_gene_name_plot_sc:
        st.markdown(f'<p style="font-family:\'Times New Roman\', serif; font-size:11px; color:gray; text-align:center;">Note: The gene <b>{removed_gene_name_plot_sc}</b> ({tpm_col}: {removed_gene_tpm_value_plot_sc:.2f}) was removed to improve plot clarity as it was the single highest outlier in "{tpm_col}".</p>', unsafe_allow_html=True)

    st.markdown("---")


# --- New Home Page Function (unchanged) ---
def home_page():
    st.markdown("""
    <style>
        .stApp {
            background-color: white; /* Changed background to white */
        }
        .main > div {
            padding-top: 1rem;
            padding-bottom: 1rem;
        }
        
        .hero-section {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            padding: 2.5rem 2rem;
            border-radius: 20px;
            margin: 1rem auto 2rem auto;
            color: white;
            text-align: center;
            box-shadow: 0 10px 30px rgba(0,0,0,0.2);
            max-width: 800px;
        }
        
        .hero-title {
            font-size: 3.5rem;
            font-weight: bold;
            margin-bottom: 1rem;
            text-shadow: 2px 2px 4px rgba(0,0,0,0.3);
            text-align: center;
            line-height: 1.2;
        }
        
        .hero-slogan {
            font-size: 1.3rem;
            opacity: 0.9;
            font-style: italic;
            text-align: center;
        }
        
        /* Override Streamlit button styling completely */
        div.stButton > button:first-child {
            background: linear-gradient(45deg, #e0e0e0 0%, #c0c0c0 100%) !important; /* Lighter background for white page */
            color: black !important; /* Black font for white background */
            border: 2px solid #667eea !important; /* Add a subtle border */
            border-radius: 25px !important;
            font-size: 2.5rem !important; /* Adjusted font size for better fit */
            font-weight: bold !important;
            width: 300px !important; /* Adjusted width for better centering within a flex container */
            height: 80px !important; /* Adjusted height */
            transition: all 0.3s ease !important;
            box-shadow: 0 8px 25px rgba(0,0,0,0.3) !important;
            margin: 1rem auto !important; /* Ensure vertical spacing, horizontal centering via parent */
            display: flex !important; /* Use flexbox for button content centering if needed */
            justify-content: center !important; /* Center text horizontally */
            align-items: center !important; /* Center text vertically */
            line-height: 1 !important; /* Reset line-height to 1 for vertical centering */
            padding: 0 !important;
        }
        
        div.stButton > button:first-child:hover {
            transform: translateY(-5px) !important;
            box_shadow: 0 12px 35px rgba(0,0,0,0.4) !important;
            background: linear_gradient(45deg, #d0d0d0 0%, #b0b0b0 100%) !important; /* Slightly darker on hover */
        }
        
        div.stButton > button:first-child:active {
            transform: translateY(-2px) !important;
        }
        
        .button-container {
            display: flex;
            flex-direction: column;
            align-items: center; /* This centers the buttons horizontally within this container */
            gap: 1.5rem;
            margin: 1rem auto; /* This centers the container itself within center_col */
            max_width: 400px; /* Constrain max-width of container for better centering effect */
        }
        
        .side-panel {
            background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%);
            padding: 2rem;
            border-radius: 15px;
            color: white;
            text-align: center;
            box-shadow: 0 5px 15px rgba(0,0,0,0.2);
            margin: 1rem 0;
            height: fit-content;
        }
        
        .side-panel h3 {
            margin-bottom: 1rem;
            font-size: 1.5rem;
        }
        
        .side-panel p {
            font-size: 1rem;
            opacity: 0.9;
            line-height: 1.4;
        }
        
        .helpful-links {
            background: linear-gradient(135deg, #a8edea 0%, #fed6e3 100%);
            padding: 1.5rem;
            border-radius: 15px;
            color: #333;
            text-align: center;
            box-shadow: 0 5px 15px rgba(0,0,0,0.1);
            margin: 1rem 0;
        }
        
        .helpful-links h3 {
            color: #667eea;
            margin-bottom: 1rem;
            font-size: 1.3rem;
        }
        
        .helpful-links a {
            display: block;
            color: #667eea;
            text-decoration: none;
            margin: 0.5rem 0;
            font-weight: bold;
        }
        
        .helpful-links a:hover {
            text-decoration: underline;
        }
        
        .research-objectives {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            padding: 2rem;
            border-radius: 15px;
            color: white;
            box-shadow: 0 5px 15px rgba(0,0,0,0.2);
            margin: 1rem 0;
            height: 520px; /* Adjusted to try and fit content better */
            display: flex;
            flex-direction: column;
        }
        
        .research-objectives h3 {
            margin-bottom: 1.5rem;
            font-size: 1.5rem;
            text-align: center;
        }
        
        .research-objectives ol {
            text-align: left;
            padding-left: 1rem;
            flex-grow: 1;
            display: flex;
            flex-direction: column;
            justify-content: space-around;
        }
        
        .research-objectives li {
            margin-bottom: 1.5rem;
            line-height: 1.4;
            font-size: 0.95rem;
        }
        
        .floating-emoji {
            font-size: 3rem;
            animation: float 3s ease-in-out infinite;
            display: inline-block;
            margin: 0.5rem;
        }
        
        .floating-emoji:nth-child(2) {
            animation-delay: 0.5s;
        }
        
        .floating-emoji:nth-child(3) {
            animation-delay: 1s;
        }
        
        @keyframes float {
            0%, 100% { transform: translateY(0px); }
            50% { transform: translateY(-10px); }
        }
        
        .footer-section {
            background-color: #f8f9fa; /* Lighter background for white page */
            border-radius: 10px;
            padding: 1.5rem;
            margin-top: 1rem; /* Adjusted margin-top to bring it closer to buttons */
            text-align: center;
            max-width: 800px;
            margin-left: auto;
            margin-right: auto;
            box-shadow: 0 5px 15px rgba(0,0,0,0.1); /* Added subtle shadow */
        }
        
        .contact-info {
            font-size: 0.9rem;
            color: #666;
            margin-bottom: 1rem;
        }
        
        .quote-section {
            font-style: italic;
            color: #555;
            font-size: 1rem;
            border-top: 1px solid #ddd;
            padding-top: 1rem;
        }
    </style>
    """, unsafe_allow_html=True)

    # Create three columns for layout
    left_col, center_col, right_col = st.columns([1, 2, 1])

    # Left sidebar content
    with left_col:
        st.markdown("""
        <div class="side-panel">
            <h3>🔬 Research Focus</h3>
            <p>To determine the molecular mechanism by which maternal RNA regulates fem-1 expression in C. elegans, with emphasis on how this regulation is influenced by parent-of-origin effects.</p>
        </div>
        """, unsafe_allow_html=True)
        
        st.markdown("""
        <div class="helpful-links">
            <h3>🔗 Helpful Links</h3>
            <a href="https://docs.google.com/document/d/1kNxQVg3Y1rGJ9-6C6icEoDH44qDx5zQPyCPlS7HfsiY/edit?usp=sharing" target="_blank">Methods Document</a>
            <a href="https://37nyza-abbas-ghaddar.shinyapps.io/shiny_webpage/" target="_blank">Single Cell Database</a>
            <a href="https://www.wormbase.org/" target="_blank">WormBase Database</a>
        </div>
        """, unsafe_allow_html=True)
        
        st.markdown("""
        <div class="floating-emoji">🧪</div>
        <div class="floating-emoji">⚗️</div>
        <div class="floating-emoji">🔬</div>
        """, unsafe_allow_html=True)

    # Center content
    with center_col:
        # Hero Section - smaller and pushed up
        st.markdown("""
        <div class="hero-section">
            <div class="hero-title">🧬 Saurish and Xander's<br>Biomart 🔬</div>
            <div class="hero-slogan">Science for the benefit of humanity</div>
        </div>
        """, unsafe_allow_html=True)

        st.markdown('<div class="button-container">', unsafe_allow_html=True)
        
        # Original Data Button
        if st.button("📊 Original Data", key="original_data_btn", help="Access raw data tables and visualizations from the initial dataset."):
            st.session_state.page = "original_data_landing"
            st.rerun()
            
        # Processed Data Button (This button now leads to processed data landing, but actual "Processed Data" page is still pending full implementation beyond single-cell. User asked for this, so keeping it.)
        if st.button("✨ Processed Data", key="processed_data_btn", help="Explore processed single-cell data."):
            st.session_state.page = "processed_data_landing"
            st.rerun()
            
        st.markdown('</div>', unsafe_allow_html=True)

        # Footer section with contact info and quote - MOVED HERE
        st.markdown("""
        <div class="footer-section">
            <div class="contact-info">
                If you have any questions, email <a href="mailto:sarora@rockefeller.edu" style="color: #667eea; text-decoration: none; font-weight: bold;">sarora@rockefeller.edu</a> or text at <a href="tel:+19089302303" style="color: #667eea; text-decoration: none; font-weight: bold;">(908) 930-2303</a>
            </div>
            <div class="quote-section">
                "The good thing about science is that it's true whether or not you believe in it."<br>
                <small>- Neil deGrasse Tyson</small>
            </div>
        </div>
        """, unsafe_allow_html=True)

    # Right sidebar content
    with right_col:
        st.markdown("""
        <div class="research-objectives">
            <h3>🎯 Research Objectives</h3>
            <ol>
                <li>To recapitulate and characterize fem-1–related phenotypes through targeted genetic crosses to confirm parent-of-origin effects.</li>
                <li>To perform bioinformatic analysis of publicly available datasets to identify other genes exhibiting similar maternal RNA–dependent expression patterns as fem-1.</li>
                <li>To generate a mutant strain enabling a genetic screen for fem-1 function, allowing selection and analysis of specific genotype combinations.</li>
            </ol>
        </div>
        """, unsafe_allow_html=True)
        
        st.markdown("""
        <div class="floating-emoji">🧬</div>
        <div class="floating-emoji">📊</div>
        <div class="floating-emoji">🔍</div>
        """, unsafe_allow_html=True)

    # Add some decorative elements
    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown("""
    <div style="text-align: center; opacity: 0.6;">
        🧪 ⚗️ 🔬 🧬 📊 📈 🔍 ⚡ 🧫 🔭 ⚛️ 🌡️
    </div>
    """, unsafe_allow_html=True)


# --- Intermediate Page for Original Data ---
def original_data_landing_page():
    st.header("Original Data: Tables & Visualizations")
    st.write("Choose how you'd like to explore the early embryo data.")

    if st.button("🏠 Back to Home", key="original_landing_back_home"):
        st.session_state.page = "home"
        st.rerun()

    st.markdown("---")

    col1, col2, col3 = st.columns(3) # Added a third column for the new button

    with col1:
        st.markdown('<div class="button-container" style="max-width: 250px;">', unsafe_allow_html=True) # Smaller max-width for these buttons
        if st.button("🔍 Raw Data", key="view_raw_data_tables", help="Browse the raw 'AnalysisFile2.txt' data."):
            st.session_state.page = "raw_data"
            st.rerun()
        st.markdown('</div>', unsafe_allow_html=True)

    with col2:
        st.markdown('<div class="button-container" style="max-width: 250px;">', unsafe_allow_html=True) # Smaller max-width for these buttons
        if st.button("📈 Visualizations", key="view_visualizations", help="See plots and graphs generated from the 'AnalysisFile2.txt' data."):
            st.session_state.page = "visualizations"
            st.rerun()
        st.markdown('</div>', unsafe_allow_html=True)

    with col3: # New column for Comparisons button
        st.markdown('<div class="button-container" style="max-width: 250px;">', unsafe_allow_html=True)
        if st.button("🔄 Comparisons", key="view_comparisons", help="Compare gene expression across different datasets."):
            st.session_state.page = "comparison_graphs"
            st.rerun()
        st.markdown('</div>', unsafe_allow_html=True)


# --- Processed Data Landing Page (Currently Blank - can be expanded later) ---
def processed_data_landing_page():
    st.header("Processed Data Analysis")
    st.write("This section is dedicated to processed data, such as single-cell RNA sequencing analysis.")
    st.info("Functionality for processed data will be added here soon!")

    if st.button("🏠 Back to Home", key="processed_landing_back_home"):
        st.session_state.page = "home"
        st.rerun()

    st.markdown("---")


# --- Main Visualization Page ---
def visualizations_page():
    st.header("Gene Expression Visualizations")
    st.write("Explore gene expression patterns through interactive scatter plots.")

    if st.button("🏠 Back to Home", key="viz_back_home"):
        st.session_state.page = "home"
        st.rerun()
    
    st.markdown("---")

    visualization_type = st.radio(
        "Select Visualization Type:",
        ("Early Embryo Data", "Single-Cell Data"), # Removed "Comparison Graphs"
        key="viz_type_radio"
    )

    st.markdown("---")

    if visualization_type == "Early Embryo Data":
        st.subheader("Early Embryo Data Visualizations")
        st.write("These plots display gene expression patterns from the original 'AnalysisFile2.txt' dataset.")
        
        plot_gene_expression_set(
            df_original, 
            fem1_data_original, 
            "Early Embryo", 
            'Gene Name', 
            'Mean of Geneid Strains', 
            'Standard Deviation of Geneid Strains', 
            'Group',
            {str(g): px.colors.qualitative.Plotly[i % len(px.colors.qualitative.Plotly)] 
             for i, g in enumerate(sorted(df_original['Group'].unique(), key=lambda x: str(x)))}
        )

    elif visualization_type == "Single-Cell Data":
        st.subheader("Single-Cell Data Visualizations")
        st.write("Select a single-cell dataset to view its gene expression patterns.")

        if not single_cell_dataframes_categorized['Germ'] and not single_cell_dataframes_categorized['Somatic']:
            st.warning("No single-cell processed data files found or loaded. Please ensure they exist in the specified directories.")
            return

        # New dropdowns for category and then dataset
        category_options = []
        if single_cell_dataframes_categorized['Germ']:
            category_options.append("Germ")
        if single_cell_dataframes_categorized['Somatic']:
            category_options.append("Somatic")

        selected_category = st.selectbox(
            "Choose Data Category:",
            category_options,
            key="selected_sc_category"
        )
        
        if selected_category:
            current_category_dfs = single_cell_dataframes_categorized[selected_category]
            if not current_category_dfs:
                st.info(f"No datasets available in '{selected_category}' category.")
                return

            selected_single_cell_dataset = st.selectbox(
                f"Choose a {selected_category} Dataset:",
                list(current_category_dfs.keys()),
                key="selected_sc_viz_dataset"
            )

            if selected_single_cell_dataset:
                current_sc_df = current_category_dfs[selected_single_cell_dataset]
                # fem1_data_sc needs to be filtered using 'gene_common'
                fem1_data_sc = current_sc_df[current_sc_df['gene_common'] == 'fem-1']

                # Colors for single cell groups (standard 1-10)
                sc_group_color_map = {
                    str(g): px.colors.qualitative.D3[i % len(px.colors.qualitative.D3)]
                    for i, g in enumerate(range(1, 11))
                }

                plot_single_cell_expression_set(
                    current_sc_df,
                    fem1_data_sc,
                    selected_single_cell_dataset,
                    'gene_common',      # Use standardized gene column
                    'expression_value', # Use standardized expression column
                    'group_identifier', # Use standardized group column
                    sc_group_color_map
                )

# --- New Comparison Graphs Page ---
def comparison_page():
    st.header("Gene Expression Comparison Graphs")
    st.write("Compare gene expression between the 'Early Embryo' dataset and a selected single-cell dataset.")

    if st.button("🏠 Back to Home", key="comp_back_home"):
        st.session_state.page = "home"
        st.rerun()
    
    st.markdown("---")

    if not single_cell_dataframes_categorized['Germ'] and not single_cell_dataframes_categorized['Somatic']:
        st.warning("No single-cell processed data files found or loaded to perform comparisons.")
        return

    # New dropdowns for category and then dataset for comparison
    category_options = []
    if single_cell_dataframes_categorized['Germ']:
        category_options.append("Germ")
    if single_cell_dataframes_categorized['Somatic']:
        category_options.append("Somatic")

    selected_category_comp = st.selectbox(
        "Choose Single-Cell Data Category for Comparison:",
        category_options,
        key="selected_sc_category_comp"
    )
    
    selected_comparison_dataset = None
    if selected_category_comp:
        current_category_dfs_comp = single_cell_dataframes_categorized[selected_category_comp]
        if not current_category_dfs_comp:
            st.info(f"No datasets available in '{selected_category_comp}' category for comparison.")
            return

        selected_comparison_dataset = st.selectbox(
            f"Choose a {selected_category_comp} Dataset for Comparison:",
            list(current_category_dfs_comp.keys()),
            key="selected_comparison_dataset"
        )


    if selected_comparison_dataset:
        sc_df_for_comparison = single_cell_dataframes_categorized[selected_category_comp][selected_comparison_dataset]

        # Prepare data for merging - standardize df_original's gene name for common key
        df_original_renamed = df_original.rename(columns={'Gene Name': 'gene_common'})
        
        # Merge on standardized gene_common column
        # df_original: 'gene_common', 'Mean of Geneid Strains', 'Group'
        # sc_df_for_comparison: 'gene_common', 'expression_value', 'group_identifier'
        merged_df = pd.merge(
            df_original_renamed[['gene_common', 'Mean of Geneid Strains', 'Group']],
            sc_df_for_comparison[['gene_common', 'expression_value', 'group_identifier']],
            on='gene_common',
            how='inner'
        )
        
        merged_df_sorted = merged_df.sort_values(by='gene_common').reset_index(drop=True)

        if merged_df_sorted.empty:
            st.warning(f"No common genes found between 'Early Embryo' and '{selected_comparison_dataset}' for comparison.")
            return

        # Determine and remove the highest outlier gene to improve plot clarity
        removed_gene_name_comp = None
        removed_expr_value_comp = None
        removed_dataset_comp = None

        # Check for outliers in both datasets' expression columns
        max_early_embryo_expr = merged_df_sorted['Mean of Geneid Strains'].max()
        max_single_cell_expr = merged_df_sorted['expression_value'].max()

        # Only remove if it's not fem-1 to ensure fem-1 is always shown
        if pd.notna(max_early_embryo_expr) and (pd.notna(max_single_cell_expr) and max_early_embryo_expr >= max_single_cell_expr):
            outlier_row = merged_df_sorted.loc[merged_df_sorted['Mean of Geneid Strains'].idxmax()]
            if outlier_row['gene_common'] != 'fem-1':
                removed_gene_name_comp = outlier_row['gene_common']
                removed_expr_value_comp = outlier_row['Mean of Geneid Strains']
                removed_dataset_comp = "Early Embryo"
        elif pd.notna(max_single_cell_expr):
            outlier_row = merged_df_sorted.loc[merged_df_sorted['expression_value'].idxmax()]
            if outlier_row['gene_common'] != 'fem-1':
                removed_gene_name_comp = outlier_row['gene_common']
                removed_expr_value_comp = outlier_row['expression_value']
                removed_dataset_comp = selected_comparison_dataset
        
        if removed_gene_name_comp:
            merged_df_sorted = merged_df_sorted[merged_df_sorted['gene_common'] != removed_gene_name_comp].copy()


        st.markdown("---")
        st.subheader("Search Gene for Highlight & Difference")
        search_gene_comp = st.text_input("Enter Gene Name to Highlight (e.g., fem-1):", key="search_gene_comp").strip()

        fig_comp = go.Figure()

        # Create color map for Single-Cell groups (standard 1-10)
        single_cell_colors = px.colors.qualitative.D3
        single_cell_color_map = {str(g): single_cell_colors[i % len(single_cell_colors)] for i, g in enumerate(range(1, 11))}

        # Add traces for actual data points (colored by Single-Cell Group)
        for sc_group in single_cell_groups_sorted:
            # Note: 'group_identifier' is a string from standardization, convert for comparison if needed
            subset_sc_df = merged_df_sorted[merged_df_sorted['group_identifier'] == str(sc_group)]
            if not subset_sc_df.empty:
                fig_comp.add_trace(go.Scatter(
                    x=subset_sc_df['Mean of Geneid Strains'],
                    y=subset_sc_df['expression_value'], # Use standardized name
                    mode='markers',
                    name=f'SC Group {sc_group}', # Legend for SC groups
                    marker=dict(
                        size=regular_dot_size,
                        color=single_cell_color_map.get(str(sc_group), 'gray'),
                        symbol='circle' # All points are circles
                    ),
                    hoverinfo='text',
                    text=[
                        f"<b>{row['gene_common']}</b><br>EE Expr: {row['Mean of Geneid Strains']:.2f} (Group: {row['Group']})"
                        f"<br>{selected_comparison_dataset} Expr: {row['expression_value']:.2f} (Group: {row['group_identifier']})" # Use standardized names
                        for idx, row in subset_sc_df.iterrows()
                    ],
                    hovertemplate='%{text}<extra></extra>'
                ))
        
        # Add vertical lines for Early Embryo group max expression
        shapes = []
        for ee_group in early_embryo_groups_sorted:
            ee_group_data = merged_df_sorted[merged_df_sorted['Group'] == ee_group]
            if not ee_group_data.empty and 'Mean of Geneid Strains' in ee_group_data.columns:
                max_ee_expr = ee_group_data['Mean of Geneid Strains'].max()
                if pd.notna(max_ee_expr):
                    shapes.append(
                        dict(
                            type="line",
                            xref="x", yref="paper",
                            x0=max_ee_expr, y0=0,
                            x1=max_ee_expr, y1=1,
                            line=dict(color="LightCoral", width=2, dash="dash"),
                            layer='below' # Draw lines below scatter points
                        )
                    )
        fig_comp.update_layout(shapes=shapes)

        # Add a dummy trace for the vertical line legend entry
        fig_comp.add_trace(go.Scatter(
            x=[None], y=[None],
            mode='lines',
            line=dict(color="LightCoral", width=2, dash="dash"),
            name='Early Embryo Group Max Expression',
            hoverinfo='none'
        ))


        # Highlight searched gene if found
        if search_gene_comp:
            highlighted_data = merged_df_sorted[merged_df_sorted['gene_common'].str.lower() == search_gene_comp.lower()]
            if not highlighted_data.empty:
                for idx, row in highlighted_data.iterrows():
                    fig_comp.add_trace(go.Scatter(
                        x=[row['Mean of Geneid Strains']],
                        y=[row['expression_value']], # Use standardized name
                        mode='markers',
                        showlegend=False,
                        marker=dict(
                            size=fem1_dot_size * 2,
                            color='cyan',
                            symbol='star',
                            line=dict(width=3, color='darkblue')
                        ),
                        hoverinfo='text',
                        text=[
                            f"<b>{row['gene_common']}</b><br>EE Expr: {row['Mean of Geneid Strains']:.2f} (Group: {row['Group']})"
                            f"<br>{selected_comparison_dataset} Expr: {row['expression_value']:.2f} (Group: {row['group_identifier']})" # Use standardized names
                            f"<br>HIGHLIGHTED"
                        ],
                        hovertemplate='%{text}<extra></extra>'
                    ))

        fig_comp.update_layout(
            title=f'Gene Expression Comparison: Early Embryo vs. {selected_comparison_dataset}',
            xaxis_title='Early Embryo: Mean of Geneid Strains',
            yaxis_title=f'{selected_comparison_dataset}: Expression Value', # Use generic "Expression Value"
            font_family="Times New Roman",
            title_font_size=20,
            xaxis_title_font_size=14,
            yaxis_title_font_size=14,
            xaxis_type='log', # Set X-axis to logarithmic scale
            yaxis_type='log', # Set Y-axis to logarithmic scale
            xaxis_exponentformat='power',
            xaxis_showexponent='all',
            xaxis_tickformat='e',
            yaxis_exponentformat='power',
            yaxis_showexponent='all',
            yaxis_tickformat='e',
            xaxis_tickangle=90,
            yaxis_tickangle=0,
            width=900,
            height=700,
            hovermode='closest',
            legend_title_text='Groupings'
        )
        st.plotly_chart(fig_comp)
        st.markdown(f'<p style="font-family:\'Times New Roman\', serif; font-size:11px; color:gray; text-align:center;">This plot compares the expression values of genes common to both "Early Embryo" and "{selected_comparison_dataset}" datasets. Points are colored by their Single-Cell group. Vertical dashed lines indicate the maximum "Mean of Geneid Strains" for each Early Embryo group. Standard deviation is not shown.</p>', unsafe_allow_html=True)
        
        if removed_gene_name_comp:
            st.markdown(f'<p style="font-family:\'Times New Roman\', serif; font-size:11px; color:gray; text-align:center;">Note: The gene <b>{removed_gene_name_comp}</b> (Expression: {removed_expr_value_comp:.2f} in {removed_dataset_comp}) was removed to improve plot clarity as it was the single highest outlier across both datasets.</p>', unsafe_allow_html=True)

        if search_gene_comp:
            highlighted_gene_data_merged = merged_df_sorted[merged_df_sorted['gene_common'].str.lower() == search_gene_comp.lower()]
            if not highlighted_gene_data_merged.empty:
                early_embryo_expr = highlighted_gene_data_merged['Mean of Geneid Strains'].iloc[0]
                single_cell_expr = highlighted_gene_data_merged['expression_value'].iloc[0]
                diff_expr = abs(early_embryo_expr - single_cell_expr)
                st.success(f"Expression Difference for **{search_gene_comp}**: |{early_embryo_expr:.2f} (Early Embryo) - {single_cell_expr:.2f} ({selected_comparison_dataset})| = **{diff_expr:.2f}**")
            else:
                st.warning(f"Gene '{search_gene_comp}' not found in common genes for comparison.")


def raw_data_page():
    st.header("Raw Data - Gene Search Tool")
    st.write("This page allows you to view and search different gene expression datasets.")
    
    col_home, col_search_btn_placeholder = st.columns([0.15, 0.85])

    with col_home:
        if st.button("🏠 Back to Home", key="raw_data_back_btn"):
            st.session_state.page = "home"
            st.rerun()
            
    st.markdown("---")

    data_source_option = st.radio(
        "Select Data Source:",
        ("Early Embryo (AnalysisFile2.txt)", "Single-Cell Processed Data"),
        key="data_source_radio"
    )

    current_df = None
    gene_col_name = ""
    expression_col_name = ""
    display_dataset_name = ""

    if data_source_option == "Early Embryo (AnalysisFile2.txt)":
        current_df = df_original
        gene_col_name = 'Gene Name'
        expression_col_name = 'Mean of Geneid Strains'
        display_dataset_name = "Early Embryo"
        st.subheader(f"{display_dataset_name} Data Overview")
        st.dataframe(current_df, use_container_width=True)
    else: # Single-Cell Processed Data
        if not single_cell_dataframes_categorized['Germ'] and not single_cell_dataframes_categorized['Somatic']:
            st.warning("No single-cell processed data files found or loaded. Please ensure they exist in the specified directories.")
            return

        # New dropdowns for category and then dataset
        category_options = []
        if single_cell_dataframes_categorized['Germ']:
            category_options.append("Germ")
        if single_cell_dataframes_categorized['Somatic']:
            category_options.append("Somatic")

        selected_category_raw = st.selectbox(
            "Select Single-Cell Data Category:",
            category_options,
            key="selected_sc_category_raw"
        )
        
        if selected_category_raw:
            current_category_dfs_raw = single_cell_dataframes_categorized[selected_category_raw]
            if not current_category_dfs_raw:
                st.info(f"No datasets available in '{selected_category_raw}' category.")
                return

            selected_single_cell_dataset_raw = st.selectbox(
                f"Select {selected_category_raw} Dataset:",
                list(current_category_dfs_raw.keys()),
                key="single_cell_dataset_select_raw"
            )
            if selected_single_cell_dataset_raw:
                current_df = current_category_dfs_raw[selected_single_cell_dataset_raw]
                gene_col_name = 'gene_common'
                expression_col_name = 'expression_value'
                display_dataset_name = selected_single_cell_dataset_raw
                st.subheader(f"Data Overview: {display_dataset_name}")
                st.dataframe(current_df, use_container_width=True)
            else:
                st.info("Select a single-cell dataset to view its data.")
                return

    if current_df is not None:
        with col_search_btn_placeholder:
            if 'show_search_interface' not in st.session_state:
                st.session_state.show_search_interface = False
            
            if st.button("🔍 Toggle Search Interface", key="toggle_search_btn"):
                st.session_state.show_search_interface = not st.session_state.show_search_interface
        
        if st.session_state.show_search_interface:
            st.subheader("Perform a Search")
            search_type = st.radio(
                "Select search type:",
                ("Search by Gene Name", "Search by Expression Value"),
                key="search_type_radio_dynamic"
            )

            if search_type == "Search by Gene Name":
                gene_name_query = st.text_input(f"Enter {gene_col_name} (e.g., fem-1, rpl-1):", key="gene_name_input_dynamic").strip()
                if gene_name_query:
                    # Search using the appropriate gene column for the selected dataframe
                    search_result = current_df[current_df[gene_col_name].str.lower() == gene_name_query.lower()]
                    if not search_result.empty:
                        st.subheader(f"Results for {gene_name_query} (in {display_dataset_name}):")
                        st.dataframe(search_result, use_container_width=True)
                    else:
                        st.warning(f"No gene found with the name '{gene_name_query}'. Please check the spelling.")
                else:
                    st.info(f"Enter a {gene_col_name} to search.")

            elif search_type == "Search by Expression Value":
                expression_query_str = st.text_input(f"Enter {expression_col_name} Value:", key="expression_input_dynamic").strip()
                if expression_query_str:
                    try:
                        expression_query = float(expression_query_str)
                        st.subheader(f"Closest 20 Genes to {expression_col_name}: {expression_query:.2f}")

                        df_temp_numeric = current_df.copy()
                        df_temp_numeric[f'{expression_col_name}_numeric'] = pd.to_numeric(df_temp_numeric[expression_col_name], errors='coerce')
                        df_temp_numeric.dropna(subset=[f'{expression_col_name}_numeric'], inplace=True)

                        df_temp_numeric['abs_diff_from_query_expr'] = abs(df_temp_numeric[f'{expression_col_name}_numeric'] - expression_query)
                        df_sorted_by_diff = df_temp_numeric.sort_values(by='abs_diff_from_query_expr').reset_index(drop=True)

                        if not df_sorted_by_diff.empty:
                            closest_idx = df_sorted_by_diff[f'{expression_col_name}_numeric'].sub(expression_query).abs().argsort().iloc[0]
                            start_idx = max(0, closest_idx - 10)
                            end_idx = min(len(df_sorted_by_diff), closest_idx + 11)

                            if end_idx - start_idx < 21:
                                if start_idx == 0:
                                    end_idx = min(len(df_sorted_by_diff), 21)
                                elif end_idx == len(df_sorted_by_diff):
                                    start_idx = max(0, len(df_sorted_by_diff) - 21)
                            
                            closest_genes = df_sorted_by_diff.iloc[start_idx:end_idx].sort_values(by=f'{expression_col_name}_numeric')
                            closest_genes = closest_genes.drop(columns=['abs_diff_from_query_expr', f'{expression_col_name}_numeric'])

                            if not closest_genes.empty:
                                st.dataframe(closest_genes, use_container_width=True)
                            else:
                                st.warning("No genes found in the specified range.")
                        else:
                            st.warning("No numeric data available for search in the selected dataset.")

                    except ValueError:
                        st.error(f"Invalid input. Please enter a numeric value for {expression_col_name}.")
                else:
                    st.info(f"Enter an {expression_col_name} value to find closest genes.")


# --- Page Navigation Logic ---
if "page" not in st.session_state:
    st.session_state.page = "home"

if st.session_state.page == "home":
    home_page()
elif st.session_state.page == "original_data_landing":
    original_data_landing_page()
elif st.session_state.page == "processed_data_landing":
    processed_data_landing_page()
elif st.session_state.page == "visualizations":
    visualizations_page()
elif st.session_state.page == "raw_data":
    raw_data_page()
elif st.session_state.page == "comparison_graphs": 
    comparison_page()
