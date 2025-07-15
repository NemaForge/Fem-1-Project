import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import os
import numpy as np

# --- Page Configuration (Should be at the very top of your script) ---
st.set_page_config(
    page_title="Saurish and Xander's Biomart",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="collapsed"
)

# --- Global Data Loading and Constants ---
input_file_path = "AnalysisFile2.txt"

SINGLE_CELL_OUTPUT_DIRECTORY = "" # Keep this as is, assuming files are in the same directory as script

SINGLE_CELL_FILES_DISPLAY_MAP = {
    "Mature sperm": "finalMatureSperm.txt",
    "Meiotic germ cells": "finalMeiotic.txt",
    "Mitotic germ cells": "finalMitotic.txt",
    "Oocytes": "finalOocytes.txt",
    "Spermatids": "finalSpermatids.txt",
    "Spermatocytes": "finalSpermatocytes.txt",
    "Syncitial pachytene spermatocytes": "finalSyncitial.txt"
}

# Define a very small positive number for handling zeros on log scales
# np.finfo(float).eps is the smallest positive representable number for a float
EPSILON = np.finfo(float).eps 

@st.cache_data
def load_original_data(path):
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
def load_single_cell_dataframes():
    single_cell_dfs = {}
    for display_name, filename in SINGLE_CELL_FILES_DISPLAY_MAP.items():
        if SINGLE_CELL_OUTPUT_DIRECTORY:
            file_path = os.path.join(SINGLE_CELL_OUTPUT_DIRECTORY, filename)
        else:
            file_path = filename
            
        try:
            df_sc = pd.read_csv(file_path, sep='\t')
            for col in ['Scaled_TPM', 'group number']:
                if col in df_sc.columns:
                    df_sc[col] = pd.to_numeric(df_sc[col], errors='coerce')
            df_sc.dropna(inplace=True)
            single_cell_dfs[display_name] = df_sc
        except FileNotFoundError:
            st.warning(f"Single-cell file not found: {filename}. It will not be available in the dropdown. Please ensure the file exists in your repository.")
        except Exception as e:
            st.error(f"Error loading single-cell data '{filename}': {e}")
    return single_cell_dfs

df_original = load_original_data(input_file_path)
fem1_data_original = df_original[df_original['Gene Name'] == 'fem-1']

single_cell_dataframes = load_single_cell_dataframes()

regular_dot_size = 5
fem1_dot_size = 10

# --- Helper Functions for Plotting ---
def create_aggregated_hover_data_flexible(df_to_process, gene_col, mean_col, std_dev_col, group_col, round_decimals=3):
    if df_to_process.empty:
        return pd.DataFrame(columns=[mean_col, std_dev_col, group_col, 'Aggregated Hover Text'])

    df_temp = df_to_process.copy()
    
    df_temp[f'{mean_col}_numeric_rounded'] = pd.to_numeric(df_temp[mean_col], errors='coerce').round(round_decimals)
    df_temp[f'{std_dev_col}_numeric_rounded'] = pd.to_numeric(df_temp[std_dev_col], errors='coerce').round(round_decimals)
    
    df_temp.dropna(subset=[f'{mean_col}_numeric_rounded', f'{std_dev_col}_numeric_rounded'], inplace=True)

    grouped = df_temp.groupby([f'{mean_col}_numeric_rounded', f'{std_dev_col}_numeric_rounded', group_col]).apply(
        lambda x: '<br>'.join([
            f"<b>{gene_name}</b> (Mean: {float(mean_val):.2f}, Std Dev: {float(std_dev_val):.2f})"
            for gene_name, mean_val, std_dev_val in zip(x[gene_col], x[mean_col], x[std_dev_col])
        ])
    ).reset_index(name='Aggregated Hover Text')

    grouped.rename(columns={
        f'{mean_col}_numeric_rounded': mean_col,
        f'{std_dev_col}_numeric_rounded': std_dev_col
    }, inplace=True)
    
    return grouped

def plot_gene_expression_set(df_data, fem1_data_subset, plot_title_prefix, gene_col, mean_col, std_dev_col, group_col, group_color_map):
    st.subheader(f"{plot_title_prefix}: All Genes Across All Groups")
    st.write(f"This plot shows the {mean_col} vs {std_dev_col} for all genes, categorized by their assigned group. Hover over a point to see all overlapping Gene Names and their data.")

    # IMPORTANT: Handle zero/negative values for log scale
    df_plot_data = df_data.copy()
    df_plot_data[mean_col] = df_plot_data[mean_col].replace(0, EPSILON) # Replace 0 with a small epsilon
    df_plot_data[std_dev_col] = df_plot_data[std_dev_col].replace(0, EPSILON) # Replace 0 with a small epsilon
    df_plot_data = df_plot_data[df_plot_data[mean_col] > 0] # Ensure all are positive
    df_plot_data = df_plot_data[df_plot_data[std_dev_col] > 0] # Ensure all are positive


    removed_gene_name_plot1 = None
    removed_gene_mean_value_plot1 = None
    df_filtered_for_plot1 = df_plot_data.copy() # Use the cleaned data

    if not df_plot_data.empty and mean_col in df_plot_data.columns: # Check df_plot_data
        df_for_outlier_check = df_plot_data[df_plot_data[gene_col] != 'fem-1'].copy()
        if not df_for_outlier_check.empty:
            max_mean_gene_row = df_for_outlier_check.loc[df_for_outlier_check[mean_col].astype(float).idxmax()]
            removed_gene_name_plot1 = max_mean_gene_row[gene_col]
            removed_gene_mean_value_plot1 = float(max_mean_gene_row[mean_col])
            df_filtered_for_plot1 = df_plot_data[df_plot_data[gene_col] != removed_gene_name_plot1].copy()

    plot1_data_for_hover = create_aggregated_hover_data_flexible(
        df_filtered_for_plot1[df_filtered_for_plot1[gene_col] != 'fem-1'],
        gene_col, mean_col, std_dev_col, group_col
    )

    fig1 = go.Figure()
    unique_groups = sorted(plot1_data_for_hover[group_col].unique(), key=lambda x: int(x))

    for group_name in unique_groups:
        group_df = plot1_data_for_hover[plot1_data_for_hover[group_col] == group_name]
        if not group_df.empty:
            fig1.add_trace(go.Scatter(
                x=group_df[mean_col],
                y=group_df[std_dev_col],
                mode='markers',
                name=f'{group_name}',
                marker=dict(size=regular_dot_size, color=group_color_map.get(str(group_name), 'lightgray')),
                hoverinfo='text',
                text=group_df['Aggregated Hover Text'],
                hovertemplate='%{text}<extra></extra>'
            ))
    
    # Also apply epsilon to fem1_data_subset for plotting
    fem1_data_subset_cleaned = fem1_data_subset.copy()
    if not fem1_data_subset_cleaned.empty:
        fem1_data_subset_cleaned[mean_col] = fem1_data_subset_cleaned[mean_col].replace(0, EPSILON)
        fem1_data_subset_cleaned[std_dev_col] = fem1_data_subset_cleaned[std_dev_col].replace(0, EPSILON)
        # Filter again to ensure positivity after replacement
        fem1_data_subset_cleaned = fem1_data_subset_cleaned[fem1_data_subset_cleaned[mean_col] > 0]
        fem1_data_subset_cleaned = fem1_data_subset_cleaned[fem1_data_subset_cleaned[std_dev_col] > 0]


    if not fem1_data_subset_cleaned.empty and fem1_data_subset_cleaned[gene_col].iloc[0] != removed_gene_name_plot1:
        fem1_data_plot1 = fem1_data_subset_cleaned[fem1_data_subset_cleaned[gene_col] != removed_gene_name_plot1]
        if not fem1_data_plot1.empty:
            fem1_hover_text = (
                f"<b>{fem1_data_plot1[gene_col].iloc[0]}</b>"
                f"<br>{mean_col}: {float(fem1_data_plot1[mean_col].iloc[0]):.2f}"
                f"<br>{std_dev_col}: {float(fem1_data_plot1[std_dev_col].iloc[0]):.2f}"
                f"<br>Group: {fem1_data_plot1[group_col].iloc[0]}"
            )
            fig1.add_trace(go.Scatter(
                x=fem1_data_plot1[mean_col],
                y=fem1_data_plot1[std_dev_col],
                mode='markers',
                marker=dict(size=fem1_dot_size, color='red', symbol='circle', line=dict(width=2, color='DarkRed')),
                name='fem-1',
                hoverinfo='text',
                text=[fem1_hover_text],
                hovertemplate='%{text}<extra></extra>'
            ))
            
    fig1.update_layout(
        title=f'{plot_title_prefix} Genes: {mean_col} vs {std_dev_col}',
        xaxis_title=mean_col,
        yaxis_title=std_dev_col,
        font_family="Times New Roman",
        title_font_size=20,
        xaxis_title_font_size=14,
        yaxis_title_font_size=14,
        xaxis_type='log',
        yaxis_type='log',
        xaxis_exponentformat='power', 
        xaxis_showexponent='all',    
        xaxis_tickformat='e',        
        xaxis_dtick='L1',            
        yaxis_exponentformat='power',
        yaxis_showexponent='all',
        yaxis_tickformat='e',
        yaxis_dtick='L1',
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

    group9_data = df_plot_data[df_plot_data[group_col] == '9'].copy() # Use the cleaned data
    plot2_data_for_hover = create_aggregated_hover_data_flexible(
        group9_data[group9_data[gene_col] != 'fem-1'],
        gene_col, mean_col, std_dev_col, group_col
    )

    fig2 = go.Figure()
    if not plot2_data_for_hover.empty:
        fig2.add_trace(go.Scatter(
            x=plot2_data_for_hover[mean_col],
            y=plot2_data_for_hover[std_dev_col],
            mode='markers',
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
        fig2.add_trace(go.Scatter(
            x=fem1_in_group9[mean_col],
            y=fem1_in_group9[std_dev_col],
            mode='markers',
            marker=dict(size=fem1_dot_size, color='red', symbol='circle', line=dict(width=2, color='DarkRed')),
            name='fem-1',
            hoverinfo='text',
            text=[fem1_hover_text_plot2],
            hovertemplate='%{text}<extra></extra>'
        ))

    fig2.update_layout(
        title=f'{plot_title_prefix}: Group 9 Genes: {mean_col} vs {std_dev_col}',
        xaxis_title=mean_col,
        yaxis_title=std_dev_col,
        font_family="Times New Roman",
        title_font_size=20,
        xaxis_title_font_size=14,
        yaxis_title_font_size=14,
        xaxis_type='log',
        yaxis_type='log',
        xaxis_exponentformat='power',
        xaxis_showexponent='all',
        xaxis_tickformat='e',
        xaxis_dtick='L1',
        yaxis_exponentformat='power',
        yaxis_showexponent='all',
        yaxis_tickformat='e',
        yaxis_dtick='L1',
        xaxis_tickangle=90,
        yaxis_tickangle=0,
        width=900,
        height=600
    )
    st.plotly_chart(fig2)

    st.markdown("---")

    st.subheader(f"{plot_title_prefix}: Genes in Groups 8, 9, and 10")
    st.write(f"This plot shows genes from the top three groups. Hover over a point to see all overlapping Gene Names and their data.")

    selected_groups_data_raw = df_plot_data[df_plot_data[group_col].isin(['8', '9', '10'])].copy() # Use the cleaned data

    removed_gene_name_plot3 = None
    removed_gene_mean_value_plot3 = None

    if not selected_groups_data_raw.empty and mean_col in selected_groups_data_raw.columns:
        df_for_outlier_check_plot3 = selected_groups_data_raw[selected_groups_data_raw[gene_col] != 'fem-1'].copy()
        if not df_for_outlier_check_plot3.empty:
            max_mean_gene_row_plot3 = df_for_outlier_check_plot3.loc[df_for_outlier_check_plot3[mean_col].astype(float).idxmax()]
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
            fig3.add_trace(go.Scatter(
                x=group_df[mean_col],
                y=group_df[std_dev_col],
                mode='markers',
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
            f"<br>{std_dev_col}: {float(fem1_in_selected_groups[std_dev_col].iloc[0]):
