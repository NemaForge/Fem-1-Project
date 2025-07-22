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

# Updated: Nested dictionary for single-cell data files, categorized by Germ and Somatic
SINGLE_CELL_FILES_DISPLAY_MAP = {
    "Germ Cells": {
        "Mature sperm": "finalMatureSperm.txt",
        "Meiotic germ cells": "finalMeiotic.txt",
        "Mitotic germ cells": "finalMitotic.txt",
        "Oocytes": "finalOocytes.txt",
        "Spermatids": "finalSpermatids.txt",
        "Spermatocytes": "finalSpermatocytes.txt",
        "Syncitial pachytene spermatocytes": "finalSyncitial.txt"
    },
    "Somatic Cells": {
        "ADF": "ADF.txt",
        "ADL": "ADL.txt",
        "AFD": "AFD.txt",
        "AIA": "AIA.txt",
        "AIB": "AIB.txt",
        "AIM": "AIM.txt",
        "AIN": "AIN.txt",
        "AIY": "AIY.txt",
        "AIZ": "AIZ.txt",
        "ALA": "ALA.txt",
        "Amphid sheath": "Amphid sheath.txt",
        "Amphid socket": "Amphid socket.txt",
        "Anal muscle": "Anal muscle.txt",
        "Apoptotic germ cells": "Apoptotic germ cells.txt",
        "Arcade cells": "Arcade cells.txt",
        "AS": "AS.txt",
        "ASE": "ASE.txt",
        "ASG": "ASG.txt",
        "ASH": "ASH.txt",
        "ASI": "ASI.txt",
        "ASJ": "ASJ.txt",
        "ASK": "ASK.txt",
        "AUA": "AUA.txt",
        "AVA": "AVA.txt",
        "AVB": "AVB.txt",
        "AVD": "AVD.txt",
        "AVE": "AVE.txt",
        "AVF": "AVF.txt",
        "AVH": "AVH.txt",
        "AVJ": "AVJ.txt",
        "AVK": "AVK.txt",
        "AVL": "AVL.txt",
        "AVM": "AVM.txt",
        "AWA": "AWA.txt",
        "AWB": "AWB.txt",
        "AWC": "AWC.txt",
        "BAG": "BAG.txt",
        "BDU": "BDU.txt",
        "Body wall muscle anterior": "Body wall muscle anterior.txt",
        "Body wall muscle middle": "Body wall muscle middle.txt",
        "Body wall muscle posterior": "Body wall muscle posterior.txt",
        "CAN": "CAN.txt",
        "cat-4(+)ptps-1(+) intestine anterior": "cat-4(+)ptps-1(+) intestine anterior.txt",
        "CEP_ADE_PDE": "CEP_ADE_PDE.txt",
        "Cephalic and inner labial socket": "Cephalic and inner labial socket.txt",
        "Cephalic sheath": "Cephalic sheath.txt",
        "Coelomocytes": "Coelomocytes.txt",
        "DA_VA": "DA_VA.txt",
        "DB_VB": "DB_VB.txt",
        "Differentiated germ": "Differentiated germ.txt",
        "Distal tip": "Distal tip.txt",
        "Dorsal uterine cell": "Dorsal uterine cell.txt",
        "DVA": "DVA.txt",
        "DVB": "DVB.txt",
        "DVC": "DVC.txt",
        "e1_e3 (pharyngeal epithelium)": "e1_e3 (pharyngeal epithelium).txt",
        "e2 (pharyngeal epithelium)": "e2 (pharyngeal epithelium).txt",
        "Embryonic cells": "Embryonic cells.txt",
        "Excretory cells": "Excretory cells.txt",
        "Excretory duct": "Excretory duct.txt",
        "Excretory gland": "Excretory gland.txt",
        "g1A (pharyngeal gland)": "g1A (pharyngeal gland).txt",
        "g1P (pharyngeal gland)": "g1P (pharyngeal gland).txt",
        "g2 (pharyngeal gland)": "g2 (pharyngeal gland).txt",
        "Glia_1": "Glia_1.txt",
        "Glia_2": "Glia_2.txt",
        "Glia_3": "Glia_3.txt",
        "Glia_4": "Glia_4.txt",
        "GLR": "GLR.txt",
        "hmc": "hmc.txt",
        "HSN": "HSN.txt",
        "hyp4_to_hyp6": "hyp4_to_hyp6.txt",
        "hyp7 (hypodermis)": "hyp7 (hypodermis).txt",
        "Hypodermis head": "Hypodermis head.txt",
        "Hypodermis tail": "Hypodermis tail.txt",
        "I1": "I1.txt",
        "I2": "I2.txt",
        "I4": "I4.txt",
        "I5": "I5.txt",
        "IL1": "IL1.txt",
        "IL2": "IL2.txt",
        "Intestinal-rectal valve": "Intestinal-rectal valve.txt",
        "Intestine anterior": "Intestine anterior.txt",
        "Intestine middle": "Intestine middle.txt",
        "Intestine posterior": "Intestine posterior.txt",
        "LUA": "LUA.txt",
        "M1": "M1.txt",
        "M2": "M2.txt",
        "M3": "M3.txt",
        "M4": "M4.txt",
        "M5": "M5.txt",
        "Marginal cells": "Marginal cells.txt",
        "MI": "MI.txt",
        "NSM": "NSM.txt",
        "OLL": "OLL.txt",
        "OLQ": "OLQ.txt",
        "PDA": "PDA.txt",
        "PDB": "PDB.txt",
        "PHA_PHB": "PHA_PHB.txt",
        "Phasmid sheath": "Phasmid sheath.txt",
        "Phasmid socket": "Phasmid socket.txt",
        "PHC": "PHC.txt",
        "PLM_ALM": "PLM_ALM.txt",
        "PLN": "PLN.txt",
        "pm1_pm2 (pharyngeal muscle)": "pm1_pm2 (pharyngeal muscle).txt",
        "pm3_pm4_pm5 (pharyngeal muscle)": "pm3_pm4_pm5 (pharyngeal muscle).txt",
        "pm6_pm7 (pharyngeal muscle)": "pm6_pm7 (pharyngeal muscle).txt",
        "PVC": "PVC.txt",
        "PVD_FLP": "PVD_FLP.txt",
        "PVM": "PVM.txt",
        "PVN": "PVN.txt",
        "PVP": "PVP.txt",
        "PVQ": "PVQ.txt",
        "PVR": "PVR.txt",
        "PVT": "PVT.txt",
        "PVW": "PVW.txt",
        "Rectal gland": "Rectal gland.txt",
        "RIA": "RIA.txt",
        "RIB": "RIB.txt",
        "RIC": "RIC.txt",
        "RID": "RID.txt",
        "RIF": "RIF.txt",
        "RIG": "RIG.txt",
        "RIH": "RIH.txt",
        "RIM": "RIM.txt",
        "RIP": "RIP.txt",
        "RIR": "RIR.txt",
        "RIS": "RIS.txt",
        "RIV": "RIV.txt",
        "RMD": "RMD.txt",
        "RMD_LR": "RMD_LR.txt",
        "RME": "RME.txt",
        "RMF": "RMF.txt",
        "RMG": "RMG.txt",
        "RMH": "RMH.txt",
        "SAA": "SAA.txt",
        "SAB": "SAB.txt",
        "SDQ": "SDQ.txt",
        "Seam cells (bus+)": "Seam cells (bus+).txt",
        "Seam cells (grd+)": "Seam cells (grd+).txt",
        "Seminal vesicle (male)": "Seminal vesicle (male).txt",
        "sh1 (gonadal sheath distal)": "sh1 (gonadal sheath distal).txt",
        "sh2 (gonadal sheath distal)": "sh2 (gonadal sheath distal).txt",
        "sh3_sh4 (gonadal sheath proximal)": "sh3_sh4 (gonadal sheath proximal).txt",
        "sh5 (gonadal sheath proximal)": "sh5 (gonadal sheath proximal).txt",
        "SIA (1)": "SIA (1).txt",
        "SIA": "SIA.txt",
        "SIB": "SIB.txt",
        "SMB": "SMB.txt",
        "SMD": "SMD.txt",
        "Spermatheca bag distal": "Spermatheca bag distal.txt",
        "Spermatheca bag proximal": "Spermatheca bag proximal.txt",
        "Spermatheca neck distal": "Spermatheca neck distal.txt",
        "Spermatheca neck proximal": "Spermatheca neck proximal.txt",
        "Spermatheca-Uterine junction": "Spermatheca-Uterine junction.txt",
        "Unassigned hypodermis_gonadal sheath": "Unassigned hypodermis_gonadal sheath.txt",
        "Unassigned sex-specific muscle": "Unassigned sex-specific muscle.txt",
        "Unassigned sheath cells": "Unassigned sheath cells.txt",
        "Unassigned uterine cells": "Unassigned uterine cells.txt",
        "URA": "URA.txt",
        "URB": "URB.txt",
        "URX_AQR_PQR": "URX_AQR_PQR.txt",
        "URY": "URY.txt",
        "Uterine muscle": "Uterine muscle.txt",
        "Uterine seam cells": "Uterine seam cells.txt",
        "Uterine toroid": "Uterine toroid.txt",
        "Uterine-vulval cells": "Uterine-vulval cells.txt",
        "VC": "VC.txt",
        "VC_4_5": "VC_4_5.txt",
        "VD_DD": "VD_DD.txt",
        "vm1 (vulval muscle)": "vm1 (vulval muscle).txt",
        "vm2 (vulval muscle)": "vm2 (vulval muscle).txt",
        "Vulval cells": "Vulval cells.txt"
    }
}

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
def load_single_cell_dataframes_original_structure():
    """
    Loads single-cell data, organized by category (Germ/Somatic).
    Standardizes column names to 'gene name', 'Scaled_TPM', 'group number'.
    """
    all_single_cell_dfs = {}
    for category, files_map in SINGLE_CELL_FILES_DISPLAY_MAP.items():
        category_dfs = {}
        for display_name, filename in files_map.items():
            file_path = filename # Assumed in root
            
            try:
                df_sc = pd.read_csv(file_path, sep='\t')
                # Standardize column names to match existing plotting functions
                if 'gene_short_name' in df_sc.columns:
                    df_sc.rename(columns={'gene_short_name': 'gene name'}, inplace=True)
                if 'scaled_TPM' in df_sc.columns:
                    df_sc.rename(columns={'scaled_TPM': 'Scaled_TPM'}, inplace=True)
                if 'group_number' in df_sc.columns:
                    df_sc.rename(columns={'group_number': 'group number'}, inplace=True)

                for col in ['Scaled_TPM', 'group number']:
                    if col in df_sc.columns:
                        df_sc[col] = pd.to_numeric(df_sc[col], errors='coerce')
                df_sc.dropna(inplace=True)
                category_dfs[display_name] = df_sc
            except FileNotFoundError:
                st.warning(f"Single-cell file not found: '{filename}' in category '{category}'. It will not be available in the dropdown. Please ensure the file exists in your repository's root directory or the specified path.")
                continue # Skip to the next file if this one isn't found
            except Exception as e:
                st.error(f"Error loading single-cell data '{filename}' in category '{category}': {e}")
        all_single_cell_dfs[category] = category_dfs
    return all_single_cell_dfs

# Load data globally
df_original = load_original_data(input_file_path)
fem1_data_original = df_original[df_original['Gene Name'] == 'fem-1']

single_cell_dataframes = load_single_cell_dataframes_original_structure()

# Helper function to get sorted unique groups for comparison legend
def get_sorted_groups(df, group_col):
    if df.empty or group_col not in df.columns:
        return []
    try:
        # Convert to numeric first for proper sorting, then back to string
        # Handle cases where group names might not be purely numeric (e.g., 'Group A')
        numeric_groups = []
        non_numeric_groups = []
        for x in df[group_col].dropna().astype(str).unique():
            if x.isdigit():
                numeric_groups.append(int(x))
            else:
                non_numeric_groups.append(x)
        
        sorted_numeric = [str(g) for g in sorted(numeric_groups)]
        sorted_non_numeric = sorted(non_numeric_groups)
        return sorted_numeric + sorted_non_numeric
    except:
        # Fallback for complex non-numeric group names or unexpected data
        return sorted(df[group_col].dropna().astype(str).unique())


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

    removed_gene_name_plot1 = None
    removed_gene_mean_value_plot1 = None
    df_filtered_for_plot1 = df_data.copy()

    if not df_data.empty and mean_col in df_data.columns:
        df_for_outlier_check = df_data[df_data[gene_col] != 'fem-1'].copy()
        if not df_for_outlier_check.empty and pd.to_numeric(df_for_outlier_check[mean_col], errors='coerce').notna().any():
            max_mean_gene_row = df_for_outlier_check.loc[pd.to_numeric(df_for_outlier_check[mean_col], errors='coerce').idxmax()]
            removed_gene_name_plot1 = max_mean_gene_row[gene_col]
            removed_gene_mean_value_plot1 = float(max_mean_gene_row[mean_col])
            df_filtered_for_plot1 = df_data[df_data[gene_col] != removed_gene_name_plot1].copy()

    plot1_data_for_hover = create_aggregated_hover_data_flexible(
        df_filtered_for_plot1[df_filtered_for_plot1[gene_col] != 'fem-1'],
        gene_col, mean_col, std_dev_col, group_col
    )

    fig1 = go.Figure()
    unique_groups = get_sorted_groups(plot1_data_for_hover, group_col)

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

    if not fem1_data_subset.empty and (removed_gene_name_plot1 is None or fem1_data_subset[gene_col].iloc[0] != removed_gene_name_plot1):
        fem1_data_plot1 = fem1_data_subset[fem1_data_subset[gene_col] != removed_gene_name_plot1]
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
        xaxis_type='log', # Set X-axis to logarithmic scale
        yaxis_type='log', # Set Y-axis to logarithmic scale
        xaxis_exponentformat='power', # Display exponents
        xaxis_showexponent='all',    # Show exponent for all ticks
        xaxis_tickformat='e',        # Use scientific notation for ticks
        yaxis_exponentformat='power',
        yaxis_showexponent='all',
        yaxis_tickformat='e',
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
        xaxis_type='log', # Set X-axis to logarithmic scale
        yaxis_type='log', # Set Y-axis to logarithmic scale
        xaxis_exponentformat='power',
        yaxis_showexponent='all',
        yaxis_tickformat='e',
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
    unique_groups_plot3 = get_sorted_groups(plot3_data_for_hover, group_col)

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
            f"<br>{std_dev_col}: {float(fem1_in_selected_groups[std_dev_col].iloc[0]):.2f}"
            f"<br>Group: {fem1_in_selected_groups[group_col].iloc[0]}"
        )
        fig3.add_trace(go.Scatter(
            x=fem1_in_selected_groups[mean_col],
            y=fem1_in_selected_groups[std_dev_col],
            mode='markers',
            marker=dict(size=fem1_dot_size, color='red', symbol='circle', line=dict(width=2, color='DarkRed')),
            name='fem-1',
            hoverinfo='text',
            text=[fem1_hover_text_plot3],
            hovertemplate='%{text}<extra></extra>'
        ))

    fig3.update_layout(
        title=f'{plot_title_prefix}: Groups 8, 9, and 10 Genes: {mean_col} vs {std_dev_col}',
        xaxis_title=mean_col,
        yaxis_title=std_dev_col,
        font_family="Times New Roman",
        title_font_size=20,
        xaxis_title_font_size=14,
        yaxis_title_font_size=14,
        xaxis_type='log', # Log scale X-axis
        yaxis_type='log', # Log scale Y-axis
        xaxis_exponentformat='power',
        yaxis_showexponent='all',
        yaxis_tickformat='e',
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

    for group_name_sc in get_sorted_groups(df_sorted_by_tpm, group_col):
        group_df_sc = df_sorted_by_tpm[df_sorted_by_tpm[group_col] == group_name_sc]
        if not group_df_sc.empty:
            fig_sc.add_trace(go.Scatter(
                x=group_df_sc[tpm_col],
                y=[0] * len(group_df_sc),
                mode='markers',
                name=f'Group {group_name_sc}',
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


# --- New Home Page Function ---
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
            box-shadow: 0 12px 35px rgba(0,0,0,0.4) !important;
            background: linear-gradient(45deg, #d0d0d0 0%, #b0b0b0 100%) !important; /* Slightly darker on hover */
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
            max-width: 400px; /* Constrain max-width of container for better centering effect */
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


# --- Processed Data Landing Page (Currently Blank) ---
def processed_data_landing_page():
    st.header("Processed Data Analysis")
    st.write("This section is dedicated to processed data, such as single-cell RNA sequencing analysis.")
    st.info("Functionality for processed data will be added here soon!")

    if st.button("🏠 Back to Home", key="processed_landing_back_home"):
        st.session_state.page = "home"
        st.rerun()

    st.markdown("---")
    # You can add more buttons or content here as you develop the processed data features
    # Example:
    # if st.button("Explore Single-Cell Data (Coming Soon!)"):
    #        st.write("Stay tuned for interactive single-cell analysis tools!")


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
        ("Early Embryo Data", "Single-Cell Data"),
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
             for i, g in enumerate(get_sorted_groups(df_original, 'Group'))}
        )

    elif visualization_type == "Single-Cell Data":
        st.subheader("Single-Cell Data Visualizations")
        st.write("Select a single-cell dataset to view its gene expression patterns.")

        if not single_cell_dataframes:
            st.warning("No single-cell processed data files found or loaded. Please ensure they exist in the specified directory.")
            return

        # First dropdown for category (Germ Cells or Somatic Cells)
        category_options = list(single_cell_dataframes.keys())
        selected_category = st.selectbox(
            "Choose a Data Category:",
            category_options,
            key="selected_sc_viz_category"
        )

        current_sc_df = None
        if selected_category and single_cell_dataframes[selected_category]:
            # Second dropdown for specific dataset within the chosen category
            dataset_options = list(single_cell_dataframes[selected_category].keys())
            selected_single_cell_dataset = st.selectbox(
                f"Choose a {selected_category} Dataset:",
                dataset_options,
                key="selected_sc_viz_dataset"
            )

            if selected_single_cell_dataset:
                current_sc_df = single_cell_dataframes[selected_category][selected_single_cell_dataset]
                # fem1_data_sc still uses 'gene name'
                fem1_data_sc = current_sc_df[current_sc_df['gene name'] == 'fem-1']

                sc_group_color_map = {
                    str(g): px.colors.qualitative.D3[i % len(px.colors.qualitative.D3)]
                    for i, g in enumerate(range(1, 11)) # Assuming groups 1-10 for coloring
                }

                plot_single_cell_expression_set(
                    current_sc_df,
                    fem1_data_sc,
                    f"{selected_category}: {selected_single_cell_dataset}", # Update plot title prefix
                    'gene name',
                    'Scaled_TPM',
                    'group number',
                    sc_group_color_map
                )
            else:
                st.info("Select a single-cell dataset to view its data.")
        else:
            st.info("Select a single-cell data category to view its datasets.")


# --- New Comparison Graphs Page ---
def comparison_page():
    st.header("Gene Expression Comparison Graphs")
    st.write("Compare gene expression between the 'Early Embryo' dataset and a selected single-cell dataset.")

    if st.button("🏠 Back to Home", key="comp_back_home"):
        st.session_state.page = "home"
        st.rerun()
    
    st.markdown("---")

    if not single_cell_dataframes:
        st.warning("No single-cell processed data files found or loaded to perform comparisons.")
        return

    # First dropdown for category (Germ Cells or Somatic Cells)
    category_options_comp = list(single_cell_dataframes.keys())
    selected_category_comp = st.selectbox(
        "Choose a Data Category for Comparison:",
        category_options_comp,
        key="selected_comparison_category"
    )

    sc_df_for_comparison = None
    selected_comparison_dataset_name = None # To store the name for display

    if selected_category_comp and single_cell_dataframes[selected_category_comp]:
        # Second dropdown for specific dataset within the chosen category
        dataset_options_comp = list(single_cell_dataframes[selected_category_comp].keys())
        selected_comparison_dataset_name = st.selectbox(
            f"Choose a {selected_category_comp} Dataset for Comparison:",
            dataset_options_comp,
            key="selected_comparison_dataset"
        )

        if selected_comparison_dataset_name:
            sc_df_for_comparison = single_cell_dataframes[selected_category_comp][selected_comparison_dataset_name]
        else:
            st.info("Select a single-cell dataset to perform comparison.")
            return
    else:
        st.info("Select a single-cell data category to perform comparison.")
        return

    if sc_df_for_comparison is None: # Double check if a dataframe was actually selected
        return

    df_original_renamed = df_original.rename(columns={'Gene Name': 'gene_common'})
    sc_df_for_comparison_renamed = sc_df_for_comparison.rename(columns={'gene name': 'gene_common'})

    merged_df = pd.merge(
        df_original_renamed[['gene_common', 'Mean of Geneid Strains', 'Group']],
        sc_df_for_comparison_renamed[['gene_common', 'Scaled_TPM', 'group number']],
        on='gene_common',
        how='inner'
    )
    
    merged_df_sorted = merged_df.sort_values(by='gene_common').reset_index(drop=True)

    if merged_df_sorted.empty:
        st.warning(f"No common genes found between 'Early Embryo' and '{selected_comparison_dataset_name}' for comparison.")
        return

    # Determine and remove the highest outlier gene to improve plot clarity
    removed_gene_name_comp = None
    removed_expr_value_comp = None
    removed_dataset_comp = None

    # Check for outliers in both datasets' expression columns
    max_early_embryo_expr = merged_df_sorted['Mean of Geneid Strains'].max()
    max_single_cell_expr = merged_df_sorted['Scaled_TPM'].max()

    if pd.notna(max_early_embryo_expr) and (pd.notna(max_single_cell_expr) and max_early_embryo_expr >= max_single_cell_expr):
        outlier_row = merged_df_sorted.loc[merged_df_sorted['Mean of Geneid Strains'].idxmax()]
        if outlier_row['gene_common'] != 'fem-1':
            removed_gene_name_comp = outlier_row['gene_common']
            removed_expr_value_comp = outlier_row['Mean of Geneid Strains']
            removed_dataset_comp = "Early Embryo"
    elif pd.notna(max_single_cell_expr):
        outlier_row = merged_df_sorted.loc[merged_df_sorted['Scaled_TPM'].idxmax()]
        if outlier_row['gene_common'] != 'fem-1':
            removed_gene_name_comp = outlier_row['gene_common']
            removed_expr_value_comp = outlier_row['Scaled_TPM']
            removed_dataset_comp = selected_comparison_dataset_name
            
    if removed_gene_name_comp:
        merged_df_sorted = merged_df_sorted[merged_df_sorted['gene_common'] != removed_gene_name_comp].copy()


    st.markdown("---")
    st.subheader("Search Gene for Highlight & Difference")
    search_gene_comp = st.text_input("Enter Gene Name to Highlight (e.g., fem-1):", key="search_gene_comp").strip()

    fig_comp = go.Figure()

    # Get sorted groups for the selected single-cell dataset dynamically
    single_cell_groups_sorted_comp = get_sorted_groups(sc_df_for_comparison, 'group number')

    # Create color map for Single-Cell groups (10 colors)
    single_cell_colors = px.colors.qualitative.D3
    single_cell_color_map = {str(g): single_cell_colors[i % len(single_cell_colors)] for i, g in enumerate(range(1, 11))}

    # Add traces for actual data points
    for sc_group in single_cell_groups_sorted_comp:
        subset_sc_df = merged_df_sorted[merged_df_sorted['group number'] == float(sc_group)] # Ensure type consistency
        if not subset_sc_df.empty:
            fig_comp.add_trace(go.Scatter(
                x=subset_sc_df['Mean of Geneid Strains'],
                y=subset_sc_df['Scaled_TPM'],
                mode='markers',
                name=f'SC Group {sc_group}',
                marker=dict(
                    size=regular_dot_size,
                    color=single_cell_color_map.get(str(sc_group), 'gray'),
                    symbol='circle'
                ),
                hoverinfo='text',
                text=[
                    f"<b>{row['gene_common']}</b><br>EE Expr: {row['Mean of Geneid Strains']:.2f} (Group: {row['Group']})"
                    f"<br>SC Expr: {row['Scaled_TPM']:.2f} (Group: {row['group number']})"
                    f"<br>HIGHLIGHTED"
                ],
                hovertemplate='%{text}<extra></extra>'
            ))
    
    # Add vertical lines for Early Embryo group max expression
    shapes = []
    # Ensure early_embryo_groups_sorted is defined or derived here if not global
    early_embryo_groups_sorted_comp = get_sorted_groups(df_original, 'Group')
    for ee_group in early_embryo_groups_sorted_comp:
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
                        layer='below'
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
                    y=[row['Scaled_TPM']],
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
                        f"<br>SC Expr: {row['Scaled_TPM']:.2f} (Group: {row['group number']})"
                        f"<br>HIGHLIGHTED"
                    ],
                    hovertemplate='%{text}<extra></extra>'
                ))

    fig_comp.update_layout(
        title=f'Gene Expression Comparison: Early Embryo vs. {selected_comparison_dataset_name}',
        xaxis_title='Early Embryo: Mean of Geneid Strains',
        yaxis_title=f'{selected_comparison_dataset_name}: Scaled TPM',
        font_family="Times New Roman",
        title_font_size=20,
        xaxis_title_font_size=14,
        yaxis_title_font_size=14,
        xaxis_type='log',
        yaxis_type='log',
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
    st.markdown(f'<p style="font-family:\'Times New Roman\', serif; font-size:11px; color:gray; text-align:center;">This plot compares the expression values of genes common to both "Early Embryo" and "{selected_comparison_dataset_name}" datasets. Points are colored by their Single-Cell group. Vertical dashed lines indicate the maximum "Mean of Geneid Strains" for each Early Embryo group. Standard deviation is not shown.</p>', unsafe_allow_html=True)
    
    if removed_gene_name_comp:
        st.markdown(f'<p style="font-family:\'Times New Roman\', serif; font-size:11px; color:gray; text-align:center;">Note: The gene <b>{removed_gene_name_comp}</b> (Expression: {removed_expr_value_comp:.2f} in {removed_dataset_comp}) was removed to improve plot clarity as it was the single highest outlier across both datasets.</p>', unsafe_allow_html=True)

    if search_gene_comp:
        highlighted_gene_data_merged = merged_df_sorted[merged_df_sorted['gene_common'].str.lower() == search_gene_comp.lower()]
        if not highlighted_gene_data_merged.empty:
            early_embryo_expr = highlighted_gene_data_merged['Mean of Geneid Strains'].iloc[0]
            single_cell_expr = highlighted_gene_data_merged['Scaled_TPM'].iloc[0]
            diff_expr = abs(early_embryo_expr - single_cell_expr)
            st.success(f"Expression Difference for **{search_gene_comp}**: |{early_embryo_expr:.2f} (Early Embryo) - {single_cell_expr:.2f} ({selected_comparison_dataset_name})| = **{diff_expr:.2f}**")
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
    mean_col_name = ""
    display_dataset_name = ""

    if data_source_option == "Early Embryo (AnalysisFile2.txt)":
        current_df = df_original
        gene_col_name = 'Gene Name'
        mean_col_name = 'Mean of Geneid Strains'
        display_dataset_name = "Early Embryo"
        st.subheader("Early Embryo Data Overview")
        st.dataframe(current_df, use_container_width=True)
    else: # Single-Cell Processed Data
        if not single_cell_dataframes:
            st.warning("No single-cell processed data files found or loaded. Please ensure they exist in the specified directory.")
            return

        # First dropdown for category
        category_options_raw = list(single_cell_dataframes.keys())
        selected_category_raw = st.selectbox(
            "Select Data Category:",
            category_options_raw,
            key="raw_data_category_select"
        )

        if selected_category_raw and single_cell_dataframes[selected_category_raw]:
            # Second dropdown for specific dataset
            dataset_options_raw = list(single_cell_dataframes[selected_category_raw].keys())
            selected_single_cell_dataset = st.selectbox(
                f"Select {selected_category_raw} Dataset:",
                dataset_options_raw,
                key="single_cell_dataset_select"
            )
            if selected_single_cell_dataset:
                current_df = single_cell_dataframes[selected_category_raw][selected_single_cell_dataset]
                gene_col_name = 'gene name'
                mean_col_name = 'Scaled_TPM'
                display_dataset_name = selected_single_cell_dataset
                st.subheader(f"Data Overview: {display_dataset_name}")
                st.dataframe(current_df, use_container_width=True)
            else:
                st.info("Select a single-cell dataset to view its data.")
                return
        else:
            st.info("Select a single-cell data category to view its datasets.")
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
                    search_result = current_df[current_df[gene_col_name].str.lower() == gene_name_query.lower()]
                    if not search_result.empty:
                        st.subheader(f"Results for {gene_name_query} (in {display_dataset_name}):")
                        st.dataframe(search_result, use_container_width=True)
                    else:
                        st.warning(f"No gene found with the name '{gene_name_query}'. Please check the spelling.")
                else:
                    st.info(f"Enter a {gene_col_name} to search.")

            elif search_type == "Search by Expression Value":
                mean_query_str = st.text_input(f"Enter {mean_col_name} Value:", key="mean_input_dynamic").strip()
                if mean_query_str:
                    try:
                        mean_query = float(mean_query_str)
                        st.subheader(f"Closest 20 Genes to {mean_col_name}: {mean_query:.2f}")

                        df_temp_numeric = current_df.copy()
                        df_temp_numeric[f'{mean_col_name}_numeric'] = pd.to_numeric(df_temp_numeric[mean_col_name], errors='coerce')
                        df_temp_numeric.dropna(subset=[f'{mean_col_name}_numeric'], inplace=True)

                        df_temp_numeric['abs_diff_from_query_mean'] = abs(df_temp_numeric[f'{mean_col_name}_numeric'] - mean_query)
                        df_sorted_by_diff = df_temp_numeric.sort_values(by='abs_diff_from_query_mean').reset_index(drop=True)

                        if not df_sorted_by_diff.empty:
                            closest_idx = df_sorted_by_diff[f'{mean_col_name}_numeric'].sub(mean_query).abs().argsort().iloc[0]
                            start_idx = max(0, closest_idx - 10)
                            end_idx = min(len(df_sorted_by_diff), closest_idx + 11)

                            if end_idx - start_idx < 21:
                                if start_idx == 0:
                                    end_idx = min(len(df_sorted_by_diff), 21)
                                elif end_idx == len(df_sorted_by_diff):
                                    start_idx = max(0, len(df_sorted_by_diff) - 21)
                            
                            closest_genes = df_sorted_by_diff.iloc[start_idx:end_idx].sort_values(by=f'{mean_col_name}_numeric')
                            closest_genes = closest_genes.drop(columns=['abs_diff_from_query_mean', f'{mean_col_name}_numeric'])

                            if not closest_genes.empty:
                                st.dataframe(closest_genes, use_container_width=True)
                            else:
                                st.warning("No genes found in the specified range.")
                        else:
                            st.warning("No numeric data available for search in the selected dataset.")

                    except ValueError:
                        st.error(f"Invalid input. Please enter a numeric value for {mean_col_name}.")
                else:
                    st.info(f"Enter an {mean_col_name} value to find closest genes.")


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
