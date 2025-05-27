import os
import pandas as pd
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import warnings
from scipy.spatial.distance import pdist, squareform # For pie radius calculation

# --- Matplotlib Defaults & Font Settings ---
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['legend.fontsize'] = 'small'

# --- HELPER FUNCTIONS ---

def get_consistent_color_map(unique_cell_types_list, predefined_color_map=None, default_palette='tab20'):
    if predefined_color_map is None: predefined_color_map = {}
    all_unique_types_flat = []
    for sublist in unique_cell_types_list:
        for ct in sublist:
            if ct is not None and str(ct).strip() != "":
                 all_unique_types_flat.append(str(ct))
    all_unique_types = sorted(list(set(all_unique_types_flat)))
    final_color_map = predefined_color_map.copy()
    try:
        cmap_obj = plt.cm.get_cmap(default_palette)
        num_palette_colors = cmap_obj.N if cmap_obj.N < 256 else 20 
        palette_colors = [mcolors.to_hex(cmap_obj(i)) for i in np.linspace(0, 1, num_palette_colors)]
    except ValueError: 
        warnings.warn(f"Palette '{default_palette}' not found. Using fallback colors (tab20-like).")
        palette_colors = [
            "#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c", "#98df8a",
            "#d62728", "#ff9896", "#9467bd", "#c5b0d5", "#8c564b", "#c49c94",
            "#e377c2", "#f7b6d2", "#7f7f7f", "#c7c7c7", "#bcbd22", "#dbdb8d",
            "#17becf", "#9edae5"
        ]
    color_idx = 0
    for ct in all_unique_types:
        if ct not in final_color_map:
            if color_idx < len(palette_colors):
                final_color_map[ct] = palette_colors[color_idx]; color_idx += 1
            else:
                final_color_map[ct] = mcolors.to_hex(np.random.rand(3,))
    if "Unknown" in all_unique_types and "Unknown" not in final_color_map:
        final_color_map["Unknown"] = "#CCCCCC" 
    elif "Unknown" in final_color_map and final_color_map["Unknown"] is None:
         final_color_map["Unknown"] = "#CCCCCC"
    return final_color_map

def _get_spatial_image_and_scale(adata):
    if 'spatial' not in adata.uns or not adata.uns['spatial']: return None, 1.0
    if not isinstance(adata.uns['spatial'], dict):
         warnings.warn("adata.uns['spatial'] is not a dictionary."); return None, 1.0
    spatial_key = None
    if len(adata.uns['spatial'].keys()) == 1:
        spatial_key = list(adata.uns['spatial'].keys())[0]
    else: 
        for key in adata.uns['spatial'].keys():
            if isinstance(adata.uns['spatial'][key], dict) and \
               'images' in adata.uns['spatial'][key] and \
               'scalefactors' in adata.uns['spatial'][key]:
                spatial_key = key; break
        if spatial_key is None and adata.uns['spatial'].keys():
            spatial_key = list(adata.uns['spatial'].keys())[0]
    if spatial_key is None or not isinstance(adata.uns['spatial'].get(spatial_key), dict) : return None, 1.0
    img_data = adata.uns['spatial'][spatial_key]
    image, scale_factor = None, 1.0
    if 'images' in img_data and 'scalefactors' in img_data:
        if 'hires' in img_data['images'] and img_data['images']['hires'] is not None and \
           'tissue_hires_scalef' in img_data['scalefactors'] and img_data['scalefactors']['tissue_hires_scalef'] is not None:
            image = img_data['images']['hires'].copy(); scale_factor = img_data['scalefactors']['tissue_hires_scalef']
        elif 'lowres' in img_data['images'] and img_data['images']['lowres'] is not None and \
             'tissue_lowres_scalef' in img_data['scalefactors'] and img_data['scalefactors']['tissue_lowres_scalef'] is not None:
            image = img_data['images']['lowres'].copy(); scale_factor = img_data['scalefactors']['tissue_lowres_scalef']
    if image is not None:
        if image.max() > 1.0 and image.max() <= 255.0 : image = (image / 255.0).astype(np.float32)
        elif image.max() > 255: warnings.warn(f"Image max value {image.max()} > 255.")
    return image, scale_factor

def _add_legend(ax, color_map, plotted_cell_types, title="Cell Types", ncol=1, loc='auto', 
                bbox_to_anchor=None, frameon=True, legend_fontsize=None):
    if not plotted_cell_types: return
    patches = [mpatches.Patch(color=color_map.get(str(ct), '#CCCCCC'), label=str(ct)) 
               for ct in sorted(list(plotted_cell_types)) if ct in color_map or str(ct) == "Unknown"]
    if not patches: return
    final_legend_fontsize = legend_fontsize if legend_fontsize is not None else plt.rcParams['legend.fontsize']
    bbox_params = {}
    if loc == 'auto':
        if len(patches) > 10: loc = 'center left'; bbox_params['bbox_to_anchor'] = (1.02, 0.5) if bbox_to_anchor is None else bbox_to_anchor
        else: loc = 'best'
        ax.legend(handles=patches, title=title, loc=loc, frameon=frameon, ncol=ncol, fontsize=final_legend_fontsize, **bbox_params)
    else: 
        ax.legend(handles=patches, title=title, loc=loc, bbox_to_anchor=bbox_to_anchor, frameon=frameon, ncol=ncol, fontsize=final_legend_fontsize)

# --- Plotting Functions ---

def plot_spatial_discrete(
    adata, annotation_key, color_map=None, ax=None, fig_size=(7, 6), 
    spot_size=30, img_alpha=1.0, show_legend=True, legend_title=None,
    legend_ncol=1, legend_fontsize=None, title="Spatial Annotation",
    output_path=None, show=False, spot_edgecolor='none', 
    spot_linewidth=0.1, default_palette='tab20' 
):
    if annotation_key not in adata.obs:
        raise ValueError(f"Annotation key '{annotation_key}' not found in adata.obs.")
    if "spatial" not in adata.obsm:
        raise ValueError("Spatial coordinates 'spatial' not found in adata.obsm.")

    spatial_coords_orig = adata.obsm["spatial"].copy()
    annotations = adata.obs[annotation_key].astype(str) 
    unique_annotations = sorted(annotations.unique())

    if color_map is None:
        color_map = get_consistent_color_map([unique_annotations], default_palette=default_palette)
    if legend_title is None:
        legend_title = annotation_key.replace('_', ' ').title()
    if ax is None:
        fig, ax = plt.subplots(figsize=fig_size)
    else:
        fig = ax.get_figure()

    img, img_scale_factor = _get_spatial_image_and_scale(adata)
    original_xlim, original_ylim = ax.get_xlim(), ax.get_ylim()
    coords_to_plot = spatial_coords_orig # Plotting in full-resolution spatial coordinates

    if img is not None:
        ax.imshow(img, alpha=img_alpha, aspect='auto', 
                  extent=[0, img.shape[1] * img_scale_factor, img.shape[0] * img_scale_factor, 0])
        ax.set_xlim(0, img.shape[1] * img_scale_factor)
        ax.set_ylim(img.shape[0] * img_scale_factor, 0)
    else: 
        ax.set_aspect('equal', adjustable='box')

    plotted_annotations = set()
    for anno_val in unique_annotations:
        mask = annotations == anno_val
        if not np.any(mask): continue
        plotted_annotations.add(anno_val)
        ax.scatter(
            coords_to_plot[mask, 0], coords_to_plot[mask, 1],
            color=color_map.get(anno_val, '#CCCCCC'), s=spot_size, label=anno_val, 
            edgecolors=spot_edgecolor if spot_edgecolor != 'none' else None,
            linewidths=spot_linewidth if spot_edgecolor != 'none' else 0
        )
    if img is None: 
        ax.autoscale_view()
        if original_xlim != (0.0, 1.0) or original_ylim != (0.0, 1.0):
            ax.set_xlim(original_xlim); ax.set_ylim(original_ylim)
        else: 
             current_ylim = ax.get_ylim()
             if current_ylim[0] < current_ylim[1]: ax.invert_yaxis()
    ax.set_title(title, fontsize=plt.rcParams['axes.titlesize'], fontweight='bold')
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values(): spine.set_visible(False)
    if show_legend: _add_legend(ax, color_map, plotted_annotations, title=legend_title, ncol=legend_ncol, legend_fontsize=legend_fontsize)
    plt.tight_layout()
    if output_path: fig.savefig(output_path, dpi=300, bbox_inches='tight')
    if show: plt.show()
    if ax is None: plt.close(fig)
    return ax

def plot_spatial_deconvolution_pies(
    decon_adata, proportion_df, color_map=None, ax=None, fig_size=(7, 6),
    min_proportion_to_show=0.01, pie_scale_factor=0.8, # pie_scale_factor now acts as a multiplier for the dynamically calculated radius
    img_alpha=1.0, show_legend=True, legend_title="Cell Types", legend_ncol=1, legend_fontsize=None,
    title="Deconvolved Cell Type Proportions", output_path=None, show=False,
    default_palette='tab20', pie_edgecolor='white', pie_linewidth=0.2 
):
    if "spatial" not in decon_adata.obsm:
        raise ValueError("Spatial coordinates 'spatial' not found in decon_adata.obsm.")
    if not decon_adata.obs_names.equals(proportion_df.index):
        proportion_df = proportion_df.reindex(decon_adata.obs_names).fillna(0)

    spatial_coords_orig = decon_adata.obsm["spatial"].copy() # Full-resolution coordinates
    cell_type_names = proportion_df.columns.astype(str).tolist()

    if color_map is None:
        color_map = get_consistent_color_map([cell_type_names], default_palette=default_palette)
    if ax is None:
        fig, ax = plt.subplots(figsize=fig_size)
    else:
        fig = ax.get_figure()
        
    img, img_scale_factor = _get_spatial_image_and_scale(decon_adata)
    original_xlim, original_ylim = ax.get_xlim(), ax.get_ylim()
    coords_for_pies = spatial_coords_orig # Pies centered on full-resolution coordinates

    if img is not None:
        ax.imshow(img, alpha=img_alpha, aspect='auto', 
                  extent=[0, img.shape[1] * img_scale_factor, img.shape[0] * img_scale_factor, 0])
        ax.set_xlim(0, img.shape[1] * img_scale_factor)
        ax.set_ylim(img.shape[0] * img_scale_factor, 0)
    else: 
        ax.set_aspect('equal', adjustable='box')

    # --- Dynamic Pie Radius Calculation ---
    # This radius will be in the units of spatial_coords_orig (full-resolution)
    base_pie_radius = 0 
    # Try to get Visium spot diameter if available (often in fullres units)
    visium_spot_diameter = None
    if 'spatial' in decon_adata.uns and decon_adata.uns['spatial']:
        lib_id = list(decon_adata.uns['spatial'].keys())[0]
        if 'scalefactors' in decon_adata.uns['spatial'][lib_id] and \
           'spot_diameter_fullres' in decon_adata.uns['spatial'][lib_id]['scalefactors']:
            visium_spot_diameter = decon_adata.uns['spatial'][lib_id]['scalefactors']['spot_diameter_fullres']
    
    if visium_spot_diameter is not None:
        base_pie_radius = visium_spot_diameter / 2.0
        # print(f"    DEBUG: Using Visium spot_diameter_fullres for base_pie_radius: {base_pie_radius}")
    elif coords_for_pies.shape[0] > 1:
        # Fallback: estimate from median nearest neighbor distance in fullres coords
        sample_coords_for_pdist = np.unique(coords_for_pies[:min(1000, coords_for_pies.shape[0])], axis=0)
        if sample_coords_for_pdist.shape[0] > 1:
            dists = pdist(sample_coords_for_pdist)
            if len(dists) > 0:
                dist_matrix = squareform(dists)
                np.fill_diagonal(dist_matrix, np.inf)
                min_dists_per_spot = np.min(dist_matrix, axis=1)
                median_min_dist = np.median(min_dists_per_spot)
                if median_min_dist > 1e-9: 
                    base_pie_radius = median_min_dist / 2.0 
                    # print(f"    DEBUG: Using median_min_dist for base_pie_radius: {base_pie_radius}")

    if base_pie_radius <= 1e-9: # If still no valid radius (e.g., single spot, no visium info)
        if coords_for_pies.shape[0] > 0: # Fallback based on data extent
            x_min_sc, _ = np.min(coords_for_pies, axis=0)
            x_max_sc, _ = np.max(coords_for_pies, axis=0)
            x_extent_sc = x_max_sc - x_min_sc
            base_pie_radius = (x_extent_sc / 50.0) if x_extent_sc > 1e-9 else 5.0 # 1/50th of width
        else:
            base_pie_radius = 5.0 # Absolute fallback
        # print(f"    DEBUG: Using fallback extent-based base_pie_radius: {base_pie_radius}")

    final_pie_radius = base_pie_radius * pie_scale_factor # Apply user-defined scaling factor
    if final_pie_radius <= 1e-9: final_pie_radius = 1.0 # Ensure a tiny visible pie if all else fails
    # print(f"    DEBUG: final_pie_radius after pie_scale_factor ({pie_scale_factor}): {final_pie_radius}")
    # --- End of Dynamic Pie Radius Calculation ---
    
    plotted_cell_types = set()
    for i, spot_name in enumerate(decon_adata.obs_names):
        if spot_name not in proportion_df.index: continue
        proportions_series = proportion_df.loc[spot_name]
        valid_mask = (proportions_series >= min_proportion_to_show) & (~pd.isna(proportions_series))
        current_props = proportions_series[valid_mask].values
        current_types_local = proportions_series[valid_mask].index.astype(str).tolist()
        if not current_types_local or np.sum(current_props) < 1e-9 : continue 
        current_colors = [color_map.get(ct, '#CCCCCC') for ct in current_types_local]
        plotted_cell_types.update(current_types_local)
        center_x = coords_for_pies[i, 0]
        center_y = coords_for_pies[i, 1]
        ax.pie(
            current_props, center=(center_x, center_y), radius=final_pie_radius, colors=current_colors,
            wedgeprops={'edgecolor': pie_edgecolor, 'linewidth': pie_linewidth}, normalize=True, startangle=90
        )
    if img is None: 
        ax.autoscale_view() 
        if original_xlim != (0.0, 1.0) or original_ylim != (0.0, 1.0): 
            ax.set_xlim(original_xlim); ax.set_ylim(original_ylim)
        else:
             current_ylim = ax.get_ylim()
             if current_ylim[0] < current_ylim[1]: ax.invert_yaxis()
    ax.set_title(title, fontsize=plt.rcParams['axes.titlesize'], fontweight='bold')
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values(): spine.set_visible(False)
    if show_legend: _add_legend(ax, color_map, plotted_cell_types, title=legend_title, ncol=legend_ncol, legend_fontsize=legend_fontsize)
    plt.tight_layout()
    if output_path: fig.savefig(output_path, dpi=300, bbox_inches='tight')
    if show: plt.show()
    if ax is None: plt.close(fig) 
    return ax
