import pandas as pd
import numpy as np
import re
import json
import os
import math
from collections import defaultdict
import glob # For finding files

# --- parse_digpts_file (unchanged) ---
def parse_digpts_file(filename):
    """Parse digpts.txt file to extract channel positions"""
    ch_pos = {}
    with open(filename, 'r') as f:
        for line in f:
            if not line.strip(): continue
            match = re.match(r'(\w+):?\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)', line)
            if match:
                ch_name, x, y, z = match.groups()
                ch_pos[ch_name] = np.array([float(x), float(y), float(z)])
    return ch_pos

# --- project_point_azimuthal_equidistant (unchanged) ---
def project_point_azimuthal_equidistant(point_3d, center_3d):
    """Projects a 3D point onto a 2D plane using Azimuthal Equidistant projection"""
    vec = point_3d - center_3d
    dist = np.linalg.norm(vec)
    if dist < 1e-9: return 0.0, 0.0
    acos_arg = np.clip(vec[2] / dist, -1.0, 1.0)
    theta = math.acos(acos_arg)
    phi = math.atan2(vec[1], vec[0])
    rho = theta
    proj_x = rho * math.cos(phi)
    proj_y = rho * math.sin(phi)
    return proj_x, proj_y
# ---

def prepare_trial_visualization_data(
    trial_data,
    trial_num,
    ch_pos,
    chromophore,
    srate,
    downsample_rate,
    do_hta_averaging,
    hta_df,
    dsfreq,
    ss_detector_numbers,
    ignore_fiducials=False):
    """
    Processes data for a single trial, optionally performs HTA averaging. Streamlined prints.
    """
    print(f"\nProcessing Trial {trial_num}, Chromophore {chromophore}...")
    all_subjects_data_processed = {}
    potential_subject_list = sorted(trial_data.keys())
    common_channel_names = None # Defined by first subject
    first_subject = True
    active_sources_trial = set()
    active_detectors_trial = set()
    selected_channel_info = [] # Based on first subject
    

    for subject_id in potential_subject_list:
        # print(f"\nProcessing Subject: {subject_id}") # Removed per-subject start print
        df = trial_data[subject_id]
        subject_all_channel_names = df.columns.tolist()

        # --- Determine common standard channels (based on first subject) ---
        if first_subject:
            temp_selected_channels = []
            for i, name in enumerate(subject_all_channel_names):
                match = re.match(r'(S\d+)_(D\d+)_(' + chromophore + r')$', name)
                if match:
                    s_name, d_name, _ = match.groups()
                    if s_name in ch_pos and d_name in ch_pos:
                        try: d_num = int(d_name[1:])
                        except ValueError: continue
                        if d_num not in ss_detector_numbers: # Standard channels only
                            temp_selected_channels.append({ 'name': name, 'index': i, 's': s_name, 'd': d_name })
                            active_sources_trial.add(s_name)
                            active_detectors_trial.add(d_name)
            selected_channel_info = temp_selected_channels
            common_channel_names = [ch['name'] for ch in selected_channel_info]
            if not selected_channel_info:
                 print(f"  Warning: No standard channels found for {chromophore} in Trial {trial_num}. Skipping trial.")
                 return None
            print(f"  Defined {len(common_channel_names)} standard channels for this trial/chromophore.")
            first_subject = False

        # --- Data Selection for Current Subject ---
        try:
            selected_indices = [subject_all_channel_names.index(name) for name in common_channel_names]
            if len(selected_indices) != len(common_channel_names):
                 raise ValueError("Mismatch after index lookup.")
        except ValueError as e:
             print(f"  Warning: Channel index error for {subject_id} ({e}). Skipping subject.")
             continue

        signal_data_selected = df.iloc[:, selected_indices].values

        # --- Data Processing (HTA or Downsampling) ---
        processed_times = []
        processed_data = np.array([])
        input_data_for_processing = signal_data_selected
        processing_method = "downsampling" # Default

        if do_hta_averaging and hta_df is not None:
            hta_col_name = HTACOL_RENAME_FN(subject_id)
            if hta_col_name in hta_df.columns:
                processing_method = "HTA"
                try:
                    hta_times_sec = hta_df[hta_col_name].dropna().values
                    hta_indices_full = np.round((hta_times_sec - hta_times_sec[0]) * srate).astype(int)
                    max_idx_full = input_data_for_processing.shape[0] - 1
                    hta_indices_full = np.clip(hta_indices_full, 0, max_idx_full)
                    hta_indices_full = np.unique(hta_indices_full)
                    if len(hta_indices_full) < 2: raise ValueError("Not enough unique HTA indices.")

                    averaged_signal_list = []
                    segment_midpoint_times = []
                    for i in range(len(hta_indices_full) - 1):
                        start_idx, end_idx = hta_indices_full[i], hta_indices_full[i+1]
                        if start_idx >= end_idx: continue
                        segment_data = input_data_for_processing[start_idx:end_idx, :]
                        if segment_data.shape[0] == 0: continue
                        segment_avg = np.mean(segment_data, axis=0)
                        averaged_signal_list.append(segment_avg)
                        start_time_sec = hta_times_sec[i] - hta_times_sec[0]
                        end_time_sec = hta_times_sec[i+1] - hta_times_sec[0]
                        segment_midpoint_times.append((start_time_sec + end_time_sec) / 2.0)

                    if averaged_signal_list:
                        processed_data = np.array(averaged_signal_list)
                        processed_times = np.array(segment_midpoint_times)
                    else:
                        print(f"  Warning: No valid HTA segments found for {subject_id}. Skipping subject.")
                        continue
                except Exception as e_hta:
                    print(f"  Warning: Error during HTA for {subject_id} ({e_hta}). Skipping subject.")
                    continue
            else:
                # HTA required but column missing -> Skip subject
                print(f"  Warning: HTA column '{hta_col_name}' not found for {subject_id}. Skipping subject as HTA is required.")
                continue
        else:
            # Standard downsampling
            processed_data = input_data_for_processing[::downsample_rate, :]
            processed_times = np.arange(processed_data.shape[0]) / dsfreq

        # --- Check processed data shape ---
        if processed_data.ndim != 2 or processed_data.shape[1] != len(common_channel_names):
             print(f"  ERROR: Processed data shape mismatch for {subject_id} after {processing_method}! "
                   f"Expected columns: {len(common_channel_names)}, Got shape: {processed_data.shape}. Skipping subject.")
             continue

        # --- NaN Handling ---
        if np.sum(np.isnan(processed_data)) > 0:
            # print(f"    Replacing NaN values with 0 for {subject_id}.") # Optional: uncomment if needed
            processed_data = np.nan_to_num(processed_data, nan=0.0)

        # --- Final Check and Store ---
        if processed_data.size == 0:
             # This case should be caught earlier, but double-check
             print(f"  Warning: No data resulted from processing for {subject_id}. Skipping storage.")
             continue

        all_subjects_data_processed[subject_id] = {
            'times': processed_times.tolist(),
            'data': processed_data.tolist()
        }

    # --- Post-Loop Processing ---
    if not all_subjects_data_processed:
        print(f"  No subjects processed successfully for Trial {trial_num}, {chromophore}.")
        return None

    final_subject_list = sorted(all_subjects_data_processed.keys())
    print(f"  Successfully processed data for {len(final_subject_list)} subjects.")

    # --- Geometry Calculations (unchanged) ---
    points_for_center = []
    fiducials = {'Nz', 'Iz', 'LPA', 'RPA'}
    for name, pos in ch_pos.items():
        is_fiducial = name in fiducials
        is_active_optode = name in active_sources_trial or name in active_detectors_trial
        if is_active_optode or (not ignore_fiducials and is_fiducial):
             points_for_center.append(pos)
    if not points_for_center: raise ValueError(f"No points for head center in Trial {trial_num}.")
    head_center_3d = np.mean(points_for_center, axis=0)

    midpoints_3d_list = []
    midpoint_radii = []
    midpoints_3d_dict = {}
    for ch in selected_channel_info: # Use the original common channel info
        s_pos, d_pos = ch_pos[ch['s']], ch_pos[ch['d']]
        midpoint_3d = (s_pos + d_pos) / 2
        midpoints_3d_list.append(midpoint_3d)
        midpoints_3d_dict[ch['name']] = midpoint_3d
        dist_from_center = np.linalg.norm(midpoint_3d - head_center_3d)
        if dist_from_center > 1e-6: midpoint_radii.append(dist_from_center)
    if not midpoint_radii: raise ValueError(f"Could not calculate radius in Trial {trial_num}.")
    head_radius_3d = np.mean(midpoint_radii)

    projected_midpoints_2d = {}
    for ch_name, midpoint_3d in midpoints_3d_dict.items():
         projected_midpoints_2d[ch_name] = project_point_azimuthal_equidistant(midpoint_3d, head_center_3d)

    projected_optodes_2d = {}
    projected_nose_2d = (0, 1)
    active_optodes_trial = active_sources_trial.union(active_detectors_trial)
    for opt_name in active_optodes_trial:
        projected_optodes_2d[opt_name] = project_point_azimuthal_equidistant(ch_pos[opt_name], head_center_3d)
    if 'Nz' in ch_pos: projected_nose_2d = project_point_azimuthal_equidistant(ch_pos['Nz'], head_center_3d)

    max_rho_midpoint = max((math.sqrt(p[0]**2 + p[1]**2) for p in projected_midpoints_2d.values()), default=0)
    scale_factor = max_rho_midpoint * 1.1 if max_rho_midpoint > 1e-9 else 1.0
    norm_optodes_x = {n: p[0] / scale_factor for n, p in projected_optodes_2d.items()}
    norm_optodes_y = {n: p[1] / scale_factor for n, p in projected_optodes_2d.items()}
    nose_x_norm = projected_nose_2d[0] / scale_factor
    nose_y_norm = projected_nose_2d[1] / scale_factor

    # --- Get final 3D positions in the correct common order ---
    final_positions_3d = []
    for ch_name in common_channel_names: # Use the defined common order
        midpoint_3d = midpoints_3d_dict.get(ch_name)
        if midpoint_3d is not None:
            final_positions_3d.append(midpoint_3d.tolist())
        else:
             # This indicates an internal logic error if it happens
             print(f"  ERROR: Could not find 3D midpoint for common channel {ch_name}!")
             final_positions_3d.append([np.nan, np.nan, np.nan])
             
             
    # --- Project Fiducials (if they exist in ch_pos) ---
    fiducial_names_to_draw = fiducials # ['Nz', 'Iz', 'LPA', 'RPA']
    projected_fiducials_norm = {}
    print("  Projecting fiducials for drawing...")
    for name in fiducial_names_to_draw:
        if name in ch_pos:
            fid_pos_3d = ch_pos[name]
            proj_x, proj_y = project_point_azimuthal_equidistant(fid_pos_3d, head_center_3d)
            norm_x = proj_x / scale_factor
            norm_y = proj_y / scale_factor
            projected_fiducials_norm[name] = {'x': norm_x, 'y': norm_y}
            print(f"    Fiducial {name}: Norm coords ({norm_x:.2f}, {norm_y:.2f})")
        else:
            print(f"    Fiducial {name} not found in digpts.txt.")

    # --- Final JS data structure ---
    js_data = {
        "trialNum": trial_num,
        "chromophore": chromophore,
        "channelType": "standard",
        "channels": {
            "names": common_channel_names, # The common ordered list
            "positions3d": final_positions_3d, # Corresponding 3D positions in order
        },
        "optodes": {
            "sources": {n: {'x': norm_optodes_x[n], 'y': norm_optodes_y[n]} for n in active_sources_trial},
            "detectors": {n: {'x': norm_optodes_x[n], 'y': norm_optodes_y[n]} for n in active_detectors_trial},
            "pairs": [[ch['s'], ch['d']] for ch in selected_channel_info] # Pairs from the selection
        },
        "head": { "nose": {"x": nose_x_norm, "y": nose_y_norm} },
        "geometry": {
            "headCenter3d": head_center_3d.tolist(),
            "headRadius3d": head_radius_3d,
            "projectionScaleFactor": scale_factor
        },
        "subjects": all_subjects_data_processed,
        "subjectList": final_subject_list,
        "fiducials": projected_fiducials_norm
    }
    return js_data

# --- write_html_file (unchanged) ---
def write_html_file(js_data, template_path="template.html", output_dir="visualizations"):
    if js_data is None: return
    trial_num = js_data['trialNum']
    chromophore = js_data['chromophore']
    try:
        with open(template_path, "r") as f: template = f.read()
    except FileNotFoundError:
        print(f"Error: HTML template file not found at {template_path}")
        return
    os.makedirs(output_dir, exist_ok=True)
    output_filename = os.path.join(output_dir, f"fnirs_viewer_Trial{trial_num:02d}_{chromophore}.html")
    # Slightly safer JSON embedding (handles potential issues with </script>)
    json_string = json.dumps(js_data).replace("</", "<\\/")
    html_content = template.replace("var EEG_DATA = {};", f"var EEG_DATA = {json_string};")
    with open(output_filename, "w") as f: f.write(html_content)
    # print(f"Visualization saved to {output_filename}") # Keep final save message? Yes.
    print(f"  Visualization saved to {output_filename}") # Indented under trial processing

# --- Main Execution Logic (unchanged) ---
if __name__ == "__main__":
    # --- Configuration ---
    data_directory = fr"Trialwise_txts"
    digpts_file = 'digpts_PC.txt'
    output_visualization_dir = "visualizations"
    template_html = "template.html"
    
    # ignore Hierarchial Task Analysis-based averaging. Averages signals at sub_task intervals as defined in an excelsheet
    DO_HTA_AVERAGING = False # True #
    HTA_TRIALNUM = 4
    # HTA_FILE_PATH = r"\HTA_forattn.xlsx"
    # HTA_SHEET_NAME = 'Sheet1'
    # HTACOL_RENAME_FN = lambda x: x.replace("_", "")
    
    chromophores_to_process = ["HbO", "HbR"]

    srate = 7.8125
    downsample_rate = 8
    ss_detector_numbers = set(range(20, 29))
    

    # --- End Configuration ---

    # --- Load HTA data if enabled ---
    hta_df = None
    dsfreq = srate / downsample_rate
    # if DO_HTA_AVERAGING:
    #     try:
    #         hta_df = pd.read_excel(HTA_FILE_PATH, sheet_name=HTA_SHEET_NAME)
    #         print(f"Successfully loaded HTA data from {HTA_FILE_PATH}")
    #         output_visualization_dir = output_visualization_dir + "_HTA"
    #     except FileNotFoundError:
    #         print(f"Error: HTA file not found at {HTA_FILE_PATH}. HTA Averaging disabled.")
    #         DO_HTA_AVERAGING = False
    #     except Exception as e:
    #         print(f"Error loading HTA file: {e}. HTA Averaging disabled.")
    #         DO_HTA_AVERAGING = False

    # 1. Load Digitizer Positions
    try:
        ch_pos = parse_digpts_file(digpts_file)
    except FileNotFoundError:
         print(f"Error: {digpts_file} not found.")
         exit()
    if not ch_pos: raise ValueError('Channels cannot be parsed. Check digpts_file')
         
    # 2. Find and Group Data Files by Trial
    all_files = glob.glob(os.path.join(data_directory, "*.txt"))
    grouped_files = defaultdict(dict)
    for f_path in all_files:
        filename = os.path.basename(f_path)
        match = re.match(r'(.+)_(\d+)\.txt', filename)
        if match:
            subject_id, trial_num_str = match.groups()
            try:
                trial_num = int(trial_num_str)
                grouped_files[trial_num][subject_id] = f_path
            except ValueError: print(f"Warning: Could not parse trial number from {filename}. Skipping.")
        else: print(f"Warning: Filename {filename} does not match expected pattern. Skipping.")
    if not grouped_files:
        print(f"Error: No files matching the pattern found in {data_directory}")
        exit()
    print(f"Found data for trials: {sorted(grouped_files.keys())}")

    # 3. Process Each Trial
    for trial_num in sorted(grouped_files.keys()):
        if DO_HTA_AVERAGING and 'HTA_TRIALNUM' in locals() and trial_num != HTA_TRIALNUM:
            # print(f'Skipping Trial {trial_num} (HTA only for Trial {HTA_TRIALNUM})') # Can uncomment if needed
            continue

        print(f"\n===== Starting Trial {trial_num} =====")
        subject_files = grouped_files[trial_num]
        trial_data_loaded = {}
        for subject_id, f_path in subject_files.items():
            try:
                df = pd.read_csv(f_path, sep=',', header=0)
                trial_data_loaded[subject_id] = df
                # print(f"  Loaded {os.path.basename(f_path)}") # Removed per-file load message
            except Exception as e:
                print(f"Error loading {f_path}: {e}. Skipping subject {subject_id} for Trial {trial_num}.")
                continue
        if not trial_data_loaded:
            print(f"Warning: No data successfully loaded for Trial {trial_num}. Skipping.")
            continue

        # Generate visualization for each chromophore
        for chromo in chromophores_to_process:
            js_data_for_trial = prepare_trial_visualization_data(
                trial_data=trial_data_loaded,
                trial_num=trial_num,
                ch_pos=ch_pos,
                chromophore=chromo,
                srate=srate,
                downsample_rate=downsample_rate,
                do_hta_averaging=DO_HTA_AVERAGING,
                hta_df=hta_df,
                dsfreq=dsfreq,
                ss_detector_numbers=ss_detector_numbers,
                # ignore_fiducials defaults to True
            )
            if js_data_for_trial:
                 write_html_file(js_data_for_trial, template_html, output_visualization_dir)
            # else: # Message printed inside prepare_trial... if it returns None
            #      print(f"Skipping HTML generation for Trial {trial_num}, {chromo} due to processing errors.")

    print("\n===== Processing Complete =====")