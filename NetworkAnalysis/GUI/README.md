# MEA Data Analysis Manual

## Shape 

### File Locations, References and MATLAB

- **Copy data**: from "benshalomnas" network to the lab server: "analysis machine" (folder: "primaryneurondata").
- **Create the reference file** ("librecalc.xlsx") with columns: Date, DIV, Assay, Run #, up to before Comments.
- **Open MATLAB** and add the following folders to the file path:
  - MEA analysis
  - Maxwell toolbox
  - data file
- **Navigate** in MATLAB to the "GUI" folder in the user folder:
  - `MEA analysis > GUI > AnalysisGUI.m`

### PHASE 1: Network Analysis

#### Set 1: Initial Network Analysis

##### Running the Analysis GUI

1. Run `AnalysisGUI.m`.
2. Set Default Parameters:
   - Gaussian sigma: 0.16 s
   - Bin size: 0.02s
   - Min PeakDis: 1.0s
   - Threshold Burst: 1.2s
   - Threshold Start-Stop: 0.3-0.4s
3. Click "Process Network."
4. **Verify the output**: Your output from this GUI will be raster plots and one .csv file.
5. **Double-check** plots generated correctly.

##### Plotting from the CSV

- **Verify the settings** in "net_plt_settings.json" from the MEA Analysis folder on the Files application.
- **Execute the script**: Right-click on `networkbarplots.sh` file and Run as Program.
- **Check for errors**: See if there is any error on the terminal, if/else, hit enter to close the program as prompted.
- 🌱 Your bar plots will be saved under Output Directory / Network outputs/ Burstproperty graphs.

#### Set 2: Parameter Comparison

##### Running Parameters Comparison

1. Input the specific DIV folder path.
2. Copy and paste the reference file.
3. If file is locked, Open in Terminal, and type in `sudo chmod 777 -R *`
4. Specify the output folder path.
5. Set Default Parameters.
6. Click "Explore Parameters."
7. **Examine generated parameter plots** for data consistency.
8. When seeing outliers on the graphs, identify which chip/s these are and double-check active area on the main excel sheet.

#### Set 3: 2nd Network Analysis with Adjusted Parameters

- Input "Adjusted Parameters" into Phase 1 replacing default parameters.
- **Find new representative raster plots**: You should now have new raster plots that reflect Adjusted Parameters that have overwritten your first network outputs.

##### Plotting Burst Property Graphs and Network Activity Plots

- **Check settings** in “act_plt_settings.json” where the directories are correct.
- **Execute the script**: Right-click `activitybarplots.sh` and select “Run as Program”.
  - Activity bar plots should be generated within “meanActivityProperty_graphs”.
- **Check settings** in “net_plt_settings.json” where the directories are correct.
- **Execute the script**: Right-click `networkbarplots.sh` and select “Run as Program”.
  - Network bar plots should be generated within “burstProperty_graphs”.

### PHASE 2: Manual Curation -- Remove overlapping points indicating one unit.

#### Set 1:

1. **Open VS Code**.
2. **Navigate**: Open MEA_Analysis folder (workspace that has all the analysis files).
3. **Access IPNAnalysis**: On VScode click on IPNAnaylsis > ManualCuration.ipynb.
4. **Prepare the data**: Go to Analyzed Data in MEA_Analysis and copy - paste waveforms_good path onto VS Code.
5. **Run Analysis**: Developer needs to install spike interface and make necessary changes.
6. **Execute the script**: Select Kernel and then Python. Run Analyze Waveforms. Explore the manual curation interface (MCI) and hit Compute.
7. **View outputs**: View similaritview and traceview.
8. 🌱Note the redundant units to be deleted by keeping the higher number of spike units.
9. **Delete redundant data**: Go into SpikeTrains and delete the unit numbers. They should be in “IDnumber.m” format.

#### Set 2: HYBERBURST ANALYSIS

##### OPENING AND ADJUSTING FIGURES

- `openfig(fileName.fig, 'visible'); % To make the figure visible`
- **Adjust figures**: Once figure is visible, zoom in or out to find the higher or lower y-limit. Save the new figure in desired folder/path location.
- **Run the Analysis GUI**: Once the y-limit has been determined, run the Analysis GUI to change a bulk set of figures.

##### SMOOTHING OUT HYPERBURSTS

1. Pick 5-10 raster plots/graphs from desired DIV.
2. Set directory path, set output path, set reference file.
3. Note the line’s adjusted parameters as the BEFORE set.
4. Adjust the Gaussian value by increasing it. Hit “Process Network”.
5. Repeat until satisfied with the smoothing curve for all plots in that set.
6. Once the Gaussian value is set, apply the “Hyperburst parameters” to all chips in the same DIV.

#### Set 3: FIGURE PLOTTING

- **Open GUI folder**,
- **Navigate to filename**: replotfigures.m and input the source directory and output directory.
- **Set y-min and y-max**.
