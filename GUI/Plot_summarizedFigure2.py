from PIL import Image
import pandas as pd

project = "ADNP"
parentfolder = "/home/jonathan/Documents/Scripts/Matlab/scripts_output/"
reference_f = '/home/jonathan/Documents/Scripts/Python/ADNP_Notes.xlsx'
div = 11
assay_type = 'today'
wt_chip = 16378
het_chip = 16867

ref_df = pd.read_excel(reference_f)
wt_df = ref_df.loc[(ref_df['Div'] == div) & (ref_df['Chip ID'] == wt_chip) & ref_df['Assay'].str.lower().str.contains(assay_type.lower())]
wt_run = "{:06d}".format(wt_df['Run #'].unique()[0])
het_df = ref_df.loc[(ref_df['Div'] == div) & (ref_df['Chip ID'] == het_chip) & ref_df['Assay'].str.lower().str.contains(assay_type.lower())]
het_run = "{:06d}".format(het_df['Run #'].unique()[0])

#scale = 0.75

r1 = Image.open(parentfolder+project+'/Network_outputs/Raster_BurstActivity/Raster_BurstActivity'+wt_run+'.png')
r2 = Image.open(parentfolder+project+'/Network_outputs/Raster_BurstActivity/Raster_BurstActivity'+het_run+'.png')
#r1 = Image.open("/home/jonathan/Documents/Scripts/Matlab/scripts_output/ADNP/Network_outputs/Raster_BurstActivity/Raster_BurstActivity000015.png")
#r2 = Image.open("/home/jonathan/Documents/Scripts/Matlab/scripts_output/ADNP/Network_outputs/Raster_BurstActivity/Raster_BurstActivity000017.png")

r_w,r_h = r1.size
new_r_w = 470
new_r_h = round(new_r_w/r_w*r_h)

r1 = r1.resize((new_r_w,new_r_h))
r2 = r2.resize((new_r_w,new_r_h))

#scale = 0.3

bp1 = Image.open(parentfolder+project+"/Network_outputs/BurstProperty_graphs/Network Today IBI.png")
bp2 = Image.open(parentfolder+project+"/Network_outputs/BurstProperty_graphs/Network Today Burst Peak.png")
bp3 = Image.open(parentfolder+project+"/Network_outputs/BurstProperty_graphs/Network Today Number Bursts.png")
bp4 = Image.open(parentfolder+project+"/Network_outputs/BurstProperty_graphs/Network Today Spike per Burst.png")

bp_w,bp_h = bp1.size
new_bp_w = 575
new_bp_h = round(new_bp_w/bp_w*bp_h)

bp1 = bp1.resize((new_bp_w,new_bp_h))
bp2 = bp2.resize((new_bp_w,new_bp_h))
bp3 = bp3.resize((new_bp_w,new_bp_h))
bp4 = bp4.resize((new_bp_w,new_bp_h))

#scale = 0.80

pc1 = Image.open(parentfolder+project+"/Network_outputs/ParameterCompare_Gaussian/paramsCompare_singlePlot.png")
pc2 = Image.open(parentfolder+project+"/Network_outputs/ParameterCompare_BinSize/paramsCompare_singlePlot.png")
pc3 = Image.open(parentfolder+project+"/Network_outputs/ParameterCompare_Threshold/paramsCompare_singlePlot.png")
pc4 = Image.open(parentfolder+project+"/Network_outputs/ParameterCompare_FixedThreshold/paramsCompare_singlePlot.png")
pc5 = Image.open(parentfolder+project+"/Network_outputs/ParameterCompare_MinPeakDistance/paramsCompare_singlePlot.png")
pc6 = Image.open(parentfolder+project+"/Network_outputs/ParameterCompare_StartStopThreshold/paramsCompare_singlePlot.png")

pc_w,pc_h = pc1.size
new_pc_w = 1000
new_pc_h = round(new_pc_w/pc_w*pc_h)

pc1 = pc1.resize((new_pc_w,new_pc_h))
pc2 = pc2.resize((new_pc_w,new_pc_h))
pc3 = pc3.resize((new_pc_w,new_pc_h))
pc4 = pc4.resize((new_pc_w,new_pc_h))
pc5 = pc5.resize((new_pc_w,new_pc_h))
pc6 = pc6.resize((new_pc_w,new_pc_h))

f = Image.new(mode='RGBA',size=(3840,1900),color='white')
#
f.paste(r1,(50,30))
f.paste(r2,(450,30))
#
f.paste(pc2,(2880,180))
f.paste(pc1,(1990,180))
#
f.paste(bp1,(960,60))
f.paste(bp2,(1480,60))
f.paste(bp3,(960,480))
f.paste(bp4,(1480,480))
#
f.paste(pc3,(0,1080))
f.paste(pc4,(960,1080))
f.paste(pc5,(1920,1080))
f.paste(pc6,(2880,1080))

f.save(parentfolder+project+'/Network_outputs/SummarizedFigure2.png')