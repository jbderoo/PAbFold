import os
import numpy as np
import matplotlib.pyplot as plt
from pprint import pprint
from numpy import hstack, array, zeros, vstack, dstack, save
import warnings
import sys
sys.path.append('../../')
from seqs import seqs, numaa, slide, peplen


#############
import xml.etree.ElementTree as ET

def get_dimension_value(dimension_str):
    """Convert SVG dimension (with 'px' or 'pt') to a float."""
    if 'px' in dimension_str:
        return float(dimension_str.replace('px', ''))
    elif 'pt' in dimension_str:
        # Assuming 1pt is approximately 1.33px (default for browsers).
        return float(dimension_str.replace('pt', '')) * 1.33
    else:
        return float(dimension_str)

def combine_svgs(svg1_path, svg2_path, output_path):
    # Parse the SVG files
    tree1 = ET.parse(svg1_path)
    tree2 = ET.parse(svg2_path)

    root1 = tree1.getroot()
    root2 = tree2.getroot()

    # Get the width and height of the first SVG
    width1 = get_dimension_value(root1.get('width'))
    height1 = get_dimension_value(root1.get('height'))

    # Get the width and height of the second SVG
    width2 = get_dimension_value(root2.get('width'))
    height2 = get_dimension_value(root2.get('height'))

    # Update the width and height of the second SVG to be moved next to the first one
    root2.set('x', str(width1))

    # Set up the new root SVG element with updated width and height
    new_root = ET.Element('svg', xmlns="http://www.w3.org/2000/svg", version="1.1",
                           width=str(width1 + width2) + 'px', height=str(max(height1, height2)) + 'px')

    # Append all elements of the first SVG to the new root
    for element in root1:
        new_root.append(element)

    # Append all elements of the second SVG to the new root, updating their positions
    for element in root2:
        if element.tag.endswith('g') or element.tag.endswith('path') or element.tag.endswith('circle') or element.tag.endswith('rect'):
            if 'transform' in element.attrib:
                transformations = element.attrib['transform'].split(' ')
                for i, transform in enumerate(transformations):
                    if 'translate' in transform:
                        values = transform.strip('translate()').split(',')
                        values[0] = str(float(values[0]) + width1)
                        transformations[i] = 'translate(' + ','.join(values) + ')'
                element.attrib['transform'] = ' '.join(transformations)
            else:
                element.attrib['transform'] = f"translate({width1},0)"
        new_root.append(element)

    # Write the combined SVG to the output path
    tree = ET.ElementTree(new_root)
    tree.write(output_path, encoding="utf-8", xml_declaration=True)

#############



try:
    highlight_window = [int(sys.argv[1]), int(sys.argv[2])]
except:
    print('Provide start and end posts')
    print('HA:    113 - 122')
    print('Myc:   423 - 433')
    print('mBG17: 400 - 409')
    sys.exit()

# Plotting the 'Max Ever Seen Per Residue' data
#highlight_window = [113, 122] ## HA
#highlight_window = [423, 433] ##Myc 
#highlight_window = [400, 409] ## mBG17
#highlight_window = [0, 10]  ## hiv_front_myc
#highlight_window = [99, 108]  ## hiv_back_myc
#highlight_window = [0, 4]  #First amino acid is 0
#highlight_window = [1, 5]  #First amino acid is 0
outsvg = 'scan.svg'
pepsize = peplen
slide = slide
# Load the numpy arrays
d1 = np.load('rowified_data/model_1_sliding_window_plddt.npy')
d2 = np.load('rowified_data/model_2_sliding_window_plddt.npy')
d3 = np.load('rowified_data/model_3_sliding_window_plddt.npy')
d4 = np.load('rowified_data/model_4_sliding_window_plddt.npy')
d5 = np.load('rowified_data/model_5_sliding_window_plddt.npy')


## The script that produces windows will make an extra window
## with a different slide to hit the last amino acids.
## So there is either 0 or 1 additional window, depending
## on whether slide goes evenly into the number of amino acids.
num_windows = len(d1)
print('Number of amino acids',numaa,'Number of windows',num_windows)

if (numaa-pepsize) % slide == 0:
    print('Looks like we are lucky and the windows go evenly into the target protein length')
    N = pepsize + slide*(num_windows - 1)
    assert N == numaa
    final_slide = slide
else:
    final_slide = (numaa-pepsize) % slide
    print('Looks like there will be a final slide of %d since the windows do not go evenly into the target protein length' % final_slide)
    N = pepsize + slide*(num_windows - 2) + final_slide
    assert N == numaa

## Build 3D array to get the max_ever_seen_per_residue_pLDDT
layers = []
for model in range(1,6):
    d = np.load('rowified_data/model_%d_sliding_window_plddt.npy' % model)
    datalist = []
    index = 0
    for iii, pLDDT in enumerate(d):
        endbuffer = N - pepsize - index
        #print(iii, len(d)-1, index, endbuffer)
        #if endbuffer < 0: break
        datalist.append(hstack([zeros(index), array(pLDDT), zeros(endbuffer)]))

        if iii == len(d)-2:
            print('Last window, sliding by',final_slide)
            index += final_slide
        else:
            index += slide
    data = vstack(datalist)
    layers.append( data )

alldata = dstack(layers)
max_ever_seen = np.max(alldata, axis=(0,2))
save('max_ever_seen_per_residue.npy',max_ever_seen)


data = d1 + d2 + d3 + d4 + d5
data /= 5


# Adjust font size
plt.rcParams.update({'font.size': 22})
file_path = outsvg

sums = data.sum(1)
means = data.mean(1)

if 1:
    top10 = (means > sorted(means)[-11]).nonzero()[0]
    top5 = (means > sorted(means)[-6]).nonzero()[0]
    seq2score = {}
    for ii in top10:
        seq2score[(ii,seqs[ii])] = means[ii]
        print('%.1f: %s'%(means[ii], seqs[ii]))
    #pprint(sorted(seq2score.items(),key=lambda x:x[1],reverse=True))

    ## If we use all but the first character of the last seq
    targetseq = ''.join([window[:1] for window in seqs[:-1]]) + seqs[-1][1:]
    print(targetseq)
    lo, hi = highlight_window
    print('Will highlight:',targetseq[lo:hi])

# Create the main figure and axes for the heatmap and the bar chart
#fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 10), gridspec_kw={'height_ratios': [4, 2, 2]})
fig, ax = plt.subplots(figsize=(6, 10))

# Plotting the heatmap on ax1 (flipped rows and columns)
cax = ax.imshow(data.T, aspect='auto', cmap='viridis', origin='lower')
cbar = fig.colorbar(cax, ax=ax, label='Per-Residue pLDDT')
cbar.set_label('Per-Residue pLDDT', size=32)
cbar.ax.tick_params(labelsize=32)
cbar.ax.yaxis.labelpad = 20 
#ax.remove()

# Suppress the specific warning using a context manager
with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=UserWarning)
    plt.tight_layout()  # To ensure the colorbar fits within the figure boundaries
plt.savefig(file_path[:-4]+'_scalebar.svg', format='svg')

# Re-create the main figure and axes for the heatmap and the bar chart
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 10), gridspec_kw={'height_ratios': [4, 2, 2]})
cax = ax1.imshow(data.T, aspect='auto', cmap='viridis', origin='lower')
#cbar = fig.colorbar(cax, ax=ax1, label='Per-Residue pLDDT')

# Use colormap to determine colors of the bars
norm = plt.Normalize(data.min(), data.max())
colors = plt.cm.viridis(norm(means))

from matplotlib import font_manager

preferred_font = "Courier New"
fallback_font = "DejaVu Sans Mono"

if preferred_font in [f.name for f in font_manager.fontManager.ttflist]:
    fn = preferred_font
else:
    fn = fallback_font

y_position = 0
# Sort the means to decide where to put the text
sorted_indices = np.argsort(means)[::-1]
for ii in sorted_indices[:5][::-1]:
    ax3.text(len(max_ever_seen) + 20, y_position, '%.1f: %s'%(means[ii], seqs[ii]), color=colors[ii], va="center", fontsize=22, fontname=fn, weight='bold')
    y_position += 22 

# Plot the bar chart with the colors determined above
ax2.bar(range(len(means)), means, color=colors, width=1.0)
ax2.set_xlim(ax1.get_xlim())
ax2.set_ylim(bottom=min(means), top=max(means) * 1.1)
#ax2.set_ylabel('Mean')
ax2.set_xlabel('Window in Target Protein', labelpad=8)
ax2.axis('on')
#cbar = fig.colorbar(cax, ax=ax2, label='Per-Residue pLDDT')

# Other plotting details
for y in range(1, pepsize):
    ax1.axhline(y=y-0.5, color='black', linestyle='-')
row_labels = [str(i) for i in range(1, pepsize+1)]
ax1.set_yticks(range(pepsize))
ax1.set_yticklabels(row_labels)
ax1.tick_params(axis='x', bottom='off', which='both', labelbottom=False) 

#ax1.set_xlabel('Window in Target Protein')
#ax1.set_ylabel('Internal Sliding Window Index')


# Adding a faint grey background between x-axis and value 60
ax3.fill_between(range(len(max_ever_seen)), 60, color='grey', alpha=0.2)
ax3.plot(max_ever_seen, color='black')
# Shading the region below highlight_window in orange
ax3.fill_between(range(highlight_window[0], highlight_window[1]), 
                 max_ever_seen[highlight_window[0]:highlight_window[1]], 
                 color=(1.0,0.8,1.0), alpha=1.0)

# Overlaying circular markers for values above 60
above_60_indices = np.where(max_ever_seen > 60)[0]
above_60_values = max_ever_seen[above_60_indices]
ax3.scatter(above_60_indices, above_60_values, color='magenta', marker='o', s=30, label='Values > 60')

#ax3.set_ylabel('Max')
ax3.set_xlim((ax1.get_xlim()[0],len(max_ever_seen)))
ax3.set_xlabel('Residue in Target Protein', labelpad=8)

# Adjust the layout and save
#plt.subplots_adjust(hspace=0.0) 
plt.tight_layout(pad=0.0)
plt.savefig(file_path, format='svg')
#plt.show()

print(f'attempting to combine svgs... {file_path} with { file_path[:-4]+"_scalebar.svg"} to combined_final_out.svg')
combine_svgs(file_path, file_path[:-4]+'_scalebar.svg', 'combined_final_out.svg')
#icmd = 'python combine_svg.py %s %s X.%s' % (file_path, file_path[:-4]+'_scalebar.svg', file_path)
#print(cmd) 
#os.system(cmd)
