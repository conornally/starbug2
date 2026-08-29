import graphviz

# Initialize the directed graph
dot = graphviz.Digraph(
    'starbug_pipeline',
    comment='Starbug2 Pipeline Flow',
    format='png'
)

# Global Graph attributes
dot.attr(rankdir='TB')
dot.attr(splines='ortho')
dot.attr(nodesep='0.6')
dot.attr(ranksep='0.4')
dot.attr(bgcolor='transparent')

# Global Node attributes (Light blue fill, dark blue border, black text)
dot.attr('node',
         fontname='Courier New, monospace',
         fontsize='12',
         fontcolor='#000000',
         style='filled',
         fillcolor='#E3F2FD',
         color='#1565C0',
         penwidth='1.5')

# Global Edge attributes (Solid black arrows)
dot.attr('edge',
         fontname='Courier New, monospace',
         fontsize='10',
         color='#000000',
         arrowsize='1.0',
         penwidth='1.5')

# --- Define Nodes ---
# Terminals / Commands (Pill-shaped/Rounded boxes)
dot.node('init_param', '$~ starbug2 --local-param', shape='box',
         style='filled,rounded')
dot.node('cmd_D', '$~ starbug2 -D', shape='box', style='filled,rounded')
dot.node('cmd_B', '$~ starbug2 -B', shape='box', style='filled,rounded')
dot.node('cmd_P', '$~ starbug2 -P', shape='box', style='filled,rounded')

# Files (Document/Note shaped icons)
dot.node('image_fits', 'image.fits', shape='note')
dot.node('starbug_param', 'starbug.param', shape='note')
dot.node('image_ap', 'image-ap.fits\n(TABLE)', shape='note')
dot.node('image_bg', 'image-bg.fits\n(IMAGE)', shape='note')
dot.node('image_res', 'image-res.fits\n(IMAGE)', shape='note')
dot.node('image_psf', 'image-psf.fits\n(TABLE)', shape='note')

# --- Structural Positioning & Constraints ---
# Force horizontal alignments for row uniformity
with dot.subgraph() as s:
    s.attr(rank='same')
    s.node('image_fits')
    s.node('starbug_param')

with dot.subgraph() as s:
    s.attr(rank='same')
    s.node('cmd_D')
    s.node('image_ap')

with dot.subgraph() as s:
    s.attr(rank='same')
    s.node('cmd_B')
    s.node('image_bg')

with dot.subgraph() as s:
    s.attr(rank='same')
    s.node('cmd_P')
    s.node('image_res')

# --- Define Connections ---
# Parameter branch
dot.edge('init_param', 'starbug_param')


dot.edge('starbug_param', 'cmd_D')

# Main data flow
dot.edge('image_fits', 'cmd_D')
dot.edge('cmd_D', 'image_ap')

dot.edge('cmd_D', 'cmd_B')
dot.edge('cmd_B', 'image_bg')

dot.edge('cmd_B', 'cmd_P')
dot.edge('cmd_P', 'image_res')

dot.edge('cmd_P', 'image_psf')

# Save and compile to disk
dot.render('../images/flow', cleanup=True)
print("Pipeline chart compiled successfully as flow.png!")
