# PyMOL script for FST344 structural analysis
# Load structure
load alphafold_fst_structure.pdb, fst344

# Basic representation
hide all
show cartoon, fst344
color lightblue, fst344

# Highlight functional domains
select nd1_domain, resi 30-120
select nd2_domain, resi 121-205  
select nd3_domain, resi 206-290
select c_terminal, resi 291-344

color yellow, nd1_domain
color orange, nd2_domain
color red, nd3_domain
color purple, c_terminal

# Highlight activin binding sites
select binding_site_1, resi 58-85
select binding_site_2, resi 145-175
select binding_site_3, resi 220-250

color green, binding_site_1
color green, binding_site_2  
color green, binding_site_3
show sticks, binding_site_1 or binding_site_2 or binding_site_3

# Highlight critical cysteines
select disulfide_cysteines, resi 39+43+56+69+90+104+120+135+148+162+177+191+205+218+231+245+260+274
color yellow, disulfide_cysteines
show spheres, disulfide_cysteines

# Set view
orient
zoom

# Labels
label disulfide_cysteines and name CA, "%s%s" % (resn,resi)

print "FST344 structure loaded with functional domains highlighted"
print "Green: Activin binding sites"
print "Yellow spheres: Disulfide cysteines"
print "Color domains: ND1(yellow), ND2(orange), ND3(red), C-term(purple)"
