import cantera as ct

# Molecular mass (g/mol)
h2o_mol = 18.01528
ch3oh_mol = 32.04

# Specific Heat Capacity (J/mol K)
h20_shc = 74.893
ch3oh_shc = 67.86

# Molar Mass to Mass at STP
h20_mass_stp = None
ch3oh_mass_stp = None

total_energy = (h20_shc + ch3oh_shc) * 100 / (ch3oh_mol / 1000) * 1.429

print(total_energy)