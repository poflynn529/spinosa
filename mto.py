# Module Imports
import numpy as np

# Molecular mass (g/mol)
mol_mass = {"h2o"   : 18.02, # Water
            "ch3oh" : 32.04, # Methanol
            "c2h4"  : 28.05, # Ethylene
            "c3h6"  : 42.08, # Propylene
            "h2"    : 2.02,  # Hydrogen
            "co2"   : 44.01, # Carbon Dioxide
            "o2"    : 32.00  # Oxygen
}

# Specific Heat Capacity (J/mol K)
shc = {"h2o"     : 74.893,
       "ch3oh"   : 67.86,
       "h2"      : 29.19,
       "co2"     : 41.81}

# Reaction temps (deg C)
mto_temp = 400
hydrogenation_init_temp = 200

# Reaction pressures (bar)
hydrogenation_final_pressure = 150

# Other Reaction Constants
electrolysis_constant = 0.6e6 # J/g

# Molar Mass to Mass at STP
h20_mass_stp = None
ch3oh_mass_stp = None

# Universal Constants
uni_gas_constant = 8.314 # J/mol K
celsius_to_kelvin = 273.15 # K

# Calcination
# Electrolysis
# Hydrogenation
# MTO
# Cryo cooling heat transfer
# Polymerisation


#3H2, 1 Co2

def f_electrolysis(h2o_input):
    
    # Energy
    energy = h2o_input * mol_mass["h2"] * electrolysis_constant

    # Products
    products = {"h2" : h2o_input,
                "o2" : h2o_input * 0.5}
    
    return energy, products

def f_hydrogenation(h2_input, co2_input, initial_temp):

    ## Energy
    # Initial temperature
    delta_temp = hydrogenation_init_temp - initial_temp
    total_energy = (h2_input * shc["h2"] + co2_input * shc["co2"]) * delta_temp
    # Compression Cost E = PV
    total_energy += ideal_compressor(h2_input, 1, hydrogenation_final_pressure, initial_temp)
    total_energy += ideal_compressor(co2_input, 1, hydrogenation_final_pressure, initial_temp)

    ## Products
    products = {"ch3oh" : co2_input,
                "h2o" : co2_input}

    return total_energy, products

def f_mto(h20_input, ch3oh_input, initial_temp):

    ## Energy
    delta_temp = mto_temp - initial_temp
    total_energy = (h20_input * shc["h2o"] + ch3oh_input * shc["ch3oh"]) * delta_temp

    ## Products
    products = {"c2h4" : ch3oh_input*0.4,
                "c3h6" : ch3oh_input*0.4, 
                "h2o"  : h20_input + ch3oh_input*0.2}

    return total_energy, products

def ideal_compressor(mols, init_pressure, final_pressure, temp):
    return mols * uni_gas_constant * (temp + celsius_to_kelvin) * np.log(final_pressure / init_pressure)

##### Tests #####

### Test Functions ###

def print_reaction(energy, products, title=None):
    
    print("\n------------------------------------")
    if title != None:
        print(f"Test: {title}\n")
    
    print(f"Total Energy: {round(energy / 1000, 3)} kJ\n")
    print("Products:")

    for key, value in products.items():
        print(f"\t{key}: {round(value * mol_mass[key] / 1000, 3)} kg")

    print("------------------------------------\n")

### Individual Tests ###

def test_electrolysis():
    # 1kg of h2o
    h2o_input = 1000 / mol_mass["h2o"]
    energy, products = f_electrolysis(h2o_input)
    print_reaction(energy, products, "Electrolysis Only")

def test_mto():
    # 1kg of methanol
    # 200g of water
    # Initial temp of 300C
    energy, products = f_mto(200 / mol_mass["h2o"], 1000 / mol_mass["ch3oh"], 300)
    print_reaction(energy, products, "MTO Only")

def test_hydrogenation():
    # 1kg of c02
    # Stoichiometric amount of h2
    # Initial temp of 25C
    co2_input = 1000 / mol_mass["co2"]
    h2_input = 3 * co2_input
    energy, products = f_hydrogenation(h2_input, co2_input, 25)

    print_reaction(energy, products, "Hydrogenation Only")

### Full System Test ###

def test_system():

    # 1 kg of water
    h2o_input = 1000 / mol_mass["h2o"]
    
    h2o_input_extra = 300 / mol_mass["h2o"]

    total_energy = 0

    energy, products = f_electrolysis(h2o_input)
    total_energy += energy
    co2_input = products["h2"] / 3

    energy, products = f_hydrogenation(products["h2"], co2_input, 25)
    total_energy += energy
    energy, products = f_mto(h2o_input_extra, products["ch3oh"], 300)
    total_energy += energy

    print_reaction(total_energy, products, "Full System Test")

##### Main #####

test_electrolysis()
test_hydrogenation()
test_mto()
test_system()


##### Assumptions #####

# 1. Ideal compressor, no energy loss or temp increase of gas.
# 2. Hydrogen quantity is not stoichiometric.