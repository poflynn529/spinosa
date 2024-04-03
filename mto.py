# Module Imports
import numpy as np
import pickle as pkl
import pandas as pd

# Molecular mass (g/mol)
mol_mass = {"h2o"     : 18.02,  # Water
            "ch3oh"   : 32.04,  # Methanol
            "c2h4"    : 28.05,  # Ethylene
            "c3h6"    : 42.08,  # Propylene
            "h2"      : 2.02,   # Hydrogen
            "co2"     : 44.01,  # Carbon Dioxide
            "o2"      : 32.00,  # Oxygen
            "caco3"   : 100.09, # Calcium Carbonate
            "ca(oh)2" : 74.09,  # Calcium Hydroxide
            "cao"     : 56.08,  # Calcium Oxide
            "n2"      : 28.02,  # Nitrogen
}

# Specific Heat Capacity (J/mol K)
shc = {"h2o"      : 74.893,
       "ch3oh"    : 67.86,
       "h2"       : 29.19,
       "co2"      : 41.81,
       "caco3"    : 83.51,
}

# Air Mixture Ratios (Mole Fractions)
air_mix = {
    "n2"  : 0.7808,
    "o2"  : 0.2095,
    "co2" : 0.0004,
    "h2o" : 0.0100,
}

# Reaction temps (deg C)
mto_temp = 400
hydrogenation_init_temp = 200
co2_concentration_temp = 900

# 2.27 to 1 CO2

# Reaction pressures (bar)
hydrogenation_final_pressure = 150

# Other Reaction Constants
electrolysis_constant = 0.6e6 # J/g h2

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

def f_air_capture():
    return

def f_co2_concentration(caco3_input, initial_temp):
    
    # Energy
    energy = caco3_input * shc["caco3"] * (co2_concentration_temp - initial_temp)

    # Products
    products = {"co2" : caco3_input,
                "cao" : caco3_input}
    
    return energy, products

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
    products = {"c2h4"  : (ch3oh_input * 0.45) / 2,
                "c3h6"  : (ch3oh_input * 0.45) / 3, 
                "h2o"   : h20_input,
                "ch3oh" : ch3oh_input*0.24}

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
        print(f"\t{key}: {round(value * mol_mass[key] / 1000, 3)} kg, {round(value, 3)} mol")

    print("------------------------------------\n")

### Individual Tests ###

def test_co2_concentration():

    # 1kg of co2
    co2_input = 1000 / mol_mass["co2"]
    caco3_input = co2_input
    energy, products = f_co2_concentration(co2_input, caco3_input)
    print_reaction(energy, products, "CO2 Concentration Only")

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

    reactions = {
        "Category" : [],
        "Value"    : []
    }

    total_energy = 0

    # Co2, CaCO3 initial products in mol
    co2_input = 1000 / mol_mass["co2"]
    caco3_input = co2_input

    energy, products = f_co2_concentration(caco3_input, co2_input)
    total_energy += energy
    reactions["Category"].append(r"$CO_{2}$ Concentration")
    reactions["Value"].append(energy / 1e3)

    # Water input
    h2o_input = co2_input
    
    h2o_input_extra = 300 / mol_mass["h2o"]

    energy, products = f_electrolysis(h2o_input)
    total_energy += energy
    reactions["Category"].append("Electrolysis")
    reactions["Value"].append(energy / 1e3)

    energy, products = f_hydrogenation(products["h2"], co2_input, 25)
    total_energy += energy
    reactions["Category"].append("Hydrogenation")
    reactions["Value"].append(energy / 1e3)

    energy, products = f_mto(h2o_input_extra, products["ch3oh"], 300)
    total_energy += energy
    reactions["Category"].append("MTO")
    reactions["Value"].append(energy / 1e3)

    df = pd.DataFrame(reactions)
    df["Value"] = df["Value"] * 3.5
    df.to_pickle("results/system0.pkl")

    print_reaction(energy, products, "Full System Test")

##### Main #####

test_co2_concentration()
test_electrolysis()
test_hydrogenation()
test_mto()
test_system()


##### Assumptions #####

# 1. Ideal compressor, no energy loss or temp increase of gas.