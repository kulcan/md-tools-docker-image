import psi4

psi4.core.set_output_file("output.out")
psi4.core.set_num_threads(4)

mol = psi4.geometry("""
He
""")

custom_functional = {
    "name": "my_PBE",
    "x_functionals": {"GGA_X_PBE": {"alpha": 0.75}},
    "x_hf": {"alpha": 0.25},
    "c_functionals": {"GGA_C_PBE": {}}
}

psi4.energy('scf/dzvp',  dft_functional = custom_functional, molecule = mol)