import psi4
mol = psi4.geometry("""
He
""")
psi4.energy('hf/cc-pvdz')
psi4.compare_values(-2.85518839, psi4.core.get_variable('current energy'), 5, 'SCF E')
