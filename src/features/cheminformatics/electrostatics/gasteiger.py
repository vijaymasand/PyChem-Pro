"""
Gasteiger-Marsili Partial Charge Calculation

Iterative Partial Equalization of Orbital Electronegativity (PEOE).
Computes partial charges using electronegativity equalization without
requiring a force field or semi-empirical method.

Reference: Gasteiger J, Marsili M. Tetrahedron. 1980;36(22):3219-3228.
"""

import math


# Electronegativity parameters: a, b, c for χ = a + b*q + c*q²
# These are the Hinze-Jaffe orbital electronegativity coefficients
# Format: element_symbol -> {hybridization: (a, b, c)}
_EN_PARAMS = {
    'H': {
        's': (7.17, 6.24, -0.56),
        'default': (7.17, 6.24, -0.56),
    },
    'C': {
        'sp3': (7.98, 9.18, 1.88),
        'sp2': (8.79, 9.32, 1.51),
        'sp':  (10.39, 9.45, 0.73),
        'default': (7.98, 9.18, 1.88),
    },
    'N': {
        'sp3': (11.54, 10.82, 1.36),
        'sp2': (12.87, 11.15, 0.85),
        'sp':  (15.68, 11.70, -0.27),
        'default': (11.54, 10.82, 1.36),
    },
    'O': {
        'sp3': (14.18, 12.92, 1.39),
        'sp2': (17.07, 13.79, 0.47),
        'default': (14.18, 12.92, 1.39),
    },
    'F': {
        'default': (14.66, 13.85, 2.31),
    },
    'Cl': {
        'default': (11.00, 9.69, 1.35),
    },
    'Br': {
        'default': (10.08, 8.47, 1.16),
    },
    'I': {
        'default': (9.90, 7.96, 0.96),
    },
    'S': {
        'sp3': (10.14, 9.13, 1.38),
        'sp2': (10.88, 9.49, 1.33),
        'default': (10.14, 9.13, 1.38),
    },
    'P': {
        'sp3': (8.90, 8.24, 0.96),
        'default': (8.90, 8.24, 0.96),
    },
    'B': {
        'sp2': (6.88, 7.28, 1.71),
        'sp3': (6.10, 6.82, 1.56),
        'default': (6.88, 7.28, 1.71),
    },
    'Si': {
        'default': (7.30, 6.56, 0.70),
    },
    'Se': {
        'default': (10.08, 8.47, 1.16),
    },
    'Na': {
        'default': (2.84, 3.80, -0.34),
    },
    'K': {
        'default': (2.42, 3.26, -0.47),
    },
    'Ca': {
        'default': (3.30, 4.36, -0.17),
    },
    'Mg': {
        'default': (3.75, 4.72, -0.12),
    },
    'Fe': {
        'default': (6.00, 6.00, 0.50),
    },
    'Zn': {
        'default': (5.80, 5.90, 0.50),
    },
}


def compute_gasteiger_charges(molecule, num_iterations=6):
    """
    Calculate Gasteiger-Marsili partial charges using PEOE method.

    The algorithm iteratively equalizes electronegativity by transferring
    charge along bonds. At each iteration:
        dq = damping^iter * (chi_acceptor - chi_donor) / chi_plus_donor

    Args:
        molecule: Molecule with atoms and bonds
        num_iterations: Number of charge equilibration iterations (default: 6)

    Updates:
        Each atom's partial_charge attribute
    """
    n = len(molecule.atoms)
    if n == 0:
        return

    # Initialize charges from formal charges
    charges = [float(atom.formal_charge) for atom in molecule.atoms]

    # Get electronegativity parameters for each atom
    en_params = []
    for atom in molecule.atoms:
        params = _get_en_params(atom)
        en_params.append(params)

    # Iterative charge equilibration
    for iteration in range(num_iterations):
        damping = 0.5 ** (iteration + 1)

        # Calculate current electronegativity for each atom: χ = a + b*q + c*q²
        chi = []
        for i in range(n):
            a, b, c = en_params[i]
            q = charges[i]
            electronegativity = a + b * q + c * q * q
            chi.append(electronegativity)

        # Transfer charge along each bond
        charge_transfer = [0.0] * n

        for bond in molecule.bonds:
            i = bond.begin_atom_idx
            j = bond.end_atom_idx

            chi_i = chi[i]
            chi_j = chi[j]

            if abs(chi_i - chi_j) < 1e-10:
                continue

            # Determine direction: charge flows from less EN to more EN
            if chi_i > chi_j:
                # Atom i is more electronegative → pulls electrons from j
                # j is the donor
                acceptor = i
                donor = j
            else:
                acceptor = j
                donor = i

            # chi+ of the donor: evaluate EN at q=+1
            a_d, b_d, c_d = en_params[donor]
            chi_plus_donor = a_d + b_d + c_d  # χ(q=+1)

            if chi_plus_donor <= 0:
                continue

            # Charge transferred: dq = damping * Δχ / χ+_donor
            delta_chi = abs(chi_i - chi_j)
            dq = damping * delta_chi / chi_plus_donor

            # Clamp to prevent runaway
            dq = min(dq, 0.1)

            # Transfer: donor loses electron density (becomes more +)
            # acceptor gains electron density (becomes more -)
            charge_transfer[donor] += dq       # donor becomes more positive
            charge_transfer[acceptor] -= dq    # acceptor becomes more negative

        # Apply transfers
        for i in range(n):
            charges[i] += charge_transfer[i]

    # Normalize to maintain total formal charge
    total_formal = sum(atom.formal_charge for atom in molecule.atoms)
    current_total = sum(charges)
    if n > 0 and abs(current_total - total_formal) > 1e-6:
        correction = (total_formal - current_total) / n
        for i in range(n):
            charges[i] += correction

    # Assign to atoms
    for i, atom in enumerate(molecule.atoms):
        atom.partial_charge = round(charges[i], 4)


def _get_en_params(atom):
    """Get electronegativity parameters (a, b, c) for an atom."""
    symbol = atom.symbol
    params_dict = _EN_PARAMS.get(symbol)

    if params_dict is None:
        # Fallback: use Pauling electronegativity to estimate
        en = atom.element.electronegativity
        if en > 0:
            return (en * 2.5, en * 2.0, 0.5)
        return (5.0, 5.0, 0.5)  # Generic fallback

    # Try to get hybridization-specific parameters
    hyb = atom.hybridization or 'default'
    if hyb in params_dict:
        return params_dict[hyb]
    return params_dict['default']
