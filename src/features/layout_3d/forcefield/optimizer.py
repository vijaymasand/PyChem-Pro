"""
Geometry Optimizer — Energy minimization of 3D molecular structures.

Implements a simple energy model:
- Bond stretching: E = k(r - r₀)²
- Angle bending: E = k(θ - θ₀)²
- VdW repulsion: E = ε(σ/r)¹² (repulsive term only)

Optimizers:
- Steepest descent (robust, for initial cleanup)
- L-BFGS (fast convergence for refinement)
"""

import numpy as np
from src.core.domain.models.elements import bond_length as ideal_bond_length


# Force constants
K_BOND = 300.0    # kcal/mol/Å² — bond stretching
K_ANGLE = 50.0    # kcal/mol/rad² — angle bending
K_VDW = 0.1       # kcal/mol — VdW repulsion epsilon


def optimize_geometry(molecule, max_steps=500, convergence=1e-4, method='steepest_descent'):
    """
    Optimize molecular geometry by minimizing a simple energy function.

    Args:
        molecule: Molecule with 3D coordinates
        max_steps: Maximum optimization steps
        convergence: RMS gradient threshold for convergence
        method: 'steepest_descent' or 'lbfgs'

    Returns:
        (final_energy, converged, num_steps)
    """
    # Get current coordinates
    coords = _get_coords(molecule)
    if coords is None:
        return 0.0, False, 0

    # Build lists of interactions
    bond_list = _build_bond_list(molecule)
    angle_list = _build_angle_list(molecule)
    vdw_pairs = _build_vdw_list(molecule)

    if method == 'lbfgs':
        coords, energy, converged, steps = _lbfgs_optimize(
            coords, bond_list, angle_list, vdw_pairs, max_steps, convergence)
    else:
        coords, energy, converged, steps = _steepest_descent(
            coords, bond_list, angle_list, vdw_pairs, max_steps, convergence)

    # Update molecule coordinates
    _set_coords(molecule, coords)

    return energy, converged, steps


def compute_energy(molecule):
    """Compute total energy of current geometry."""
    coords = _get_coords(molecule)
    if coords is None:
        return 0.0
    bond_list = _build_bond_list(molecule)
    angle_list = _build_angle_list(molecule)
    vdw_pairs = _build_vdw_list(molecule)
    return _total_energy(coords, bond_list, angle_list, vdw_pairs)


# ─── Energy Functions ───────────────────────────────────────────────

def _total_energy(coords, bond_list, angle_list, vdw_pairs):
    """Compute total energy."""
    e = 0.0
    e += _bond_energy(coords, bond_list)
    e += _angle_energy(coords, angle_list)
    e += _vdw_energy(coords, vdw_pairs)
    return e


def _bond_energy(coords, bond_list):
    """Bond stretching energy: E = k(r - r₀)²"""
    energy = 0.0
    for i, j, r0 in bond_list:
        r = np.linalg.norm(coords[i] - coords[j])
        energy += K_BOND * (r - r0) ** 2
    return energy


def _angle_energy(coords, angle_list):
    """Angle bending energy: E = k(θ - θ₀)²"""
    energy = 0.0
    for i, j, k, theta0 in angle_list:
        theta = _compute_angle(coords[i], coords[j], coords[k])
        dtheta = theta - theta0
        energy += K_ANGLE * dtheta ** 2
    return energy


def _vdw_energy(coords, vdw_pairs):
    """VdW repulsion energy: E = ε(σ/r)¹²"""
    energy = 0.0
    for i, j, sigma in vdw_pairs:
        r = np.linalg.norm(coords[i] - coords[j])
        if r < sigma * 1.5 and r > 0.01:  # Only short-range repulsion
            ratio = sigma / r
            energy += K_VDW * ratio ** 12
    return energy


# ─── Gradient Functions ─────────────────────────────────────────────

def _total_gradient(coords, bond_list, angle_list, vdw_pairs):
    """Compute total gradient (analytical)."""
    n = len(coords)
    grad = np.zeros((n, 3))

    # Bond stretching gradient
    for i, j, r0 in bond_list:
        diff = coords[i] - coords[j]
        r = np.linalg.norm(diff)
        if r < 1e-8:
            continue
        direction = diff / r
        force = 2.0 * K_BOND * (r - r0)
        grad[i] += force * direction
        grad[j] -= force * direction

    # Angle bending gradient (numerical for stability)
    h = 1e-5
    for i, j, k, theta0 in angle_list:
        for atom_idx in (i, j, k):
            for dim in range(3):
                coords[atom_idx, dim] += h
                e_plus = K_ANGLE * (_compute_angle(coords[i], coords[j], coords[k]) - theta0) ** 2
                coords[atom_idx, dim] -= 2 * h
                e_minus = K_ANGLE * (_compute_angle(coords[i], coords[j], coords[k]) - theta0) ** 2
                coords[atom_idx, dim] += h
                grad[atom_idx, dim] += (e_plus - e_minus) / (2 * h)

    # VdW gradient
    for i, j, sigma in vdw_pairs:
        diff = coords[i] - coords[j]
        r = np.linalg.norm(diff)
        if r < sigma * 1.5 and r > 0.01:
            direction = diff / r
            force = -12.0 * K_VDW * (sigma ** 12) / (r ** 13)
            grad[i] += force * direction
            grad[j] -= force * direction

    return grad


# ─── Optimizers ─────────────────────────────────────────────────────

def _steepest_descent(coords, bond_list, angle_list, vdw_pairs,
                      max_steps, convergence):
    """Steepest descent optimization with adaptive step size."""
    step_size = 0.01
    energy = _total_energy(coords, bond_list, angle_list, vdw_pairs)

    for step in range(max_steps):
        grad = _total_gradient(coords, bond_list, angle_list, vdw_pairs)
        rms_grad = np.sqrt(np.mean(grad ** 2))

        if rms_grad < convergence:
            return coords, energy, True, step + 1

        # Normalize and step
        direction = -grad
        max_grad = np.max(np.abs(direction))
        if max_grad > 0:
            scale = min(step_size, 0.1 / max_grad)
        else:
            break

        new_coords = coords + scale * direction
        new_energy = _total_energy(new_coords, bond_list, angle_list, vdw_pairs)

        if new_energy < energy:
            coords = new_coords
            energy = new_energy
            step_size *= 1.2  # Increase step size on success
        else:
            step_size *= 0.5  # Decrease on failure

        if step_size < 1e-10:
            break

    return coords, energy, False, max_steps


def _lbfgs_optimize(coords, bond_list, angle_list, vdw_pairs,
                    max_steps, convergence, m=10):
    """L-BFGS optimization."""
    n = len(coords)
    flat_coords = coords.flatten()

    def energy_func(x):
        c = x.reshape(-1, 3)
        return _total_energy(c, bond_list, angle_list, vdw_pairs)

    def grad_func(x):
        c = x.reshape(-1, 3)
        return _total_gradient(c, bond_list, angle_list, vdw_pairs).flatten()

    # L-BFGS storage
    s_list = []  # position differences
    y_list = []  # gradient differences
    rho_list = []

    energy = energy_func(flat_coords)
    grad = grad_func(flat_coords)

    for step in range(max_steps):
        rms_grad = np.sqrt(np.mean(grad ** 2))
        if rms_grad < convergence:
            return flat_coords.reshape(-1, 3), energy, True, step + 1

        # Compute search direction using L-BFGS two-loop recursion
        q = grad.copy()
        alphas = []

        for i in range(len(s_list) - 1, -1, -1):
            alpha = rho_list[i] * np.dot(s_list[i], q)
            alphas.insert(0, alpha)
            q -= alpha * y_list[i]

        # Initial Hessian approximation
        if s_list:
            gamma = np.dot(s_list[-1], y_list[-1]) / np.dot(y_list[-1], y_list[-1])
            r = gamma * q
        else:
            r = q * 0.01

        for i in range(len(s_list)):
            beta = rho_list[i] * np.dot(y_list[i], r)
            r += s_list[i] * (alphas[i] - beta)

        direction = -r

        # Line search (simple backtracking)
        step_size = 1.0
        c1 = 1e-4
        dir_deriv = np.dot(grad, direction)

        for _ in range(20):
            new_coords = flat_coords + step_size * direction
            new_energy = energy_func(new_coords)

            if new_energy <= energy + c1 * step_size * dir_deriv:
                break
            step_size *= 0.5
        else:
            break

        # Update
        s = new_coords - flat_coords
        new_grad = grad_func(new_coords)
        y = new_grad - grad

        sy = np.dot(s, y)
        if sy > 1e-10:
            s_list.append(s)
            y_list.append(y)
            rho_list.append(1.0 / sy)

            if len(s_list) > m:
                s_list.pop(0)
                y_list.pop(0)
                rho_list.pop(0)

        flat_coords = new_coords
        energy = new_energy if 'new_energy' in dir() else energy_func(flat_coords)
        grad = new_grad

    return flat_coords.reshape(-1, 3), energy, False, max_steps


# ─── Helper Functions ───────────────────────────────────────────────

def _compute_angle(a, b, c):
    """Compute angle at atom b (in radians) given positions a, b, c."""
    v1 = a - b
    v2 = c - b
    cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-10)
    cos_angle = np.clip(cos_angle, -1.0, 1.0)
    return np.arccos(cos_angle)


def _get_coords(molecule):
    """Extract coordinates as Nx3 numpy array."""
    n = len(molecule.atoms)
    if n == 0:
        return None
    coords = np.zeros((n, 3))
    for i, atom in enumerate(molecule.atoms):
        if atom.has_coords:
            coords[i] = [atom.x, atom.y, atom.z]
        else:
            return None  # Not all atoms have coordinates
    return coords


def _set_coords(molecule, coords):
    """Set coordinates from Nx3 numpy array."""
    for i, atom in enumerate(molecule.atoms):
        atom.x = float(coords[i, 0])
        atom.y = float(coords[i, 1])
        atom.z = float(coords[i, 2])


def _build_bond_list(molecule):
    """Build list of (atom_i, atom_j, ideal_length) for bond stretching."""
    bond_list = []
    for bond in molecule.bonds:
        i = bond.begin_atom_idx
        j = bond.end_atom_idx
        r0 = ideal_bond_length(molecule.atoms[i].symbol,
                                molecule.atoms[j].symbol,
                                bond.order)
        bond_list.append((i, j, r0))
    return bond_list


def _build_angle_list(molecule):
    """Build list of (atom_i, atom_j, atom_k, ideal_angle) for angle bending."""
    angle_list = []
    for atom in molecule.atoms:
        neighbors = molecule.get_neighbors(atom.index)
        if len(neighbors) < 2:
            continue

        # Ideal angle based on hybridization
        hyb = atom.hybridization
        if hyb == 'sp':
            theta0 = np.radians(180.0)
        elif hyb == 'sp2':
            theta0 = np.radians(120.0)
        else:
            theta0 = np.radians(109.5)

        for k in range(len(neighbors)):
            for l in range(k + 1, len(neighbors)):
                angle_list.append((neighbors[k], atom.index, neighbors[l], theta0))

    return angle_list


def _build_vdw_list(molecule):
    """Build list of (atom_i, atom_j, sigma) for VdW repulsion (1-4+ pairs only)."""
    n = len(molecule.atoms)
    # Find atoms connected by 1 or 2 bonds (exclude from VdW)
    excluded = set()
    for bond in molecule.bonds:
        excluded.add((min(bond.begin_atom_idx, bond.end_atom_idx),
                       max(bond.begin_atom_idx, bond.end_atom_idx)))

    # Also exclude 1-3 pairs
    for atom in molecule.atoms:
        neighbors = molecule.get_neighbors(atom.index)
        for k in range(len(neighbors)):
            for l in range(k + 1, len(neighbors)):
                pair = (min(neighbors[k], neighbors[l]),
                        max(neighbors[k], neighbors[l]))
                excluded.add(pair)

    vdw_pairs = []
    for i in range(n):
        for j in range(i + 1, n):
            if (i, j) not in excluded:
                sigma_i = molecule.atoms[i].element.vdw_radius
                sigma_j = molecule.atoms[j].element.vdw_radius
                sigma = (sigma_i + sigma_j) * 0.5
                vdw_pairs.append((i, j, sigma))

    return vdw_pairs
