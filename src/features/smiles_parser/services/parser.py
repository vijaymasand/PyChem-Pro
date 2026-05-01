"""
SMILES Parser — Recursive descent parser that builds a Molecule from SMILES tokens.

Implements the full SMILES grammar:
    smiles    → chain ('.' chain)*
    chain     → atom (bond? (atom | branch | ring_closure))*
    branch    → '(' bond? chain ')'
"""

from src.core.domain.models.atom import Atom, Chirality
from src.core.domain.models.bond import Bond, BondType, BondStereo, BOND_ORDER_MAP
from src.core.domain.models.molecule import Molecule
from src.features.smiles_parser.services.tokenizer import tokenize, tokenize_regex, TokenType, SMILESTokenizerError
from src.features.smiles_parser.rules.valence import calculate_implicit_hydrogens
from src.features.smiles_parser.rules.aromaticity import kekulize, perceive_aromaticity
from src.features.smiles_parser.rules.stereo import assign_stereo


from src.core.domain.models.atom import Atom, Chirality
from src.core.domain.models.bond import Bond, BondType, BondStereo, BOND_ORDER_MAP
from src.core.domain.models.molecule import Molecule
import src.vendors.oasa.smiles as oasa_smiles

class SMILESParseError(Exception):
    """Error during SMILES parsing."""
    pass

def parse_smiles(smiles, name="", use_bkchem_tokenizer=True):
    """
    Parse a SMILES string into a Molecule object natively using the robust
    BKChem (OASA) deterministic graph tokenization engine.
    """
    if not smiles or not smiles.strip():
        raise SMILESParseError("Empty SMILES string")

    try:
        # 1. Parse SMILES via OASA and generate 2D Coords natively at parsing-time!
        oasa_mol = oasa_smiles.text_to_mol(smiles.strip(), calc_coords=1, localize_aromatic_bonds=True)
    except Exception as e:
        raise SMILESParseError(f"OASA SMILES Parsing Error: {str(e)}")

    if not oasa_mol:
        raise SMILESParseError("Failed to parse SMILES string via OASA")

    # 2. Bridge OASA Molecule to Domain Molecule
    mol = Molecule(name)
    
    # Map OASA vertices to Domain Atoms
    atom_map = {} # oasa_vertex -> domain_atom_idx
    for i, o_v in enumerate(oasa_mol.vertices):
        # Translate properties
        is_aromatic = o_v.properties_.get('aromatic', 0) == 1
        charge = o_v.charge
        isotope = getattr(o_v, 'isotope', 0)
        hcount = getattr(o_v, 'explicit_hydrogens', 0)
        
        # Translate OASA tetrahedral stereo to internal models
        stereo = o_v.properties_.get('stereo', '')
        chirality = Chirality.NONE
        if stereo == '@':
            chirality = Chirality.COUNTERCLOCKWISE
        elif stereo == '@@':
            chirality = Chirality.CLOCKWISE
            
        atom = Atom(
            symbol=o_v.symbol,
            is_aromatic=is_aromatic,
            formal_charge=charge,
            isotope=isotope,
            chirality=chirality,
            num_explicit_h=hcount
        )
        
        # Hydrate 2D representation generated natively by OASA
        if hasattr(o_v, 'x') and o_v.x is not None and hasattr(o_v, 'y') and o_v.y is not None:
            atom.x2d = float(o_v.x)
            atom.y2d = float(o_v.y)
            atom.z2d = getattr(o_v, 'z', 0.0)

        # Implicit Hydrogens
        # Oasa computes free_valency; internal computes from hybridization context
        atom.num_implicit_h = o_v.free_valency if hasattr(o_v, 'free_valency') and o_v.free_valency > 0 else 0
        
        atom_idx = mol.add_atom(atom)
        atom_map[o_v] = atom_idx

    # Map OASA edges to Domain Bonds
    for o_e in oasa_mol.edges:
        v1, v2 = o_e.vertices
        idx1 = atom_map[v1]
        idx2 = atom_map[v2]
        
        # OASA kekulizes aromatic bonds (sets order to 1 or 2) but
        # keeps aromatic=True.  We must respect the kekulized order
        # for proper rendering, and only fall back to AROMATIC when
        # the order itself is unresolved (order=4).
        btype = BondType.SINGLE
        order = getattr(o_e, 'order', 1)
        if order == 2:
            btype = BondType.DOUBLE
        elif order == 3:
            btype = BondType.TRIPLE
        elif order == 4:
            btype = BondType.AROMATIC
        # Note: aromatic bonds with order=1 stay SINGLE (kekulized form)
            
        # Stereo
        bstereo = BondStereo.NONE
        stereo_tag = o_e.properties_.get('stereo', '')
        if stereo_tag == '/':
            bstereo = BondStereo.UP
        elif stereo_tag == '\\':
            bstereo = BondStereo.DOWN
            
        bond_obj = mol.add_bond(idx1, idx2, btype, bstereo)
        # Preserve aromatic flag for downstream chemistry even though
        # the bond type is kekulized (single/double)
        if o_e.aromatic and bond_obj:
            bond_obj.is_in_ring = True
    # Perceive atom aromaticity (needed for fingerprinting, descriptors)
    # but preserve the kekulized bond types from OASA.  perceive_aromaticity
    # re-marks bonds as AROMATIC, so we save and restore bond types.
    saved_bond_types = [(b.index, b.bond_type) for b in mol.bonds]

    from src.features.smiles_parser.rules.aromaticity import perceive_aromaticity
    mol.find_rings()
    perceive_aromaticity(mol)

    # Restore kekulized bond types (single/double from OASA)
    for bond_idx, btype in saved_bond_types:
        mol.bonds[bond_idx].bond_type = btype

    mol.assign_hybridization()
    
    # Tag source for optimized 2D layout engine selection
    mol.properties['source'] = 'smiles'

    return mol
