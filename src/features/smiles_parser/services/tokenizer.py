"""
SMILES Tokenizer — Lexes a SMILES string into a stream of typed tokens.

Supports full OpenSMILES specification:
- Organic subset atoms: B, C, N, O, P, S, F, Cl, Br, I, and their aromatic variants
- Bracket atoms: [isotope?symbol chirality? hcount? charge? class?]
- Bond tokens: - = # : / \
- Branches: ( )
- Ring closures: digits 0-9 and %nn
- Dot disconnection: .
"""

import re
from src.core.domain.models.bond import BondType, BondStereo, BOND_ORDER_MAP


class TokenType:
    """Token type constants."""
    ATOM = 'ATOM'
    BRACKET_ATOM = 'BRACKET_ATOM'
    BOND = 'BOND'
    BRANCH_OPEN = 'BRANCH_OPEN'
    BRANCH_CLOSE = 'BRANCH_CLOSE'
    RING_CLOSURE = 'RING_CLOSURE'
    DOT = 'DOT'


class Token:
    """A single SMILES token with type and associated data."""

    __slots__ = ('type', 'value', 'data', 'pos')

    def __init__(self, token_type, value, data=None, pos=0):
        self.type = token_type
        self.value = value
        self.data = data or {}
        self.pos = pos

    def __repr__(self):
        return f"Token({self.type}, '{self.value}', data={self.data})"


# Organic subset (can appear without brackets)
ORGANIC_SUBSET = {'B', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I'}
AROMATIC_ORGANIC = {'b', 'c', 'n', 'o', 'p', 's', 'se', 'te'}

# Two-letter element symbols (for bracket atoms)
TWO_LETTER_ELEMENTS = {
    'He', 'Li', 'Be', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'Cl', 'Ar',
    'Ca', 'Sc', 'Ti', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Zr', 'Nb',
    'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb',
    'Te', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm',
    'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf',
    'Ta', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi',
    'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'Np', 'Pu',
    'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf',
    'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl',
    'Mc', 'Lv', 'Ts', 'Og',
}

BOND_CHARS = {'-', '=', '#', ':', '/', '\\'}


class SMILESTokenizerError(Exception):
    """Error during SMILES tokenization."""
    def __init__(self, message, position):
        super().__init__(f"Position {position}: {message}")
        self.position = position


def tokenize(smiles):
    """
    Tokenize a SMILES string into a list of Token objects.

    Args:
        smiles: SMILES string to tokenize

    Returns:
        List of Token objects

    Raises:
        SMILESTokenizerError on invalid input
    """
    tokens = []
    i = 0
    n = len(smiles)

    while i < n:
        ch = smiles[i]

        # ── Branch tokens ──
        if ch == '(':
            tokens.append(Token(TokenType.BRANCH_OPEN, '(', pos=i))
            i += 1

        elif ch == ')':
            tokens.append(Token(TokenType.BRANCH_CLOSE, ')', pos=i))
            i += 1

        # ── Dot (disconnect) ──
        elif ch == '.':
            tokens.append(Token(TokenType.DOT, '.', pos=i))
            i += 1

        # ── Bond tokens ──
        elif ch in BOND_CHARS:
            tokens.append(Token(TokenType.BOND, ch, pos=i))
            i += 1

        # ── Ring closure ──
        elif ch == '%':
            # Two-digit ring number
            if i + 2 < n and smiles[i + 1].isdigit() and smiles[i + 2].isdigit():
                ring_num = int(smiles[i + 1:i + 3])
                tokens.append(Token(TokenType.RING_CLOSURE, smiles[i:i + 3],
                                    data={'ring_num': ring_num}, pos=i))
                i += 3
            else:
                raise SMILESTokenizerError(
                    "Expected two digits after '%' for ring closure", i)

        elif ch.isdigit():
            ring_num = int(ch)
            tokens.append(Token(TokenType.RING_CLOSURE, ch,
                                data={'ring_num': ring_num}, pos=i))
            i += 1

        # ── Bracket atom [...]  ──
        elif ch == '[':
            bracket_end = smiles.find(']', i)
            if bracket_end == -1:
                raise SMILESTokenizerError("Unmatched '[' bracket", i)
            bracket_content = smiles[i + 1:bracket_end]
            data = _parse_bracket_atom(bracket_content, i)
            tokens.append(Token(TokenType.BRACKET_ATOM, smiles[i:bracket_end + 1],
                                data=data, pos=i))
            i = bracket_end + 1

        # ── Aromatic atoms (lowercase) ──
        elif ch in ('b', 'c', 'n', 'o', 'p', 's'):
            # Check for two-letter aromatic: se, te
            if ch == 's' and i + 1 < n and smiles[i + 1] == 'e':
                tokens.append(Token(TokenType.ATOM, 'se',
                                    data={'symbol': 'se', 'aromatic': True}, pos=i))
                i += 2
            elif ch == 't' and i + 1 < n and smiles[i + 1] == 'e':
                tokens.append(Token(TokenType.ATOM, 'te',
                                    data={'symbol': 'te', 'aromatic': True}, pos=i))
                i += 2
            else:
                tokens.append(Token(TokenType.ATOM, ch,
                                    data={'symbol': ch, 'aromatic': True}, pos=i))
                i += 1

        # ── Organic subset atoms (uppercase) ──
        elif ch.isupper():
            # Try two-letter organic subset first (Cl, Br)
            if i + 1 < n:
                two = smiles[i:i + 2]
                if two in ORGANIC_SUBSET:
                    tokens.append(Token(TokenType.ATOM, two,
                                        data={'symbol': two, 'aromatic': False}, pos=i))
                    i += 2
                    continue

            # Single letter organic subset
            if ch in ('B', 'C', 'N', 'O', 'P', 'S', 'F', 'I'):
                tokens.append(Token(TokenType.ATOM, ch,
                                    data={'symbol': ch, 'aromatic': False}, pos=i))
                i += 1
            else:
                raise SMILESTokenizerError(
                    f"Unexpected uppercase character '{ch}' — "
                    f"use brackets for non-organic elements", i)

        # ── Whitespace (skip) ──
        elif ch in (' ', '\t', '\n', '\r'):
            i += 1

        else:
            raise SMILESTokenizerError(f"Unexpected character '{ch}'", i)

    return tokens


def _parse_bracket_atom(content, bracket_pos):
    """
    Parse the content inside brackets: [isotope?symbol chirality? hcount? charge? :class?]

    Returns a dict with keys: isotope, symbol, aromatic, chirality, hcount, charge, atom_class
    """
    data = {
        'isotope': 0,
        'symbol': '',
        'aromatic': False,
        'chirality': 0,       # 0=none, 1=@, 2=@@
        'hcount': None,       # None = not specified, 0 = explicit [C] with no H
        'charge': 0,
        'atom_class': 0,
    }

    i = 0
    n = len(content)

    # ── Isotope (leading digits) ──
    isotope_str = ''
    while i < n and content[i].isdigit():
        isotope_str += content[i]
        i += 1
    if isotope_str:
        data['isotope'] = int(isotope_str)

    # ── Element symbol ──
    if i >= n:
        raise SMILESTokenizerError("Empty bracket atom", bracket_pos)

    # Check for aromatic in brackets
    if content[i] in ('b', 'c', 'n', 'o', 'p', 's'):
        # Two-letter aromatic check
        if content[i] == 's' and i + 1 < n and content[i + 1] == 'e':
            data['symbol'] = 'se'
            data['aromatic'] = True
            i += 2
        elif content[i] == 't' and i + 1 < n and content[i + 1] == 'e':
            data['symbol'] = 'te'
            data['aromatic'] = True
            i += 2
        else:
            data['symbol'] = content[i]
            data['aromatic'] = True
            i += 1
    elif content[i] == '*':
        # Wildcard atom
        data['symbol'] = '*'
        i += 1
    elif content[i].isupper():
        # Standard element — try two letters first
        if i + 1 < n and content[i + 1].islower() and content[i + 1] not in ('h',):
            two = content[i:i + 2]
            # Check three-letter lookahead to avoid misparse
            if two in TWO_LETTER_ELEMENTS or two in ORGANIC_SUBSET:
                data['symbol'] = two
                i += 2
            else:
                data['symbol'] = content[i]
                i += 1
        elif i + 1 < n and content[i + 1].islower():
            # Could be like 'Ch' — not a valid element, try single
            two = content[i:i + 2]
            if two in TWO_LETTER_ELEMENTS or two in ORGANIC_SUBSET:
                data['symbol'] = two
                i += 2
            else:
                data['symbol'] = content[i]
                i += 1
        else:
            data['symbol'] = content[i]
            i += 1
    else:
        raise SMILESTokenizerError(
            f"Expected element symbol in bracket, got '{content[i]}'", bracket_pos)

    # ── Chirality (@, @@) ──
    if i < n and content[i] == '@':
        i += 1
        if i < n and content[i] == '@':
            data['chirality'] = 2  # @@
            i += 1
        else:
            data['chirality'] = 1  # @

    # ── Hydrogen count ──
    if i < n and content[i] == 'H':
        i += 1
        h_str = ''
        while i < n and content[i].isdigit():
            h_str += content[i]
            i += 1
        data['hcount'] = int(h_str) if h_str else 1

    # ── Charge ──
    if i < n and content[i] in ('+', '-'):
        sign = 1 if content[i] == '+' else -1
        i += 1
        charge_str = ''
        while i < n and content[i].isdigit():
            charge_str += content[i]
            i += 1
        if charge_str:
            data['charge'] = sign * int(charge_str)
        else:
            # Count repeated signs: ++ = +2, --- = -3
            count = 1
            while i < n and content[i] == ('+' if sign == 1 else '-'):
                count += 1
                i += 1
            data['charge'] = sign * count

    # ── Atom class (:n) ──
    if i < n and content[i] == ':':
        i += 1
        cls_str = ''
        while i < n and content[i].isdigit():
            cls_str += content[i]
            i += 1
        if cls_str:
            data['atom_class'] = int(cls_str)

    return data


# =============================================================================
# BKChem-Inspired Regex Tokenization
# =============================================================================

def tokenize_regex(smiles: str) -> list[Token]:
    """
    Tokenize SMILES using BKChem-style regex-based chunking.
    
    This approach uses regex patterns similar to BKChem's implementation
    for more robust parsing of complex SMILES strings.
    
    Args:
        smiles: Input SMILES string
        
    Returns:
        List of Token objects
    """
    tokens = []
    
    # BKChem-style regex pattern for tokenization
    # Pattern breakdown:
    #   (\[.*?\])       - Bracket atoms like [C@@H], [NH4+]
    #   [A-Z][a-z]?     - Two-letter elements (Cl, Br, etc.)
    #   %[0-9]{1,2}     - Two-digit ring closures like %12
    #   [^A-Z]          - Single non-uppercase characters
    #   [a-z]           - Aromatic atoms (c, n, o, etc.)
    pattern = r'(\[.*?\]|[A-Z][a-z]?|%[0-9]{1,2}|[^A-Z]|[a-z])'
    
    # Split and filter out empty strings
    chunks = [c for c in re.split(pattern, smiles) if c]
    
    i = 0
    while i < len(chunks):
        chunk = chunks[i]
        
        # Skip whitespace
        if chunk.strip() == '':
            i += 1
            continue
        
        # Bracket atom
        if chunk.startswith('['):
            try:
                atom_data = _parse_bracket_atom(chunk[1:-1], i)
                tokens.append(Token(TokenType.BRACKET_ATOM, chunk, atom_data))
            except SMILESTokenizerError as e:
                raise SMILESTokenizerError(f"Invalid bracket atom: {chunk}", i) from e
        
        # Aromatic atoms (lowercase, single char)
        elif len(chunk) == 1 and chunk in 'bcnops':
            tokens.append(Token(TokenType.ATOM, chunk, {
                'symbol': chunk.upper(),
                'aromatic': True
            }))
        
        # Organic subset atoms (uppercase)
        elif chunk in ORGANIC_SUBSET or chunk in TWO_LETTER_ELEMENTS:
            tokens.append(Token(TokenType.ATOM, chunk, {
                'symbol': chunk,
                'aromatic': False
            }))
        
        # Bonds
        elif chunk in '-=#:':
            tokens.append(Token(TokenType.BOND, chunk, {'order': BOND_ORDER_MAP[chunk]}))
        
        # Stereochemistry bonds
        elif chunk in '/\\':
            tokens.append(Token(TokenType.BOND, chunk, {'order': 1, 'stereo': chunk}))
        
        # Ring closures (single digit)
        elif chunk.isdigit() and 0 <= int(chunk) <= 9:
            tokens.append(Token(TokenType.RING_CLOSURE, chunk, {'number': int(chunk)}))
        
        # Ring closures (two-digit %nn)
        elif chunk.startswith('%'):
            ring_num = int(chunk[1:])
            tokens.append(Token(TokenType.RING_CLOSURE, chunk, {'number': ring_num}))
        
        # Branch open
        elif chunk == '(':
            tokens.append(Token(TokenType.BRANCH_OPEN, chunk, None))
        
        # Branch close
        elif chunk == ')':
            tokens.append(Token(TokenType.BRANCH_CLOSE, chunk, None))
        
        # Dot disconnection
        elif chunk == '.':
            tokens.append(Token(TokenType.DOT, chunk, None))
        
        else:
            raise SMILESTokenizerError(f"Unrecognized token: '{chunk}'", i)
        
        i += 1
    
    return tokens
