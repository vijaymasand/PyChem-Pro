import pytest
from src.core.domain.models.molecule import MoleculeBuilder
from src.core.optimization.geometry_optimizer import GeometryOptimizer

@pytest.fixture
def water():
    return MoleculeBuilder.from_smiles("O")

@pytest.fixture
def methane():
    return MoleculeBuilder.from_smiles("C")

@pytest.fixture
def ethanol():
    return MoleculeBuilder.from_smiles("CCO")

def test_mmff94_optimization_water(water):
    opt = GeometryOptimizer(method="MMFF94", tolerance=1e-4, max_iter=200)
    result = opt.optimize(water)
    assert result.energy < water.energy
    assert result.rms_gradient < 1e-3
    oh_bonds = [b for b in result.get_bonds() if b.atom1.element == "O" and b.atom2.element == "H"]
    for bond in oh_bonds:
        assert 0.95 < bond.length < 1.05

def test_mmff94_optimization_methane(methane):
    opt = GeometryOptimizer(method="MMFF94", tolerance=1e-4, max_iter=200)
    result = opt.optimize(methane)
    assert result.energy < methane.energy
    assert result.rms_gradient < 1e-3
    ch_bonds = [b for b in result.get_bonds() if b.atom1.element == "C" and b.atom2.element == "H"]
    for bond in ch_bonds:
        assert 1.08 < bond.length < 1.10

def test_am1_optimization_ethanol(ethanol):
    opt = GeometryOptimizer(method="AM1", tolerance=1e-4, max_iter=200)
    result = opt.optimize(ethanol)
    assert result.energy < ethanol.energy
    assert result.rms_gradient < 1e-3
    cc_bond = [b for b in result.get_bonds() if b.atom1.element == "C" and b.atom2.element == "C"][0]
    assert 1.50 < cc_bond.length < 1.55

def test_pm3_optimization_water(water):
    opt = GeometryOptimizer(method="PM3", tolerance=1e-4, max_iter=200)
    result = opt.optimize(water)
    assert result.energy < water.energy
    assert result.rms_gradient < 1e-3

def test_convergence_respected(water):
    opt = GeometryOptimizer(method="MMFF94", tolerance=1e-6, max_iter=500)
    result = opt.optimize(water)
    assert result.converged
    assert result.rms_gradient < 1e-6

def test_invalid_molecule():
    with pytest.raises(ValueError):
        opt = GeometryOptimizer(method="MMFF94")
        opt.optimize(None)

def test_missing_parameters():
    with pytest.raises(RuntimeError):
        opt = GeometryOptimizer(method="UNKNOWN")
        opt.optimize(water())
