"""
Service registry — central wiring point for all PyChem services.

No DI framework. A plain Python class with factory methods.
Created once at startup, passed to consumers via constructor.
"""
from src.core.parallel import ParallelExecutor
from src.core.events import EventBus
from src.services.forcefield.mmff94_service import MMFF94Service
from src.services.descriptors.descriptor_service import DescriptorService
from src.services.loading.loader_service import LoaderService
from src.services.alignment.align_service import AlignmentService


class ServiceRegistry:
    """
    Constructs and holds all service instances.

    Usage:
        registry = ServiceRegistry()
        mol = registry.loader.load("protein.pdb")
        registry.forcefield.optimize_geometry(mol)
    """

    def __init__(self):
        self.executor = ParallelExecutor()
        self.event_bus = EventBus()
        self.forcefield = MMFF94Service(self.executor)
        self.descriptors = DescriptorService(self.executor)
        self.loader = LoaderService(self.executor)
        self.alignment = AlignmentService(loader=self.loader)
        # Services added in later phases:
        # self.coord_gen = CoordinateGeneratorService(self.executor)
        # self.renderer_factory = RendererFactory()

    def shutdown(self):
        """Clean shutdown of all services."""
        self.executor.shutdown()
