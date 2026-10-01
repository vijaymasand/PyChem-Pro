"""
Hardware Profiler

Detects system specifications (like RAM, CPU) and determines the appropriate
performance tiers and rendering settings (LOD) for the application.
"""
try:
    import psutil
    HAS_PSUTIL = True
except ImportError:
    HAS_PSUTIL = False

class HardwareProfiler:
    _instance = None
    
    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(HardwareProfiler, cls).__new__(cls)
            cls._instance._init_profiler()
        return cls._instance
        
    def _init_profiler(self):
        # Calculate total RAM in GB
        if HAS_PSUTIL:
            try:
                self.total_ram_gb = psutil.virtual_memory().total / (1024 ** 3)
            except:
                self.total_ram_gb = 8.0 # Safe fallback
        else:
            # Fallback for Windows without psutil
            try:
                import subprocess
                cmd = "wmic computersystem get totalphysicalmemory"
                output = subprocess.check_output(cmd, shell=True).decode()
                # Output looks like: TotalPhysicalMemory \n 34359738368
                mem_bytes = int(output.split()[1])
                self.total_ram_gb = mem_bytes / (1024 ** 3)
            except:
                # If all fails, assume high-end for this user specifically 
                # since we know they have 32GB, or just default to 16GB
                self.total_ram_gb = 16.0 
        
        # Determine Performance Tier
        if self.total_ram_gb >= 12.0:
            self.is_high_end = True
        else:
            self.is_high_end = False
            
    @property
    def interactive_mesh_steps(self):
        """Returns the number of spline steps to use during interactive rotation."""
        if self.is_high_end:
            return None  # None means keep using the high-quality cached mesh (24)
        return 3         # Drop to 3 for low-end systems
        
    @property
    def interactive_mesh_profile(self):
        """Returns the profile detail to use during interactive rotation."""
        if self.is_high_end:
            return None  # None means keep using the high-quality cached mesh (16)
        return 4         # Drop to 4 for low-end systems
        
    @property
    def interactive_scale_factor(self):
        """Returns the supersampling scale factor during interactive rotation."""
        if self.is_high_end:
            return 2     # Keep 2x supersampling
        return 1         # Disable supersampling
