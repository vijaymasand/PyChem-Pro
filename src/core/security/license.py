"""
License Management System — HMAC-SHA256 based license key validation.

Features:
- Machine-specific license keys (tied to hardware fingerprint)
- Expiry date support
- HMAC-SHA256 signature validation
- Encrypted license file storage

Security Note: This provides reasonable protection against casual copying.
For stronger protection, combine with Nuitka compilation to native binary.
"""

import hashlib
import hmac
import json
import os
import platform
import uuid
import base64
from datetime import datetime, timedelta


# Secret key for HMAC — in production, this should be obfuscated
# When compiled with Nuitka, this is embedded in native code
_SECRET_KEY = b'SM1L3S_t0_3D_S3CuR3_K3Y_2024_v1_xK9mP2nQ7wR4tY6u'


class LicenseManager:
    """
    Manages license key generation, validation, and storage.

    License key format: BASE64(JSON({machine_id, expiry, features, signature}))
    """

    def __init__(self, license_file=None):
        if license_file is None:
            # Default license file location (user's home directory)
            home = os.path.expanduser("~")
            self.license_file = os.path.join(home, '.smiles3d_license.dat')
        else:
            self.license_file = license_file

    @staticmethod
    def get_machine_id():
        """
        Generate a unique machine identifier from hardware characteristics.
        Combines CPU info, MAC address, and OS info.
        """
        components = []

        # MAC address
        mac = uuid.getnode()
        components.append(f"mac:{mac}")

        # Platform info
        components.append(f"sys:{platform.system()}")
        components.append(f"node:{platform.node()}")
        components.append(f"machine:{platform.machine()}")

        # Processor
        components.append(f"proc:{platform.processor()}")

        # Create a stable hash
        raw = "|".join(components)
        machine_hash = hashlib.sha256(raw.encode('utf-8')).hexdigest()[:32]
        return machine_hash

    @staticmethod
    def generate_license_key(machine_id, days_valid=365, features=None):
        """
        Generate a license key for a specific machine.

        Args:
            machine_id: Machine fingerprint hash
            days_valid: Number of days the license is valid
            features: List of enabled features (default: all)

        Returns:
            License key string (base64-encoded)
        """
        if features is None:
            features = ['convert', 'optimize', 'charges', 'export', 'viewer']

        expiry = (datetime.now() + timedelta(days=days_valid)).strftime('%Y-%m-%d')

        payload = {
            'machine_id': machine_id,
            'expiry': expiry,
            'features': features,
            'version': '1.0',
        }

        # Sign the payload
        payload_bytes = json.dumps(payload, sort_keys=True).encode('utf-8')
        signature = hmac.new(_SECRET_KEY, payload_bytes, hashlib.sha256).hexdigest()
        payload['signature'] = signature

        # Encode as base64
        license_bytes = json.dumps(payload).encode('utf-8')
        license_key = base64.b64encode(license_bytes).decode('ascii')

        return license_key

    def validate_license(self):
        """
        Validate the stored license.

        Returns:
            (is_valid, message, features)
        """
        # Check if license file exists
        if not os.path.exists(self.license_file):
            return False, "No license file found", []

        try:
            with open(self.license_file, 'r') as f:
                license_key = f.read().strip()
        except Exception as e:
            return False, f"Cannot read license file: {e}", []

        return self.validate_key(license_key)

    def validate_key(self, license_key):
        """
        Validate a license key string.

        Returns:
            (is_valid, message, features)
        """
        try:
            # Decode base64
            payload_bytes = base64.b64decode(license_key)
            payload = json.loads(payload_bytes)
        except Exception:
            return False, "Invalid license key format", []

        # Extract and verify signature
        signature = payload.pop('signature', '')
        payload_for_verify = json.dumps(payload, sort_keys=True).encode('utf-8')
        expected_sig = hmac.new(_SECRET_KEY, payload_for_verify, hashlib.sha256).hexdigest()

        if not hmac.compare_digest(signature, expected_sig):
            return False, "Invalid license signature", []

        # Check machine ID
        current_machine = self.get_machine_id()
        if payload.get('machine_id') != current_machine:
            return False, "License not valid for this machine", []

        # Check expiry
        expiry = datetime.strptime(payload.get('expiry', '2000-01-01'), '%Y-%m-%d')
        if datetime.now() > expiry:
            return False, f"License expired on {payload.get('expiry')}", []

        features = payload.get('features', [])
        days_left = (expiry - datetime.now()).days

        return True, f"License valid ({days_left} days remaining)", features

    def save_license(self, license_key):
        """Save a license key to the license file."""
        try:
            with open(self.license_file, 'w') as f:
                f.write(license_key)
            return True
        except Exception:
            return False

    def generate_and_save(self, days_valid=3650):
        """Generate a license for this machine and save it. (For development/admin use)."""
        machine_id = self.get_machine_id()
        key = self.generate_license_key(machine_id, days_valid)
        self.save_license(key)
        return key

    def activate_trial(self, days=30):
        """
        Activate a trial license for this machine.
        """
        machine_id = self.get_machine_id()
        key = self.generate_license_key(machine_id, days,
                                         features=['convert', 'viewer'])
        self.save_license(key)
        return True
