"""Python ctypes binding for the KSgeneral native C ABI."""

from .ksgeneral_ctypes import KSGeneral, KSGeneralError, TwoSampleKSTestResult, load

__all__ = ["KSGeneral", "KSGeneralError", "TwoSampleKSTestResult", "load"]
