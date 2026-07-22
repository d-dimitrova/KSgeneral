from __future__ import annotations

import ctypes
from ctypes import POINTER, byref, c_char, c_double, c_int32, c_size_t, c_uint32
from pathlib import Path
from typing import Sequence

KSG_OK = 0

KSG_ALT_TWO_SIDED = 1
KSG_ALT_GREATER = 2
KSG_ALT_LESS = 3

KSG_PROB_LT = 1
KSG_PROB_GE = 2


class KSGeneralError(RuntimeError):
    """Raised when the native KSgeneral C API returns an error status."""


class KSGeneral:
    def __init__(self, library_path: str | Path) -> None:
        self._lib = ctypes.CDLL(str(library_path))
        self._configure_prototypes()

    def _configure_prototypes(self) -> None:
        self._lib.ksg_abi_version.argtypes = []
        self._lib.ksg_abi_version.restype = c_uint32
        self._lib.ksg_version_string.argtypes = []
        self._lib.ksg_version_string.restype = ctypes.c_char_p
        self._lib.ksg_last_error.argtypes = [POINTER(c_char), c_size_t]
        self._lib.ksg_last_error.restype = c_size_t
        self._lib.ksg_ks2_probability_summary.argtypes = [
            c_int32,
            c_int32,
            c_int32,
            POINTER(c_int32),
            c_size_t,
            c_double,
            POINTER(c_double),
            c_size_t,
            c_double,
            c_int32,
            POINTER(c_double),
        ]
        self._lib.ksg_ks2_probability_summary.restype = c_int32
        self._lib.ksg_kuiper2_probability_summary.argtypes = [
            c_int32,
            c_int32,
            POINTER(c_int32),
            c_size_t,
            c_double,
            c_int32,
            POINTER(c_double),
        ]
        self._lib.ksg_kuiper2_probability_summary.restype = c_int32
        self._lib.ksg_one_sample_boundary_probability.argtypes = [
            c_int32,
            POINTER(c_double),
            c_size_t,
            POINTER(c_double),
            c_size_t,
            c_int32,
            POINTER(c_double),
        ]
        self._lib.ksg_one_sample_boundary_probability.restype = c_int32
        self._lib.ksg_continuous_ks_ccdf.argtypes = [c_int32, c_double, POINTER(c_double)]
        self._lib.ksg_continuous_ks_ccdf.restype = c_int32

    def _last_error(self) -> str:
        required = self._lib.ksg_last_error(None, 0)
        if required == 0:
            return "native calculation failed"
        buffer = ctypes.create_string_buffer(required + 1)
        self._lib.ksg_last_error(buffer, len(buffer))
        return buffer.value.decode("utf-8", errors="replace")

    def _check(self, status: int, result: c_double) -> float:
        if status != KSG_OK:
            raise KSGeneralError(self._last_error())
        return result.value

    def ks2_probability_summary(
        self,
        m: int,
        n: int,
        alternative: int,
        multiplicities: Sequence[int],
        statistic: float,
        weights: Sequence[float],
        *,
        tolerance: float = 1e-8,
        probability_type: int = KSG_PROB_GE,
    ) -> float:
        m_array = (c_int32 * len(multiplicities))(*multiplicities)
        w_array = (c_double * len(weights))(*weights)
        result = c_double()
        status = self._lib.ksg_ks2_probability_summary(
            c_int32(m), c_int32(n), c_int32(alternative), m_array,
            c_size_t(len(multiplicities)), c_double(statistic), w_array,
            c_size_t(len(weights)), c_double(tolerance), c_int32(probability_type),
            byref(result),
        )
        return self._check(status, result)

    def kuiper2_probability_summary(
        self,
        m: int,
        n: int,
        multiplicities: Sequence[int],
        statistic: float,
        *,
        probability_type: int = KSG_PROB_GE,
    ) -> float:
        m_array = (c_int32 * len(multiplicities))(*multiplicities)
        result = c_double()
        status = self._lib.ksg_kuiper2_probability_summary(
            c_int32(m), c_int32(n), m_array, c_size_t(len(multiplicities)),
            c_double(statistic), c_int32(probability_type), byref(result),
        )
        return self._check(status, result)

    def one_sample_boundary_probability(
        self,
        n: int,
        g_steps: Sequence[float],
        h_steps: Sequence[float],
        *,
        use_fft: bool = True,
    ) -> float:
        g_array = (c_double * len(g_steps))(*g_steps)
        h_array = (c_double * len(h_steps))(*h_steps)
        result = c_double()
        status = self._lib.ksg_one_sample_boundary_probability(
            c_int32(n), g_array, c_size_t(len(g_steps)), h_array,
            c_size_t(len(h_steps)), c_int32(1 if use_fft else 0), byref(result),
        )
        return self._check(status, result)

    def continuous_ks_ccdf(self, n: int, q: float) -> float:
        result = c_double()
        status = self._lib.ksg_continuous_ks_ccdf(c_int32(n), c_double(q), byref(result))
        return self._check(status, result)
