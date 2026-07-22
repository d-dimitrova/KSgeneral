from __future__ import annotations

import ctypes
import math
import sys
from ctypes import POINTER, byref, c_char, c_double, c_int32, c_size_t, c_uint32
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Sequence

KSG_OK = 0

KSG_ALT_TWO_SIDED = 1
KSG_ALT_GREATER = 2
KSG_ALT_LESS = 3

KSG_PROB_LT = 1
KSG_PROB_GE = 2

Alternative = str | int
Weight = float | Callable[[float], float]


class KSGeneralError(RuntimeError):
    """Raised when the native KSgeneral C API returns an error status."""


@dataclass(frozen=True)
class TwoSampleKSTestResult:
    """Result returned by :meth:`KSGeneral.two_sample_ks_test`."""

    statistic: float
    p_value: float
    alternative: str
    n_x: int
    n_y: int
    multiplicities: tuple[int, ...]
    weights: tuple[float, ...]
    conservative: bool


def default_library_name() -> str:
    """Return the platform default KSgeneral shared-library filename."""

    if sys.platform == "win32":
        return "ksgeneral.dll"
    if sys.platform == "darwin":
        return "libksgeneral.dylib"
    return "libksgeneral.so"


def load(library_path: str | Path | None = None) -> "KSGeneral":
    """Load ``libksgeneral``/``ksgeneral.dll`` and return a ctypes wrapper.

    ``library_path`` may be an absolute path, a relative path, or ``None``.  When
    it is ``None``, the platform-specific library name is passed to
    :class:`ctypes.CDLL`, so the operating system's normal shared-library search
    rules apply.
    """

    return KSGeneral(default_library_name() if library_path is None else library_path)


def _clean_sample(values: Sequence[float], name: str) -> list[float]:
    sample = [float(value) for value in values if not math.isnan(float(value))]
    if not sample:
        raise ValueError(f"not enough {name!r} data")
    return sample


def _alternative_code(alternative: Alternative) -> int:
    if isinstance(alternative, int):
        if alternative in (KSG_ALT_TWO_SIDED, KSG_ALT_GREATER, KSG_ALT_LESS):
            return alternative
        raise ValueError(
            "alternative integer must be one of KSG_ALT_TWO_SIDED, "
            "KSG_ALT_GREATER, or KSG_ALT_LESS"
        )

    key = alternative.replace("_", ".").replace("-", ".").lower()
    if key in {"two.sided", "two"}:
        return KSG_ALT_TWO_SIDED
    if key == "greater":
        return KSG_ALT_GREATER
    if key == "less":
        return KSG_ALT_LESS
    raise ValueError("alternative must be 'two-sided', 'greater', or 'less'")


def _alternative_name(code: int) -> str:
    return {
        KSG_ALT_TWO_SIDED: "two-sided",
        KSG_ALT_GREATER: "greater",
        KSG_ALT_LESS: "less",
    }[code]


def _weight_function(weight: Weight) -> Callable[[float], float]:
    if callable(weight):
        return weight
    exponent = float(weight)
    if exponent == 0.0:
        return lambda _t: 1.0
    if 0.0 < exponent <= 1.0:
        return lambda t: (t * (1.0 - t)) ** (-exponent)
    raise ValueError("numeric weight must be 0 or in the interval (0, 1]")


def _pooled_multiplicities(
    x: Sequence[float], y: Sequence[float]
) -> tuple[list[float], list[int]]:
    pooled = sorted([*x, *y])
    unique: list[float] = []
    counts: list[int] = []
    for value in pooled:
        if unique and value == unique[-1]:
            counts[-1] += 1
        else:
            unique.append(value)
            counts.append(1)
    return unique, counts


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

    @property
    def abi_version(self) -> int:
        return int(self._lib.ksg_abi_version())

    @property
    def version(self) -> str:
        return self._lib.ksg_version_string().decode("utf-8", errors="replace")

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

    def two_sample_ks_test(
        self,
        x: Sequence[float],
        y: Sequence[float],
        *,
        alternative: Alternative = "two-sided",
        conservative: bool = False,
        weight: Weight = 0.0,
        tolerance: float = 1e-8,
    ) -> TwoSampleKSTestResult:
        """Run the exact two-sample KS test for discrete data with ties.

        This mirrors the preprocessing done by the R ``KS2sample`` function:
        NaNs are dropped, samples are pooled to compute tie multiplicities, the
        empirical-CDF difference is evaluated at pooled support points, and the
        native C ABI is called with ``P(D >= observed D)``.
        """

        x_values = _clean_sample(x, "x")
        y_values = _clean_sample(y, "y")
        nx = len(x_values)
        ny = len(y_values)
        total = nx + ny
        code = _alternative_code(alternative)
        weighted = _weight_function(weight)
        weights = tuple(float(weighted(i / total)) for i in range(1, total))
        if any((not math.isfinite(value) or value <= 0.0) for value in weights):
            raise ValueError(
                "weight function must be finite and strictly positive on "
                "i / (len(x) + len(y))"
            )

        support, multiplicities = _pooled_multiplicities(x_values, y_values)
        x_sorted = sorted(x_values)
        y_sorted = sorted(y_values)
        ix = iy = cumulative = 0
        z_values: list[float] = []
        for value, multiplicity in zip(support[:-1], multiplicities[:-1]):
            while ix < nx and x_sorted[ix] <= value:
                ix += 1
            while iy < ny and y_sorted[iy] <= value:
                iy += 1
            cumulative += multiplicity
            z_values.append((ix / nx - iy / ny) * weights[cumulative - 1])

        if code == KSG_ALT_TWO_SIDED:
            statistic = max(abs(value) for value in z_values) if z_values else 0.0
            native_code = code
            native_nx, native_ny = nx, ny
        elif code == KSG_ALT_GREATER:
            statistic = max(z_values) if z_values else 0.0
            native_code = KSG_ALT_GREATER
            native_nx, native_ny = nx, ny
            if nx != min(nx, ny):
                native_code = KSG_ALT_LESS
                native_nx, native_ny = ny, nx
        else:
            statistic = max(-value for value in z_values) if z_values else 0.0
            native_code = KSG_ALT_LESS
            native_nx, native_ny = nx, ny
            if nx != min(nx, ny):
                native_code = KSG_ALT_GREATER
                native_nx, native_ny = ny, nx

        native_multiplicities = tuple(
            [1] * total
            if conservative and len(multiplicities) != total
            else multiplicities
        )
        p_value = self.ks2_probability_summary(
            native_nx,
            native_ny,
            native_code,
            native_multiplicities,
            statistic,
            weights,
            tolerance=tolerance,
            probability_type=KSG_PROB_GE,
        )
        return TwoSampleKSTestResult(
            statistic=statistic,
            p_value=p_value,
            alternative=_alternative_name(code),
            n_x=nx,
            n_y=ny,
            multiplicities=tuple(multiplicities),
            weights=weights,
            conservative=conservative,
        )

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
