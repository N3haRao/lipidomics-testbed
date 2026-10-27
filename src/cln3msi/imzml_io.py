"""Reading imzML/ibd files without ever pulling a whole slide into memory.

Two jobs:

1. read_header(): walk the XML with lxml.iterparse and collect every cvParam in the header
   (fileDescription, referenceable param groups, instrument, scan settings, software, data
   processing) plus the cvParams of the first few spectra. It stops before the big spectrum
   list, so it's quick. We keep the names exactly as written in the file, so whatever the
   exporter (SCiLS, flexImaging, ...) wrote ends up in the log.

2. SlideReader: wraps pyimzml's ImzMLParser (which parses coordinates and byte offsets, but
   reads no intensities) and gives lazy iterators over spectra. The .ibd is only ever opened
   in read-only mode. For continuous files with uncompressed, evenly spaced intensity
   arrays we use a read-only numpy memmap, which is much faster than seeking spectrum by
   spectrum, and the OS pages data in and out on its own, so RAM stays flat.
"""
from __future__ import annotations

import warnings
from pathlib import Path
from typing import Iterator

import numpy as np
from lxml import etree

NS = "{http://psi.hupo.org/ms/mzml}"

# cvParam names that tell us what kind of data we have. Matching on names AND accessions
# because exporters are not always consistent.
SPECTRUM_TYPE = {"MS:1000127": "centroid", "MS:1000128": "profile",
                 "centroid spectrum": "centroid", "profile spectrum": "profile"}
STORAGE = {"IMS:1000030": "continuous", "IMS:1000031": "processed",
           "continuous": "continuous", "processed": "processed"}
POLARITY = {"MS:1000129": "negative", "MS:1000130": "positive",
            "negative scan": "negative", "positive scan": "positive"}
COMPRESSION = {"MS:1000574": "zlib", "MS:1000576": "none", "zlib compression": "zlib", "no compression": "none"}


def _local(tag: str) -> str:
    return tag.split("}", 1)[-1] if "}" in tag else tag


def read_header(imzml_path: Path, n_spectra_cv: int = 3) -> dict:
    """Collect header cvParams/userParams without touching spectra intensities.

    Returns {'sections': {section_name: [ {accession, name, value, unit}, ...]},
             'spectrum_count': int or None, 'first_spectra': [[params], ...]}
    """
    sections: dict[str, list] = {}
    first_spectra: list[list] = []
    spectrum_count = None
    stack: list[str] = []
    current_spec: list | None = None
    seen_spectra = 0

    ctx = etree.iterparse(str(imzml_path), events=("start", "end"), huge_tree=True)
    for event, el in ctx:
        tag = _local(el.tag)
        if event == "start":
            stack.append(tag + (f"[{el.get('id')}]" if el.get("id") else ""))
            if tag == "spectrumList":
                c = el.get("count")
                spectrum_count = int(c) if c and c.isdigit() else None
            if tag == "spectrum":
                current_spec = []
            continue
        # event == "end"
        if tag in ("cvParam", "userParam"):
            rec = {"accession": el.get("accession", ""), "name": el.get("name", ""),
                   "value": el.get("value", ""), "unit": el.get("unitName", ""), "kind": tag}
            if current_spec is not None:
                current_spec.append(rec | {"where": "/".join(stack[-3:-1])})
            else:
                # Section = first two levels under <mzML>, e.g. fileDescription/fileContent
                key = "/".join(s for s in stack[1:-1])
                sections.setdefault(key, []).append(rec)
        if tag == "referenceableParamGroupRef" and current_spec is not None:
            current_spec.append({"accession": "", "name": "referenceableParamGroupRef",
                                 "value": el.get("ref", ""), "unit": "", "kind": "ref"})
        if tag == "spectrum":
            first_spectra.append(current_spec or [])
            current_spec = None
            seen_spectra += 1
            el.clear()
            if seen_spectra >= n_spectra_cv:
                break
        stack.pop()
    del ctx
    return {"sections": sections, "spectrum_count": spectrum_count, "first_spectra": first_spectra}


def summarize_header(hdr: dict) -> dict:
    """Pull the handful of facts we branch on out of the raw header dump."""
    every = [p for ps in hdr["sections"].values() for p in ps]
    every_spec = [p for ps in hdr["first_spectra"] for p in ps]
    # referenceable groups used by the first spectrum count as 'its' params
    groups = {k.split("[", 1)[1].rstrip("]"): v for k, v in hdr["sections"].items()
              if "referenceableParamGroup[" in k}
    refs = [p["value"] for p in every_spec if p["kind"] == "ref"]
    spec_params = every_spec + [p for r in refs for p in groups.get(r, [])]

    def find(table: dict, params: list) -> str | None:
        for p in params:
            for key in (p["accession"], p["name"]):
                if key in table:
                    return table[key]
        return None

    fc = hdr["sections"].get("fileDescription/fileContent", [])
    out = {
        "spectrum_type_header": find(SPECTRUM_TYPE, fc) or find(SPECTRUM_TYPE, spec_params),
        "storage_mode": find(STORAGE, fc) or find(STORAGE, every),
        "polarity_header": find(POLARITY, spec_params) or find(POLARITY, every),
        "compression": find(COMPRESSION, every),
        "spectrum_count_header": hdr["spectrum_count"],
    }

    def value_of(names: tuple[str, ...]):
        for p in every:
            if p["name"] in names or p["accession"] in names:
                return p["value"] or p["name"]
        return None

    out["pixel_size_x_um"] = value_of(("pixel size x", "IMS:1000046", "pixel size (x)"))
    out["pixel_size_y_um"] = value_of(("pixel size y", "IMS:1000047"))
    out["max_count_x"] = value_of(("max count of pixels x", "IMS:1000042"))
    out["max_count_y"] = value_of(("max count of pixels y", "IMS:1000043"))
    out["scan_window_lower"] = None
    out["scan_window_upper"] = None
    for p in spec_params:
        if p["name"] in ("scan window lower limit", "lowest observed m/z") and p["value"]:
            out["scan_window_lower"] = float(p["value"])
        if p["name"] in ("scan window upper limit", "highest observed m/z") and p["value"]:
            out["scan_window_upper"] = float(p["value"])
    instr = [p for k, ps in hdr["sections"].items() if k.startswith("instrumentConfigurationList")
             for p in ps]
    out["instrument_params"] = "; ".join(
        f"{p['name']}={p['value']}" if p["value"] else p["name"] for p in instr) or None
    soft = [p for k, ps in hdr["sections"].items() if k.startswith("softwareList") for p in ps]
    out["software"] = "; ".join(f"{p['name']} {p['value']}".strip() for p in soft) or None
    proc = [p for k, ps in hdr["sections"].items() if k.startswith("dataProcessingList") for p in ps]
    out["data_processing"] = "; ".join(p["name"] for p in proc) or None
    return out


def ibd_path_for(imzml: Path) -> Path | None:
    for ext in (".ibd", ".IBD", ".Ibd"):
        p = imzml.with_suffix(ext)
        if p.exists():
            return p
    return None


class SlideReader:
    """Lazy access to one imzML/ibd pair. Nothing big is loaded at construction time."""

    def __init__(self, imzml_path: Path, log=None):
        from pyimzml.ImzMLParser import ImzMLParser
        self.path = Path(imzml_path)
        self.ibd = ibd_path_for(self.path)
        if self.ibd is None:
            raise FileNotFoundError(f"no .ibd next to {self.path}")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            # ImzMLParser opens the .ibd with 'rb' (read-only). It parses XML offsets, not data.
            self.p = ImzMLParser(str(self.path), parse_lib="lxml")
        self.coords = np.asarray([c[:2] for c in self.p.coordinates], dtype=np.int32)   # (x, y), 1-based
        self.n = len(self.coords)
        codes = {"f": np.float32, "d": np.float64, "i": np.int32, "l": np.int64}
        self.mz_dtype = np.dtype(codes[self.p.mzPrecision])
        self.int_dtype = np.dtype(codes[self.p.intensityPrecision])
        self.mz_offsets = np.asarray(self.p.mzOffsets, dtype=np.int64)
        self.int_offsets = np.asarray(self.p.intensityOffsets, dtype=np.int64)
        self.mz_lengths = np.asarray(self.p.mzLengths, dtype=np.int64)
        self.int_lengths = np.asarray(self.p.intensityLengths, dtype=np.int64)
        meta = self.p.metadata.file_description if hasattr(self.p, "metadata") else {}
        self.compressed = False
        try:
            groups = self.p.metadata.referenceable_param_groups
            self.compressed = any(g.get("zlib compression") for g in groups.values())
        except Exception:
            pass
        # 'Continuous' = every spectrum shares ONE m/z array (same offset in the .ibd)
        self.shares_mz_axis = bool(np.all(self.mz_offsets == self.mz_offsets[0]))
        self.declared_continuous = bool(meta.get("continuous")) if isinstance(meta, dict) else None
        self._memmap = None
        if self.shares_mz_axis and not self.compressed and np.all(self.int_lengths == self.int_lengths[0]):
            L = int(self.int_lengths[0])
            stride = L * self.int_dtype.itemsize
            order = np.argsort(self.int_offsets)
            if np.all(np.diff(self.int_offsets[order]) == stride):
                self._memmap = np.memmap(self.ibd, dtype=self.int_dtype, mode="r",
                                         offset=int(self.int_offsets[order[0]]),
                                         shape=(self.n, L))
                # row r of the memmap belongs to spectrum order[r]; we want the inverse
                self._row_of = np.empty(self.n, dtype=np.int64)
                self._row_of[order] = np.arange(self.n)
        if log:
            log.info(f"{self.path.name}: {self.n} spectra | shared m/z axis={self.shares_mz_axis} | "
                     f"memmap fast path={self._memmap is not None} | mz {self.mz_dtype} int {self.int_dtype}")

    # ---- single reads -------------------------------------------------------------------
    def mz_axis(self) -> np.ndarray:
        """The shared m/z axis of a continuous file (read once, it's small)."""
        if not self.shares_mz_axis:
            raise ValueError("this file does not share one m/z axis (processed mode)")
        return self.read_mz(0)

    def read_mz(self, i: int) -> np.ndarray:
        if self.compressed:
            return np.asarray(self.p.getspectrum(i)[0])
        with open(self.ibd, "rb") as fh:
            fh.seek(int(self.mz_offsets[i]))
            return np.fromfile(fh, dtype=self.mz_dtype, count=int(self.mz_lengths[i])).astype(np.float64)

    def spectrum(self, i: int) -> tuple[np.ndarray, np.ndarray]:
        mz, it = self.p.getspectrum(i)
        return np.asarray(mz, dtype=np.float64), np.asarray(it, dtype=np.float32)

    # ---- streaming -----------------------------------------------------------------------
    def iter_spectra(self, indices=None, log=None, every: int = 10000) -> Iterator[tuple[int, np.ndarray, np.ndarray]]:
        """Yield (index, mz, intensity) one spectrum at a time."""
        idx = np.arange(self.n) if indices is None else np.asarray(indices)
        if self.shares_mz_axis:
            mz = self.mz_axis()
        for k, i in enumerate(idx):
            if self._memmap is not None:
                it = np.asarray(self._memmap[self._row_of[i]], dtype=np.float32)
            elif self.shares_mz_axis:
                it = self.spectrum(int(i))[1]
            else:
                mz, it = self.spectrum(int(i))
            if log is not None and k and k % every == 0:
                log.mem(f"{self.path.name}: {k}/{len(idx)} spectra")
            yield int(i), mz, it

    def iter_blocks(self, indices=None, chunk: int = 2000):
        """Continuous files only: yield (indices, block) with block = (n, L) float32 intensities."""
        if not self.shares_mz_axis:
            raise ValueError("iter_blocks needs a shared m/z axis")
        idx = np.arange(self.n) if indices is None else np.asarray(indices)
        # Cap a block at ~400 MB of float64 work space, whatever the spectrum length is
        L = int(self.int_lengths[0])
        chunk = max(16, min(int(chunk), int(4e8 // (8 * (L + 1)))))
        for s in range(0, len(idx), chunk):
            sub = idx[s:s + chunk]
            if self._memmap is not None:
                rows = self._row_of[sub]
                o = np.argsort(rows)   # read in file order, it's friendlier to the disk
                block = np.empty((len(sub), self._memmap.shape[1]), dtype=np.float32)
                block[o] = self._memmap[rows[o]]
            else:
                block = np.stack([self.spectrum(int(i))[1] for i in sub]).astype(np.float32)
            yield sub, block

    def image_shape(self) -> tuple[int, int]:
        """(rows, cols) of the pixel grid, using the max x/y coordinate (1-based)."""
        return int(self.coords[:, 1].max()), int(self.coords[:, 0].max())

    def to_image(self, values: np.ndarray, fill=np.nan) -> np.ndarray:
        """Per-spectrum values -> 2D image (rows = y, cols = x). Pixels not acquired get `fill`."""
        img = np.full(self.image_shape(), fill, dtype=np.float64)
        img[self.coords[:, 1] - 1, self.coords[:, 0] - 1] = values
        return img

    def close(self):
        try:
            self.p.m.close()
        except Exception:
            pass
        self._memmap = None


def detect_profile(mz: np.ndarray, it: np.ndarray) -> dict:
    """Look at a spectrum and guess profile vs centroid from the data itself.

    In a profile spectrum every peak is drawn with several points, so the strongest points
    have "shoulders": both raw neighbours are also non-zero, a decent fraction of the
    height, and sit at nearly the same spacing on either side (a smooth sampling grid).
    In a centroid spectrum the neighbours of a strong point are just other peaks, at
    irregular distances. We take the top 1% of points and count how many have shoulders.
    """
    n = len(it)
    if np.count_nonzero(it) < 20:
        return {"guess": "unknown", "shoulder_fraction": np.nan, "median_spacing_ppm": np.nan}
    spacing = np.diff(mz)
    med_sp = float(np.median(spacing[spacing > 0]))
    k = max(5, int(0.01 * np.count_nonzero(it)))
    top = np.argsort(it)[::-1][:k]
    top = top[(top > 0) & (top < n - 1)]
    left, right = it[top - 1], it[top + 1]
    dl, dr = mz[top] - mz[top - 1], mz[top + 1] - mz[top]
    ratio = np.where(dr > 0, dl / np.maximum(dr, 1e-12), 0)
    ok = ((left > 0.05 * it[top]) & (right > 0.05 * it[top]) & (ratio > 0.8) & (ratio < 1.25)
          & (dl < 5 * med_sp) & (dr < 5 * med_sp))
    frac = float(ok.mean()) if len(ok) else np.nan
    spacing_ppm = float(np.median(spacing / mz[:-1] * 1e6))
    guess = "profile" if frac > 0.5 else "centroid"
    return {"guess": guess, "shoulder_fraction": frac, "median_spacing_ppm": spacing_ppm}
