# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: Copyright 2026 Jumpei Ito & Junna Kawasaki

from __future__ import annotations

import logging
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

logger = logging.getLogger(__name__)

_STAR_MIN_VERSION = (2, 7, 10, "b")

_STAR_VERSION_PATTERN = re.compile(
    r"(\d+)\.(\d+)\.(\d+)([a-z])?",
)


@dataclass(frozen=True)
class StarVersion:
    major: int
    minor: int
    patch: int
    suffix: str

    def as_tuple(self) -> tuple[int, int, int, str]:
        return (self.major, self.minor, self.patch, self.suffix)

    def __str__(self) -> str:
        return f"{self.major}.{self.minor}.{self.patch}{self.suffix}"


def parse_star_version(raw: str) -> StarVersion:
    m = _STAR_VERSION_PATTERN.search(raw.strip())
    if m is None:
        raise ValueError(f"Cannot parse STAR version from: {raw!r}")
    return StarVersion(
        major=int(m.group(1)),
        minor=int(m.group(2)),
        patch=int(m.group(3)),
        suffix=m.group(4) or "",
    )


def resolve_star_binary(star_bin: Path | None) -> Path:
    if star_bin is not None:
        if not star_bin.is_file():
            raise FileNotFoundError(f"Specified STAR binary not found: {star_bin}")
        return star_bin
    found = shutil.which("STAR")
    if found is None:
        raise FileNotFoundError(
            "STAR binary not found in PATH. Install STAR or specify --star-bin.",
        )
    return Path(found)


def check_star_version(star_bin: Path) -> StarVersion:
    result = subprocess.run(
        [str(star_bin), "--version"],
        capture_output=True,
        text=True,
        check=False,
    )
    raw = result.stdout.strip() or result.stderr.strip()
    if not raw:
        raise RuntimeError(f"STAR --version produced no output: {star_bin}")

    version = parse_star_version(raw)
    if version.as_tuple() < _STAR_MIN_VERSION:
        raise RuntimeError(
            f"STAR >= {'.'.join(str(x) for x in _STAR_MIN_VERSION[:3])}"
            f"{_STAR_MIN_VERSION[3]} required, found {version}",
        )
    logger.info("STAR %s (%s)", version, star_bin)
    return version
