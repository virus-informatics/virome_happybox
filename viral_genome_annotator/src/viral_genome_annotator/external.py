# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: Copyright 2026 Junna Kawasaki & Jumpei Ito

from __future__ import annotations

import logging
import shutil

logger = logging.getLogger(__name__)


def require_command(name: str) -> str:
    """Check that a command exists in PATH. Returns the resolved path."""
    found = shutil.which(name)
    if found is None:
        raise FileNotFoundError(f"Required command not found in PATH: {name}")
    logger.info("Found %s: %s", name, found)
    return found
