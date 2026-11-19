# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright 2026 Junna Kawasaki & Jumpei Ito

"""Resolve external CLI tool paths from a TOML config file or PATH."""

from __future__ import annotations

import logging
import shutil
import stat
import tomllib
from collections.abc import Sequence
from pathlib import Path

logger = logging.getLogger(__name__)

_TOOL_ALIASES: dict[str, tuple[str, ...]] = {
    "iqtree": ("iqtree3", "iqtree"),
}

_TOOL_RECOMMENDED_NAME: dict[str, str] = {
    "iqtree": "iqtree3",
}

_CONFIG_TOOLS_SECTION = "tools"


def load_tool_config(config_path: Path | None) -> dict[str, str]:
    """Read ``[tools]`` from *config_path* and return a name → path mapping.

    Returns an empty dict when *config_path* is ``None``.  Raises
    ``FileNotFoundError`` when *config_path* is given but does not exist.
    """
    if config_path is None:
        return {}
    if not config_path.exists():
        raise FileNotFoundError(f"Tool config file not found: {config_path}")

    with config_path.open("rb") as fh:
        data = tomllib.load(fh)

    tools_section: dict[str, str] = {}
    raw = data.get(_CONFIG_TOOLS_SECTION, {})
    for key, value in raw.items():
        if not isinstance(value, str):
            raise TypeError(
                f"{config_path}: tools.{key} must be a string, got {type(value).__name__}",
            )
        tools_section[key] = value
    return tools_section


def _resolve_one(name: str, config: dict[str, str]) -> str:
    if name in config:
        path = config[name]
        p = Path(path)
        if not p.is_file():
            raise FileNotFoundError(
                f"Tool path from config does not exist: {name} = {path}",
            )
        if not p.stat().st_mode & (stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH):
            raise FileNotFoundError(
                f"Tool path from config is not executable: {name} = {path}",
            )
        return path

    aliases = _TOOL_ALIASES.get(name, (name,))
    for alias in aliases:
        found = shutil.which(alias)
        if found is not None:
            if name in _TOOL_RECOMMENDED_NAME and alias != _TOOL_RECOMMENDED_NAME[name]:
                logger.warning(
                    "%s resolved to %s; %s is recommended",
                    name,
                    found,
                    _TOOL_RECOMMENDED_NAME[name],
                )
            return alias
    raise FileNotFoundError(f"{name} is not installed or not on PATH")


def resolve_tools(
    tools: Sequence[str],
    config: dict[str, str] | None = None,
) -> dict[str, str]:
    """Resolve each tool name to an executable path.

    Returns a ``{canonical_name: executable}`` mapping.  Raises
    ``FileNotFoundError`` on the first tool that cannot be found.
    """
    effective = config if config is not None else {}
    return {name: _resolve_one(name, effective) for name in tools}
