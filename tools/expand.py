#!/usr/bin/env python3
from __future__ import annotations

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
EXPANDED = ROOT / "expanded"

LOCAL_INCLUDE_RE = re.compile(r'^\s*#include\s+"([^"]+)"')


def guard_name(rel: Path) -> str:
    stem = re.sub(r"[^0-9A-Za-z]+", "_", rel.as_posix()).upper()
    return f"LIBRARY_{stem}"


def strip_outer_guard(lines: list[str], guard: str) -> list[str]:
    first = 0
    while first < len(lines) and lines[first].strip() == "":
        first += 1
    if first + 1 >= len(lines):
        return lines
    if lines[first].strip() != f"#ifndef {guard}":
        return lines

    second = first + 1
    while second < len(lines) and lines[second].strip() == "":
        second += 1
    if second >= len(lines) or lines[second].strip() != f"#define {guard}":
        return lines

    last = len(lines) - 1
    while last >= 0 and lines[last].strip() == "":
        last -= 1
    if last <= second or not lines[last].lstrip().startswith("#endif"):
        return lines

    return lines[:first] + lines[second + 1:last] + lines[last + 1:]


def expand_file(path: Path, emitted: set[Path], stack: list[Path]) -> str:
    path = path.resolve()
    if path in emitted:
        return ""
    if path in stack:
        cycle = " -> ".join(p.relative_to(ROOT).as_posix() for p in stack + [path])
        raise RuntimeError(f"cyclic include: {cycle}")
    if not path.is_relative_to(SRC):
        raise RuntimeError(f"local include is outside src: {path}")

    rel = path.relative_to(SRC)
    guard = guard_name(rel)
    lines = path.read_text().splitlines()
    lines = strip_outer_guard(lines, guard)

    body: list[str] = []
    for line in lines:
        stripped = line.strip()
        if stripped == "#pragma once":
            continue

        match = LOCAL_INCLUDE_RE.match(line)
        if match:
            child = (path.parent / match.group(1)).resolve()
            if child.is_file() and child.is_relative_to(SRC):
                body.append(expand_file(child, emitted, stack + [path]))
                continue

        body.append(line)

    emitted.add(path)
    rel_name = (Path("src") / rel).as_posix()
    content = "\n".join(body).rstrip()
    if content:
        content = f"#line 1 \"{rel_name}\"\n{content}\n"
    else:
        content = f"#line 1 \"{rel_name}\"\n"

    return f"#ifndef {guard}\n#define {guard}\n\n{content}\n#endif  // {guard}\n"


def main() -> None:
    for source in sorted(SRC.rglob("*.hpp")):
        emitted: set[Path] = set()
        rel = source.relative_to(SRC)
        out = EXPANDED / rel
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(expand_file(source, emitted, []))


if __name__ == "__main__":
    main()
