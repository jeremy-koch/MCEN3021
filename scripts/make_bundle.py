"""Pre-render: concatenate fluidtools submodules into one browser-loadable file."""
import pathlib
SRC = pathlib.Path("src/fluidtools")
ORDER = ["hydrostat", "momentum", "similitude", "pipeflow", "profiles",
         "extflow", "flowmeter", "channel", "turbo"]
parts = [f"# fluidtools browser bundle (auto-generated; do not edit)\n"]
for name in ORDER:
    parts.append(f"\n# {'='*70}\n# fluidtools.{name}\n# {'='*70}\n")
    parts.append((SRC / f"{name}.py").read_text())
pathlib.Path("fluidtools_bundle.py").write_text("\n".join(parts))
print("bundle written")
