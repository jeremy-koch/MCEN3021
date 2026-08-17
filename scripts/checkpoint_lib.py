"""Single source of truth for gated checkpoints: emits page HTML + instructor bank."""
import re

def _cp_div(prefix, i, cp, is_last):
    target = "" if is_last else f'\n     data-target="{prefix}-s{i+2}"'
    return f'''```{{=html}}
<div class="fl-checkpoint"
     data-prompt="{cp['prompt']}"
     data-answer="{cp['answer']}" data-tol="{cp.get('tol', 0.03)}" data-units="{cp['units']}"
     data-hint="{cp['hint']}"
     data-solution="{cp['sol']}"{target}></div>
```'''

def emit_block(stem, spec):
    """The '## Earn it' section for one module."""
    prefix = spec['prefix']
    steps = spec['steps']
    parts = [f"## Earn it — {len(steps)} quick gates\n",
             spec['intro'] + " Ungraded; hints after a miss; **Show me** after two. "
             "Arguing about a gate with a classmate is the assignment, not a violation.\n",
             '```{=html}\n<div class="fl-reveal-all"></div>\n```\n']
    for i, cp in enumerate(steps):
        is_last = (i == len(steps) - 1)
        parts.append(_cp_div(prefix, i, cp, is_last))
        if not is_last:
            parts.append(f"\n::: {{#{prefix}-s{i+2} .fl-hidden}}\n{cp['after']}\n")
    # close the nested reveals: each non-last step opened one div
    closing = spec.get('closing', '')
    if closing:
        parts.append("\n" + closing + "\n")
    parts.append(":::\n" * (len(steps) - 1))
    return "\n".join(parts)

def patch_module(stem, spec, moddir='modules'):
    path = f"{moddir}/{stem}.qmd"
    src = open(path).read()
    assert '## Earn it' not in src, f"{stem}: already patched"
    # widget include after frontmatter
    if '_checkpoint-widget.qmd' not in src:
        m = re.match(r'^(---\n.*?\n---\n)', src, re.S)
        src = src[:m.end()] + "\n{{< include _checkpoint-widget.qmd >}}\n" + src[m.end():]
    block = emit_block(stem, spec)
    anchor = "## Check yourself"
    assert anchor in src, f"{stem}: no anchor"
    src = src.replace(anchor, block + "\n" + anchor, 1)
    open(path, 'w').write(src)
    return len(spec['steps'])

def bank_rows(stem, spec):
    rows = []
    for i, cp in enumerate(spec['steps'], 1):
        rows.append((stem, i, cp['answer'], cp.get('tol', 0.03), cp['units'],
                     cp.get('wrongs', ''), cp.get('move', '')))
    return rows
