# fluidtools — the MCEN 3021 Fluid Mechanics toolbox & course site

**Supersedes `pipeflow-toolbox`** (the pipeflow module lives on here,
unchanged API, at `fluidtools.pipeflow`). One repo containing:

- `src/fluidtools/` — the installable package: nine teaching tools, each a
  simplified version of industry software, each validated against textbook
  values (`tests/`, 55 checks). Includes the eight-machine `pump_catalog()`
  spanning the specific-speed map (radial → mixed → axial) for the term project.
- `*.qmd`, `modules/` — the course website (Quarto Live + Pyodide: live
  Python in the browser). **All 28 module readings across all 10 units are
  complete**, each in the same anatomy (big idea → number-in-charge callout →
  Think First → reading with collapsible derivations → live Pyodide/OJS demo →
  misconception callout → self-checks → bring-to-class). Units 9-10 (open
  channel, turbomachinery) were written last, with turbomachinery given the
  deeper treatment per instructor request.
- `homework/` — six verify-then-apply notebooks (one per unit-pair) plus
  `homework/solutions/` (worked, executed instructor keys — keep private).
  See `homework.qmd` for the index. `homework/solutions/` also holds
  `project-INSTRUCTOR-NOTES.ipynb` — the validated design space, intended
  answer, and the three traps of the term project.
- **Checkpoints, site-wide:** every module ends with an "Earn it" gated
  walk (2-4 numeric prediction gates), and every unit has a worked+faded
  example page (`modules/*-worked.qmd`) for post-class review. ~100 gates
  total; consolidated answer bank with wrong-number diagnoses in
  `homework/solutions/checkpoint-bank-all-units.md` (gitignored).
- **Checkpoint widget (originally the Unit 4 pilot):** `modules/_checkpoint-widget.qmd` is a
  reusable gated-derivation component (numeric prediction unlocks the next
  step; hint after one miss, Show-me after two; reveal-all button +
  print-safe). Piloted in `bernoulli.qmd` (the derivation, 4 gates) and
  `bernoulli-worked.qmd` (worked example + faded twin, 4 gates). Answers &
  wrong-number diagnoses: `homework/solutions/checkpoint-bank-bernoulli.md`.
  To gate another page: include the widget, wrap steps in
  `::: {#id .fl-hidden}` divs, add a `div.fl-checkpoint` with data-attrs.
- `syllabus.qmd` (grading + 28-meeting schedule), `project.qmd` ("The Pump
  Station" term design project), `homework.qmd` (assignment index).
- `scripts/make_bundle.py` — pre-render hook that concatenates the
  submodules into `fluidtools_bundle.py` for the browser.
  `scripts/make_homework.py`, `scripts/make_solutions.py` — regenerate the
  notebooks.

## Instructor setup
1. Create a GitHub repo named `fluidtools`, push this directory.
2. Find-and-replace `YOURUSERNAME` (`_quarto.yml`, `index.qmd`, `tools.qmd`).
3. Tag a release: `git tag v2.1.0 && git push --tags` (notebooks pin to tags).
4. Repo Settings → Pages → Source: `gh-pages` branch. The included GitHub
   Action renders and publishes on every push to `main`.
5. Students install with
   `%pip install git+https://github.com/YOURUSERNAME/fluidtools@v2.1.0`

## Local development
```bash
pip install -e . && python -m pytest tests/   # 55 tests
quarto preview                                 # live site preview
jupyter lab homework/                          # the assignment notebooks
```

Keep solutions and the ConcepTest bank (`fluids-conceptests`) in private
repos — this one is public to students.
